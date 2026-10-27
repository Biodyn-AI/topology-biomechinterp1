"""E4a — H3 edge relation (pre-registered primary, CONVENTION3 §4) + edge sensitivity analyses, one context.

PRIMARY (exactly per CONVENTION3 §4 H3):
  pairs      cache/<ctx>/h3_pairs.npz (TRRUST edges within V vs target-matched TRRUST-negative candidate pairs)
  scorer     signed cosine of the V-centred full 512-d embedding of TF and target, per layer (12 layers)
  statistic  within-TF stratified AUROC (stats3.StratifiedPairs.auroc)
  detection  one within-TF label permutation applied jointly to all 12 layers, N = 10,000, max-T over the
             12 layers (stats3.max_t_adjust, tail 'greater': H3 asks whether the relation is recoverable,
             i.e. AUROC > 0.5; the two-sided max-T p is reported alongside as a sensitivity column)
  intervals  two-way TF x target pigeonhole cluster bootstrap, B = 2,000, per layer: percentile 95% CI and
             one-sided 95% upper bound (stats3.one_sided_upper)
  comparator |co-expression| (|Pearson r|) on the same pairs, on expr_extraction.npy (primary comparator, the
             exact 256/192 cells behind the embedding) and on expr.npy (3,000-cell tier); paired delta
             (layer L minus |co-expression|) on the SAME bootstrap resamples.
SENSITIVITY (no pass/fail): |cos|, rank-k cosine, SV bands, macro per-TF AUROC, pooled AUROC on the legacy
all-V non-edges, shortcut rows, strata, ChIP-only independent positives, TF-detection restriction,
PAX5->CD79A worked example. Label-noise and the decision rule are applied in e4a_aggregate.py.

Checkpoints in work/<ctx>/ (scores per layer, permutation, bootstrap per layer, sensitivity blocks); a
re-run skips every finished block. Usage: python e4a_main.py <ctx> [<ctx> ...]
"""
from __future__ import annotations
import json, os, sys, time, zlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3)); sys.path.insert(0, str(HERE))
import stats3 as S          # noqa: E402
import contexts as X        # noqa: E402
import e4a_lib as L         # noqa: E402

EMB_LOCAL = Path(os.path.expanduser("~/biomechinterp-local/subproject41_r3/emb"))
CACHE = R3 / "cache"
N_PERM = 10_000
B = 2_000
B_SECONDARY = 1_000            # rank-k / band scorers, strata on large contexts, ChIP outside lung pair
RANKS = [1, 2, 3, 5, 7, 10, 16, 32, 64, 128, 256]
BANDS = {"SV2-SV4": (1, 4), "SV5-SV7": (4, 7), "SV2-SV7": (1, 7)}
NL = 12


def log(ctx, *a):
    print(f"[{time.strftime('%H:%M:%S')}] [{ctx}]", *a, flush=True)


def seed_of(ctx):
    return 20261026 + zlib.crc32(ctx.encode()) % 100000


def jdump(obj, p):
    tmp = Path(str(p) + ".tmp"); tmp.write_text(json.dumps(obj, indent=1, default=float)); tmp.replace(p)


def load_ctx(ctx):
    d = CACHE / ctx
    D = dict(genes=json.loads((d / "obs_genes.json").read_text()), obs_idx=np.load(d / "obs_idx.npy"),
             ncells=np.load(d / "ncells.npy"), h3=dict(np.load(d / "h3_pairs.npz")),
             leg=dict(np.load(d / "h3_legacy_pairs.npz")),
             chip=dict(np.load(d / "h3_chip_pairs.npz")) if (d / "h3_chip_pairs.npz").exists() else None,
             edges_trrust=np.load(d / "edges_trrust.npy"))
    D["emb"] = np.load(EMB_LOCAL / f"{ctx}.npy", mmap_mode="r")
    assert D["emb"].shape[0] == NL
    return D


# ===================================================================================== stage A: scores
def stage_scores(ctx, D, W):
    sp = D["sp"]; leg = D["sp_leg"]; ch = D["sp_chip"]
    for l in range(NL):
        f = W / f"scores_L{l:02d}.npz"
        if f.exists(): continue
        Xl = np.asarray(D["emb"][l], np.float64)[D["obs_idx"]]
        assert np.all(np.linalg.norm(Xl, axis=1) > 0), "never-observed row in V"
        Xc = Xl - Xl.mean(0, keepdims=True); En = S.l2norm_rows(Xc)
        out = dict(prim=L.rowcos(En, sp.tf, sp.tg), leg=L.rowcos(En, leg.tf, leg.tg))
        if ch is not None: out["chip"] = L.rowcos(En, ch.tf, ch.tg)
        US, Sv, _ = S.svd_observed(D["emb"], l, D["obs_idx"])
        # consistency: cosines of the full US rows equal the centred-embedding cosines (rotation)
        chk = L.rowcos(S.l2norm_rows(US), sp.tf[:2000], sp.tg[:2000])
        assert np.max(np.abs(chk - out["prim"][:2000])) < 1e-8
        for k in RANKS:
            out[f"rk{k}"] = np.round(L.rowcos(S.l2norm_rows(US[:, :k]), sp.tf, sp.tg), 10)
        for nm, (a, b) in BANDS.items():
            out[f"band_{nm}"] = np.round(L.rowcos(S.l2norm_rows(US[:, a:b]), sp.tf, sp.tg), 10)
        out["sv"] = Sv
        np.savez(f, **out)
        log(ctx, f"scores L{l} done")
    f = W / "coexp.npz"
    if not f.exists():
        d = CACHE / ctx
        Ex = np.load(d / "expr_extraction.npy"); E3 = np.load(d / "expr.npy")
        out = {}
        for tag, M in (("cx_ext", Ex), ("cx_3k", E3)):
            out[f"{tag}"] = L.abs_corr_pairs(M, sp.tf, sp.tg)
            out[f"{tag}_leg"] = L.abs_corr_pairs(M, leg.tf, leg.tg)
            if ch is not None: out[f"{tag}_chip"] = L.abs_corr_pairs(M, ch.tf, ch.tg)
        out["tf_detfrac_ext"] = (Ex > 0).mean(0)          # per observed gene: fraction of extraction cells
        # shortcut scorers (target-level only)
        indeg = np.bincount(D["edges_trrust"][:, 1], minlength=len(D["genes"])).astype(float)
        out["indeg"] = indeg
        np.savez(f, **out)
        log(ctx, "co-expression done")


def scores_layer(W, l):
    return dict(np.load(W / f"scores_L{l:02d}.npz"))


# ===================================================================================== stage B: permutation
def stage_perm(ctx, D, W):
    f = W / "perm.json"
    if f.exists(): return json.loads(f.read_text())
    sp = D["sp"]; fs = L.FastStrat(sp)
    P = np.column_stack([scores_layer(W, l)["prim"] for l in range(NL)])
    obs = np.array([sp.auroc(P[:, l]) for l in range(NL)])            # stats3 directly
    t0 = time.time()
    null, labs = fs.perm_null(P, N_PERM, np.random.default_rng([seed_of(ctx), 1]), return_labels_first=3)
    # verify fast null against stats3 on the first permuted label vectors, every layer
    vmax = max(abs(S.StratifiedPairs(sp.tf, sp.tg, yy).auroc(P[:, l]) - null[r, l])
               for r, yy in enumerate(labs) for l in (0, 6, 11))
    assert vmax < 1e-12, vmax
    floor = S.assert_resolution(N_PERM, 0.05, NL)
    z, p_adj, p_raw = S.max_t_adjust(obs, null, tail="greater")
    z2, p_adj2, p_raw2 = S.max_t_adjust(obs, null, tail="two-sided")
    res = dict(N=N_PERM, seed=[seed_of(ctx), 1], floor=floor, fast_vs_stats3_max_absdiff=float(vmax),
               obs=obs.tolist(), z=z.tolist(), p_maxT=p_adj.tolist(), p_raw=p_raw.tolist(),
               p_maxT_two_sided=p_adj2.tolist(), p_raw_two_sided=p_raw2.tolist(),
               null_summary=[L.summarise_null(null[:, l]) for l in range(NL)],
               null_max_z_summary=L.summarise_null(((null - null.mean(0)) / null.std(0, ddof=1)).max(1)),
               seconds=time.time() - t0)
    jdump(res, f); log(ctx, f"perm done: min p_maxT={min(p_adj):.4g}")
    return res


# ===================================================================================== stage C: bootstrap
def boot_block(fs, s, wts, Bn=None):
    if Bn is not None:
        wts = (wts[0][:Bn], None if wts[1] is None else wts[1][:Bn])
    return fs.auc(s, wts)


def rowdict(name, layer, est, boot, ref_est=None, ref_boot=None, ref2_est=None, ref2_boot=None):
    r = dict(scorer=name, layer=layer, auroc=float(est))
    c = L.ci_row(est, boot); r.update(ci_lo=c["ci_lo"], ci_hi=c["ci_hi"], ub95=c["ub95"], boot_sd=c["boot_sd"],
                                      B=int(np.isfinite(boot).sum()))
    for tag, re_, rb in (("cxext", ref_est, ref_boot), ("cx3k", ref2_est, ref2_boot)):
        if re_ is None: continue
        dboot = boot - rb[:boot.size]
        c = L.ci_row(est - re_, dboot)
        r.update({f"delta_{tag}": c["est"], f"delta_{tag}_ci_lo": c["ci_lo"], f"delta_{tag}_ci_hi": c["ci_hi"],
                  f"delta_{tag}_ub95": c["ub95"], f"delta_{tag}_sd": c["boot_sd"]})
    return r


def stage_boot(ctx, D, W):
    sp = D["sp"]; fs = L.FastStrat(sp)
    wts = fs.boot_weights(B, np.random.default_rng([seed_of(ctx), 2]), two_way=True)
    wts1 = fs.boot_weights(B, np.random.default_rng([seed_of(ctx), 3]), two_way=False)
    cx = dict(np.load(W / "coexp.npz"))
    fc = W / "boot_coexp.npz"
    if fc.exists():
        bc = dict(np.load(fc))
    else:
        bc = {}
        for tag in ("cx_ext", "cx_3k"):
            bc[tag] = fs.auc(cx[tag], wts); bc[tag + "_tfonly"] = fs.auc(cx[tag], wts1)
            bc[tag + "_est"] = np.array([sp.auroc(cx[tag])])
        for tag in ("indeg",):
            s = cx["indeg"][sp.tg]; bc["indeg"] = fs.auc(s, wts); bc["indeg_est"] = np.array([sp.auroc(s)])
        s = np.ones(sp.tg.size)          # 'target is a TRRUST target of anyone in V': all h3 targets are
        is_t = np.zeros(len(D["genes"]), bool); is_t[D["edges_trrust"][:, 1]] = True
        s = is_t[sp.tg].astype(float); bc["anytarget"] = fs.auc(s, wts); bc["anytarget_est"] = np.array([sp.auroc(s)])
        bc["anytarget_frac_pos"] = np.array([s[sp.y].mean()]); bc["anytarget_frac_neg"] = np.array([s[~sp.y].mean()])
        np.savez(fc, **bc)
    ce, c3 = float(bc["cx_ext_est"][0]), float(bc["cx_3k_est"][0])
    rows_cx = []
    for tag, nm in (("cx_ext", "|co-expression| extraction cells"), ("cx_3k", "|co-expression| 3,000-cell tier")):
        r = rowdict(nm, -1, bc[tag + "_est"][0], bc[tag])
        c1 = L.ci_row(bc[tag + "_est"][0], bc[tag + "_tfonly"]); r.update(tfonly_ci_lo=c1["ci_lo"], tfonly_ci_hi=c1["ci_hi"],
                                                                         tfonly_ub95=c1["ub95"])
        rows_cx.append(r)
    rows_cx.append(rowdict("shortcut: target TRRUST in-degree", -1, bc["indeg_est"][0], bc["indeg"]))
    r = rowdict("shortcut: target is a TRRUST target of any TF", -1, bc["anytarget_est"][0], bc["anytarget"])
    r.update(frac_pos=float(bc["anytarget_frac_pos"][0]), frac_neg=float(bc["anytarget_frac_neg"][0])); rows_cx.append(r)
    jdump(rows_cx, W / "rows_coexp.json")
    for l in range(NL):
        f = W / f"rows_L{l:02d}.json"
        if f.exists(): continue
        sc = scores_layer(W, l); rows = []
        s = sc["prim"]; est = sp.auroc(s); bt = fs.auc(s, wts); bt1 = fs.auc(s, wts1)
        r = rowdict("signed cosine full-512 (PRIMARY)", l, est, bt, ce, bc["cx_ext"], c3, bc["cx_3k"])
        c1 = L.ci_row(est, bt1); r.update(tfonly_ci_lo=c1["ci_lo"], tfonly_ci_hi=c1["ci_hi"], tfonly_ub95=c1["ub95"])
        r["macro_auroc"] = sp.macro_auroc(s)
        rows.append(r)
        s = np.abs(sc["prim"]); rows.append(rowdict("|cos| full-512", l, sp.auroc(s), fs.auc(s, wts), ce, bc["cx_ext"], c3, bc["cx_3k"]))
        for k in RANKS:
            s = sc[f"rk{k}"]
            rows.append(rowdict(f"rank-{k} cosine", l, sp.auroc(s), boot_block(fs, s, wts, B_SECONDARY), ce, bc["cx_ext"], c3, bc["cx_3k"]))
        for nm in BANDS:
            s = sc[f"band_{nm}"]
            rows.append(rowdict(f"band {nm} (nominal) cosine", l, sp.auroc(s), boot_block(fs, s, wts, B_SECONDARY), ce, bc["cx_ext"], c3, bc["cx_3k"]))
        jdump(rows, f); log(ctx, f"bootstrap L{l} done (primary AUROC {est:.4f})")


# ===================================================================================== stage D: sensitivity
def macro_block(sp, s, rng, Bm=2000):
    """macro per-TF AUROC (stats3.macro_auroc) with a TF bootstrap of the per-TF AUROCs."""
    vals = []
    for c in range(sp.tfs.size):
        m = sp.tf_code == c
        if (m & sp.y).sum() >= 3: vals.append(S.auroc(s[m & sp.y], s[m & ~sp.y]))
    vals = np.array(vals); est = sp.macro_auroc(s)
    assert abs(est - vals.mean()) < 1e-12
    bt = np.array([vals[rng.integers(0, vals.size, vals.size)].mean() for _ in range(Bm)])
    c = L.ci_row(est, bt)
    return dict(macro_auroc=est, ci_lo=c["ci_lo"], ci_hi=c["ci_hi"], ub95=c["ub95"], n_tf=int(vals.size),
                frac_tf_above_half=float(np.mean(vals > 0.5)))


def stage_sens(ctx, D, W):
    sp = D["sp"]; cx = dict(np.load(W / "coexp.npz")); seed = seed_of(ctx)
    # ---- macro
    f = W / "sens_macro.json"
    if not f.exists():
        rows = []
        for l in range(NL):
            r = macro_block(sp, scores_layer(W, l)["prim"], np.random.default_rng([seed, 10, l])); r.update(scorer="signed cosine full-512", layer=l); rows.append(r)
        for tag in ("cx_ext", "cx_3k"):
            r = macro_block(sp, cx[tag], np.random.default_rng([seed, 11])); r.update(scorer=tag, layer=-1); rows.append(r)
        jdump(rows, f); log(ctx, "macro done")
    # ---- legacy all-V non-edges: pooled (legacy definition) and within-TF stratified
    f = W / "sens_legacy_pooled.json"
    if not f.exists():
        lg = D["sp_leg"]; fp = L.FastStrat(lg, "pooled"); ft = L.FastStrat(lg)
        wts = fp.boot_weights(B, np.random.default_rng([seed, 20]))
        rows = []
        is_t = np.zeros(len(D["genes"]), bool); is_t[D["edges_trrust"][:, 1]] = True
        series = [(f"signed cosine full-512", l, scores_layer(W, l)["leg"]) for l in range(NL)]
        series += [("|co-expression| extraction cells", -1, cx["cx_ext_leg"]), ("|co-expression| 3,000-cell tier", -1, cx["cx_3k_leg"]),
                   ("shortcut: target is a TRRUST target of any TF", -1, is_t[lg.tg].astype(float)),
                   ("shortcut: target TRRUST in-degree", -1, cx["indeg"][lg.tg])]
        ce_p = S.auroc(cx["cx_ext_leg"][lg.y], cx["cx_ext_leg"][~lg.y]); ce_b = fp.auc(cx["cx_ext_leg"], wts)
        for nm, l, s in series:
            est = S.auroc(s[lg.y], s[~lg.y]); bt = fp.auc(s, wts)
            assert abs(est - fp.auc(s)[0]) < 1e-12
            r = rowdict(nm, l, est, bt, ce_p, ce_b); r["metric"] = "pooled AUROC, all-V non-edges (legacy)"
            r["stratified_auroc_same_pairs"] = lg.auroc(s); rows.append(r)
        jdump(dict(rows=rows, n_pos=int(lg.y.sum()), n_neg=int((~lg.y).sum()), n_tf=int(lg.tfs.size)), f); log(ctx, "legacy pooled done")
    # ---- strata (positives restricted, every negative of the TF kept) and TF-detection restriction
    f = W / "sens_strata.json"
    if not f.exists():
        h3 = D["h3"]; base_keep = D["h3_keep"]
        y = h3["y"][base_keep].astype(bool); npm = h3["n_pmid"][base_keep]; mode = h3["mode"][base_keep]
        det = cx["tf_detfrac_ext"][sp.tf]
        defs = {"PMID count = 1": (~y) | (npm == 1), "PMID count >= 2": (~y) | (npm >= 2),
                "mode Activation": (~y) | (mode == "Activation"), "mode Repression": (~y) | (mode == "Repression"),
                "TF detected in >= 10% of extraction cells": det >= 0.10,
                "TF detected in < 10% of extraction cells": det < 0.10}
        out = []
        for sname, m in defs.items():
            if m.sum() == 0: out.append(dict(stratum=sname, n_pos=0)); continue
            sub = S.StratifiedPairs(sp.tf[m], sp.tg[m], sp.y[m])
            if sub.tf.size == 0 or sub.y.sum() == 0: out.append(dict(stratum=sname, n_pos=0)); continue
            # map subset rows back to the full arrays (StratifiedPairs drops TFs without both labels)
            idx_m = np.flatnonzero(m); keep2 = np.isin(sp.tf[m], sub.tfs)
            ii = idx_m[keep2]; assert np.array_equal(sp.tg[ii], sub.tg)
            fs = L.FastStrat(sub); Bs = B if sub.tf.size < 20000 else B_SECONDARY
            wts = fs.boot_weights(Bs, np.random.default_rng([seed, 30, zlib.crc32(sname.encode())]))
            ce = sub.auroc(cx["cx_ext"][ii]); cb = fs.auc(cx["cx_ext"][ii], wts)
            c3 = sub.auroc(cx["cx_3k"][ii]); c3b = fs.auc(cx["cx_3k"][ii], wts)
            rws = [rowdict("|co-expression| extraction cells", -1, ce, cb)]
            for l in range(NL):
                s = scores_layer(W, l)["prim"][ii]
                rws.append(rowdict("signed cosine full-512", l, sub.auroc(s), fs.auc(s, wts), ce, cb, c3, c3b))
            for r in rws: r.update(stratum=sname, n_pos=int(sub.y.sum()), n_neg=int((~sub.y).sum()), n_tf=int(sub.tfs.size))
            out += rws
        jdump(out, f); log(ctx, "strata done")
    # ---- ChIP-only independent positives
    f = W / "sens_chip.json"
    if not f.exists() and D["sp_chip"] is not None:
        ch = D["sp_chip"]; fs = L.FastStrat(ch)
        Bc = B if ctx in ("lung", "lung_external") or ch.tf.size < 50000 else B_SECONDARY
        wts = fs.boot_weights(Bc, np.random.default_rng([seed, 40]))
        ce = ch.auroc(cx["cx_ext_chip"]); cb = fs.auc(cx["cx_ext_chip"], wts)
        c3 = ch.auroc(cx["cx_3k_chip"]); c3b = fs.auc(cx["cx_3k_chip"], wts)
        rows = [rowdict("|co-expression| extraction cells", -1, ce, cb), rowdict("|co-expression| 3,000-cell tier", -1, c3, c3b)]
        for l in range(NL):
            s = scores_layer(W, l)["chip"]
            rows.append(rowdict("signed cosine full-512", l, ch.auroc(s), fs.auc(s, wts), ce, cb, c3, c3b))
            log(ctx, f"chip L{l} {rows[-1]['auroc']:.4f}")
        jdump(dict(rows=rows, n_pos=int(ch.y.sum()), n_neg=int((~ch.y).sum()), n_tf=int(ch.tfs.size), B=Bc), f)
        log(ctx, "chip done")
    # ---- PAX5 -> CD79A worked example
    f = W / "sens_example.json"
    if not f.exists():
        g = D["genes"]; pos = {x: i for i, x in enumerate(g)}; ex = {}
        if "PAX5" in pos and "CD79A" in pos:
            a, b = pos["PAX5"], pos["CD79A"]
            m = sp.tf == a; hit = np.flatnonzero(m & (sp.tg == b))
            ex = dict(PAX5_in_benchmark=bool(m.any()), pair_in_benchmark=bool(hit.size), PAX5_det_frac=float(cx["tf_detfrac_ext"][a]),
                      CD79A_det_frac=float(cx["tf_detfrac_ext"][b]))
            if hit.size:
                h = hit[0]; ex.update(label=int(sp.y[h]), n_candidates=int(m.sum()), n_pos=int((m & sp.y).sum()), n_neg=int((m & ~sp.y).sum()))
                per = []
                for l in range(NL):
                    s = scores_layer(W, l)["prim"]
                    sm = s[m]; rank = int((sm > s[h]).sum() + 1); negs = s[m & ~sp.y]
                    per.append(dict(layer=l, cosine=float(s[h]), rank_among_all_candidates=rank,
                                    frac_negatives_below=float(np.mean(negs < s[h])), pax5_tf_auroc=S.auroc(s[m & sp.y], negs)))
                for tag in ("cx_ext", "cx_3k"):
                    s = cx[tag]; sm = s[m]; negs = s[m & ~sp.y]
                    per.append(dict(layer=tag, cosine=float(s[h]), rank_among_all_candidates=int((sm > s[h]).sum() + 1),
                                    frac_negatives_below=float(np.mean(negs < s[h])), pax5_tf_auroc=S.auroc(s[m & sp.y], negs)))
                ex["per_layer"] = per
        jdump(ex, f)


# ===================================================================================== driver
def run(ctx):
    W = HERE / "work" / ctx; W.mkdir(parents=True, exist_ok=True)
    if (W / "DONE").exists(): log(ctx, "already done"); return
    t0 = time.time(); D = load_ctx(ctx); h3 = D["h3"]
    D["sp"] = S.StratifiedPairs(h3["tf"], h3["tg"], h3["y"])
    # rows kept by StratifiedPairs (for strata labels)
    keep = np.isin(h3["tf"], D["sp"].tfs); D["h3_keep"] = keep
    assert np.array_equal(h3["tg"][keep], D["sp"].tg)
    D["sp_leg"] = S.StratifiedPairs(D["leg"]["tf"], D["leg"]["tg"], D["leg"]["y"])
    D["sp_chip"] = S.StratifiedPairs(D["chip"]["tf"], D["chip"]["tg"], D["chip"]["y"]) if D["chip"] is not None else None
    meta_f = W / "meta.json"
    if not meta_f.exists():
        sp = D["sp"]; s0 = None
        stage_scores(ctx, D, W)
        s0 = scores_layer(W, 0)["prim"]
        ver = L.verify_against_stats3(sp, s0, seed_of(ctx), n_boot=10, n_perm=5)
        jdump(dict(context=ctx, n_pairs=int(sp.tf.size), n_pos=int(sp.y.sum()), n_neg=int((~sp.y).sum()),
                   n_tf=int(sp.tfs.size), n_targets=int(sp.tgs.size), n_pairs_file=int(h3["tf"].size),
                   frac_neg_fully_matched=float(np.mean(h3["match_level"][keep][~sp.y] == 0)),
                   n_cells_extraction=int(np.load(CACHE / ctx / "expr_extraction.npy", mmap_mode="r").shape[0]),
                   seed=seed_of(ctx), fast_vs_stats3=ver,
                   emb=str(EMB_LOCAL / f"{ctx}.npy"), emb_source=str(X.CONTEXTS[ctx]["emb"])), meta_f)
        log(ctx, f"verification {ver}")
    stage_scores(ctx, D, W)
    stage_perm(ctx, D, W)
    stage_boot(ctx, D, W)
    stage_sens(ctx, D, W)
    (W / "DONE").write_text(f"{time.time() - t0:.0f}s\n"); log(ctx, f"DONE in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    for c in sys.argv[1:]:
        run(c)
