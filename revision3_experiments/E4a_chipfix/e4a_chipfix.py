"""E4a_chipfix -- corrected ChIP-seq-only independent-positive benchmark and its E4a sensitivity re-run.

AUDIT FINDING. prep_pairs.main() built cache/<ctx>/h3_chip_pairs.npz with
    build(ch_v, genes, pos, anyreg - set(ch_v), indeg_c, det_dec, rng, other_c)
so every ChIP-only positive pair (and its reverse) was removed from the exclusion set: a factor's own
ChIP-only targets stayed eligible as its "candidate non-edges" (in lung 52,908 of 557,447 candidate rows are
also positive rows of the same TF). CONVENTION3 section 3-4 requires every pair of any regulatory resource
(anyreg = TRRUST U DoRothEA main U DoRothEA ChIP, all levels) to be excluded from every candidate set.

WHAT THIS SCRIPT DOES (nothing outside E4a_chipfix/ is written; prep_pairs.py, stats3.py, contexts.py and
../cache/ are only read):
  build      For every context: replays prep_pairs.main()'s RNG stream exactly (rng = default_rng(SEED) per
             context; the h3_pairs build and the legacy draws consume it first) and verifies the replay
             reproduces ../cache/<ctx>/h3_pairs.npz, h3_legacy_pairs.npz and h3_chip_pairs.npz bit for bit.
             From the SAME generator state as the old ChIP draw it then calls prep_pairs.build with the full
             anyreg (the only change). Writes cache/<ctx>/h3_chip_pairs_fixed.npz (tf, tg, y, match_level; same
             layout as the old file) and build_audit.csv / build_audit.json.
  sens CTX.. E4a's ChIP sensitivity code path on the fixed pairs: signed cosine of the V-centred full 512-d
             embedding (stats3.l2norm_rows, e4a_lib.rowcos) per layer; |Pearson r| on expr_extraction.npy and
             expr.npy rounded to 12 decimals (the verifier's tie fix, e4a verify/fix_ties.py); within-TF
             stratified AUROC (stats3.StratifiedPairs); two-way TF x target pigeonhole bootstrap via
             e4a_lib.FastStrat with seed [seed_of(ctx), 40] (E4a's ChIP stream), B = 2,000 in EVERY context
             (E4a used 1,000 for the four large kidney contexts); paired deltas vs |co-expression| on the same
             resamples (e4a_main.rowdict). Second block: candidates with match_level 0 only (+ all positives),
             seed [seed_of(ctx), 41], as e4a_extra.py does for h3. Third: old pairs with the contaminated
             candidate rows simply deleted (point estimates only, decomposition of old -> new).
  repro CTX  Code-path check: the same functions on the OLD pairs with the OLD settings (unrounded comparator
             outside the lung pair, E4a's B) must reproduce E4a_edges/sens_chip.csv.
  aggregate  sens_chip_fixed.csv (columns identical to E4a_edges/sens_chip.csv), sens_chip_fixed_matchlevel0.csv,
             chip_old_vs_new.csv, chip_match_levels.csv, results.json.
"""
from __future__ import annotations
import copy, json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; R3 = HERE.parent; E4 = R3 / "E4a_edges"
sys.path.insert(0, str(R3)); sys.path.insert(0, str(E4))
import stats3 as S            # noqa: E402
import contexts as X          # noqa: E402
import prep_pairs as P        # noqa: E402
import e4a_lib as L           # noqa: E402
import e4a_main as M          # noqa: E402

OLD_CACHE = R3 / "cache"; NEW_CACHE = HERE / "cache"; W0 = HERE / "work"
B = 2000; NL = 12; DEC = 12
BASIS = "tie-corrected"


def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)
def jl(p): return json.loads(Path(p).read_text())


# ===================================================================================== build
def replay_context(ctx, refs):
    """prep_pairs.main() for one context, verbatim up to the ChIP draw; returns old and fixed ChIP builds."""
    trrust, anyreg, chip_indep, lambert = refs
    d = OLD_CACHE / ctx; rng = np.random.default_rng(P.SEED)
    genes = json.loads((d / "obs_genes.json").read_text()); pos = {g: k for k, g in enumerate(genes)}
    V = set(genes); nobs = np.load(d / "ncells.npy").astype(float)
    log_nobs = np.log(np.maximum(nobs, 1)); det_dec = S.quantile_strata(log_nobs, 10)
    tr_v = sorted(p for p in trrust if p[0] in V and p[1] in V)
    tfs_tr = {a for a, _ in tr_v}
    indeg = {}
    for a, b in tr_v: indeg[b] = indeg.get(b, 0) + 1
    other = {tf: {b for a, b in tr_v if a != tf} for tf in tfs_tr}
    tf, tg, y, lev = P.build(tr_v, genes, pos, anyreg, indeg, det_dec, rng, other)
    old_h3 = np.load(d / "h3_pairs.npz")
    h3_ok = all(np.array_equal(a, old_h3[k]) for a, k in ((tf, "tf"), (tg, "tg"), (y, "y"), (lev, "match_level")))
    leg = []
    for t in sorted(tfs_tr):
        tg_t = [b for a, b in tr_v if a == t]
        if len(tg_t) < 3: continue
        cand = np.array([pos[g] for g in genes if g != t and (t, g) not in anyreg and (g, t) not in anyreg])
        take = rng.choice(cand, min(P.NEG_PER_POS * len(tg_t), cand.size), replace=False)
        leg += [(pos[t], pos[b], 1) for b in tg_t] + [(pos[t], int(k), 0) for k in take]
    leg = np.array(leg, dtype=np.int64)
    old_leg = np.load(d / "h3_legacy_pairs.npz")
    leg_ok = np.array_equal(leg[:, 0], old_leg["tf"]) and np.array_equal(leg[:, 1], old_leg["tg"])
    ch_v = sorted(p for p in chip_indep if p[0] in V and p[1] in V)
    indeg_c = {}
    for a, b in ch_v: indeg_c[b] = indeg_c.get(b, 0) + 1
    other_c = {tf_: {b for a, b in ch_v if a != tf_} | {b for a, b in tr_v if a != tf_} for tf_ in {a for a, _ in ch_v}}
    rng_fixed = copy.deepcopy(rng)                       # identical generator state for both ChIP draws
    state_equal = rng_fixed.bit_generator.state == rng.bit_generator.state
    rc_old = P.build(ch_v, genes, pos, anyreg - set(ch_v), indeg_c, det_dec, rng, other_c)      # as prep_pairs
    rc_new = P.build(ch_v, genes, pos, anyreg, indeg_c, det_dec, rng_fixed, other_c)            # the fix
    old_file = dict(np.load(d / "h3_chip_pairs.npz")) if (d / "h3_chip_pairs.npz").exists() else None
    chip_ok = old_file is not None and rc_old is not None and all(
        np.array_equal(rc_old[i], old_file[k]) for i, k in enumerate(("tf", "tg", "y", "match_level")))
    return dict(genes=genes, ch_v=ch_v, anyreg=anyreg, rc_old=rc_old, rc_new=rc_new, old_file=old_file,
                h3_ok=bool(h3_ok), leg_ok=bool(leg_ok), chip_ok=bool(chip_ok), state_equal=bool(state_equal))


def audit_context(ctx, R):
    genes = R["genes"]; chs = set(R["ch_v"]); anyreg = R["anyreg"]
    o = R["rc_old"]; n = R["rc_new"]
    def pairs(rc, lab): return [(genes[a], genes[b]) for a, b, yy in zip(rc[0], rc[1], rc[2]) if yy == lab]
    po, pn = pairs(o, 1), pairs(n, 1); no, nn = pairs(o, 0), pairs(n, 0)
    so, sn = set(no), set(nn)
    contam_fwd = sum(p in chs for p in no)                       # candidate == a ChIP-only positive of that TF
    contam_rev = sum((b, a) in chs for a, b in no)               # reverse is a ChIP-only positive
    contam_any = sum((p in anyreg) or ((p[1], p[0]) in anyreg) for p in no)
    tfs_o = sorted({a for a, _ in po}); tfs_n = sorted({a for a, _ in pn})
    dropped = sorted(set(tfs_o) - set(tfs_n))
    first_div = None                                             # first TF (sorted order) whose candidates differ
    for t in tfs_o:
        a = [p for p in no if p[0] == t]; b = [p for p in nn if p[0] == t]
        if a != b: first_div = t; break
    lev_o = o[3][o[2] == 0]; lev_n = n[3][n[2] == 0]
    row = dict(context=ctx, replay_h3_pairs_identical=R["h3_ok"], replay_legacy_identical=R["leg_ok"],
               replay_old_chip_identical=R["chip_ok"], rng_state_identical_at_chip_draw=R["state_equal"],
               n_chip_indep_within_V=len(R["ch_v"]),
               old_n_pos=len(po), old_n_neg=len(no), old_n_tf=len(tfs_o),
               new_n_pos=len(pn), new_n_neg=len(nn), new_n_tf=len(tfs_n),
               positives_identical_rows=bool(po == pn), positives_new_subset_of_old=bool(set(pn) <= set(po)),
               n_pos_dropped=len(set(po) - set(pn)), tfs_dropped=";".join(dropped), n_tfs_dropped=len(dropped),
               old_neg_is_own_chip_positive=contam_fwd, old_neg_reverse_is_chip_positive=contam_rev,
               old_neg_in_anyreg_either_direction=contam_any,
               new_neg_in_anyreg_either_direction=sum((p in anyreg) or ((p[1], p[0]) in anyreg) for p in nn),
               neg_common=len(so & sn), neg_only_old=len(so - sn), neg_only_new=len(sn - so),
               neg_only_old_contaminated=sum(((p in anyreg) or ((p[1], p[0]) in anyreg)) for p in so - sn),
               frac_old_neg_retained=len(so & sn) / max(len(so), 1),
               first_tf_with_changed_draw=first_div)
    for k in (0, 1, 2):
        row[f"old_share_level{k}"] = float(np.mean(lev_o == k)); row[f"new_share_level{k}"] = float(np.mean(lev_n == k))
    return row


def cmd_build(ctxs):
    refs = P.load_refs(); rows = []
    for ctx in ctxs:
        t0 = time.time(); R = replay_context(ctx, refs)
        assert R["h3_ok"] and R["leg_ok"] and R["chip_ok"] and R["state_equal"], (ctx, R["h3_ok"], R["leg_ok"], R["chip_ok"])
        n = R["rc_new"]; (NEW_CACHE / ctx).mkdir(parents=True, exist_ok=True)
        np.savez(NEW_CACHE / ctx / "h3_chip_pairs_fixed.npz", tf=n[0], tg=n[1], y=n[2], match_level=n[3])
        row = audit_context(ctx, R); rows.append(row)
        log(ctx, {k: row[k] for k in ("old_n_pos", "new_n_pos", "old_n_neg", "new_n_neg", "old_n_tf", "new_n_tf",
                                      "old_neg_is_own_chip_positive", "neg_common", "positives_identical_rows")}, f"{time.time()-t0:.0f}s")
    df = pd.DataFrame(rows); f = HERE / "build_audit.csv"
    if f.exists():
        old = pd.read_csv(f); df = pd.concat([old[~old.context.isin(df.context)], df])
    df.to_csv(f, index=False)
    (HERE / "build_audit.json").write_text(df.to_json(orient="records", indent=1))


# ===================================================================================== scoring (E4a path)
def chip_scores(ctx, tf, tg, W, tag):
    """signed cosine of the V-centred full 512-d embedding, per layer (e4a_main.stage_scores)."""
    f = W / f"scores_{tag}.npz"
    if f.exists(): return dict(np.load(f))
    emb = np.load(M.EMB_LOCAL / f"{ctx}.npy", mmap_mode="r"); obs = np.load(OLD_CACHE / ctx / "obs_idx.npy")
    assert emb.shape[0] == NL
    out = {}
    for l in range(NL):
        Xl = np.asarray(emb[l], np.float64)[obs]
        assert np.all(np.linalg.norm(Xl, axis=1) > 0), "never-observed row in V"
        Xc = Xl - Xl.mean(0, keepdims=True); En = S.l2norm_rows(Xc)
        out[f"L{l:02d}"] = L.rowcos(En, tf, tg)
    d = OLD_CACHE / ctx
    Ex = np.load(d / "expr_extraction.npy"); E3 = np.load(d / "expr.npy")
    for t, Mx in (("cx_ext", Ex), ("cx_3k", E3)):
        raw = L.abs_corr_pairs(Mx, tf, tg); out[t + "_raw"] = raw; out[t] = np.round(raw, DEC)
    np.savez(f, **out); log(ctx, f"scores {tag} done")
    return out


def sens_rows(ch, fs, sc, wts, cx_key="cx_ext", cx3_key="cx_3k", basis=BASIS):
    ce = ch.auroc(sc[cx_key]); cb = fs.auc(sc[cx_key], wts); c3 = ch.auroc(sc[cx3_key]); c3b = fs.auc(sc[cx3_key], wts)
    rows = [M.rowdict("|co-expression| extraction cells", -1, ce, cb), M.rowdict("|co-expression| 3,000-cell tier", -1, c3, c3b)]
    for l in range(NL):
        s = sc[f"L{l:02d}"]
        rows.append(M.rowdict("signed cosine full-512", l, ch.auroc(s), fs.auc(s, wts), ce, cb, c3, c3b))
    for r in rows: r["delta_cx_basis"] = basis
    return rows


def cmd_sens(ctxs):
    for ctx in ctxs:
        W = W0 / ctx; W.mkdir(parents=True, exist_ok=True); t0 = time.time(); seed = M.seed_of(ctx)
        fx = dict(np.load(NEW_CACHE / ctx / "h3_chip_pairs_fixed.npz"))
        ch = S.StratifiedPairs(fx["tf"], fx["tg"], fx["y"]); keep = np.isin(fx["tf"], ch.tfs)
        assert np.array_equal(fx["tg"][keep], ch.tg)
        sc_all = chip_scores(ctx, fx["tf"], fx["tg"], W, "fixed")
        sc = {k: v[keep] for k, v in sc_all.items()}
        fs = L.FastStrat(ch)
        # ---- main block (E4a stream [seed, 40], B = 2,000 everywhere)
        f = W / "sens_chip_fixed.json"
        if not f.exists():
            wts = fs.boot_weights(B, np.random.default_rng([seed, 40]))
            rows = sens_rows(ch, fs, sc, wts)
            M.jdump(dict(rows=rows, n_pos=int(ch.y.sum()), n_neg=int((~ch.y).sum()), n_tf=int(ch.tfs.size), B=B), f)
            log(ctx, "fixed main", f"L0 {rows[2]['auroc']:.4f} cx {rows[0]['auroc']:.4f}", f"{time.time()-t0:.0f}s")
        # ---- match_level 0 only (+ all positives), as e4a_extra block (b)
        f = W / "sens_chip_fixed_ml0.json"
        if not f.exists():
            lev = fx["match_level"][keep]; m = ch.y | (lev == 0)
            sub = S.StratifiedPairs(ch.tf[m], ch.tg[m], ch.y[m]); ii = np.flatnonzero(m)[np.isin(ch.tf[m], sub.tfs)]
            assert np.array_equal(ch.tg[ii], sub.tg)
            fsb = L.FastStrat(sub); wb = fsb.boot_weights(B, np.random.default_rng([seed, 41]))
            rows = sens_rows(sub, fsb, {k: v[ii] for k, v in sc.items()}, wb)
            M.jdump(dict(rows=rows, n_pos=int(sub.y.sum()), n_neg=int((~sub.y).sum()), n_tf=int(sub.tfs.size), B=B), f)
            log(ctx, "fixed ml0", f"{time.time()-t0:.0f}s")
        # ---- decomposition: old pairs with contaminated candidate rows deleted (points only)
        f = W / "old_minus_contaminated.json"
        if not f.exists():
            od = dict(np.load(OLD_CACHE / ctx / "h3_chip_pairs.npz"))
            osc = chip_scores(ctx, od["tf"], od["tg"], W, "old")
            posset = set(zip(od["tf"][od["y"] == 1].tolist(), od["tg"][od["y"] == 1].tolist()))
            bad = np.array([(yy == 0) and (((a, b) in posset) or ((b, a) in posset)) for a, b, yy in zip(od["tf"].tolist(), od["tg"].tolist(), od["y"].tolist())])
            out = dict(n_rows_removed=int(bad.sum()))
            for tag, msk in (("old_points_recomputed", np.ones(bad.size, bool)), ("old_minus_contaminated", ~bad)):
                sp2 = S.StratifiedPairs(od["tf"][msk], od["tg"][msk], od["y"][msk]); k2 = np.isin(od["tf"][msk], sp2.tfs)
                f2 = L.FastStrat(sp2); res = dict(n_pos=int(sp2.y.sum()), n_neg=int((~sp2.y).sum()), n_tf=int(sp2.tfs.size))
                for key in ["cx_ext", "cx_3k"] + [f"L{l:02d}" for l in range(NL)]:
                    res[key] = float(f2.auc(osc[key][msk][k2])[0])
                out[tag] = res
            M.jdump(out, f); log(ctx, "old-minus-contaminated", f"{time.time()-t0:.0f}s")
        (W / "DONE").write_text(f"{time.time()-t0:.0f}s\n"); log(ctx, "DONE", f"{time.time()-t0:.0f}s")


def cmd_repro(ctxs):
    """Old pairs, old settings -> must equal E4a_edges/sens_chip.csv (cosine rows always; comparator rows with
    the comparator basis E4a used in that context)."""
    old = pd.read_csv(E4 / "sens_chip.csv"); out = []
    for ctx in ctxs:
        W = W0 / ctx; W.mkdir(parents=True, exist_ok=True); seed = M.seed_of(ctx)
        od = dict(np.load(OLD_CACHE / ctx / "h3_chip_pairs.npz")); ch = S.StratifiedPairs(od["tf"], od["tg"], od["y"])
        keep = np.isin(od["tf"], ch.tfs); sc = {k: v[keep] for k, v in chip_scores(ctx, od["tf"], od["tg"], W, "old").items()}
        fs = L.FastStrat(ch); Bc = M.B if ctx in ("lung", "lung_external") or ch.tf.size < 50000 else M.B_SECONDARY
        wts = fs.boot_weights(Bc, np.random.default_rng([seed, 40]))
        tie = ctx in ("lung", "lung_external")
        rows = sens_rows(ch, fs, sc, wts, "cx_ext" if tie else "cx_ext_raw", "cx_3k" if tie else "cx_3k_raw")
        o = old[old.context == ctx].reset_index(drop=True)
        cols = ["auroc", "ci_lo", "ci_hi", "ub95", "boot_sd", "delta_cxext", "delta_cxext_ub95", "delta_cx3k", "delta_cx3k_ub95"]
        mx = 0.0
        for i, r in enumerate(rows):
            assert o.scorer[i] == r["scorer"] and int(o.layer[i]) == r["layer"]
            for c in cols:
                if c in r and np.isfinite(o[c][i]): mx = max(mx, abs(o[c][i] - r[c]))
        out.append(dict(context=ctx, B=Bc, comparator="tie-corrected" if tie else "raw float (as E4a)", max_absdiff_vs_sens_chip_csv=mx,
                        n_pos=int(ch.y.sum()), n_neg=int((~ch.y).sum())))
        log(ctx, "repro max |diff|", mx)
    f = HERE / "repro_old_sens_chip.csv"; df = pd.DataFrame(out)
    if f.exists():
        prev = pd.read_csv(f); df = pd.concat([prev[~prev.context.isin(df.context)], df])
    df.to_csv(f, index=False)


# ===================================================================================== aggregate
def cmd_aggregate():
    old = pd.read_csv(E4 / "sens_chip.csv"); cols = list(old.columns)
    aud = pd.read_csv(HERE / "build_audit.csv").set_index("context")
    main_rows, ml0_rows, cmp_rows = [], [], []
    ctxs = [c for c in X.CONTEXTS if (W0 / c / "DONE").exists()]
    for c in ctxs:
        j = jl(W0 / c / "sens_chip_fixed.json")
        for r in j["rows"]: main_rows.append(dict(context=c, n_pos=j["n_pos"], n_neg=j["n_neg"], n_tf=j["n_tf"], **r))
        j2 = jl(W0 / c / "sens_chip_fixed_ml0.json")
        for r in j2["rows"]: ml0_rows.append(dict(context=c, block="match_level 0 candidates only (+ all positives)", n_pos=j2["n_pos"], n_neg=j2["n_neg"], n_tf=j2["n_tf"], **r))
        dec = jl(W0 / c / "old_minus_contaminated.json")
        o = old[old.context == c].reset_index(drop=True); nw = pd.DataFrame(j["rows"]); m0 = pd.DataFrame(j2["rows"])
        for i in range(len(nw)):
            key = "cx_ext" if nw.layer[i] == -1 and "extraction" in nw.scorer[i] else ("cx_3k" if nw.layer[i] == -1 else f"L{int(nw.layer[i]):02d}")
            cmp_rows.append(dict(context=c, scorer=nw.scorer[i], layer=int(nw.layer[i]),
                                 old_n_pos=int(o.n_pos[0]), old_n_neg=int(o.n_neg[0]), old_n_tf=int(o.n_tf[0]),
                                 new_n_pos=j["n_pos"], new_n_neg=j["n_neg"], new_n_tf=j["n_tf"],
                                 old_B=int(o.B[i]), new_B=int(nw.B[i]), old_delta_cx_basis=o.delta_cx_basis[i],
                                 old_auroc=o.auroc[i], new_auroc=nw.auroc[i], diff_auroc=nw.auroc[i] - o.auroc[i],
                                 old_ci_lo=o.ci_lo[i], old_ci_hi=o.ci_hi[i], new_ci_lo=nw.ci_lo[i], new_ci_hi=nw.ci_hi[i],
                                 old_ub95=o.ub95[i], new_ub95=nw.ub95[i],
                                 old_ci_excludes_05=bool(o.ci_lo[i] > 0.5 or o.ci_hi[i] < 0.5), new_ci_excludes_05=bool(nw.ci_lo[i] > 0.5 or nw.ci_hi[i] < 0.5),
                                 old_delta_cxext=o.delta_cxext[i], new_delta_cxext=nw.get("delta_cxext", pd.Series([np.nan]*len(nw)))[i],
                                 old_delta_cxext_ci_lo=o.delta_cxext_ci_lo[i], new_delta_cxext_ci_lo=nw.get("delta_cxext_ci_lo", pd.Series([np.nan]*len(nw)))[i],
                                 old_delta_cxext_ub95=o.delta_cxext_ub95[i], new_delta_cxext_ub95=nw.get("delta_cxext_ub95", pd.Series([np.nan]*len(nw)))[i],
                                 old_delta_cx3k=o.delta_cx3k[i], new_delta_cx3k=nw.get("delta_cx3k", pd.Series([np.nan]*len(nw)))[i],
                                 old_delta_cx3k_ub95=o.delta_cx3k_ub95[i], new_delta_cx3k_ub95=nw.get("delta_cx3k_ub95", pd.Series([np.nan]*len(nw)))[i],
                                 old_points_recomputed_tiecorrected=dec["old_points_recomputed"][key],
                                 old_minus_contaminated_auroc=dec["old_minus_contaminated"][key],
                                 n_contaminated_rows_removed=dec["n_rows_removed"],
                                 new_ml0_auroc=m0.auroc[i], new_ml0_ci_lo=m0.ci_lo[i], new_ml0_ci_hi=m0.ci_hi[i],
                                 new_ml0_delta_cxext=m0.get("delta_cxext", pd.Series([np.nan]*len(m0)))[i],
                                 new_ml0_delta_cxext_ub95=m0.get("delta_cxext_ub95", pd.Series([np.nan]*len(m0)))[i]))
    sf = pd.DataFrame(main_rows); assert set(sf.columns) == set(cols), set(sf.columns) ^ set(cols)
    sf[cols].to_csv(HERE / "sens_chip_fixed.csv", index=False)
    pd.DataFrame(ml0_rows).to_csv(HERE / "sens_chip_fixed_matchlevel0.csv", index=False)
    cmp = pd.DataFrame(cmp_rows); cmp.to_csv(HERE / "chip_old_vs_new.csv", index=False)
    ml = aud.reset_index()[["context", "old_n_neg", "new_n_neg", "old_share_level0", "old_share_level1", "old_share_level2",
                            "new_share_level0", "new_share_level1", "new_share_level2"]]
    ml.to_csv(HERE / "chip_match_levels.csv", index=False)
    cv = (R3 / "FROZEN.txt").read_text().split()[2]
    res = dict(experiment="E4a_chipfix", convention3_sha256=cv, contexts_done=ctxs,
               contexts_missing=[c for c in X.CONTEXTS if c not in ctxs],
               change="candidate non-edges exclude anyreg in full (prep_pairs used anyreg - set(ch_v))",
               bootstrap="two-way TF x target, B=2000, seed [seed_of(ctx),40] (main), [seed_of(ctx),41] (match-level-0)",
               comparator="|Pearson r| rounded to 12 decimals (tie-corrected), extraction cells and 3,000-cell tier",
               build_audit=jl(HERE / "build_audit.json"),
               repro=pd.read_csv(HERE / "repro_old_sens_chip.csv").to_dict(orient="records") if (HERE / "repro_old_sens_chip.csv").exists() else None,
               files=sorted(p.name for p in HERE.glob("*.csv")))
    (HERE / "results.json").write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    log("aggregate done", ctxs)


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "build": cmd_build(args or list(X.CONTEXTS))
    elif cmd == "sens": cmd_sens(args)
    elif cmd == "repro": cmd_repro(args)
    elif cmd == "aggregate": cmd_aggregate()
