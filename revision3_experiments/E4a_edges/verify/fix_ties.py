"""Verifier correction for E4a_edges: exact tie handling in two scorers.

Bug 1 (|co-expression| comparator). e4a_lib.abs_corr_pairs returns |Pearson r| in float64 without rounding. Pairs
whose |r| is mathematically identical (e.g. targets with identical expression vectors on the 256/192 extraction cells)
get values that differ in the last bits, so the stratified AUROC breaks these ties arbitrarily instead of counting
them 1/2 (the pre-registered definition, stats3.StratifiedPairs: "Ties count 1/2"). Fix: |r| rounded to 12 decimals
(stable: rounding to 10, 11 or 13 decimals gives the same AUROC, recorded below).
Bug 2 (rank-1 cosine). stats3.l2norm_rows divides by (norm + 1e-12), so a 1-D cosine is +-(1 - eps*|x|^-1) and exact
ties are ordered by magnitude; rounding to 10 decimals left ~0.8 % of pairs off +-1. Fix: exact sign(x_tf) sign(x_tg),
which equals sign of the stored value.

What is recomputed (same seeds and weight streams as e4a_main.py, via e4a_lib.FastStrat / e4a_main.rowdict):
  all 16 contexts: coexp scores, rows_coexp (|co-expr| rows), and in rows_L*.json the PRIMARY row, the |cos| row and
  the rank-1 row (all delta columns vs the corrected comparator); sens_macro co-expression rows.
  Contextual stratum ('TF detected in >= 10% of extraction cells') in every context; ChIP-only positives in lung and
  lung_external (pass --chip).
Rows not recomputed keep their AUROC/CI (unaffected: they do not involve |co-expr|) and get
'delta_cx_basis' = 'original float-tie-broken |co-expression| (not recomputed)'; recomputed rows get 'tie-corrected'.
Originals are copied to work/<ctx>/pre_verify_backup/ before anything is overwritten. Resumable (verify/fixed/<ctx>.done).
Usage: python fix_ties.py <ctx> ... [--chip]
"""
from __future__ import annotations
import json, shutil, sys, time, zlib
from pathlib import Path
import numpy as np

V = Path(__file__).resolve().parent; E4 = V.parent; R3 = E4.parent
sys.path.insert(0, str(R3)); sys.path.insert(0, str(E4))
import stats3 as S          # noqa: E402
import e4a_lib as L         # noqa: E402
import e4a_main as M        # noqa: E402

DEC = 12
BASIS_NEW = "tie-corrected"
BASIS_OLD = "original float-tie-broken |co-expression| (not recomputed)"


def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def cxr(M_, a, b): return np.round(L.abs_corr_pairs(M_, a, b), DEC)


def backup(W, names):
    bk = W / "pre_verify_backup"; bk.mkdir(exist_ok=True)
    for n in names:
        if (W / n).exists() and not (bk / n).exists(): shutil.copy2(W / n, bk / n)


def jl(p): return json.loads(Path(p).read_text())


def run(ctx, do_chip):
    W = E4 / "work" / ctx; done = V / "fixed" / f"{ctx}{'_chip' if do_chip else ''}.done"
    if done.exists(): log(ctx, "already fixed"); return
    t0 = time.time(); seed = M.seed_of(ctx); d = M.CACHE / ctx
    h3 = dict(np.load(d / "h3_pairs.npz")); sp = S.StratifiedPairs(h3["tf"], h3["tg"], h3["y"]); keep = np.isin(h3["tf"], sp.tfs)
    Ex = np.load(d / "expr_extraction.npy"); E3 = np.load(d / "expr.npy")
    backup(W, ["coexp.npz", "boot_coexp.npz", "rows_coexp.json", "sens_macro.json", "sens_strata.json", "sens_chip.json",
               "sens_extra.json", "sens_legacy_pooled.json"] + [f"rows_L{l:02d}.json" for l in range(12)])
    bk = W / "pre_verify_backup"
    old = dict(np.load(bk / "coexp.npz"))
    audit = dict(context=ctx)
    if do_chip:
        ch = dict(np.load(d / "h3_chip_pairs.npz")); chs = S.StratifiedPairs(ch["tf"], ch["tg"], ch["y"]); fs = L.FastStrat(chs)
        Bc = M.B if ctx in ("lung", "lung_external") or chs.tf.size < 50000 else M.B_SECONDARY
        wts = fs.boot_weights(Bc, np.random.default_rng([seed, 40]))
        ce_s = cxr(Ex, chs.tf, chs.tg); c3_s = cxr(E3, chs.tf, chs.tg)
        ce = chs.auroc(ce_s); cb = fs.auc(ce_s, wts); c3 = chs.auroc(c3_s); c3b = fs.auc(c3_s, wts)
        rows = [M.rowdict("|co-expression| extraction cells", -1, ce, cb), M.rowdict("|co-expression| 3,000-cell tier", -1, c3, c3b)]
        oldj = jl(bk / "sens_chip.json")
        for l in range(12):
            s = M.scores_layer(W, l)["chip"]; bt = fs.auc(s, wts)
            r = M.rowdict("signed cosine full-512", l, chs.auroc(s), bt, ce, cb, c3, c3b)
            o = oldj["rows"][l + 2]; assert o["layer"] == l and abs(o["ci_lo"] - r["ci_lo"]) < 1e-12 and abs(o["ub95"] - r["ub95"]) < 1e-12
            rows.append(r); log(ctx, "chip", l, f"{time.time()-t0:.0f}s")
        for r in rows: r["delta_cx_basis"] = BASIS_NEW
        audit.update(chip_cx_ext_old=oldj["rows"][0]["auroc"], chip_cx_ext_new=ce, chip_cx_3k_old=oldj["rows"][1]["auroc"], chip_cx_3k_new=c3)
        M.jdump(dict(oldj, rows=rows), W / "sens_chip.json")
        (V / "fixed" / f"{ctx}_chip_audit.json").write_text(json.dumps(audit, indent=1, default=float)); done.write_text("ok\n"); return
    # ---------------- corrected comparator scores (h3, legacy; chip kept for point audit)
    new = dict(old)
    leg = dict(np.load(d / "h3_legacy_pairs.npz")); lg = S.StratifiedPairs(leg["tf"], leg["tg"], leg["y"])
    for tag, Mx in (("cx_ext", Ex), ("cx_3k", E3)):
        raw = L.abs_corr_pairs(Mx, sp.tf, sp.tg)
        assert np.max(np.abs(raw - old[tag])) == 0.0, "stored comparator != recomputed raw |r|"
        new[tag] = np.round(raw, DEC); new[tag + "_leg"] = cxr(Mx, lg.tf, lg.tg)
        audit[f"{tag}_auroc_old"] = sp.auroc(old[tag])
        for dd in (10, 11, 12, 13): audit[f"{tag}_auroc_round{dd}"] = sp.auroc(np.round(raw, dd))
        audit[f"{tag}_n_unique_raw"] = int(np.unique(raw).size); audit[f"{tag}_n_unique_round12"] = int(np.unique(new[tag]).size)
        audit[f"{tag}_leg_pooled_old"] = S.auroc(old[tag + "_leg"][lg.y], old[tag + "_leg"][~lg.y])
        audit[f"{tag}_leg_pooled_new"] = S.auroc(new[tag + "_leg"][lg.y], new[tag + "_leg"][~lg.y])
    if "cx_ext_chip" in old:
        ch = dict(np.load(d / "h3_chip_pairs.npz")); chs = S.StratifiedPairs(ch["tf"], ch["tg"], ch["y"])
        for tag, Mx in (("cx_ext", Ex), ("cx_3k", E3)):
            new[tag + "_chip"] = cxr(Mx, chs.tf, chs.tg)
            audit[f"{tag}_chip_old"] = chs.auroc(old[tag + "_chip"]); audit[f"{tag}_chip_new"] = chs.auroc(new[tag + "_chip"])
    np.savez(W / "coexp.npz", **new)
    # ---------------- primary-stream bootstrap
    fs = L.FastStrat(sp)
    wts = fs.boot_weights(M.B, np.random.default_rng([seed, 2]), two_way=True)
    wts1 = fs.boot_weights(M.B, np.random.default_rng([seed, 3]), two_way=False)
    bc = dict(np.load(bk / "boot_coexp.npz"))
    for tag in ("cx_ext", "cx_3k"):
        bc[tag] = fs.auc(new[tag], wts); bc[tag + "_tfonly"] = fs.auc(new[tag], wts1); bc[tag + "_est"] = np.array([sp.auroc(new[tag])])
    np.savez(W / "boot_coexp.npz", **bc)
    ce, c3 = float(bc["cx_ext_est"][0]), float(bc["cx_3k_est"][0])
    rcx = jl(bk / "rows_coexp.json")
    for i, (tag, nm) in enumerate((("cx_ext", "|co-expression| extraction cells"), ("cx_3k", "|co-expression| 3,000-cell tier"))):
        r = M.rowdict(nm, -1, bc[tag + "_est"][0], bc[tag]); c1 = L.ci_row(bc[tag + "_est"][0], bc[tag + "_tfonly"])
        r.update(tfonly_ci_lo=c1["ci_lo"], tfonly_ci_hi=c1["ci_hi"], tfonly_ub95=c1["ub95"], delta_cx_basis=BASIS_NEW)
        assert rcx[i]["scorer"] == nm; rcx[i] = r
    M.jdump(rcx, W / "rows_coexp.json")
    log(ctx, "comparator done", {k: v for k, v in audit.items() if "auroc" in k}, f"{time.time()-t0:.0f}s")
    for l in range(12):
        rows = jl(bk / f"rows_L{l:02d}.json"); sc = M.scores_layer(W, l)
        s = sc["prim"]; bt = fs.auc(s, wts)
        r = M.rowdict("signed cosine full-512 (PRIMARY)", l, sp.auroc(s), bt, ce, bc["cx_ext"], c3, bc["cx_3k"])
        o = rows[0]; assert o["scorer"] == r["scorer"]
        for k in ("auroc", "ci_lo", "ci_hi", "ub95", "boot_sd"): assert abs(o[k] - r[k]) < 1e-12, (k, o[k], r[k])
        r.update({k: o[k] for k in ("tfonly_ci_lo", "tfonly_ci_hi", "tfonly_ub95", "macro_auroc")}); r["delta_cx_basis"] = BASIS_NEW
        rows[0] = r
        s = np.abs(sc["prim"]); r = M.rowdict("|cos| full-512", l, sp.auroc(s), fs.auc(s, wts), ce, bc["cx_ext"], c3, bc["cx_3k"])
        assert rows[1]["scorer"] == r["scorer"] and abs(rows[1]["ub95"] - r["ub95"]) < 1e-12; r["delta_cx_basis"] = BASIS_NEW; rows[1] = r
        s1 = np.sign(sc["rk1"]); assert not np.any(s1 == 0)
        r = M.rowdict("rank-1 cosine", l, sp.auroc(s1), M.boot_block(fs, s1, wts, M.B_SECONDARY), ce, bc["cx_ext"], c3, bc["cx_3k"])
        assert rows[2]["scorer"] == r["scorer"]
        r["auroc_before_exact_sign_fix"] = rows[2]["auroc"]; r["delta_cx_basis"] = BASIS_NEW + " (exact-sign rank-1)"; rows[2] = r
        for rr in rows[3:]: rr["delta_cx_basis"] = BASIS_OLD
        M.jdump(rows, W / f"rows_L{l:02d}.json"); log(ctx, f"L{l} primary delta UB {rows[0]['delta_cxext_ub95']:.4f}", f"{time.time()-t0:.0f}s")
    # ---------------- macro: comparator rows
    mac = jl(bk / "sens_macro.json")
    for r in mac:
        if r["scorer"] in ("cx_ext", "cx_3k"):
            nr = M.macro_block(sp, new[r["scorer"]], np.random.default_rng([seed, 11])); nr.update(scorer=r["scorer"], layer=-1); r.clear(); r.update(nr)
    M.jdump(mac, W / "sens_macro.json")
    # ---------------- strata: recompute the TF-detection strata (contextual); flag the rest
    st = jl(bk / "sens_strata.json"); det = new["tf_detfrac_ext"][sp.tf]
    for sname, m in (("TF detected in >= 10% of extraction cells", det >= 0.10), ("TF detected in < 10% of extraction cells", det < 0.10)):
        if m.sum() == 0: continue
        sub = S.StratifiedPairs(sp.tf[m], sp.tg[m], sp.y[m])
        if sub.tf.size == 0 or sub.y.sum() == 0: continue
        idx_m = np.flatnonzero(m); ii = idx_m[np.isin(sp.tf[m], sub.tfs)]; assert np.array_equal(sp.tg[ii], sub.tg)
        f2 = L.FastStrat(sub); Bs = M.B if sub.tf.size < 20000 else M.B_SECONDARY
        w2 = f2.boot_weights(Bs, np.random.default_rng([seed, 30, zlib.crc32(sname.encode())]))
        ce2 = sub.auroc(new["cx_ext"][ii]); cb2 = f2.auc(new["cx_ext"][ii], w2); c32 = sub.auroc(new["cx_3k"][ii]); c3b2 = f2.auc(new["cx_3k"][ii], w2)
        rws = [M.rowdict("|co-expression| extraction cells", -1, ce2, cb2)]
        for l in range(12):
            s = M.scores_layer(W, l)["prim"][ii]
            rws.append(M.rowdict("signed cosine full-512", l, sub.auroc(s), f2.auc(s, w2), ce2, cb2, c32, c3b2))
        for r in rws: r.update(stratum=sname, n_pos=int(sub.y.sum()), n_neg=int((~sub.y).sum()), n_tf=int(sub.tfs.size), delta_cx_basis=BASIS_NEW)
        olds = [r for r in st if r.get("stratum") == sname]
        for o, r in zip(olds[1:], rws[1:]): assert abs(o["ci_lo"] - r["ci_lo"]) < 1e-12 and abs(o["ub95"] - r["ub95"]) < 1e-12
        audit[f"strata_{sname}_cx_old"] = olds[0]["auroc"]; audit[f"strata_{sname}_cx_new"] = ce2
        st = [r for r in st if r.get("stratum") != sname] + rws
        log(ctx, "stratum", sname, f"{time.time()-t0:.0f}s")
    for r in st: r.setdefault("delta_cx_basis", BASIS_OLD)
    M.jdump(st, W / "sens_strata.json")
    # ---------------- flag rows not recomputed; audit their comparator point shift
    ex = jl(bk / "sens_extra.json")
    for r in ex["rows"]: r["delta_cx_basis"] = BASIS_OLD
    for r in ex.get("mapped_bands", {}).get("rows", []): r["delta_cx_basis"] = BASIS_OLD
    M.jdump(ex, W / "sens_extra.json")
    lp = jl(bk / "sens_legacy_pooled.json")
    for r in lp["rows"]: r["delta_cx_basis"] = BASIS_OLD
    M.jdump(lp, W / "sens_legacy_pooled.json")
    if (bk / "sens_chip.json").exists() and ctx not in ("lung", "lung_external"):
        cj = jl(bk / "sens_chip.json")
        for r in cj["rows"]: r["delta_cx_basis"] = BASIS_OLD
        M.jdump(cj, W / "sens_chip.json")
    audit["seconds"] = time.time() - t0
    (V / "fixed" / f"{ctx}_audit.json").write_text(json.dumps(audit, indent=1, default=float))
    done.write_text("ok\n"); log(ctx, "FIXED", f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]; chip = "--chip" in sys.argv
    for c in args: run(c, chip)
