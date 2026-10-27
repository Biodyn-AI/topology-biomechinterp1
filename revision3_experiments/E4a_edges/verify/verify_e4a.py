"""Independent verifier for E4a_edges (does not import stats3, e4a_lib or e4a_main).

Recomputes from raw inputs (local embeddings, cache/<ctx>/h3_pairs.npz, obs_idx, expr_extraction.npy, expr.npy):
  1. primary within-TF stratified AUROC of signed V-centred full-512 cosine, all 12 layers (lung, lung_external,
     kidney_A, kidney_B, immune_Bcell) with scipy.stats.rankdata per TF;
  2. |Pearson r| co-expression AUROC on the extraction cells and on the 3k tier (same contexts);
  3. lung L0, L1, L11 two-way TF x target pigeonhole bootstrap, own loop implementation, (a) replaying the
     implementer's seed stream (should match exactly) and (b) an independent seed (should agree to MC error);
     paired delta vs |co-expr| (extraction) on the same resamples;
  4. lung within-TF label permutation, own loop (rng.permutation per TF), N = 2,000 joint over 12 layers ->
     null SD (vs analytic sqrt(sum n1 n0 (n+1)/12)/sum n1 n0) and max-T p;
  5. immune_Bcell restricted to TFs detected in >= 10% of extraction cells (L0);
  6. legacy random-split percentiles from legacy_splits_draws.csv.
Writes verify/verify_results.json and verify/verify_primary.csv.
"""
import json, os, sys, time, zlib
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent; E4 = HERE.parent; R3 = E4.parent
CACHE = R3 / "cache"; EMB = Path(os.path.expanduser("~/biomechinterp-local/subproject41_r3/emb"))


def seed_of(ctx): return 20261026 + zlib.crc32(ctx.encode()) % 100000


def strat_auc(tf, y, s, w=None):
    """sum_t sum_{p,q} w_p w_q [1(s_p>s_q)+.5 1(=)] / sum_t W_pos W_neg ; own implementation."""
    num = den = 0.0
    if w is None: w = np.ones(s.size)
    for t in np.unique(tf):
        m = tf == t
        sp, wp = s[m & y], w[m & y]; sn, wn = s[m & ~y], w[m & ~y]
        if wp.sum() == 0 or wn.sum() == 0: continue
        if w is not None and not np.all(w == 1):
            o = np.argsort(sn); sn, wn = sn[o], wn[o]; cw = np.r_[0, np.cumsum(wn)]
            lo = np.searchsorted(sn, sp, "left"); hi = np.searchsorted(sn, sp, "right")
            num += (wp * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))).sum()
        else:
            r = rankdata(np.r_[sp, sn]); n1 = sp.size
            num += r[:n1].sum() - n1 * (n1 + 1) / 2
        den += wp.sum() * wn.sum()
    return num / den


def load(ctx):
    d = CACHE / ctx
    h = np.load(d / "h3_pairs.npz"); tf, tg, y = h["tf"], h["tg"], h["y"].astype(bool)
    # keep TFs with both labels (own filter)
    ok = np.array([y[tf == t].any() and (~y[tf == t]).any() for t in np.unique(tf)])
    keepT = np.unique(tf)[ok]; k = np.isin(tf, keepT)
    return d, tf[k], tg[k], y[k], h


def cos_layers(ctx, d, tf, tg):
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r"); oi = np.load(d / "obs_idx.npy")
    genes = json.loads((d / "obs_genes.json").read_text()); assert len(genes) == oi.size
    out = []
    for l in range(12):
        X = np.asarray(emb[l], np.float64)[oi]
        assert (np.linalg.norm(X, axis=1) > 0).all()
        # never-observed rows really are zero elsewhere
        if l == 0:
            allz = np.linalg.norm(np.asarray(emb[0], np.float64), axis=1) > 0
            assert allz.sum() == oi.size and np.array_equal(np.flatnonzero(allz), oi)
        X = X - X.mean(0); X /= np.linalg.norm(X, axis=1, keepdims=True)
        out.append((X[tf] * X[tg]).sum(1))
    return np.array(out)


def abs_r(M, a, b):
    M = np.asarray(M, np.float64); out = np.empty(a.size)
    mu = M.mean(0); sd = M.std(0)
    for i in range(0, a.size, 20000):
        aa, bb = a[i:i+20000], b[i:i+20000]
        cov = ((M[:, aa] - mu[aa]) * (M[:, bb] - mu[bb])).mean(0)
        den = sd[aa] * sd[bb]
        out[i:i+20000] = np.where(den > 0, np.abs(cov / np.where(den > 0, den, 1)), 0.0)
    return out


def boot(tf, tg, y, scores, B, rng):
    tfs, tfc = np.unique(tf, return_inverse=True); tgs, tgc = np.unique(tg, return_inverse=True)
    res = np.empty((B, len(scores)))
    for r in range(B):
        wt = np.bincount(rng.integers(0, tfs.size, tfs.size), minlength=tfs.size)
        wg = np.bincount(rng.integers(0, tgs.size, tgs.size), minlength=tgs.size)
        w = (wt[tfc] * wg[tgc]).astype(float)
        for j, s in enumerate(scores): res[r, j] = strat_auc(tfc, y, s, w)
    return res


def main():
    t0 = time.time(); R = {}; prim_rows = []
    for ctx in ("lung", "lung_external", "kidney_A", "kidney_B", "immune_Bcell"):
        d, tf, tg, y, h = load(ctx)
        P = cos_layers(ctx, d, tf, tg)
        cx_e = abs_r(np.load(d / "expr_extraction.npy"), tf, tg); cx_3 = abs_r(np.load(d / "expr.npy"), tf, tg)
        au = [strat_auc(tf, y, P[l]) for l in range(12)]
        ref = pd.read_csv(E4 / f"h3_primary_{ctx}.csv")
        r = dict(context=ctx, n_pos=int(y.sum()), n_neg=int((~y).sum()), n_tf=int(np.unique(tf).size),
                 auroc_mine=au, auroc_impl=ref.auroc.tolist(), max_absdiff=float(np.max(np.abs(np.array(au) - ref.auroc.values))),
                 cx_ext_mine=strat_auc(tf, y, cx_e), cx_ext_impl=float(ref.auroc_cx_ext[0]),
                 cx_3k_mine=strat_auc(tf, y, cx_3), cx_3k_impl=float(ref.auroc_cx_3k[0]),
                 n_pos_impl=int(ref.n_pos[0]), n_neg_impl=int(ref.n_neg[0]), n_tf_impl=int(ref.n_tf[0]))
        R[ctx] = r
        for l in range(12): prim_rows.append(dict(context=ctx, layer=l, auroc_mine=au[l], auroc_impl=ref.auroc[l]))
        print(ctx, "maxdiff", r["max_absdiff"], "cx", r["cx_ext_mine"], r["cx_ext_impl"], r["cx_3k_mine"], r["cx_3k_impl"], f"{time.time()-t0:.0f}s", flush=True)
        if ctx == "immune_Bcell":
            det = (np.load(d / "expr_extraction.npy") > 0).mean(0)[tf]
            m = det >= 0.10
            sub_ok = np.array([y[m & (tf == t)].any() and (~y[m & (tf == t)]).any() for t in np.unique(tf[m])])
            keepT = np.unique(tf[m])[sub_ok]; mm = m & np.isin(tf, keepT)
            R["immune_Bcell_TFdet10_L0"] = dict(auroc_mine=strat_auc(tf[mm], y[mm], P[0][mm]), n_tf=int(keepT.size),
                                                cx_ext_mine=strat_auc(tf[mm], y[mm], cx_e[mm]))
            sc = pd.read_csv(E4 / "sens_contextual.csv")
            q = sc[(sc.context == "immune_Bcell") & (sc.layer == 0) & (sc.scorer == "signed cosine full-512")]
            R["immune_Bcell_TFdet10_L0"]["auroc_impl"] = float(q.auroc.iloc[0]); R["immune_Bcell_TFdet10_L0"]["n_tf_impl"] = int(q.n_tf.iloc[0])
            print(R["immune_Bcell_TFdet10_L0"], flush=True)
        if ctx == "lung":
            Pl, tfl, tgl, yl, cxl = P, tf, tg, y, cx_e
    # ---- lung bootstrap: replay implementer's stream, and an independent seed
    tfs, tfc = np.unique(tfl, return_inverse=True)
    scores = [Pl[0], Pl[1], Pl[11], cxl]
    ref = pd.read_csv(E4 / "h3_primary_lung.csv")
    for tag, rng in (("replay_seed", np.random.default_rng([seed_of("lung"), 2])), ("independent_seed", np.random.default_rng(987654321))):
        bs = boot(tfl, tgl, yl, scores, 2000, rng)
        rr = {}
        for j, l in enumerate((0, 1, 11)):
            dlt = bs[:, j] - bs[:, 3]
            rr[f"L{l}"] = dict(ci_lo=float(np.quantile(bs[:, j], .025)), ci_hi=float(np.quantile(bs[:, j], .975)),
                               ub95=float(np.quantile(bs[:, j], .95)), sd=float(bs[:, j].std(ddof=1)),
                               delta_ub95=float(np.quantile(dlt, .95)),
                               impl=dict(ci_lo=float(ref.ci_lo[l]), ci_hi=float(ref.ci_hi[l]), ub95=float(ref.ub95[l]),
                                         sd=float(ref.boot_sd[l]), delta_ub95=float(ref.delta_cxext_ub95[l])))
        R[f"lung_boot_{tag}"] = rr; print(tag, json.dumps(rr), f"{time.time()-t0:.0f}s", flush=True)
    # ---- lung permutation, own loop
    rng = np.random.default_rng(13579); N = 2000
    blocks = [np.flatnonzero(tfc == c) for c in range(tfs.size)]
    rk = np.empty_like(Pl)
    for b in blocks:
        for l in range(12): rk[l, b] = rankdata(Pl[l, b])
    n1 = np.array([yl[b].sum() for b in blocks]); n = np.array([b.size for b in blocks]); n0 = n - n1
    const = (n1 * (n1 + 1) / 2).sum(); den = (n1 * n0).sum()
    obs = (rk[:, yl].sum(1) - const) / den
    null = np.empty((N, 12)); yp = yl.copy()
    for i in range(N):
        for b in blocks: yp[b] = rng.permutation(yl[b])
        null[i] = (rk[:, yp].sum(1) - const) / den
    sd = null.std(0, ddof=1); mu = null.mean(0)
    z = (obs - mu) / sd; mx = ((null - mu) / sd).max(1)
    p_maxT = [(np.sum(mx >= zz) + 1) / (N + 1) for zz in z]
    sd_analytic = np.sqrt((n1 * n0 * (n + 1) / 12).sum()) / den
    perm = json.loads((E4 / "work/lung/perm.json").read_text())
    R["lung_perm_own"] = dict(N=N, obs=obs.tolist(), null_sd=sd.tolist(), sd_analytic=float(sd_analytic), z=z.tolist(), p_maxT=p_maxT,
                              impl_null_sd=[s["sd"] for s in perm["null_summary"]], impl_p_maxT=perm["p_maxT"], impl_z=perm["z"])
    print("perm", json.dumps(R["lung_perm_own"])[:800], flush=True)
    # ---- legacy split percentiles
    dr = pd.read_csv(E4 / "legacy_splits_draws.csv")
    full = dr["auroc::full-512 L0 cosine [submitted: uncentred]"]; cx = dr["auroc::|co-expression| (3,000 cells, run_trees definition)"]
    dl = dr["delta::full-512 L0 cosine [submitted: uncentred]"]
    R["legacy"] = dict(mean_full=float(full.mean()), sd_full=float(full.std()), pct_0591=float(100 * (full < 0.591).mean()),
                       pct_cx_0556=float(100 * (cx < 0.5561571567402471).mean()), pct_delta_submitted=float(100 * (dl < 0.5914202415660141 - 0.5561571567402471).mean()),
                       frac_delta_pos=float((dl > 0).mean()))
    print("legacy", R["legacy"], flush=True)
    pd.DataFrame(prim_rows).to_csv(HERE / "verify_primary.csv", index=False)
    (HERE / "verify_results.json").write_text(json.dumps(R, indent=1, default=float))
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
