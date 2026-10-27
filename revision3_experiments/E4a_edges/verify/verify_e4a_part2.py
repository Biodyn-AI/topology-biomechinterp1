"""Independent verifier, part 2 (does not import stats3 / e4a_lib / e4a_main). Raw inputs only:
local embeddings, cache/<ctx>/{h3_pairs,h3_chip_pairs}.npz, obs_idx, obs_genes, expr_extraction.npy, edges_trrust.npy.
  A. tie-corrected |co-expression| (extraction cells) AUROC, lung / immune_Bcell / kidney_ckpt_kidney / immune_CD4T
  B. ChIP-only positives, lung: L0 and L1 signed-cosine stratified AUROC and tie-corrected |co-expr| AUROC
  C. lung other-TF in-degree shortcut (counted directly from edges_trrust, not as total - y) and
     the TF x exact other-TF in-degree stratified AUROC at L0
  D. lung macro per-TF AUROC at L0 (TFs with >= 3 positives)
  E. PAX5 -> CD79A rank among PAX5's candidates, immune_Bcell and lung, every layer
  F. label-noise correction of lung max UB
  G. own power-size simulation (lung, kidney_B): TF-cluster resampling, binormal scores with target ICC and TF effect
     heterogeneity tau (tau, rho_t taken from work/<ctx>/power.json), population-mean AUROC exactly 0.5; rejection
     rate of |A-0.5| > 1.96 sd_perm (two-sided, what power.csv calls permutation size) and of A-0.5 > 1.645 sd_perm
     (one-sided, the direction of the pre-registered max-T test); plus own per-TF-AUROC heterogeneity estimate.
Writes verify/verify_part2.json.
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata, norm

HERE = Path(__file__).resolve().parent; E4 = HERE.parent; R3 = E4.parent; CACHE = R3 / "cache"
EMB = Path(os.path.expanduser("~/biomechinterp-local/subproject41_r3/emb"))


def filt(tf, tg, y, *extra):
    ok = {t for t in np.unique(tf) if y[tf == t].any() and (~y[tf == t]).any()}
    k = np.array([t in ok for t in tf])
    return (tf[k], tg[k], y[k]) + tuple(e[k] for e in extra)


def strat_auc(strata, y, s):
    num = den = 0.0
    df = pd.DataFrame(dict(g=strata, y=y, s=s))
    for _, d in df.groupby("g"):
        n1 = int(d.y.sum()); n0 = len(d) - n1
        if n1 == 0 or n0 == 0: continue
        r = rankdata(d.s.values); num += r[d.y.values].sum() - n1 * (n1 + 1) / 2; den += n1 * n0
    return num / den


def auc2(p, q):
    r = rankdata(np.r_[p, q]); return (r[:p.size].sum() - p.size * (p.size + 1) / 2) / (p.size * q.size)


def cosl(ctx, l, a, b):
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r"); oi = np.load(CACHE / ctx / "obs_idx.npy")
    X = np.asarray(emb[l], np.float64)[oi]; X = X - X.mean(0); X /= np.linalg.norm(X, axis=1, keepdims=True)
    return np.einsum("ij,ij->i", X[a], X[b])


def absr(ctx, a, b, dec=12):
    M = np.load(CACHE / ctx / "expr_extraction.npy").astype(np.float64); mu = M.mean(0); sd = M.std(0)
    out = np.empty(a.size)
    for i in range(0, a.size, 50000):
        aa, bb = a[i:i+50000], b[i:i+50000]
        out[i:i+50000] = np.abs(((M[:, aa] - mu[aa]) * (M[:, bb] - mu[bb])).mean(0) / (sd[aa] * sd[bb]))
    return np.round(out, dec)


def load(ctx, f="h3_pairs.npz"):
    h = np.load(CACHE / ctx / f); return filt(h["tf"], h["tg"], h["y"].astype(bool))


def main():
    t0 = time.time(); R = {}
    # A
    for ctx in ("lung", "immune_Bcell", "kidney_ckpt_kidney", "immune_CD4T"):
        tf, tg, y = load(ctx); c = absr(ctx, tf, tg); c10 = absr(ctx, tf, tg, 10)
        R[f"A_{ctx}_cx_ext_tiecorrected"] = strat_auc(tf, y, c); R[f"A_{ctx}_cx_ext_round10"] = strat_auc(tf, y, c10)
        print("A", ctx, R[f"A_{ctx}_cx_ext_tiecorrected"], R[f"A_{ctx}_cx_ext_round10"], f"{time.time()-t0:.0f}s", flush=True)
    # B
    tf, tg, y = load("lung", "h3_chip_pairs.npz")
    R["B_lung_chip_n"] = dict(pos=int(y.sum()), neg=int((~y).sum()), tf=int(np.unique(tf).size))
    for l in (0, 1): R[f"B_lung_chip_L{l}"] = strat_auc(tf, y, cosl("lung", l, tf, tg))
    R["B_lung_chip_cx_ext"] = strat_auc(tf, y, absr("lung", tf, tg))
    print("B", {k: v for k, v in R.items() if k.startswith("B")}, f"{time.time()-t0:.0f}s", flush=True)
    # C, D, F
    tf, tg, y = load("lung"); ed = np.load(CACHE / "lung" / "edges_trrust.npy")
    tgt_tfs = {}
    for a, b in ed: tgt_tfs.setdefault(int(b), set()).add(int(a))
    oth = np.array([len(tgt_tfs.get(int(b), set()) - {int(a)}) for a, b in zip(tf, tg)])
    R["C_lung_othertf_indeg_auroc"] = strat_auc(tf, y, oth.astype(float))
    s0 = cosl("lung", 0, tf, tg)
    key = tf.astype(np.int64) * 100000 + oth
    R["C_lung_L0_strata_TFxindeg"] = strat_auc(key, y, s0)
    vals = [auc2(s0[(tf == t) & y], s0[(tf == t) & ~y]) for t in np.unique(tf) if ((tf == t) & y).sum() >= 3]
    R["D_lung_macro_L0"] = float(np.mean(vals)); R["D_lung_frac_tf_above_half_L0"] = float(np.mean(np.array(vals) > 0.5))
    prim = pd.read_csv(E4 / "h3_primary_lung.csv"); ub = prim.ub95.max()
    R["F_lung_maxUB"] = ub; R["F_lung_maxUB_corrected"] = {pi: 0.5 + (ub - 0.5) / (1 - pi) for pi in (0.05, 0.1, 0.2)}
    print("CDF", {k: v for k, v in R.items() if k[0] in "CDF"}, f"{time.time()-t0:.0f}s", flush=True)
    # E
    for ctx in ("immune_Bcell", "lung"):
        g = json.loads((CACHE / ctx / "obs_genes.json").read_text()); p = {x: i for i, x in enumerate(g)}
        tf, tg, y = load(ctx); m = tf == p["PAX5"]; h = np.flatnonzero(m & (tg == p["CD79A"]))[0]
        ranks = []
        for l in range(12):
            s = cosl(ctx, l, tf, tg); ranks.append(dict(layer=l, rank=int((s[m] > s[h]).sum() + 1), n=int(m.sum()),
                                                         frac_neg_below=float(np.mean(s[m & ~y] < s[h]))))
        R[f"E_{ctx}_pax5_cd79a"] = ranks; print("E", ctx, [r["rank"] for r in ranks], flush=True)
    # G
    for ctx in ("lung", "kidney_B"):
        pw = json.loads((E4 / "work" / ctx / "power.json").read_text()); tau = pw["tau_tf_effect_sd"]; rho_t = pw["rho_target_icc"]
        tf, tg, y = load(ctx); tfs, tfc = np.unique(tf, return_inverse=True); tgs, tgc = np.unique(tg, return_inverse=True)
        blocks = [np.flatnonzero(tfc == c) for c in range(tfs.size)]
        # own heterogeneity estimate on observed L1 scores: var of per-TF AUROC minus mean HM variance
        s1 = cosl(ctx, 1, tf, tg); a_c, v_c = [], []
        for b in blocks:
            yy = y[b]; n1, n0 = yy.sum(), (~yy).sum()
            if n1 < 3 or n0 < 3: continue
            a = auc2(s1[b][yy], s1[b][~yy]); q1 = a / (2 - a); q2 = 2 * a * a / (1 + a)
            a_c.append(a); v_c.append((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))
        sd_auc_between = float(np.sqrt(max(0, np.var(a_c, ddof=1) - np.mean(v_c))))
        rng = np.random.default_rng(424242); nsim = 300; est = []; sdp = []
        for i in range(nsim):
            pick = rng.integers(0, tfs.size, tfs.size)
            rows = np.concatenate([blocks[c] for c in pick]); newtf = np.concatenate([np.full(blocks[c].size, j) for j, c in enumerate(pick)])
            yy = y[rows]; u = rng.normal(0, np.sqrt(rho_t), tgs.size)[tgc[rows]]
            eff = rng.normal(0, tau, tfs.size)[newtf]                       # TF effects, population mean 0 -> mean AUROC 0.5
            s = u + rng.normal(0, np.sqrt(1 - rho_t), yy.size) + eff * yy
            est.append(strat_auc(newtf, yy, s) - 0.5)
            n1 = np.bincount(newtf, weights=yy.astype(float)); n = np.bincount(newtf); n0 = n - n1
            sdp.append(np.sqrt((n1 * n0 * (n + 1) / 12).sum()) / (n1 * n0).sum())
        est = np.array(est); sdp = np.array(sdp)
        R[f"G_{ctx}"] = dict(tau=tau, rho_t=rho_t, nsim=nsim, sampling_sd_at_null=float(est.std(ddof=1)), mean_sd_perm=float(sdp.mean()),
                             reject_two_sided_1p96_sdperm=float(np.mean(np.abs(est) > 1.959964 * sdp)),
                             reject_one_sided_1p645_sdperm=float(np.mean(est > 1.644854 * sdp)),
                             own_between_tf_sd_of_auroc_L1=sd_auc_between, implied_by_tau=float(norm.pdf(0) * tau / np.sqrt(2)),
                             impl_size_permutation=pw["auroc_size_permutation"], impl_sampling_sd_at_null=pw["auroc_sampling_sd_at_null"])
        print("G", ctx, R[f"G_{ctx}"], f"{time.time()-t0:.0f}s", flush=True)
    (HERE / "verify_part2.json").write_text(json.dumps(R, indent=1, default=float)); print("done", time.time() - t0)


if __name__ == "__main__":
    main()
