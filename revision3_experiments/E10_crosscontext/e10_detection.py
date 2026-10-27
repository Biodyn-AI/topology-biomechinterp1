"""E10 step 4 (task 3): detection.

Per context (all 16):
  (a) AUROC of log n_obs ALONE: Lambert TF vs non-TF, and TRRUST TF vs target-only (stats3.auroc, ties 1/2),
      gene-bootstrap 95% CI (B = 2,000, multinomial over V, fixed scores).
  (b) per layer: logistic  is_TF ~ z(score) + z(log n_obs) + z(score) x z(log n_obs), score = the H2 LOO mean
      cosine to the Lambert TF set (full dim, V-centred); Wald 95% CIs (statsmodels Logit).
  (c) S2 (H2 primary statistic) on V restricted to n_obs >= 3, 5, 10: universe, centring, TF set and detection
      deciles are all recomputed inside the restricted V; gene-bootstrap 95% CI (Poisson, B = 500).
Also writes per-context detection descriptors (median n_obs of TFs / non-TFs; truncation of cell inputs at
max_genes = 512 from A2).
Output: detection.csv (long format; column 'part' = a / b / c / descriptor).
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.api as sm
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e10_common as K


def cci(v, bt):
    """Centred bootstrap interval: v - (quantile - mean); removes the small downward bias of the
    leave-gene-out bootstrap replicates (see SUMMARY)."""
    bt = np.asarray(bt, float); m = np.nanmean(bt)
    return float(v - (np.nanquantile(bt, .975) - m)), float(v - (np.nanquantile(bt, .025) - m))
S = K.S
B_A, B_C = 2000, 500


def main():
    rows = []
    a2 = json.loads((K.R3 / "A2_axis_verification/results.json").read_text())["contexts"]
    for ci, ctx in enumerate(K.CTXS):
        c = K.load_cache(ctx); emb = K.load_emb(ctx); nl = c["nl"]; n = len(c["genes"])
        lnobs = np.log(np.maximum(c["nobs"], 1)); is_tf = nl["is_tf_lambert"].astype(bool)
        trtf = nl["is_trrust_tf"].astype(bool); trtg = nl["is_trrust_target_only"].astype(bool)
        rng = np.random.default_rng(K.SEED + 40 + ci)
        W = rng.multinomial(n, np.full(n, 1 / n), size=B_A).T.astype(float)
        a1 = S.auroc(lnobs[is_tf], lnobs[~is_tf]); b1 = K.wauroc_fixed(lnobs, is_tf, W)
        lab = trtf | trtg
        a2_ = S.auroc(lnobs[trtf], lnobs[trtg]); b2 = K.wauroc_fixed(lnobs[lab], trtf[lab], W[lab])
        rows.append(dict(part="a", context=ctx, test="lognobs_auroc_lambertTF_vs_nonTF", layer="", value=a1,
                         lo=np.quantile(b1, .025), hi=np.quantile(b1, .975), n_pos=int(is_tf.sum()), n_neg=int((~is_tf).sum())))
        rows.append(dict(part="a", context=ctx, test="lognobs_auroc_trrustTF_vs_targetonly", layer="", value=a2_,
                         lo=np.quantile(b2, .025), hi=np.quantile(b2, .975), n_pos=int(trtf.sum()), n_neg=int(trtg.sum())))
        info = a2[ctx]
        rows.append(dict(part="descriptor", context=ctx, test="detection_descriptors", layer="", value=np.nan,
                         median_nobs_lambertTF=float(np.median(c["nobs"][is_tf])), median_nobs_nonTF=float(np.median(c["nobs"][~is_tf])),
                         median_nobs_trrustTF=float(np.median(c["nobs"][trtf])), median_nobs_targetonly=float(np.median(c["nobs"][trtg])),
                         frac_TF_nobs_ge10=float(np.mean(c["nobs"][is_tf] >= 10)), n_V=n, n_cells=info["n_cells"],
                         cells_truncated_at_512=info["cells_with_nnz_gt_max_genes"], nnz_median=info["nnz_median"]))
        for L in range(K.N_LAYERS):
            En = K.unit(K.centred(emb, L, c["obs_idx"])); sc = S.tf_scores_loo(En, np.flatnonzero(is_tf))
            zs = (sc - sc.mean()) / sc.std(); zd = (lnobs - lnobs.mean()) / lnobs.std()
            Xd = sm.add_constant(np.c_[zs, zd, zs * zd])
            try:
                f = sm.Logit(is_tf.astype(float), Xd).fit(disp=0, maxiter=200)
                ci_ = f.conf_int(); names = ["const", "score", "lognobs", "score_x_lognobs"]
                for j, nm in enumerate(names[1:], start=1):
                    rows.append(dict(part="b", context=ctx, test=f"logit_{nm}", layer=L, value=float(f.params[j]),
                                     lo=float(ci_[j, 0]), hi=float(ci_[j, 1]), p=float(f.pvalues[j]), n_pos=int(is_tf.sum()), n_neg=int((~is_tf).sum())))
            except Exception as ex:
                rows.append(dict(part="b", context=ctx, test="logit_error", layer=L, value=np.nan, note=str(ex)))
        # (c) restricted universes
        for thr in (3, 5, 10):
            keep = np.flatnonzero(c["nobs"] >= thr); m = keep.size
            tfk = np.flatnonzero(is_tf[keep]); isk = is_tf[keep]
            stk = S.quantile_strata(lnobs[keep], 10)
            Wc = rng.poisson(1.0, (m, B_C)).astype(float)
            vals, boots = [], []
            for L in range(K.N_LAYERS):
                En = K.unit(K.centred(emb, L, c["obs_idx"][keep]))
                sc = S.tf_scores_loo(En, tfk); v = S.stratified_auroc(sc, isk, stk)
                bt = K.wstrat_auroc_var(K.loo_scores_weighted(En, tfk, Wc), isk, stk, Wc)
                vals.append(v); boots.append(bt)
                lo_, hi_ = cci(v, bt)
                rows.append(dict(part="c", context=ctx, test=f"S2_nobs_ge{thr}", layer=L, value=v, lo=lo_,
                                 hi=hi_, n_pos=int(isk.sum()), n_neg=int((~isk).sum())))
            bm = np.mean(boots, 0)
            lo_, hi_ = cci(float(np.mean(vals)), bm)
            rows.append(dict(part="c", context=ctx, test=f"S2_nobs_ge{thr}", layer="mean", value=float(np.mean(vals)),
                             lo=lo_, hi=hi_, n_pos=int(isk.sum()), n_neg=int((~isk).sum())))
        print(ctx, f"lognobs AUROC TF {a1:.3f} trrust {a2_:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(K.HERE / "detection.csv", index=False)


if __name__ == "__main__":
    main()
