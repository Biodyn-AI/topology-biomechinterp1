"""E3 addendum (resumed session): baseline rows for the restricted-V sensitivity (7).

For V restricted to n_obs >= 3, 5, 10 (same subsets and re-computed detection deciles as e3_primary.py),
the detection / mean-expression baselines: unstratified and detection-stratified AUROC with gene-bootstrap
95% CI (B = 1,000). Descriptive only (no significance claim, no permutation).
writes partial/restricted_baselines.csv
"""
import numpy as np, pandas as pd
import e3_lib as L
S = L.S
rows = []
for ci, ctx in enumerate(L.CONTEXTS):
    D = L.load(ctx); y = D["is_tf"]
    lmx, lme = L.log_mean(D["mean_extr"]), L.log_mean(D["mean_expr"])
    for thr in (3, 5, 10):
        sub = D["nobs"] >= thr; ys = y[sub]; ds = S.quantile_strata(D["log_nobs"][sub], 10)
        W = L.boot_weights(int(sub.sum()), 1000, np.random.default_rng(20261026 + 1000 * ci + 900 + thr))
        for nm, sc in (("log_nobs", D["log_nobs"]), ("log_mean_expr_extraction", lmx), ("log_mean_expr_population", lme)):
            s = sc[sub]; bs = L.boot_strat_fast(s, ys, ds, W)
            rows.append(dict(context=ctx, min_nobs=thr, row_type="baseline", baseline=nm, n_genes=int(sub.sum()),
                             n_tf=int(ys.sum()), strat_auroc=S.stratified_auroc(s, ys, ds),
                             strat_lo=float(np.quantile(bs[:, 0], .025)), strat_hi=float(np.quantile(bs[:, 0], .975)),
                             auroc=S.auroc(s[ys], s[~ys]), auroc_lo=float(np.quantile(bs[:, 1], .025)),
                             auroc_hi=float(np.quantile(bs[:, 1], .975))))
    print(ctx, flush=True)
pd.DataFrame(rows).to_csv(L.HERE / "partial" / "restricted_baselines.csv", index=False)
print("DONE")
