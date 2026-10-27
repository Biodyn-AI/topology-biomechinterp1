"""Chain diagnostics: acceptance, lag-1 autocorrelation of the null r of the E2 top cell, null-max |z| 95th pct."""
import numpy as np, pandas as pd, f1_common as F
E = F.E; rows = []
for ctx in ["lung", "lung_external", "kidney_B"]:
    c = F.col_of(*F.E2_TOP[ctx][0])
    for step, key, acckey in [("ms", "null_raw", "acc"), ("cpms", "null", "acc"), ("perm", "null", None),
                              ("wb", "within_null_raw", "within_acc"), ("wb", "between_null_raw", "between_acc"),
                              ("nobs3", "null", "acc"), ("nobs5", "null", "acc"), ("nobs10", "null", "acc")]:
        z = np.load(F.CKPT / f"{ctx}_{step}.npz"); nl = z[key]
        zz = (nl - nl.mean(0)) / nl.std(0, ddof=1)
        rows.append(dict(context=ctx, null=f"{step}:{key}", N=nl.shape[0], acc_mean=float(z[acckey].mean()) if acckey else np.nan,
                         lag1_top_cell=E.lag1(nl[:, c]), lag1_median_over_cells=float(np.median([E.lag1(nl[:, j]) for j in range(nl.shape[1])])),
                         null_max_abs_z_q95=float(np.quantile(np.abs(zz).max(1), 0.95))))
df = pd.DataFrame(rows); df.to_csv(F.HERE / "chain_diagnostics.csv", index=False); print(df.round(3).to_string())
