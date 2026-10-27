"""Verifier check of the fixed-score gene bootstrap (implementer's flagged limitation).

For the best layer of lung / lung_external and for immune L2 (k ~ 113 TFs), compare the reported fixed-score
gene-bootstrap SE (strat_hi - strat_lo)/3.92 with a re-estimating delete-d (half-sampling) jackknife SE: 300 random
halves without replacement (m = n/2), the TF centroid, the LOO scores and the detection deciles' strata membership
(strata = the V deciles) re-estimated within each half; SE = sqrt(m/(n-m)) * SD(half estimates) = SD at m = n/2
(CONVENTION3_ERRATA item 3 scaling). Writes verify/verify_e3c.csv. Descriptive only.
"""
import numpy as np, pandas as pd
from pathlib import Path
from verify_common import strat, loo, ctx_data
HERE = Path(__file__).resolve().parent; E3 = HERE.parent
rows = []
for ctx, l in (("lung", 5), ("lung_external", 1), ("immune", 2), ("immune_CD4T", 2)):
    emb, obs, genes, cnt, y, det = ctx_data(ctx)
    X = np.asarray(emb[l], float)[obs]; X = X - X.mean(0); n = y.size; m = n // 2
    full = strat(loo(X, y), y, det)
    rng = np.random.default_rng(31); est = []
    for _ in range(300):
        s = np.sort(rng.choice(n, m, replace=False)); est.append(strat(loo(X[s], y[s]), y[s], det[s]))
    est = np.array(est); se_dd = np.sqrt(m / (n - m)) * est.std(ddof=1)
    r = pd.read_csv(E3 / f"h2_primary_{ctx}.csv"); r = r[(r.row_type == "H2_layer") & (r.layer == l)].iloc[0]
    rows.append(dict(context=ctx, layer=l, n=n, n_tf=int(y.sum()), strat_auroc=full, reported=float(r.strat_auroc),
                     se_fixed_score_boot=(r.strat_hi - r.strat_lo) / 3.92, se_deleted_half_reestimated=se_dd,
                     ratio=se_dd / ((r.strat_hi - r.strat_lo) / 3.92), half_mean=est.mean()))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(HERE / "verify_e3c.csv", index=False); print("DONE")
