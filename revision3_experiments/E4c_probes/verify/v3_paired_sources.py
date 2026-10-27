"""E4c verifier, part 3: paired same-fold comparison of sources (CONVENTION3 §6: comparisons between models are
paired differences on matched pairs, never 'significant here, not there').
For C3 M1 (5 repeats x 5 folds, identical folds): scGPT bilinear at layer L minus static bilinear, minus ESM2-15B
bilinear; and the same for the target-only-adjusted excess (probe - target-only of its own source).
Nadeau-Bengio corrected resampled t (as the implementer). Bookkeeping on probe_folds.csv only.
Output: verify/v3_paired_sources.csv
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import t as tdist
H = Path(__file__).resolve().parent; PF = pd.read_csv(H.parent / "probe_folds.csv")
PF = PF[(PF.model == "M1") & (PF.regime == "C3")]


def nb(d, ratio):
    J = d.size; m = d.mean(); se = np.sqrt((1 / J + ratio) * d.var(ddof=1)); q = tdist.ppf(0.975, J - 1)
    return m, m - q * se, m + q * se


def get(ctx, src, feat, layer):
    g = PF[(PF.context == ctx) & (PF.source == src) & (PF.feature == feat) & (PF.layer == layer)].sort_values(["repeat", "fold"])
    return g


rows = []
for ctx in ("lung", "lung_external", "kidney_B"):
    for other in ("static", "esm2"):
        ob = get(ctx, other, "bilinear", -1); otg = get(ctx, other, "tg_only", -1)
        for l in range(12):
            sb = get(ctx, "scgpt", "bilinear", l); stg = get(ctx, "scgpt", "tg_only", l)
            assert (sb[["repeat", "fold"]].values == ob[["repeat", "fold"]].values).all()
            ratio = sb.n_test.mean() / sb.n_train.mean()
            m, lo, hi = nb(sb.auc.values - ob.auc.values, ratio)
            ex_s = sb.auc.values - np.fmax(stg.auc.values, 0.5); ex_o = ob.auc.values - np.fmax(otg.auc.values, 0.5)
            m2, lo2, hi2 = nb(ex_s - ex_o, ratio)
            rows.append(dict(context=ctx, comparison=f"scgpt_bilinear_L{l} - {other}_bilinear", layer=l, other=other,
                             delta_auc=m, lo=lo, hi=hi, delta_excess_over_own_tg=m2, lo_excess=lo2, hi_excess=hi2, J=sb.shape[0]))
R = pd.DataFrame(rows); R.to_csv(H / "v3_paired_sources.csv", index=False)
print(R.groupby(["context", "other"]).agg(dmin=("delta_auc", "min"), dmax=("delta_auc", "max"),
                                           n_lo_gt0=("lo", lambda x: int((x > 0).sum())), n_hi_lt0=("hi", lambda x: int((x < 0).sum())),
                                           exmin=("delta_excess_over_own_tg", "min"), exmax=("delta_excess_over_own_tg", "max"),
                                           n_ex_lo_gt0=("lo_excess", lambda x: int((x > 0).sum()))).round(4).to_string())
print(R[(R.context == "lung") & (R.layer == 2)].round(4).to_string())
