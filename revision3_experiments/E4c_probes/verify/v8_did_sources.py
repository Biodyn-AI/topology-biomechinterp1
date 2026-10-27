"""E4c verifier part 8: difference-in-differences under rewiring for the layer-free sources (static, ESM2-15B) and
co-expression-free check, same construction as v1 (statistic = fold-mean AUROC(bilinear) - AUROC(concat), true labels vs the
same difference on each of the 100 rewired draws, identical folds; C3 repeat 0). Output: v8_did_sources.csv"""
from pathlib import Path
import numpy as np, pandas as pd
H = Path(__file__).resolve().parent; E = H.parent
def nm(ctx, src, feat):
    d = pd.read_csv(E / "partial" / "m1" / ctx / f"{src}_{feat}_L-1_C3_r0_null.csv")
    return d.pivot_table(index="fold", columns="label_col", values="auc").mean(0).values
rows = []
for ctx in ("lung", "lung_external", "kidney_B"):
    for src in ("static", "esm2"):
        D = nm(ctx, src, "bilinear") - nm(ctx, src, "concat")
        obs, nul = D[0], D[1:]
        rows.append(dict(context=ctx, source=src, did_obs=obs, null_mean=nul.mean(), null_sd=nul.std(ddof=1),
                         p_emp=(np.sum(nul >= obs) + 1) / (nul.size + 1), N=nul.size))
R = pd.DataFrame(rows); R.to_csv(H / "v8_did_sources.csv", index=False); print(R.round(4).to_string(index=False))
