"""Verifier check 4: per fit, count targets whose importance row is all zero although the target was fitted
(fitted == 1, non-constant expression): the model's total feature importance was 0 (no informative split).
Also: fraction of h3 pairs whose target falls in such a row, and n_iter distribution for GRNBoost2. Writes v4_zero_rows.csv."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
E = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(E))
import e4b_common as C
rows = []
for ctx in C.CONTEXTS:
    s = C.load_sets(ctx); h3t = s["tg"]
    for tier in ("T1", "T2"):
        for m in C.METHODS:
            for sd in C.SEEDS[tier]:
                z = np.load(E / "partial" / f"{ctx}_{tier}_{m}_s{sd}.npz")
                zero = z["imp"].sum(1) == 0; fitted = z["fitted"] == 1
                zr = z["targets"][zero & fitted]
                rows.append(dict(context=ctx, tier=tier, method=m, seed=sd, n_targets=int(z["targets"].size),
                                 n_unfitted_constant=int((~fitted).sum()), n_fitted_zero_importance=int((zero & fitted).sum()),
                                 frac_h3_pairs_in_zero_rows=float(np.isin(h3t, z["targets"][zero]).mean()),
                                 n_iter_min=int(z["n_iter"].min()) if m == "grnboost2" else None,
                                 n_iter_median=float(np.median(z["n_iter"])) if m == "grnboost2" else None,
                                 n_iter_le25=int((z["n_iter"] <= 26).sum()) if m == "grnboost2" else None))
df = pd.DataFrame(rows); df.to_csv(Path(__file__).parent / "v4_zero_rows.csv", index=False)
pd.set_option("display.width", 250); print(df.to_string())
