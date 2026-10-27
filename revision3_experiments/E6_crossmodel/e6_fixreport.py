"""E6 (resume agent): effect of the exact-sign 1-D-band fix on estimator (iii). Writes p1d_fix_effect.csv."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e6lib as L
J = [j for j, (a, b) in enumerate(L.BANDS) if a == b]
rows = []
for f in sorted(L.WORK.glob("*/*.p1d.npz")):
    if f.name.startswith("._"): continue
    z = np.load(f); setname = f.parent.name; emb = f.name[:-len(".p1d.npz")]
    o = np.load(f.parent / f"{emb}.npz")
    for t in ("node_lambert", "node_trrust", "ppi", "edge"):
        if f"tobs_{t}" not in z.files: continue
        to0, tn0 = o[f"tobs_{t}"].astype(float), o[f"tnull_{t}"].astype(float)
        to1, tn1 = z[f"tobs_{t}"].astype(float), z[f"tnull_{t}"].astype(float)
        _, pa0, _ = L.S.max_t_adjust(to0, tn0); _, pa1, _ = L.S.max_t_adjust(to1, tn1)
        j0, j1 = int(np.argmax(to0)), int(np.argmax(to1))
        rows.append(dict(set=setname, embedding=emb, task=t,
                         max_abs_change_1d_obs=float(np.max(np.abs(to1[J] - to0[J]))),
                         mean_change_1d_obs=float(np.mean(to1[J] - to0[J])),
                         best_band_old=L.BAND_NAMES[j0], best_band_new=L.BAND_NAMES[j1], best_value_old=float(to0[j0]),
                         best_value_new=float(to1[j1]), p_adj_best_old=float(pa0[j0]), p_adj_best_new=float(pa1[j1]),
                         n_padj05_old=int((pa0 < 0.05).sum()), n_padj05_new=int((pa1 < 0.05).sum()),
                         best_changed=bool(j0 != j1), old_best_was_1d=bool(j0 in J)))
df = pd.DataFrame(rows); df.to_csv(L.HERE / "p1d_fix_effect.csv", index=False)
print(df.groupby("task").agg(n=("best_changed", "size"), n_best_changed=("best_changed", "sum"),
                             n_old_best_1d=("old_best_was_1d", "sum"), max_abs_change=("max_abs_change_1d_obs", "max")).to_string())
x = df[df.best_changed]
print(x[["set", "embedding", "task", "best_band_old", "best_value_old", "p_adj_best_old", "best_band_new", "best_value_new", "p_adj_best_new"]].round(4).to_string())
print("rows where p_adj of best crosses 0.05:", int(((df.p_adj_best_old < 0.05) != (df.p_adj_best_new < 0.05)).sum()))
