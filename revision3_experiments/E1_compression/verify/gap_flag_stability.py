"""Monte Carlo stability of the near-degeneracy flag (basic CI lower bound <= 1) at B = 200.
SE of the bootstrap 97.5% quantile ~ boot_sd * sqrt(p(1-p)/B) / phi(z_p) (normal approximation).
A flag is called MC-unstable if |basic_lo - 1| < 2 * SE. Writes verify/gap_flag_stability.csv."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm
H = Path(__file__).resolve().parent
g = pd.read_csv(H.parent / "gaps.csv")
p, B = 0.975, 200
se = g.boot_sd * np.sqrt(p * (1 - p) / B) / norm.pdf(norm.ppf(p))
g["mc_se_q975"] = se
g["mc_unstable"] = (g.basic_lo - 1).abs() < 2 * se
rows = []
for name, sub in [("all", g), ("k<=3", g[g.k <= 3]),
                  ("core k<=3", g[(g.k <= 3) & g.context.isin(["lung", "lung_external", "kidney_B", "immune", "kidney_A"])])]:
    rows.append(dict(subset=name, n_pairs=len(sub), n_flagged=int(sub.near_degenerate.sum()),
                     n_flagged_mc_unstable=int((sub.near_degenerate & sub.mc_unstable).sum()),
                     n_unflagged_mc_unstable=int((~sub.near_degenerate & sub.mc_unstable).sum()),
                     n_flagged_stable=int((sub.near_degenerate & ~sub.mc_unstable).sum())))
out = pd.DataFrame(rows); out.to_csv(H / "gap_flag_stability.csv", index=False); print(out.to_string())
c = g[(g.k <= 3) & g.context.isin(["lung", "lung_external", "kidney_B", "immune", "kidney_A"]) & (g.near_degenerate | g.mc_unstable)]
print(c[["context", "layer", "k", "ratio", "basic_lo", "mc_se_q975", "near_degenerate", "mc_unstable"]].round(4).to_string())
c.to_csv(H / "gap_flag_stability_core_k1to3.csv", index=False)
