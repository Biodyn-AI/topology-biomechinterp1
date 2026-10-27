"""V12: does the power-matched gap to immune exceed immune's own sampling error? gap / sqrt(SE_immune^2 + SD_matched^2),
SE_immune = E10 gene-bootstrap SE of immune's layer mean (partial/stats_immune.npz), SD_matched = SD over the R = 100
power-matched subsamples (partial/powermatch_<ctx>.npz). Also with immune's S1 replaced by the 20-draw rewiring mean
(ck_v11.json). Writes verify/powermatch_z.csv."""
import json
from pathlib import Path
import numpy as np, pandas as pd
H = Path(__file__).resolve().parent; P = H.parent / "partial"
zi = np.load(P / "stats_immune.npz"); names = [str(s) for s in zi["stats"]]
im = zi["point"].mean(1); se = np.nanstd(zi["boot"].mean(1), axis=1, ddof=1)
v11 = json.loads((H / "ck_v11.json").read_text())
rows = []
for ctx in ["lung", "lung_external", "kidney_B"]:
    r = np.load(P / f"powermatch_{ctx}.npz")["res"].mean(2)
    for i, nm in enumerate(names):
        gap = im[i] - r[:, i].mean(); sdm = r[:, i].std(ddof=1); z = gap / np.sqrt(se[i] ** 2 + sdm ** 2)
        row = dict(ctx=ctx, stat=nm, immune=im[i], immune_boot_se=se[i], matched_mean=r[:, i].mean(), matched_sd=sdm, gap=gap, z=z)
        if nm == "S1":
            g2 = v11["immune"]["mean"] - r[:, i].mean()
            row.update(immune_S1_20draw_mean=v11["immune"]["mean"], gap_20draw=g2, z_20draw=g2 / np.sqrt(se[i] ** 2 + sdm ** 2 + v11["immune"]["sd"] ** 2))
        rows.append(row)
df = pd.DataFrame(rows); df.to_csv(H / "powermatch_z.csv", index=False); print(df.round(4).to_string())
