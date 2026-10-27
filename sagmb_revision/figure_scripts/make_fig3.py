"""Figure 3: H2 — transcription-factor identity across layers, bands, and leading-direction removal."""
import sys, re
import numpy as np, pandas as pd
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from figstyle import *
E3 = R3 / "E3_node"
BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
cmap = LinearSegmentedColormap.from_list("seqblue", BLUE)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=0.95, width_ratios=[1.15, 1, 1]))

ax = axs[0]
for c in CTX:
    d = pd.read_csv(E3 / f"h2_primary_{c}.csv"); d = d[d.row_type == "H2_layer"].sort_values("layer")
    ax.plot(d.layer, d.strat_auroc, STYLE[c], color=COLOR[c], label=LABEL[c])
    if c in ("lung", "lung_external"): ax.fill_between(d.layer, d.strat_lo, d.strat_hi, color=COLOR[c], alpha=0.15, lw=0)
b = pd.read_csv(E3 / "h2_primary_lung.csv"); b = b[b.row_type == "baseline"].set_index("baseline")
for key, lab in (("static_esm2", "ESM2-15B (sequence)"), ("static_scgpt", "scGPT input emb."), ("log_nobs", "detection alone")):
    v = b.loc[key, "strat_auroc"]; ax.axhline(v, color=MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.text(11, v + 0.006, lab, ha="right", va="bottom", fontsize=6, color=INK2)
ax.set_ylim(0.45, 0.92); ax.set_xticks([0, 3, 6, 9, 11]); ax.set_xlabel("layer"); ax.set_ylabel("TF AUROC (detection-stratified)")
panel(ax, "a")

# b: band scan, lung, max over layers, start i (rows) x end j (cols)
h = pd.read_csv(E3 / "partial" / "heatmap_cos_lung.csv").set_index("band")
M = np.full((12, 12), np.nan)
for band, row in h.iterrows():
    m = re.match(r"SV(\d+)(?:-(\d+))?$", band); i = int(m.group(1)); j = int(m.group(2) or i)
    M[i - 1, j - 1] = row.values.max()
ax = axs[1]; im = ax.imshow(M, cmap=cmap, vmin=0.5, vmax=0.68, origin="upper")
ax.set_xticks(range(0, 12, 2)); ax.set_xticklabels([f"{k+1}" for k in range(0, 12, 2)]); ax.set_yticks(range(0, 12, 2)); ax.set_yticklabels([f"{k+1}" for k in range(0, 12, 2)])
ax.set_xlabel("band ends at SV"); ax.set_ylabel("band starts at SV"); ax.grid(False)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.040, pad=0.02); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6)
cb.set_label("TF AUROC (stratified)", fontsize=6.5, color=INK2)
full = pd.read_csv(E3 / "h2_primary_lung.csv"); full = full[full.row_type == "H2_layer"].strat_auroc.max()
ax.set_title(f"lung, best layer per band\n(full 512-d: {full:.3f})", fontsize=6.5, color=INK2)
# letter b at the same height as a and c (the square heatmap axes sit lower than its neighbours)
fig.canvas.draw()
to_fig = fig.transFigure.inverted()
y_b = to_fig.transform(axs[0].transAxes.transform((0, 1.06)))[1]
x_b = to_fig.transform(ax.transAxes.transform((-0.32, 0)))[0]
fig.text(x_b, y_b, "b", fontsize=10, fontweight="bold", va="bottom", ha="left")

# c: removing the leading j directions (E7)
r = pd.read_csv(R3 / "E7_rank" / "remove_top_j.csv"); r = r[r.task == "NODE"]
ax = axs[2]
for c in ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]:
    d = r[r.ctx == c].copy(); d["j"] = d.rep.str[1:].astype(int); d = d.sort_values("j")
    ax.plot(np.r_[0, d.j], np.r_[0, d.diff_vs_full], STYLE[c], color=COLOR[c], marker="o", markersize=2.5, label=LABEL[c])
ax.axhline(0, color=MUTED, lw=0.8, ls=(0, (2, 2))); ax.set_xticks([0, 1, 2, 3, 5, 10])
ax.set_xlabel("directions removed"); ax.set_ylabel("change in TF AUROC"); panel(ax, "c")
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, ncol=5, loc="upper left", bbox_to_anchor=(0.06, 1.12), handlelength=2.2, columnspacing=1.4)
fig.savefig(OUT / "Figure_3.pdf"); fig.savefig(OUT / "Figure_3.png", dpi=200); print("ok", np.nanmax(M))
