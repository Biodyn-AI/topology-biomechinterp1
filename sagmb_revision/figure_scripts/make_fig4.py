"""Figure 4: H3 — which factor regulates which gene: per-layer estimate, sensitivity forest, and benchmark size vs power."""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from figstyle import *
E4A, E4B, E4C, E5, E8 = (R3 / x for x in ("E4a_edges", "E4b_grn", "E4c_probes", "E5_confounds", "E8_word2vec"))
SUM = pd.read_csv(E4A / "h3_summary.csv").set_index("context")

from matplotlib.lines import Line2D
fig = plt.figure(figsize=(7.4, 2.75))
axs = [fig.add_axes([0.065, 0.17, 0.235, 0.70]), fig.add_axes([0.535, 0.17, 0.215, 0.70]), fig.add_axes([0.835, 0.17, 0.155, 0.70])]

# a: per-layer within-TF AUROC with two-way bootstrap intervals
ax = axs[0]
for c in ("lung", "lung_external"):
    d = pd.read_csv(E4A / f"h3_primary_{c}.csv").sort_values("layer")
    ax.fill_between(d.layer, d.ci_lo, d.ci_hi, color=COLOR[c], alpha=0.15, lw=0)
    ax.plot(d.layer, d.auroc, color=COLOR[c], marker="o", markersize=2.5, label=LABEL[c])
    ax.axhline(SUM.loc[c, "auroc_cx_3k"], color=COLOR[c], lw=0.8, ls=(0, (1, 1.5)))
ax.axhline(0.5, color=MUTED, lw=0.8)
ax.axhline(0.60, color="#e34948", lw=0.9, ls=(0, (4, 2)))
ax.text(0, 0.602, "smallest effect of interest", ha="left", va="bottom", fontsize=5.8, color=INK2)
ax.axhline(SUM.loc["lung", "mde_auroc_sim_two_way_bootstrap"], color=MUTED, lw=0.8, ls=(0, (3, 2)))
ax.text(11, SUM.loc["lung", "mde_auroc_sim_two_way_bootstrap"] + 0.002, "detectable (80%)", ha="right", va="bottom", fontsize=5.8, color=INK2)
ax.text(0, SUM.loc["lung", "auroc_cx_3k"] - 0.003, "|co-expression|", ha="left", va="top", fontsize=5.8, color=INK2)
ax.set_ylim(0.47, 0.625); ax.set_xticks([0, 3, 6, 9, 11]); ax.set_xlabel("layer"); ax.set_ylabel("within-factor AUROC")
panel(ax, "a")

# b: forest of sensitivity rows (lung, layer 0 where layer-specific)
rows = []
p0 = pd.read_csv(E4A / "h3_primary_lung.csv").set_index("layer").loc[0]
rows.append(("scGPT cosine (primary)", p0.auroc, p0.ci_lo, p0.ci_hi, "scgpt"))
m = pd.read_csv(E4A / "sens_matching.csv")
r = m[(m.context == "lung") & (m.layer == 0) & (m.block == "strata TF x other-TF in-degree") & (m.scorer == "signed cosine full-512")].iloc[0]
rows.append(("  in-degree held fixed", r.auroc, r.ci_lo, r.ci_hi, "scgpt"))
r = pd.read_csv(E5 / "residual_auroc.csv"); r = r[(r.context == "lung") & (r.layer == 0)].iloc[0]
rows.append(("  residualised on covariates", r.auroc_resid, r.auroc_resid_lo, r.auroc_resid_hi, "scgpt"))
c = pd.read_csv(R3 / "E4a_chipfix" / "sens_chip_fixed.csv"); r = c[(c.context == "lung") & (c.layer == 0) & (c.scorer == "signed cosine full-512")].iloc[0]
rows.append(("  ChIP-seq positives", r.auroc, r.ci_lo, r.ci_hi, "scgpt"))
s = pd.read_csv(E4A / "sens_strata.csv"); s = s[(s.context == "lung") & (s.layer == 0) & (s.scorer == "signed cosine full-512")]
r = s[s.stratum == "PMID count >= 2"].iloc[0]; rows.append(("  $\\geq$2 publications", r.auroc, r.ci_lo, r.ci_hi, "scgpt"))
r = s[s.stratum.str.startswith("TF n_obs (scGPT input) >= 10%")].iloc[0]; rows.append(("  factors often in input", r.auroc, r.ci_lo, r.ci_hi, "scgpt"))
k = SUM.loc["lung"]; rows.append(("|co-expression|", k.auroc_cx_3k, k.cx_3k_ci_lo, k.cx_3k_ci_hi, "base"))
g = pd.read_csv(E4B / "grn_summary.csv"); g = g[(g.context == "lung") & (g.seed.astype(str) == "ens")]
for meth, lab in (("genie3", "GENIE3"), ("grnboost2", "GRNBoost2")):
    x = g[g.method == meth]
    if len(x):
        x = x[x.tier == ("T2" if (x.tier == "T2").any() else "T1")].iloc[0]
        rows.append((lab, x.auroc, x.ci_lo, x.ci_hi, "base"))
w = pd.read_csv(E8 / "w2v_results.csv")
w = w[(w.context == "lung") & (w.corpus == "full_w5_e10") & (w.task.str.startswith("H3")) & (w.seed.astype(str) == "mean")]
if len(w): rows.append(("co-occurrence (Word2Vec)", w.value.iloc[0], w.ci_lo.iloc[0], w.ci_hi.iloc[0], "base"))
pr = pd.read_csv(E4C / "probe_results.csv")
pr = pr[(pr.context == "lung") & (pr.model == "M1") & (pr.source == "scgpt") & (pr.regime == "C3") & (pr.layer == 0)].set_index("feature")
for feat, lab in (("tg_only", "probe: target only"), ("bilinear", "probe: bilinear interaction")):
    r = pr.loc[feat]; rows.append((lab, r.auc_mean, r.auc_ci_lo, r.auc_ci_hi, "probe"))
ct = pd.read_csv(E4A / "sens_contextual.csv"); r = ct[(ct.context == "immune_Bcell") & (ct.layer == 0) & (ct.scorer == "signed cosine full-512")].iloc[0]
rows.append(("B cells, active factors", r.auroc, r.ci_lo, r.ci_hi, "ctx"))
ax = axs[1]; y = np.arange(len(rows))[::-1]
col = {"scgpt": COLOR["lung"], "base": INK2, "probe": "#6250d6", "ctx": COLOR["immune"]}
for yy, (lab, v, lo, hi, kind) in zip(y, rows):
    ax.plot([lo, hi], [yy, yy], color=col[kind], lw=1.2, solid_capstyle="round")
    ax.scatter([v], [yy], s=14, color=col[kind], zorder=3, edgecolor="white", linewidth=0.4)
ax.axvline(0.5, color=MUTED, lw=0.8); ax.axvline(0.60, color="#e34948", lw=0.9, ls=(0, (4, 2)))
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=6.2); ax.grid(axis="y", visible=False)
ax.set_xlim(0.44, 0.73); ax.set_xlabel("within-factor AUROC ($L_0$)")
fig.text(0.355, 0.895, "b", fontsize=10, fontweight="bold", va="bottom")

# c: minimum detectable AUROC against the number of documented pairs
ax = axs[2]
pts = [(c, SUM.loc[c, "n_pos"], SUM.loc[c, "mde_auroc_sim_two_way_bootstrap"]) for c in ("lung", "lung_external", "kidney_B", "immune", "kidney_A")]
for c, n, mde in pts:
    if c == "lung_external":  # open marker, so that it does not hide the lung point (5,093 vs 5,266 pairs, 0.562 vs 0.564)
        ax.scatter([n], [mde], s=30, facecolor="none", edgecolor=COLOR[c], linewidth=1.0, zorder=4)
    else:
        ax.scatter([n], [mde], s=18, color=COLOR[c], zorder=3, edgecolor="white", linewidth=0.4)
# 80/20 reference split of the kidney_A pairs: pooled AUROC (positives against all held-out non-edges), two-way resampling, from power.csv
# (the raw experiment file keeps its original design/statistic labels; only the public Table_S7 CSV labels differ)
pw = pd.read_csv(E4A / "power.csv")
sp = pw[pw.design.str.startswith("kidney_A legacy") & pw.statistic.str.startswith("pooled AUROC") & (pw.test == "two_way_bootstrap")]
mde_split, n_split = float(sp.mde_auroc_sim.iloc[0]), int(sp.n_pos.iloc[0])
ax.scatter([n_split], [mde_split], s=22, marker="D", color=INK2, zorder=3, edgecolor="white", linewidth=0.4)
ax.annotate(f"80/20 split,\n{n_split} positives", (n_split, mde_split), xytext=(7, 0), textcoords="offset points", fontsize=5.8, color=INK2, va="center")
ax.axhline(0.60, color="#e34948", lw=0.9, ls=(0, (4, 2)))
ax.set_xscale("log"); ax.set_xlim(30, 9000); ax.set_ylim(0.54, 0.73); ax.minorticks_off()
ax.set_xticks([50, 200, 1000, 5000]); ax.set_xticklabels(["50", "200", "1,000", "5,000"])
ax.set_xlabel("documented pairs"); ax.set_ylabel("minimum detectable AUROC"); panel(ax, "c")
hd = [Line2D([], [], color=COLOR[c], ls=STYLE[c], marker="o", markersize=3) for c in CTX]
fig.legend(hd, [LABEL[c] for c in CTX], ncol=5, loc="upper left", bbox_to_anchor=(0.06, 1.06), handlelength=2.2, columnspacing=1.4)
fig.savefig(OUT / "Figure_4.pdf"); fig.savefig(OUT / "Figure_4.png", dpi=200)
print("ok", [(r[0].strip(), round(r[1], 3)) for r in rows], "split MDE", round(mde_split, 3))
