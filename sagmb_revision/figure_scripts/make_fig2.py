"""Figure 2: H1 — physical interaction partners share spectral coordinates."""
import sys
import numpy as np, pandas as pd
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from figstyle import *

BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
cmap = LinearSegmentedColormap.from_list("seqblue", BLUE)
fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 0.05], hspace=0.55, wspace=0.35)

for k, c in enumerate(["lung", "lung_external"]):
    d = pd.read_csv(R3 / "E2_ppi" / f"h1_primary_{c}.csv")
    M = d.pivot(index="sv", columns="layer", values="r").values
    ax = fig.add_subplot(gs[0, k])
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=0.45, aspect="auto", origin="upper")
    ax.set_xticks(range(0, 12, 2)); ax.set_xticklabels(range(0, 12, 2)); ax.set_yticks(range(10)); ax.set_yticklabels([f"SV{i}" for i in range(1, 11)])
    ax.set_xlabel("layer"); ax.grid(False); ax.set_title(LABEL[c] + f"  (120/120 significant)", fontsize=7.5, color=INK2)
    for s in ax.spines.values(): s.set_visible(False)
    panel(ax, "ab"[k])
cax = fig.add_subplot(gs[0, 2]); cb = fig.colorbar(im, cax=cax); cb.set_label("assortativity $r$", fontsize=7); cb.outline.set_visible(False)

# c: top-cell observed r vs null means under each null (lung L3 SV1; lung_external L8 SV2)
nn = pd.read_csv(R3 / "E2_ppi" / "sens_nulls_summary.csv")
f1 = pd.read_csv(R3 / "F1_h1_compartment" / "strat_null.csv")
top = {"lung": (3, 1), "lung_external": (8, 2)}
rows = [("label shuffle", "label"), ("degree strata (V)", "deg_v"), ("degree strata (genome)", "deg_g"), ("detection strata", "det"),
        ("rewiring (MS)", "maslov_sneppen"), ("compartment strata", "strat"), ("class-pair rewiring", "cpms")]
ax = fig.add_subplot(gs[1, 0])
for j, c in enumerate(["lung", "lung_external"]):
    L, k = top[c]; ff = f1[(f1.context == c) & (f1.layer == L) & (f1.sv == k)].iloc[0]
    obs = ff.r; vals = []
    for lab, key in rows:
        if key in ("strat", "cpms"): vals.append(ff[f"{key}_null_mean"])
        elif key == "maslov_sneppen": vals.append(nn[(nn.context == c) & (nn.null == key)].top_null_mean.iloc[0])
        else: vals.append(nn[(nn.context == c) & (nn.null == key)].ms_top_null_mean_this.iloc[0])
    y = np.arange(len(rows)) + (0.18 if j == 0 else -0.18)
    ax.scatter(vals, y, s=14, color=COLOR[c], zorder=3, label=LABEL[c] + " null mean")
    ax.scatter([obs] * len(rows), y, s=30, marker="|", linewidths=1.5, color=COLOR[c], zorder=3)
    for yy, v in zip(y, vals): ax.plot([v, obs], [yy, yy], color=COLOR[c], lw=0.8, alpha=0.5)
ax.set_yticks(np.arange(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=6.5); ax.invert_yaxis()
ax.set_xlim(-0.03, 0.47); ax.set_xlabel("$r$ (null mean $\\bullet$, observed $|$)"); panel(ax, "c")
ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2, fontsize=6, handletextpad=0.2, columnspacing=1.0, borderaxespad=0.1)

# d: co-pole observed/null ratio vs K/n (rewiring null; delete-d rescaled CI)
kk = pd.read_csv(R3 / "E2_ppi" / "sens_copole_K.csv"); kk = kk[(kk.null == "ms") & (kk.top_rank == 1)]
ax = fig.add_subplot(gs[1, 1])
for c in ["lung", "lung_external", "kidney_B"]:
    d = kk[kk.context == c].sort_values("pole_frac"); e = d.ratio_obs_null.values
    # the vertex-subsampling interval is for the configuration-model ratio: rescale it about its own centre
    m = d.ratio_analytic_est.values
    lo = m - np.sqrt(2) * (m - d.ratio_ci_lo.values); hi = m + np.sqrt(2) * (d.ratio_ci_hi.values - m)
    ax.fill_between(d.pole_frac, lo, hi, color=COLOR[c], alpha=0.15, lw=0)
    ax.plot(d.pole_frac, e, STYLE[c], color=COLOR[c], marker="o", markersize=3, label=LABEL[c])
ax.axhline(1, color=MUTED, lw=0.8, ls=(0, (2, 2))); ax.axvline(0.27, color=MUTED, lw=0.6, ls=":")
ax.set_xticks([0.1, 0.2, 0.27, 0.33, 0.5]); ax.set_xticklabels(["0.10", "0.20", "0.27", "0.33", "0.50"], fontsize=6)
ax.set_xlabel("pole size $K/n$"); ax.set_ylabel("co-pole rate / rewiring null"); ax.set_ylim(0.9, 2.6); panel(ax, "d")
ax.legend(loc="upper right", fontsize=6, borderaxespad=0.1)
fig.savefig(OUT / "Figure_2.pdf"); fig.savefig(OUT / "Figure_2.png", dpi=200); print("ok")
