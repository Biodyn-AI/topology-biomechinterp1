import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C
plt.rcParams.update({"font.size":8,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":200})
OUT = C.OUT/"figures"; OUT.mkdir(exist_ok=True)
reb = pd.read_csv(C.OUT/"headline_rebuilt.csv")
ALPHA = C.BONFERRONI_ALPHA

# Fig 1: co-pole significance across depth, all axes
# Layout note: the two reference lines in panel (b) are labelled in a clear right-hand
# gutter (data stop at layer 11, the axes run to 13.8) so no annotation sits over a
# series, and the series legend is shared and placed below both panels rather than
# inside panel (a).
SERIES = [("SV2_STRING_copole", "SV2 STRING", "#4C78A8"),
          ("SV3_STRING_copole", "SV3 STRING", "#F58518"),
          ("SV4_STRING_copole", "SV4 STRING", "#54A24B"),
          ("SV2_TRRUST_copole", "SV2 TRRUST", "#E45756")]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
handles = []
for tag, lbl, col in SERIES:
    d = reb[reb.claim.str.startswith(tag)].copy()
    d["L"] = d.claim.str.extract(r"_L(\d+)$").astype(int)
    d = d.sort_values("L")
    h, = ax[0].plot(d.L, d.stat, "o-", ms=3, lw=1.2, color=col, label=lbl)
    handles.append(h)
    ax[1].semilogy(d.L, d.p, "o-", ms=3, lw=1.2, color=col, label=lbl)

null_mean = float(reb[reb.claim.str.startswith("SV2_STRING_copole")].null.mean())
hnull = ax[0].axhline(null_mean, ls="--", lw=1, color="0.45", label="label-shuffle null")
handles.append(hnull)
ax[0].set_xlim(-0.5, 11.5)
ax[0].set_ylim(0.095, 0.355)
ax[0].set_xlabel("layer"); ax[0].set_ylabel("co-pole rate")
ax[0].set_title("a  Co-pole rate across depth", loc="left", fontsize=8)

ax[1].axhline(ALPHA, ls="--", lw=1, color="#B0304A")
ax[1].axhline(1e-5, ls=":", lw=1.1, color="0.45")
ax[1].set_xlim(-0.5, 13.8)          # gutter from 11.5 onward carries the line labels
ax[1].set_ylim(4.5e-6, 3.0)
ax[1].axvline(11.6, lw=0.8, color="0.85")     # visual divider before the label gutter
ax[1].text(13.6, ALPHA * 1.55, "Bonferroni\n$\\alpha=2.7\\times10^{-4}$",
           fontsize=6, color="#B0304A", ha="right", va="bottom", linespacing=1.25)
ax[1].text(13.6, 1.35e-5, "permutation floor\n$N=100{,}000$",
           fontsize=6, color="0.35", ha="right", va="bottom", linespacing=1.25)
ax[1].set_xticks([0, 2, 4, 6, 8, 10])
ax[1].set_xlabel("layer"); ax[1].set_ylabel("permutation $p$")
ax[1].set_title("b  Significance, N = 100,000", loc="left", fontsize=8)

fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=6.5,
           frameon=False, bbox_to_anchor=(0.5, -0.015), columnspacing=1.4, handlelength=1.8)
fig.tight_layout(rect=(0, 0.085, 1, 1))
fig.savefig(OUT/"Figure_1.pdf"); plt.close(fig)

# Fig 2: TF/target AUROC across depth + confound residual
# Layout: the round-1 reference line is labelled in a right-hand gutter clear of the
# series, and the legends sit in empty regions rather than over data.
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
d = reb[reb.claim.str.startswith("TF_target_AUROC_SV2_7_L")].copy()
d["L"] = d.claim.str.extract(r"_L(\d+)$").astype(int); d = d.sort_values("L")
ax[0].plot(d.L, d.stat, "o-", ms=3, lw=1.2, color="#4C78A8", label="observed (leave-one-out)")
ax[0].plot(d.L, d.null, "s--", ms=3, lw=1, color="0.5", label="label-shuffle null")
ax[0].axhline(0.744, ls=":", lw=1.1, color="#B0304A")
ax[0].set_xlim(-0.5, 14.6)
ax[0].set_ylim(0.45, 0.80)
ax[0].axvline(11.6, lw=0.8, color="0.85")
ax[0].text(14.4, 0.752, "round-1\nreported\nmean 0.744", fontsize=6, color="#B0304A",
           ha="right", va="bottom", linespacing=1.25)
ax[0].set_xticks([0, 2, 4, 6, 8, 10])
ax[0].set_xlabel("layer"); ax[0].set_ylabel("TF vs target AUROC")
# the band between the observed series (max 0.64) and the round-1 line (0.744) is empty
ax[0].legend(fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.82))
ax[0].set_title("a  TF/target separation, SV2-SV7", loc="left", fontsize=8)

for f, lbl, col in [("M1_confounds_SV2_SV7.csv", "SV2-SV7", "#4C78A8"),
                    ("M1_confounds_SV5_SV7.csv", "SV5-SV7", "#F58518")]:
    m = pd.read_csv(C.OUT/f)
    ax[1].plot(m.layer, m.r_rb, "o-", ms=3, lw=1.2, color=col, label=lbl)
ax[1].axhline(0, color="0.3", lw=1)
ax[1].set_xlim(-0.5, 11.5)
ax[1].set_xlabel("layer"); ax[1].set_ylabel(r"residualised $r_{rb}$")
# lower-left is empty: SV2-SV7 stays above +0.05, SV5-SV7 only dips below -0.04 after layer 6
ax[1].legend(fontsize=6.5, frameon=False, title="subspace", title_fontsize=6.5,
             loc="lower left", bbox_to_anchor=(0.02, 0.02))
ax[1].set_title("b  Regulatory signal after confound control", loc="left", fontsize=8)
fig.tight_layout(); fig.savefig(OUT/"Figure_2.pdf"); plt.close(fig)

# Fig 3: GRN benchmark forest plot
g = pd.read_csv(C.OUT/"GRN_benchmark_final.csv").sort_values("auroc")
fig, ax = plt.subplots(figsize=(5.0, 2.7))
y = np.arange(len(g))
ax.errorbar(g.auroc, y, xerr=[g.auroc-g.ci_lo, g.ci_hi-g.auroc], fmt="o", ms=4, lw=1.2, color="#4C78A8")
ax.axvline(0.5, ls="--", lw=1, color="0.5")
ax.text(0.5, len(g)-0.35, "chance", fontsize=6, color="0.45", ha="center", va="bottom")
ax.set_ylim(-0.7, len(g)-0.05)
ax.set_yticks(y); ax.set_yticklabels(g.method, fontsize=7)
ax.set_xlabel("AUROC on held-out TRRUST edges (49 pos / 147 neg)")
ax.set_title("Held-out regulatory-edge benchmark, 95% bootstrap CI", loc="left", fontsize=7.5)
fig.tight_layout(); fig.savefig(OUT/"Figure_3.pdf"); plt.close(fig)
print("figures written:", sorted(p.name for p in OUT.glob("*.pdf")))
