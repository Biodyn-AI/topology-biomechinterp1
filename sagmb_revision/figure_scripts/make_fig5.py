"""Figure 5: H4 — recovery as a function of the number of leading directions kept."""
import sys, json
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from figstyle import *
E7 = R3 / "E7_rank"
KS = [1, 2, 3, 5, 7, 10, 16, 32, 64, 128, 256, 512]
TASKS = [("PPI", "interaction AUROC"), ("NODE", "TF AUROC (stratified)"), ("EDGE", "within-TF edge AUROC")]
kstar = {t: v["k_star"] for t, v in json.load(open(E7 / "kstar.json"))["tasks"].items()}

fig, axs = plt.subplots(1, 4, figsize=(7.4, 2.35), gridspec_kw=dict(wspace=0.62, width_ratios=[1, 1, 1, 0.85]))
for ax, (task, ylab), letter in zip(axs[:3], TASKS, "abc"):
    for c in CTX:
        d = pd.read_csv(E7 / f"curves_{c}.csv")
        d = d[(d.family == "T") & (d.task == task)].copy(); d["k"] = d.param.astype(int)
        m = d.groupby("k").value.mean().reindex(KS)
        ax.plot(m.index, m.values, STYLE[c], color=COLOR[c], marker="o", markersize=2.2, label=LABEL[c])
    ax.axvline(kstar[task], color=MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.text(kstar[task], 1.0, f" $k^*$={kstar[task]}", transform=ax.get_xaxis_transform(), ha="right" if kstar[task] == 512 else "left",
            va="top", fontsize=6, color=INK2)
    ax.set_xscale("log", base=2); ax.set_xticks([1, 4, 16, 64, 512]); ax.set_xticklabels(["1", "4", "16", "64", "512"])
    ax.minorticks_off(); ax.set_xlabel("leading directions kept ($k$)"); ax.set_ylabel(ylab); panel(ax, letter)
axs[1].set_ylim(0.1, 0.78)

# d: static gene embeddings of other models — best truncation minus full embedding, per task
s = pd.read_csv(E7 / "static_summary.csv"); s = s[s.model != "scgpt_lung_contextual"]
ax = axs[3]
for i, (task, _) in enumerate(TASKS):
    v = s[s.task == task].best_minus_full.values
    x = i + np.linspace(-0.18, 0.18, len(v))
    ax.scatter(x, v, s=10, color="#6250d6", edgecolor="white", linewidth=0.4, zorder=3)
ax.axhline(0, color=MUTED, lw=0.8, ls=(0, (2, 2)))
ax.set_xticks(range(3)); ax.set_xticklabels(["inter-\naction", "TF", "edge"]); ax.set_xlim(-0.5, 2.5)
ax.set_ylim(-0.003, 0.03); ax.set_ylabel("best truncation $-$ full"); ax.grid(axis="x", visible=False)
ax.set_title(f"static embeddings ({s.model.nunique()})", fontsize=6.5, color=INK2); panel(ax, "d")
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, ncol=5, loc="upper left", bbox_to_anchor=(0.06, 1.13), handlelength=2.2, columnspacing=1.4)
fig.savefig(OUT / "Figure_5.pdf"); fig.savefig(OUT / "Figure_5.png", dpi=200)
print("ok", kstar, s.groupby("task").best_minus_full.max().round(4).to_dict())
