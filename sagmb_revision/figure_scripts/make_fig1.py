"""Figure 1: spectral concentration with depth, and the early-layer compartment axis."""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from figstyle import *

rd = pd.read_csv(R3 / "E1_compression" / "rank_definitions.csv")
ct = pd.read_csv(R3 / "E11_compartment" / "compartment_tests.csv")
ct = ct[(ct.universe == "V") & (ct.term == "GO:0005615")]
axis_dir = {"lung": "SV1", "lung_external": "SV1", "kidney_B": "SV1", "immune": "SV2", "kidney_A": "SV1"}

fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.15), gridspec_kw=dict(wspace=0.5))
L = np.arange(12)
for k, (col, ylab, letter) in enumerate([("pr", "participation ratio", "a"), ("er_sv", "effective rank", "b")]):
    ax = axs[k]
    for c in CTX:
        y = rd[rd.context == c].sort_values("layer")[col].values
        ax.plot(L, y, STYLE[c], color=COLOR[c], label=LABEL[c])
    ax.set_xlabel("layer"); ax.set_ylabel(ylab); ax.set_xticks([0, 3, 6, 9, 11]); panel(ax, letter)
    if col == "pr": ax.set_yscale("log"); ax.set_yticks([10, 20, 50, 100]); ax.set_yticklabels(["10", "20", "50", "100"])

ax = axs[2]
for c in CTX:
    d = ct[(ct.context == c) & (ct.direction == axis_dir[c])].sort_values("layer")
    ax.plot(d.layer, d.auroc_signfree, STYLE[c], color=COLOR[c], marker="o", markersize=2.5,
            label=LABEL[c])
null = ct.auroc_null_mean.mean()
ax.axhline(null, color=MUTED, lw=0.8, ls=(0, (2, 2)))
ax.text(0, null - 0.005, "null mean", ha="left", va="top", fontsize=6, color=MUTED)
ax.set_ylim(0.47, 0.82); ax.set_xticks([0, 3, 6, 9, 11]); ax.set_xlabel("layer")
ax.set_ylabel("extracellular AUROC"); panel(ax, "c")

# d: the lung axis replicates in the independent lung atlas (L3, SV1, matched genes)
import importlib.util
sys.path.insert(0, str(R3)); import stats3 as S
from pathlib import Path
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
g2go = pd.read_pickle(R3.parents[1] / "single_cell_mechinterp/data/perturb/gene2go_all.pkl")
coords = {}
for c in ("lung", "lung_external"):
    genes = json.loads((R3 / "cache" / c / "obs_genes.json").read_text()); idx = np.load(R3 / "cache" / c / "obs_idx.npy")
    E = np.load(EMB / f"{c}.npy", mmap_mode="r")
    US, _, _ = S.svd_observed(E, 3, idx)
    coords[c] = dict(zip(genes, US[:, 0]))
shared = sorted(set(coords["lung"]) & set(coords["lung_external"]))
x = np.array([coords["lung"][g] for g in shared]); y = np.array([coords["lung_external"][g] for g in shared])
ex = np.array(["GO:0005615" in g2go.get(g, ()) for g in shared]); cy = np.array(["GO:0005829" in g2go.get(g, ()) for g in shared]) & ~ex
# orient both so that extracellular genes have the larger mean coordinate
sx = 1 if x[ex].mean() > x[~ex].mean() else -1; sy = 1 if y[ex].mean() > y[~ex].mean() else -1
x, y = sx * x / np.abs(x).max(), sy * y / np.abs(y).max()
ax = axs[3]; other = ~ex & ~cy
ax.scatter(x[other], y[other], s=2, color="#c9c8c2", lw=0, rasterized=True, label="other")
ax.scatter(x[cy], y[cy], s=2.5, color="#e34948", lw=0, alpha=0.8, rasterized=True, label="cytosol")
ax.scatter(x[ex], y[ex], s=2.5, color="#6250d6", lw=0, alpha=0.8, rasterized=True, label="extracellular")
r = np.corrcoef(x, y)[0, 1]
ax.text(0.03, 0.97, f"r = {r:.3f}\nn = {len(shared):,}", transform=ax.transAxes, va="top", fontsize=6.5, color=INK2)
ax.set_xlabel("lung, $L_3$ SV1"); ax.set_ylabel("lung (ext.), $L_3$ SV1"); panel(ax, "d")
ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
ax.legend(loc="lower right", fontsize=6, markerscale=2.5, handletextpad=0.2, borderaxespad=0.1)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, ncol=5, loc="upper left", bbox_to_anchor=(0.06, 1.13), handlelength=2.2, columnspacing=1.4)
axs[2].text(0.98, 0.98, "immune: SV2\nothers: SV1", transform=axs[2].transAxes, ha="right", va="top", fontsize=6, color=INK2)
OUT.mkdir(exist_ok=True)
fig.savefig(OUT / "Figure_1.pdf"); fig.savefig(OUT / "Figure_1.png", dpi=200)
print("r", r, "n", len(shared), "ex", ex.sum(), "cy", cy.sum())
