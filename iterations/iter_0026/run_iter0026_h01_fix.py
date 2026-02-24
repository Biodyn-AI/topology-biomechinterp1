"""
H01 Fix: Use immune-gene family groups from the actual 209-gene vocab.
These groups serve as CORUM-proxy: within-family pairs should be closer
than cross-family pairs if the model encodes biological similarity.

Families with >=3 vocab genes:
  AP1:       FOS, FOSB, JUN, JUNB
  RUNX:      RUNX1, RUNX2, RUNX3
  KLF:       KLF2, KLF6, KLF14
  BCL2fam:   BCL2, BCL2A1, BIRC3, BIRC5
  HLA_I:     HLA-A, HLA-B, HLA-C
  HLA_II:    HLA-DMB, HLA-DOB, HLA-DPB1, HLA-DRA, HLA-DRB1
  CCL:       CCL3, CCL5, CCL20
  TNFSF:     TNF, TNFSF10, TNFSF14
  IL2_path:  IL2, IL2RA, IL2RB
  HIF:       HIF1A, HIF3A   [only 2 — skip]
  IRF:       IRF4, IRF8     [only 2 — skip]
"""

import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0026"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")

emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")  # [12, 4803, 512]
with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [l.strip() for l in f]
gene2idx = {g: i for i, g in enumerate(vocab_genes)}

rng = np.random.default_rng(42)

def auroc(pos_scores, neg_scores):
    if len(pos_scores) < 2 or len(neg_scores) < 2:
        return float("nan")
    stat, _ = mannwhitneyu(pos_scores, neg_scores, alternative="greater")
    return stat / (len(pos_scores) * len(neg_scores))

IMMUNE_FAMILIES = {
    "AP1":      ["FOS", "FOSB", "JUN", "JUNB"],
    "RUNX":     ["RUNX1", "RUNX2", "RUNX3"],
    "KLF":      ["KLF2", "KLF6", "KLF14"],
    "BCL2fam":  ["BCL2", "BCL2A1", "BIRC3", "BIRC5"],
    "HLA_I":    ["HLA-A", "HLA-B", "HLA-C"],
    "HLA_II":   ["HLA-DMB", "HLA-DOB", "HLA-DPB1", "HLA-DRA", "HLA-DRB1"],
    "CCL":      ["CCL3", "CCL5", "CCL20"],
    "TNFSF":    ["TNF", "TNFSF10", "TNFSF14"],
    "IL2_path": ["IL2", "IL2RA", "IL2RB"],
}

# Filter to genes in vocab
family_genes = {}
for fname, genes in IMMUNE_FAMILIES.items():
    valid = [g for g in genes if g in gene2idx]
    if len(valid) >= 3:
        family_genes[fname] = valid

print(f"Families with >=3 vocab genes: {len(family_genes)}")
for k, v in family_genes.items():
    print(f"  {k}: {v}")

# Build within-family pairs
within_pairs = set()
for fname, genes in family_genes.items():
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            pair = tuple(sorted([gene2idx[genes[i]], gene2idx[genes[j]]]))
            within_pairs.add(pair)

# All genes in any family
all_family_genes = []
for genes in family_genes.values():
    all_family_genes.extend(genes)
all_family_genes = list(set(all_family_genes))
all_family_idxs = [gene2idx[g] for g in all_family_genes]

# Gene -> families mapping
gene_to_fam = defaultdict(set)
for fname, genes in family_genes.items():
    for g in genes:
        gene_to_fam[g].add(fname)

# Cross-family pairs
cross_pairs = set()
for i in range(len(all_family_genes)):
    for j in range(i + 1, len(all_family_genes)):
        gi, gj = all_family_genes[i], all_family_genes[j]
        if not (gene_to_fam[gi] & gene_to_fam[gj]):
            pair = tuple(sorted([gene2idx[gi], gene2idx[gj]]))
            cross_pairs.add(pair)
cross_pairs -= within_pairs

print(f"Within-family pairs: {len(within_pairs)}")
print(f"Cross-family pairs: {len(cross_pairs)}")

# Balance
if len(cross_pairs) > len(within_pairs):
    cross_arr = list(cross_pairs)
    idx_sel = rng.choice(len(cross_arr), size=len(within_pairs), replace=False)
    cross_sampled = [cross_arr[i] for i in idx_sel]
else:
    cross_sampled = list(cross_pairs)

within_list = list(within_pairs)
print(f"Balanced: {len(within_list)} within, {len(cross_sampled)} cross")

results = {
    "hypothesis": "H01_immune_family_AUROC",
    "n_families": len(family_genes),
    "family_names": list(family_genes.keys()),
    "n_within": len(within_list),
    "n_cross": len(cross_sampled),
    "per_layer": []
}

for layer in range(emb.shape[0]):
    le = emb[layer]

    def dists(pairs):
        idx_i = np.array([p[0] for p in pairs])
        idx_j = np.array([p[1] for p in pairs])
        diff = le[idx_i] - le[idx_j]
        return np.sqrt((diff ** 2).sum(axis=1))

    d_within = dists(within_list)
    d_cross = dists(cross_sampled)
    auc = auroc(-d_within, -d_cross)

    results["per_layer"].append({
        "layer": layer,
        "AUROC": float(auc),
        "within_mean": float(d_within.mean()),
        "cross_mean": float(d_cross.mean()),
        "effect_size": float((d_cross.mean() - d_within.mean()) / d_cross.std()),
    })
    print(f"  Layer {layer:2d}: AUROC={auc:.4f} within={d_within.mean():.3f} cross={d_cross.mean():.3f}")

# Per-family AUROCs at best overall layer
best_layer = max(results["per_layer"], key=lambda x: x["AUROC"])["layer"]
print(f"\nBest layer: {best_layer}")

# Per-family breakdown at best layer
le = emb[best_layer]
per_family = {}
for fname, fgenes in family_genes.items():
    other_genes = [g for g in all_family_genes if g not in fgenes]
    if len(other_genes) < 2:
        continue
    w_pairs = [(gene2idx[fgenes[i]], gene2idx[fgenes[j]])
               for i in range(len(fgenes))
               for j in range(i + 1, len(fgenes))]
    c_pairs = [(gene2idx[fgenes[i]], gene2idx[og])
               for i in range(len(fgenes))
               for og in other_genes]

    def d(p):
        if not p: return np.array([])
        ii = np.array([x[0] for x in p])
        jj = np.array([x[1] for x in p])
        diff = le[ii] - le[jj]
        return np.sqrt((diff**2).sum(axis=1))

    dw = d(w_pairs)
    dc = d(c_pairs)
    auc = auroc(-dw, -dc)
    per_family[fname] = {"AUROC": float(auc),
                         "within_mean": float(dw.mean()),
                         "cross_mean": float(dc.mean()),
                         "n_within": len(w_pairs),
                         "n_cross": len(c_pairs)}
    print(f"  {fname}: AUROC={auc:.4f} w={dw.mean():.3f} c={dc.mean():.3f}")

results["best_layer"] = best_layer
results["per_family_best_layer"] = per_family
results["peak_AUROC"] = results["per_layer"][best_layer]["AUROC"]

# Spearman AUROC vs layer
layers = [r["layer"] for r in results["per_layer"]]
aurocs = [r["AUROC"] for r in results["per_layer"]]
from scipy.stats import spearmanr
rho, p = spearmanr(layers, aurocs)
results["spearman_auroc_vs_layer"] = {"rho": float(rho), "p": float(p)}
print(f"\nAUROC vs layer: rho={rho:.4f} p={p:.4g}")

out = ITER_DIR / "h01_immune_family_auroc.json"
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print(f"Saved: {out}")
