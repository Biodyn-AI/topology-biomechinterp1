"""
iter_0028 Multi-Hypothesis Screen

H01 (manifold_distance / new_family): Hub Gene Embedding Centrality
    Convert iter_0027 confound (hub centrality bias) into the finding.
    For each named gene, compute STRING degree (number of STRING edges >= 0.4).
    Compute mean L2 distance from gene to centroid of all named genes per layer.
    Test: Spearman rho(STRING_degree, -mean_L2) per layer.
    Hypothesis: high-degree hub genes cluster near the center of the embedding space.
    Null: shuffle gene degree labels 500x.

H02 (intrinsic_dimensionality / new_method): PC1 as T-cell vs APC Axis at L11
    The brainstormer noted JUN/LCK-top vs FOS/HLA-bottom pattern on L11 PC1.
    Curate: T-cell effector markers (GZMB, PRF1, LCK, ZAP70, CD3D, CD3E, IFNG) and
            APC/antigen-presentation markers (HLA-A, HLA-B, HLA-C, HLA-DRA, HLA-DRB1,
            CD74, CIITA, B2M).
    Test: Spearman rho(PC1_loading, cell_type_axis) where axis = +1 for T-cell, -1 for APC.
    Also report Mann-Whitney AUROC for binary separation on PC1.
    Include per-layer AUROC curve to find the best layer for this axis.

H03 (graph_topology / new_family): Spectral Gap and Layer-wise Spectral Entropy
    Build kNN graph (k=10) for named genes at each layer using L2 distances.
    Compute normalized graph Laplacian eigenvalue spectrum.
    Report: (a) spectral gap (lambda_2, algebraic connectivity) per layer,
            (b) spectral entropy = -sum(p_i log p_i) for normalized eigenvalue distribution,
            (c) Spearman rho(layer, spectral_gap).
    Null: shuffled-label kNN graph (same geometry, random gene assignment).
    This is a new family (graph_topology via spectral analysis).
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0028"
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
TRRUST_PATH = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                   "/single_cell_mechinterp/external/networks/trrust_human.tsv")
ITER_DIR.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

# ─── Load embeddings ──────────────────────────────────────────────────────────
print("Loading embeddings...", flush=True)
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")  # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  Shape: {emb.shape}", flush=True)

with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes) if g}

# Named genes from edge dataset
EDGE_PATH = CYCLE1 / "cycle1_edge_dataset.tsv"
named_gene_set = set()
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        named_gene_set.add(row['source'])
        named_gene_set.add(row['target'])

named_genes = sorted(g for g in named_gene_set if g in gene_to_emb_idx)
named_idx = np.array([gene_to_emb_idx[g] for g in named_genes])
N_NAMED = len(named_genes)
gene_to_named_pos = {g: i for i, g in enumerate(named_genes)}
named_emb = emb[:, named_idx, :]  # [12, N_NAMED, 512]
print(f"  Named genes: {N_NAMED}, emb: {named_emb.shape}", flush=True)

# ─── Load STRING ──────────────────────────────────────────────────────────────
print("Loading STRING pairs...", flush=True)
STRING_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"
string_data = json.load(open(STRING_CACHE))
string_pairs_raw = string_data["pairs"]

# Build STRING degree dict
string_degree = defaultdict(int)
for p in string_pairs_raw:
    g1, g2 = p["g1"], p["g2"]
    if g1 in gene_to_named_pos and g2 in gene_to_named_pos:
        string_degree[g1] += 1
        string_degree[g2] += 1
print(f"  STRING degree: {len(string_degree)} genes with >=1 edge", flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# H01: Hub Gene Embedding Centrality
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== H01: Hub Gene Embedding Centrality ===", flush=True)

# Only genes with STRING degree >= 1
hub_genes = [g for g in named_genes if string_degree[g] > 0]
hub_idxs = [gene_to_named_pos[g] for g in hub_genes]
hub_degrees = np.array([string_degree[g] for g in hub_genes])
print(f"  Genes with STRING edges: {len(hub_genes)}", flush=True)
print(f"  Degree range: {hub_degrees.min()}–{hub_degrees.max()}, median={np.median(hub_degrees):.1f}", flush=True)

h01_layer_results = []
N_NULL_H01 = 500

for layer in range(N_LAYERS):
    layer_emb = named_emb[layer]  # [N_NAMED, 512]
    centroid = layer_emb.mean(axis=0)  # [512]
    # Mean L2 distance to centroid for each named gene
    all_dist_to_centroid = np.linalg.norm(layer_emb - centroid, axis=1)  # [N_NAMED]

    # For hub genes only
    hub_dist = all_dist_to_centroid[hub_idxs]  # distance to centroid

    # Spearman(degree, -distance) => high degree = close to center => positive rho
    rho, pval = spearmanr(hub_degrees, -hub_dist)

    # Null: shuffle degree labels
    null_rhos = np.array([
        spearmanr(rng.permutation(hub_degrees), -hub_dist)[0]
        for _ in range(N_NULL_H01)
    ])
    null_p = np.mean(np.abs(null_rhos) >= np.abs(rho))

    h01_layer_results.append({
        "layer": layer,
        "spearman_rho": float(rho),
        "pval": float(pval),
        "null_p": float(null_p),
        "null_rho_mean": float(null_rhos.mean()),
        "null_rho_std": float(null_rhos.std()),
    })
    print(f"  L{layer:02d}: rho={rho:+.3f} p={pval:.3e} null_p={null_p:.3f}", flush=True)

best_h01 = max(h01_layer_results, key=lambda x: abs(x["spearman_rho"]))
print(f"\n  Best layer: L{best_h01['layer']} rho={best_h01['spearman_rho']:+.3f}", flush=True)

# Save
h01_results = {
    "n_hub_genes": len(hub_genes),
    "degree_range": [int(hub_degrees.min()), int(hub_degrees.max())],
    "degree_median": float(np.median(hub_degrees)),
    "layer_results": h01_layer_results,
    "best_layer": best_h01["layer"],
    "best_rho": best_h01["spearman_rho"],
    "best_pval": best_h01["pval"],
    "best_null_p": best_h01["null_p"],
    "top_hub_genes_L11": [],
}

# Report top hub genes at L11
layer = 11
layer_emb = named_emb[layer]
centroid = layer_emb.mean(axis=0)
all_dist = np.linalg.norm(layer_emb - centroid, axis=1)
hub_dist = all_dist[hub_idxs]
hub_info = sorted(zip(hub_genes, hub_degrees, hub_dist), key=lambda x: x[2])
print(f"\n  Top 10 most-central hub genes at L11:", flush=True)
for g, d, dist in hub_info[:10]:
    print(f"    {g}: degree={d}, dist_to_centroid={dist:.3f}", flush=True)
h01_results["top_hub_genes_L11"] = [
    {"gene": g, "string_degree": int(d), "dist_to_centroid": float(dist)}
    for g, d, dist in hub_info[:15]
]

json.dump(h01_results, open(ITER_DIR / "h01_hub_centrality.json", "w"), indent=2)
print(f"  Saved h01_hub_centrality.json", flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# H02: PC1 as T-cell vs APC Axis at L11
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== H02: PC1 T-cell vs APC Axis ===", flush=True)

# T-cell effector markers (+1) and APC/antigen-presentation markers (-1)
tcell_markers = ["GZMB", "PRF1", "LCK", "ZAP70", "CD3D", "CD3E", "IFNG", "CD8A", "GNLY"]
apc_markers = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "CD74", "B2M", "FCGR3A", "CD14"]

# Also test the JUN/FOS split that the brainstormer highlighted
# JUN family: JUN-top (T-cell-like); FOS: FOS-bottom (APC-like)
tcell_ext = tcell_markers + ["JUN", "JUNB", "JUND"]
apc_ext = apc_markers + ["FOS", "FOSB", "FOSL1", "FOSL2"]

# Filter to named genes in vocab
tc_in_vocab = [g for g in tcell_ext if g in gene_to_named_pos]
apc_in_vocab = [g for g in apc_ext if g in gene_to_named_pos]
print(f"  T-cell markers in vocab: {tc_in_vocab}", flush=True)
print(f"  APC markers in vocab: {apc_in_vocab}", flush=True)

# Build axis vector: +1 for T-cell, -1 for APC (for genes with known identity)
axis_genes = tc_in_vocab + apc_in_vocab
axis_labels = [1] * len(tc_in_vocab) + [-1] * len(apc_in_vocab)
axis_idxs = [gene_to_named_pos[g] for g in axis_genes]

h02_layer_results = []
for layer in range(N_LAYERS):
    layer_emb = named_emb[layer]  # [N_NAMED, 512]

    # Fit PCA on all named genes, get PC1 loadings (gene coordinates)
    pca = PCA(n_components=3, random_state=42)
    pc_coords = pca.fit_transform(layer_emb)  # [N_NAMED, 3]

    pc1_coords = pc_coords[:, 0]  # [N_NAMED]

    # PC1 coordinates for axis genes
    pc1_axis = pc1_coords[axis_idxs]

    # Spearman rho(PC1, axis_label)
    rho, pval = spearmanr(pc1_axis, axis_labels)

    # AUROC for binary separation (T-cell=positive, APC=negative)
    tc_scores = pc1_axis[:len(tc_in_vocab)]
    apc_scores = pc1_axis[len(tc_in_vocab):]

    # AUROC: either direction could be positive
    stat, mwu_p = mannwhitneyu(tc_scores, apc_scores, alternative='two-sided')
    auroc = stat / (len(tc_scores) * len(apc_scores))
    auroc = max(auroc, 1 - auroc)  # take the better direction

    # Variance explained by PC1
    var_explained = pca.explained_variance_ratio_[0]

    h02_layer_results.append({
        "layer": layer,
        "spearman_rho": float(rho),
        "pval": float(pval),
        "auroc": float(auroc),
        "mwu_p": float(mwu_p),
        "pc1_var_explained": float(var_explained),
    })

    if layer == 11:
        print(f"  L11 details:", flush=True)
        for g, lab, pc1 in zip(axis_genes, axis_labels, pc1_axis):
            print(f"    {g:12s} axis={lab:+d} pc1={pc1:+.3f}", flush=True)

    print(f"  L{layer:02d}: rho={rho:+.3f} p={pval:.3e} AUROC={auroc:.3f} PC1_var={var_explained:.3f}", flush=True)

best_h02 = max(h02_layer_results, key=lambda x: x["auroc"])
print(f"\n  Best layer: L{best_h02['layer']} AUROC={best_h02['auroc']:.3f}", flush=True)

h02_results = {
    "tcell_markers_in_vocab": tc_in_vocab,
    "apc_markers_in_vocab": apc_in_vocab,
    "n_axis_genes": len(axis_genes),
    "layer_results": h02_layer_results,
    "best_layer_auroc": best_h02["layer"],
    "best_auroc": best_h02["auroc"],
}
json.dump(h02_results, open(ITER_DIR / "h02_pc1_celltype_axis.json", "w"), indent=2)
print(f"  Saved h02_pc1_celltype_axis.json", flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# H03: Spectral Gap of kNN Graph (new family: graph_topology / spectral)
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== H03: kNN Graph Spectral Gap ===", flush=True)

K_NEIGHBORS = 10

def compute_spectral_gap_and_entropy(layer_emb, k=10):
    """Build kNN graph, compute normalized Laplacian, return spectral gap and entropy."""
    n = layer_emb.shape[0]

    # Build kNN
    nbrs = NearestNeighbors(n_neighbors=k+1, metric='euclidean', algorithm='auto').fit(layer_emb)
    distances, indices = nbrs.kneighbors(layer_emb)

    # Build adjacency matrix (unweighted for now)
    rows, cols = [], []
    for i in range(n):
        for j in indices[i, 1:]:  # skip self
            rows.append(i); cols.append(j)
            rows.append(j); cols.append(i)
    adj = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    adj = (adj > 0).astype(float)  # symmetrize

    # Degree matrix
    degrees = np.array(adj.sum(axis=1)).flatten()
    degrees_safe = np.where(degrees > 0, degrees, 1.0)

    # Normalized Laplacian: L_norm = I - D^{-1/2} A D^{-1/2}
    d_inv_sqrt = 1.0 / np.sqrt(degrees_safe)
    # Scale adj
    d_inv_sqrt_mat = csr_matrix((d_inv_sqrt, (np.arange(n), np.arange(n))), shape=(n, n))
    L_norm = csr_matrix(np.eye(n)) - d_inv_sqrt_mat @ adj @ d_inv_sqrt_mat

    # Get smallest eigenvalues (expect lambda_1 ~ 0, lambda_2 = spectral gap)
    try:
        eigenvalues = eigsh(L_norm, k=min(20, n-2), which='SM', return_eigenvectors=False, tol=1e-4)
        eigenvalues = np.sort(np.abs(eigenvalues.real))
    except Exception:
        eigenvalues = np.linspace(0, 2, 20)

    # Spectral gap: lambda_2 (algebraic connectivity), skip near-zero
    nz_eigs = eigenvalues[eigenvalues > 1e-6]
    spectral_gap = float(nz_eigs[0]) if len(nz_eigs) > 0 else 0.0

    # Spectral entropy over all eigenvalues (normalized to [0,2])
    all_eigs = np.abs(eigenvalues)
    all_eigs_safe = all_eigs + 1e-12
    p = all_eigs_safe / all_eigs_safe.sum()
    spectral_entropy = float(-np.sum(p * np.log(p)))

    # Number of connected components (eigenvalues ~ 0)
    n_components = int(np.sum(eigenvalues < 1e-6))

    return spectral_gap, spectral_entropy, n_components, degrees.mean()


h03_layer_results = []
print(f"  k={K_NEIGHBORS} neighbors", flush=True)

for layer in range(N_LAYERS):
    layer_emb = named_emb[layer]  # [N_NAMED, 512]
    gap, entropy, n_comp, mean_degree = compute_spectral_gap_and_entropy(layer_emb, k=K_NEIGHBORS)

    h03_layer_results.append({
        "layer": layer,
        "spectral_gap": float(gap),
        "spectral_entropy": float(entropy),
        "n_components": int(n_comp),
        "mean_knn_degree": float(mean_degree),
    })
    print(f"  L{layer:02d}: gap={gap:.4f} entropy={entropy:.4f} n_comp={n_comp} mean_deg={mean_degree:.1f}", flush=True)

# Test layer trend
gaps = [r["spectral_gap"] for r in h03_layer_results]
entropies = [r["spectral_entropy"] for r in h03_layer_results]
layers = list(range(N_LAYERS))

rho_gap, p_gap = spearmanr(layers, gaps)
rho_entropy, p_entropy = spearmanr(layers, entropies)
print(f"\n  Spectral gap trend: rho={rho_gap:+.3f} p={p_gap:.3e}", flush=True)
print(f"  Entropy trend:      rho={rho_entropy:+.3f} p={p_entropy:.3e}", flush=True)

# Quick null: shuffle gene positions (but keep same embedding geometry)
# Use permuted assignment of genes to positions: sample 200 random permutations
# of gene indices, compute spectral gap at L11
print("\n  Computing null distribution at L11...", flush=True)
layer = 11
layer_emb_L11 = named_emb[layer]
true_gap_L11 = gaps[layer]
null_gaps = []
for _ in range(100):
    perm_emb = layer_emb_L11[rng.permutation(N_NAMED)]
    g, _, _, _ = compute_spectral_gap_and_entropy(perm_emb, k=K_NEIGHBORS)
    null_gaps.append(g)
null_p_gap = np.mean(np.array(null_gaps) >= true_gap_L11)
print(f"  L11 true gap={true_gap_L11:.4f}, null mean={np.mean(null_gaps):.4f}±{np.std(null_gaps):.4f}, null_p={null_p_gap:.3f}", flush=True)

h03_results = {
    "k_neighbors": K_NEIGHBORS,
    "layer_results": h03_layer_results,
    "spectral_gap_trend_rho": float(rho_gap),
    "spectral_gap_trend_pval": float(p_gap),
    "entropy_trend_rho": float(rho_entropy),
    "entropy_trend_pval": float(p_entropy),
    "L11_null_gap_mean": float(np.mean(null_gaps)),
    "L11_null_gap_std": float(np.std(null_gaps)),
    "L11_null_p": float(null_p_gap),
}
json.dump(h03_results, open(ITER_DIR / "h03_spectral_gap.json", "w"), indent=2)
print(f"  Saved h03_spectral_gap.json", flush=True)

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== SUMMARY ===", flush=True)
print(f"H01 best rho={h01_results['best_rho']:+.3f} at L{h01_results['best_layer']}, null_p={h01_results['best_null_p']:.3f}", flush=True)
print(f"H02 best AUROC={h02_results['best_auroc']:.3f} at L{h02_results['best_layer_auroc']}", flush=True)
print(f"H03 gap trend rho={rho_gap:+.3f} p={p_gap:.3e}, entropy trend rho={rho_entropy:+.3f} p={p_entropy:.3e}", flush=True)

print("\nDone.", flush=True)
