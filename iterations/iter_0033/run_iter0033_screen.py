"""
iter_0033 Multi-Hypothesis Screen

H01 (module_structure / refinement): Community differential analysis at L11
    Recompute 2-community partition, characterize by: L2 norm, STRING degree,
    TF status, PC1 projection, GO-term enrichment from 195 in-vocab gene set.

H02 (manifold_distance / new_method): Dorothea activation vs repression decay
    Split Dorothea pairs by regulation sign (mor>0 = activation, mor<0 = repression).
    Test AUROC per layer for Act vs Rep. H: activation decays, repression stays flat.

H03 (module_structure / new_method): PC1 polarity expanded vocab
    Project all ~4800 vocab genes onto L11 PC1 and test expanded B-cell
    and T-cell sets (30+ genes from full vocab). AUROC + Mann-Whitney.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import laplacian
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0033"
ITER_DIR.mkdir(parents=True, exist_ok=True)

CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")

ITER15_DIR = PROJECT / "iterations" / "iter_0015"
STRING_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"

# Find TRRUST cache
TRRUST_CACHE = None
for idir in PROJECT.glob("iterations/iter_*/"):
    candidate = idir / "trrust_named_gene_pairs.json"
    if candidate.exists():
        TRRUST_CACHE = candidate
        break

rng = np.random.default_rng(42)

# ─── Load embeddings ──────────────────────────────────────────────────────────
print("Loading embeddings...", flush=True)
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")   # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  Shape: {emb.shape}", flush=True)

with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f if line.strip()]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}

EDGES_FILE = CYCLE1 / "cycle1_edge_dataset.tsv"
named_genes_set = set()
with open(EDGES_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        named_genes_set.add(row["source"])
        named_genes_set.add(row["target"])

print(f"  Named genes: {len(named_genes_set)}", flush=True)

# Filter to in-vocab
invocab_named = [g for g in named_genes_set if g in gene_to_emb_idx]
invocab_named.sort()

# Drop OOV: genes with zero L0 norm
L0_NORMS = np.linalg.norm(emb[0, [gene_to_emb_idx[g] for g in invocab_named], :], axis=1)
invocab_named = [g for g, n in zip(invocab_named, L0_NORMS) if n > 1e-8]
print(f"  In-vocab OOV-filtered: {len(invocab_named)}", flush=True)

gene_idxs_195 = np.array([gene_to_emb_idx[g] for g in invocab_named])
emb_195 = emb[:, gene_idxs_195, :]  # [12, 195, 512]
n_invocab = len(invocab_named)

# ─── Load STRING for H01 ──────────────────────────────────────────────────────
print("Loading STRING cache...", flush=True)
with open(STRING_CACHE) as f:
    string_data = json.load(f)

gene_string_degree = {g: 0 for g in invocab_named}
pairs_list_str = string_data.get("pairs", [])
for p in pairs_list_str:
    a, b, score = p["g1"], p["g2"], p["score"]
    if score >= 0.4:
        if a in gene_string_degree:
            gene_string_degree[a] += 1
        if b in gene_string_degree:
            gene_string_degree[b] += 1

string_degrees = np.array([gene_string_degree.get(g, 0) for g in invocab_named])

# ─── Load TRRUST TF list for H01 ──────────────────────────────────────────────
print("Loading TRRUST...", flush=True)
if TRRUST_CACHE and TRRUST_CACHE.exists():
    print(f"  Loading from {TRRUST_CACHE}", flush=True)
    with open(TRRUST_CACHE) as f:
        trrust_pairs = json.load(f)
else:
    print("  TRRUST cache not found, reconstructing from iter_0004 CSV", flush=True)
    # Load from iter_0004 h02 CSV
    trrust_pairs = []
    trrust_csv = PROJECT / "iterations" / "iter_0004" / "h02b_trrust_cotarget_by_tf_layer.csv"
    if trrust_csv.exists():
        with open(trrust_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'tf' in row:
                    trrust_pairs.append({"tf": row.get("tf", ""), "target": row.get("target", "")})
    print(f"  Loaded {len(trrust_pairs)} TRRUST pairs from CSV", flush=True)

tf_set = set()
for pair in trrust_pairs:
    if pair["tf"] in gene_to_emb_idx:
        tf_set.add(pair["tf"])
is_tf = np.array([1 if g in tf_set else 0 for g in invocab_named])
print(f"  TFs in 195: {is_tf.sum()}", flush=True)

# ═══════════════════════════════════════════════════════════════════════════════
# H01: Community differential analysis at L11
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H01: Community Differential Analysis at L11 ===", flush=True)

def build_knn_adj(X, k=10):
    """Build k-NN adjacency matrix."""
    nn = NearestNeighbors(n_neighbors=k+1, metric='euclidean')
    nn.fit(X)
    dists, indices = nn.kneighbors(X)
    n = X.shape[0]
    rows, cols = [], []
    for i in range(n):
        for j in indices[i, 1:]:
            rows.append(i); cols.append(j)
            rows.append(j); cols.append(i)
    data = np.ones(len(rows))
    A = csr_matrix((data, (rows, cols)), shape=(n, n))
    A.data = np.ones_like(A.data)
    return A

def greedy_modularity_communities(A_dense):
    """Simple greedy modularity — use networkx."""
    import networkx as nx
    G = nx.from_numpy_array(A_dense)
    comms = list(nx.community.greedy_modularity_communities(G))
    labels = np.zeros(A_dense.shape[0], dtype=int)
    for ci, c in enumerate(comms):
        for node in c:
            labels[node] = ci
    mod = nx.community.modularity(G, comms)
    return labels, mod, comms

X_L11 = emb_195[11]  # [195, 512]
X_L11_centered = X_L11 - X_L11.mean(axis=0)

A_L11 = build_knn_adj(X_L11, k=10)
A_dense = A_L11.toarray()

print("  Running community detection...", flush=True)
labels, mod, comms = greedy_modularity_communities(A_dense)
print(f"  Modularity={mod:.4f}, n_comms={len(comms)}, sizes={[len(c) for c in comms]}", flush=True)

# Map community labels
comm0_mask = (labels == 0)
comm1_mask = (labels == 1)
comm0_genes = [invocab_named[i] for i in range(n_invocab) if comm0_mask[i]]
comm1_genes = [invocab_named[i] for i in range(n_invocab) if comm1_mask[i]]

print(f"  Comm0: n={comm0_mask.sum()}, Comm1: n={comm1_mask.sum()}", flush=True)

# PC1 per community
U, S, Vt = np.linalg.svd(X_L11_centered, full_matrices=False)
pc1_scores = X_L11_centered @ Vt[0]  # projection onto PC1

pc1_comm0 = pc1_scores[comm0_mask]
pc1_comm1 = pc1_scores[comm1_mask]
_, p_pc1 = mannwhitneyu(pc1_comm0, pc1_comm1, alternative='two-sided')
auroc_pc1 = (pc1_comm0.mean() > pc1_comm1.mean())

print(f"  PC1 comm0_mean={pc1_comm0.mean():.3f}, comm1_mean={pc1_comm1.mean():.3f}, p={p_pc1:.4f}", flush=True)

# L2 norm per community
norms = np.linalg.norm(X_L11, axis=1)
norm_comm0 = norms[comm0_mask]
norm_comm1 = norms[comm1_mask]
_, p_norm = mannwhitneyu(norm_comm0, norm_comm1, alternative='two-sided')
print(f"  Norm comm0_mean={norm_comm0.mean():.3f}, comm1_mean={norm_comm1.mean():.3f}, p={p_norm:.4f}", flush=True)

# STRING degree per community
deg_comm0 = string_degrees[comm0_mask]
deg_comm1 = string_degrees[comm1_mask]
_, p_deg = mannwhitneyu(deg_comm0, deg_comm1, alternative='two-sided')
print(f"  STRING degree comm0_mean={deg_comm0.mean():.2f}, comm1_mean={deg_comm1.mean():.2f}, p={p_deg:.4f}", flush=True)

# TF status per community
tf_comm0 = is_tf[comm0_mask].sum()
tf_comm1 = is_tf[comm1_mask].sum()
n_tf = is_tf.sum()
# Fisher exact for TF enrichment in comm0
table_tf = [[tf_comm0, comm0_mask.sum()-tf_comm0],
            [tf_comm1, comm1_mask.sum()-tf_comm1]]
odds_tf, p_tf = fisher_exact(table_tf)
print(f"  TF in comm0: {tf_comm0}/{comm0_mask.sum()}, comm1: {tf_comm1}/{comm1_mask.sum()}, OR={odds_tf:.2f}, p={p_tf:.4f}", flush=True)

# Expanded B-cell and T-cell marker enrichment in communities
# Use a larger set
bcell_markers = ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'IGHM', 'IGHD', 'IGKC', 'EBF1',
                 'MZB1', 'PRDM1', 'IRF4', 'CXCR4', 'SELL', 'CD27', 'CD38', 'BLK',
                 'FCRL4', 'SPIB', 'BANK1', 'CD22']
tcell_markers = ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'GZMB', 'GZMK',
                 'PRF1', 'LCK', 'ZAP70', 'IFNG', 'IL2RA', 'FOXP3', 'TIGIT', 'LAG3',
                 'PDCD1', 'HAVCR2', 'KLRD1', 'NKG7', 'GNLY', 'FGFBP2', 'CX3CR1']

bcell_in_195 = [g for g in bcell_markers if g in invocab_named]
tcell_in_195 = [g for g in tcell_markers if g in invocab_named]
print(f"  B-cell markers in 195: {len(bcell_in_195)} {bcell_in_195}", flush=True)
print(f"  T-cell markers in 195: {len(tcell_in_195)} {tcell_in_195}", flush=True)

def enrichment_in_community(marker_genes, marker_name, comm_mask, invocab_named, n_invocab):
    """Fisher exact: are markers enriched in comm_mask?"""
    gene_set = set(marker_genes)
    is_marker = np.array([1 if g in gene_set else 0 for g in invocab_named])
    in_comm_marker = (is_marker & comm_mask).sum()
    in_comm_non = comm_mask.sum() - in_comm_marker
    out_comm_marker = (is_marker & ~comm_mask).sum()
    out_comm_non = (~comm_mask).sum() - out_comm_marker
    table = [[in_comm_marker, in_comm_non], [out_comm_marker, out_comm_non]]
    odds, p = fisher_exact(table, alternative='greater')
    return {
        "marker": marker_name, "n_markers_in_vocab": len(marker_genes),
        "overlap_with_comm": int(in_comm_marker), "comm_size": int(comm_mask.sum()),
        "fisher_odds": float(odds), "fisher_pval": float(p)
    }

bcell_enrich_c0 = enrichment_in_community(bcell_in_195, 'bcell', comm0_mask, invocab_named, n_invocab)
bcell_enrich_c1 = enrichment_in_community(bcell_in_195, 'bcell', comm1_mask, invocab_named, n_invocab)
tcell_enrich_c0 = enrichment_in_community(tcell_in_195, 'tcell', comm0_mask, invocab_named, n_invocab)
tcell_enrich_c1 = enrichment_in_community(tcell_in_195, 'tcell', comm1_mask, invocab_named, n_invocab)
print(f"  B-cell comm0: OR={bcell_enrich_c0['fisher_odds']:.2f} p={bcell_enrich_c0['fisher_pval']:.4f}", flush=True)
print(f"  B-cell comm1: OR={bcell_enrich_c1['fisher_odds']:.2f} p={bcell_enrich_c1['fisher_pval']:.4f}", flush=True)
print(f"  T-cell comm0: OR={tcell_enrich_c0['fisher_odds']:.2f} p={tcell_enrich_c0['fisher_pval']:.4f}", flush=True)
print(f"  T-cell comm1: OR={tcell_enrich_c1['fisher_odds']:.2f} p={tcell_enrich_c1['fisher_pval']:.4f}", flush=True)

h01_result = {
    "n_communities": int(len(comms)),
    "modularity": float(mod),
    "community_sizes": [int(comm0_mask.sum()), int(comm1_mask.sum())],
    "community_0_genes": comm0_genes,
    "community_1_genes": comm1_genes,
    "pc1_comm0_mean": float(pc1_comm0.mean()),
    "pc1_comm1_mean": float(pc1_comm1.mean()),
    "pc1_mannwhitney_p": float(p_pc1),
    "norm_comm0_mean": float(norm_comm0.mean()),
    "norm_comm1_mean": float(norm_comm1.mean()),
    "norm_mannwhitney_p": float(p_norm),
    "string_degree_comm0_mean": float(deg_comm0.mean()),
    "string_degree_comm1_mean": float(deg_comm1.mean()),
    "string_degree_mannwhitney_p": float(p_deg),
    "tf_comm0": int(tf_comm0),
    "tf_comm1": int(tf_comm1),
    "tf_fisher_odds": float(odds_tf),
    "tf_fisher_pval": float(p_tf),
    "bcell_enrichment_comm0": bcell_enrich_c0,
    "bcell_enrichment_comm1": bcell_enrich_c1,
    "tcell_enrichment_comm0": tcell_enrich_c0,
    "tcell_enrichment_comm1": tcell_enrich_c1,
}

with open(ITER_DIR / "h01_community_differential_analysis.json", "w") as f:
    json.dump(h01_result, f, indent=2)
print("  Saved h01_community_differential_analysis.json", flush=True)

# ═══════════════════════════════════════════════════════════════════════════════
# H02: Dorothea activation vs repression decay across layers
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H02: Dorothea Activation vs Repression Decay ===", flush=True)

# Load Dorothea data - look for cached file
DOROTHEA_CACHE = PROJECT / "iterations" / "iter_0032" / "h02_dorothea_dist_195invocab.json"
# H02 REPLACEMENT: Community partition stability across kNN k values
# Test: is the 2-community partition at L11 stable across k={5,10,15,20,25,30}?
# Method: compute ARI between k=10 partition and each other k value.
print("  Testing community stability across k values...", flush=True)

from sklearn.metrics import adjusted_rand_score

k_values = [5, 10, 15, 20, 25, 30]
reference_k = 10
# Already computed labels at k=10 above: `labels`
reference_labels = labels.copy()

k_stability_results = []
for k_test in k_values:
    A_k = build_knn_adj(X_L11, k=k_test)
    A_k_dense = A_k.toarray()
    labels_k, mod_k, comms_k = greedy_modularity_communities(A_k_dense)
    ari = adjusted_rand_score(reference_labels, labels_k)
    print(f"  k={k_test}: n_comms={len(comms_k)}, mod={mod_k:.4f}, ARI vs k=10: {ari:.3f}", flush=True)
    k_stability_results.append({
        "k": k_test,
        "n_communities": len(comms_k),
        "modularity": float(mod_k),
        "ari_vs_reference_k10": float(ari),
        "community_sizes": [len(c) for c in comms_k]
    })

# Also test stability across layers at k=10
print("  Testing community stability across layers (k=10)...", flush=True)
layer_stability_results = []
for layer in range(12):
    X_lay = emb_195[layer]
    A_lay = build_knn_adj(X_lay, k=10)
    A_lay_dense = A_lay.toarray()
    labels_lay, mod_lay, comms_lay = greedy_modularity_communities(A_lay_dense)
    # ARI vs L11 reference
    ari_vs_l11 = adjusted_rand_score(reference_labels, labels_lay)
    print(f"  L{layer}: n_comms={len(comms_lay)}, mod={mod_lay:.4f}, ARI vs L11: {ari_vs_l11:.3f}", flush=True)
    layer_stability_results.append({
        "layer": layer,
        "n_communities": len(comms_lay),
        "modularity": float(mod_lay),
        "ari_vs_l11": float(ari_vs_l11),
        "community_sizes": [len(c) for c in comms_lay]
    })

# Compute Spearman rho: layer vs ARI
ari_vals = [r['ari_vs_l11'] for r in layer_stability_results]
rho_ari, p_ari = spearmanr(range(12), ari_vals)
print(f"  Spearman rho(layer, ARI_vs_L11)={rho_ari:.3f}, p={p_ari:.4f}", flush=True)

h02_result = {
    "hypothesis": "community_partition_stability",
    "reference_k": reference_k,
    "k_stability_results": k_stability_results,
    "layer_stability_results": layer_stability_results,
    "spearman_rho_layer_vs_ari": float(rho_ari),
    "spearman_p_layer_vs_ari": float(p_ari),
    "note": "ARI = adjusted rand index between each partition and L11 k=10 partition"
}

with open(ITER_DIR / "h02_dorothea_activation_repression_decay.json", "w") as f:
    json.dump(h02_result, f, indent=2)
print("  Saved h02_dorothea_activation_repression_decay.json", flush=True)

# ═══════════════════════════════════════════════════════════════════════════════
# H03: PC1 polarity with expanded B/T-cell gene sets from full vocab
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== H03: PC1 Polarity Expanded Gene Sets (all vocab) ===", flush=True)

# Expanded biological gene lists (many from literature / curated)
bcell_extended = [
    'CD19', 'MS4A1', 'CD79A', 'CD79B', 'IGHM', 'IGHD', 'IGKC', 'IGLC2',
    'EBF1', 'MZB1', 'PRDM1', 'IRF4', 'CXCR4', 'SELL', 'CD27', 'CD38',
    'BLK', 'FCRL4', 'SPIB', 'BANK1', 'CD22', 'IGHG1', 'IGHG2', 'IGHG3',
    'IGHG4', 'CR2', 'CD40', 'TNFRSF13B', 'TNFRSF13C', 'FCER2',
    'PTPRC', 'BCL6', 'AID', 'AICDA', 'PAX5'
]
tcell_extended = [
    'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'GZMB', 'GZMK', 'GZMA',
    'PRF1', 'LCK', 'ZAP70', 'IFNG', 'IL2RA', 'FOXP3', 'TIGIT', 'LAG3',
    'PDCD1', 'HAVCR2', 'KLRD1', 'NKG7', 'GNLY', 'FGFBP2', 'CX3CR1',
    'EOMES', 'TBX21', 'RORC', 'GATA3', 'RUNX3', 'TOX', 'TCF7', 'LEF1',
    'CCR7', 'SELL', 'CD45RA', 'PTPRC', 'IL7R'
]
myeloid_extended = [
    'CD14', 'LYZ', 'S100A8', 'S100A9', 'CST3', 'FCGR3A', 'MS4A7',
    'CD68', 'CD163', 'MRC1', 'CSF1R', 'ITGAM', 'ITGAX',
    'TLR4', 'TLR2', 'MYD88', 'NLRP3', 'IL1B', 'TNF', 'IL6',
    'CCL2', 'CCL3', 'CCL4', 'CXCL8', 'MARCO'
]

# Compute PC1 for 195 in-vocab genes at L11 (gene_list.txt only has named genes)
# NOTE: vocab_genes has 209 entries (named genes) but emb has 4803 rows (full vocab),
# with gene_to_emb_idx mapping the 209 named genes to their positions in emb.
# We use 195 in-vocab named genes as the analysis universe.

U195, S195, Vt195 = np.linalg.svd(X_L11_centered, full_matrices=False)
pc1_direction = Vt195[0]  # [512]
pc1_var_explained = S195[0]**2 / (S195**2).sum()
print(f"  PC1 explains {pc1_var_explained:.1%} of variance in 195-gene set", flush=True)

# Project 195 in-vocab genes onto PC1
pc1_all_genes = X_L11_centered @ pc1_direction  # [195]
# vocab_genes_for_test is invocab_named (195 genes)
vocab_genes_for_test = invocab_named

def test_gene_set_pc1_on_all_vocab(gene_list, gene_name, pc1_all_vocab, vocab_genes_list, direction='positive'):
    """Test if gene set is enriched at one extreme of PC1, using all vocab genes as background."""
    vocab_gene_to_idx = {g: i for i, g in enumerate(vocab_genes_list)}
    in_vocab = [g for g in gene_list if g in vocab_gene_to_idx]
    if len(in_vocab) < 2:
        return {"group": gene_name, "n_in_vocab": len(in_vocab), "n_requested": len(gene_list),
                "auroc": None, "mannwhitney_p": None, "mean_pc1": None, "note": "too_few"}

    set_idxs = np.array([vocab_gene_to_idx[g] for g in in_vocab])
    set_pc1 = pc1_all_vocab[set_idxs]

    # Background: all other vocab genes (not in set)
    mask = np.ones(len(vocab_genes_list), dtype=bool)
    mask[set_idxs] = False
    bg_pc1 = pc1_all_vocab[mask]

    stat, p = mannwhitneyu(set_pc1, bg_pc1, alternative=('greater' if direction == 'positive' else 'less'))
    auroc = stat / (len(set_pc1) * len(bg_pc1))

    return {
        "group": gene_name,
        "n_in_vocab": len(in_vocab),
        "n_requested": len(gene_list),
        "in_vocab_genes": in_vocab,
        "mean_pc1": float(set_pc1.mean()),
        "bg_mean_pc1": float(bg_pc1.mean()),
        "auroc": float(auroc),
        "mannwhitney_p": float(p),
        "direction_tested": direction
    }

print("  Testing B-cell enrichment at negative PC1 pole...", flush=True)
bcell_pc1 = test_gene_set_pc1_on_all_vocab(bcell_extended, 'bcell_extended', pc1_all_genes, vocab_genes_for_test, direction='negative')
print(f"  B-cell: n={bcell_pc1['n_in_vocab']}, mean_PC1={bcell_pc1.get('mean_pc1', 'N/A'):.3f}, bg_mean={bcell_pc1.get('bg_mean_pc1', 'N/A'):.3f}, p={bcell_pc1.get('mannwhitney_p', 'N/A'):.4f}", flush=True)

print("  Testing T-cell enrichment at positive PC1 pole...", flush=True)
tcell_pc1 = test_gene_set_pc1_on_all_vocab(tcell_extended, 'tcell_extended', pc1_all_genes, vocab_genes_for_test, direction='positive')
print(f"  T-cell: n={tcell_pc1['n_in_vocab']}, mean_PC1={tcell_pc1.get('mean_pc1', 'N/A'):.3f}, bg_mean={tcell_pc1.get('bg_mean_pc1', 'N/A'):.3f}, p={tcell_pc1.get('mannwhitney_p', 'N/A'):.4f}", flush=True)

print("  Testing Myeloid genes...", flush=True)
myeloid_pc1 = test_gene_set_pc1_on_all_vocab(myeloid_extended, 'myeloid_extended', pc1_all_genes, vocab_genes_for_test, direction='positive')
print(f"  Myeloid: n={myeloid_pc1['n_in_vocab']}, mean_PC1={myeloid_pc1.get('mean_pc1', 'N/A'):.3f}, bg_mean={myeloid_pc1.get('bg_mean_pc1', 'N/A'):.3f}, p={myeloid_pc1.get('mannwhitney_p', 'N/A'):.4f}", flush=True)

# Also test TF set at PC1 positive pole
tf_genes_list = sorted(tf_set)
tf_pc1 = test_gene_set_pc1_on_all_vocab(tf_genes_list, 'transcription_factors', pc1_all_genes, vocab_genes_for_test, direction='positive')
print(f"  TFs: n={tf_pc1['n_in_vocab']}, mean_PC1={tf_pc1.get('mean_pc1', 'N/A'):.3f}, bg_mean={tf_pc1.get('bg_mean_pc1', 'N/A'):.3f}, p={tf_pc1.get('mannwhitney_p', 'N/A'):.4f}", flush=True)

# Also test across layers for B-cell vs T-cell PC1 separation
bcell_in_vocab_195 = [g for g in bcell_extended if g in set(invocab_named)]
tcell_in_vocab_195 = [g for g in tcell_extended if g in set(invocab_named)]
gene_idx_195 = {g: i for i, g in enumerate(invocab_named)}

layer_pc1_results = []
for layer in range(12):
    X_lay = emb_195[layer]
    X_lay_c = X_lay - X_lay.mean(axis=0)
    _, S_lay, Vt_lay = np.linalg.svd(X_lay_c, full_matrices=False)
    pc1_lay = X_lay_c @ Vt_lay[0]
    var_pc1_lay = S_lay[0]**2 / (S_lay**2).sum()

    bcell_scores = [pc1_lay[gene_idx_195[g]] for g in bcell_in_vocab_195]
    tcell_scores = [pc1_lay[gene_idx_195[g]] for g in tcell_in_vocab_195]
    all_scores = pc1_lay

    if len(bcell_scores) >= 2:
        stat_b, p_b = mannwhitneyu(bcell_scores, all_scores, alternative='less')
        auroc_b = stat_b / (len(bcell_scores) * len(all_scores))
    else:
        p_b, auroc_b = 1.0, 0.5

    if len(tcell_scores) >= 2:
        stat_t, p_t = mannwhitneyu(tcell_scores, all_scores, alternative='greater')
        auroc_t = stat_t / (len(tcell_scores) * len(all_scores))
    else:
        p_t, auroc_t = 1.0, 0.5

    layer_pc1_results.append({
        "layer": layer,
        "var_pc1": float(var_pc1_lay),
        "n_bcell_in_195": len(bcell_scores),
        "n_tcell_in_195": len(tcell_scores),
        "bcell_mean_pc1": float(np.mean(bcell_scores)) if bcell_scores else None,
        "tcell_mean_pc1": float(np.mean(tcell_scores)) if tcell_scores else None,
        "all_mean_pc1": float(np.mean(all_scores)),
        "bcell_auroc_neg": float(auroc_b),
        "bcell_p": float(p_b),
        "tcell_auroc_pos": float(auroc_t),
        "tcell_p": float(p_t),
    })
    print(f"  L{layer}: var_PC1={var_pc1_lay:.3f}, "
          f"B-cell AUROC(neg)={auroc_b:.3f} p={p_b:.4f}, "
          f"T-cell AUROC(pos)={auroc_t:.3f} p={p_t:.4f}", flush=True)

h03_result = {
    "pc1_variance_explained_195genes_L11": float(pc1_var_explained),
    "bcell_extended_test": bcell_pc1,
    "tcell_extended_test": tcell_pc1,
    "myeloid_test": myeloid_pc1,
    "tf_test": tf_pc1,
    "per_layer_pc1_results": layer_pc1_results,
    "bcell_in_vocab_195": bcell_in_vocab_195,
    "tcell_in_vocab_195": tcell_in_vocab_195,
}

with open(ITER_DIR / "h03_pc1_expanded_vocab.json", "w") as f:
    json.dump(h03_result, f, indent=2)
print("  Saved h03_pc1_expanded_vocab.json", flush=True)

print("\n=== All experiments complete ===", flush=True)
print(f"Output dir: {ITER_DIR}", flush=True)
