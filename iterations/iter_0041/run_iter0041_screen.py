"""
iter_0041 Multi-Hypothesis Screen

H01 (manifold_distance / new_family): Multi-lineage attractor screen — T-cell + Myeloid
    Test whether T-cell TFs (FOXP3, GATA3, TBX21) converge toward T-cell marker centroid
    and myeloid TFs (SPI1, CEBPA) converge toward myeloid marker centroid across layers.
    Compare attractor onset layers with B-cell/GC attractor (L3).
    Metric: rank of TF among kNN of lineage centroid, across L0..L11.

H02 (manifold_distance / refinement): BCL6 + PAX5 neighborhood characterization
    BCL6: what are its k=20 nearest neighbors at L3, L6, L11? Biological content?
    PAX5: L0 pre-wiring — top-20 neighbors at L0, L3, L11. Is L0 proximity biologically meaningful?
    Compare BCL6 and PAX5 neighbor annotation.

H03 (intrinsic_dimensionality / new_method): TwoNN intrinsic dimensionality by layer + lineage
    Estimate intrinsic dimension (TwoNN method) for:
      a) Full gene set (195 in-vocab genes)
      b) B-cell marker subset (5 genes)
      c) T-cell marker subset
      d) Myeloid marker subset
    Across all 12 layers. Does lineage-specific geometry have lower intrinsic dim than full set?
    Does L3 transition appear in ID estimates?
    (Note: PH requires ripser; fallback to TwoNN which is dependency-free)
"""

import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist, pdist

ROOT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work")
CYCLE4 = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main"
IMMUNE_H5AD = ROOT / "single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad"
ITER_DIR = ROOT / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0041"
ITER_DIR.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(42)

# Load cycle4_immune embeddings [12, 4941, 512]
# Gene list comes from h5ad var_names (aligned by construction with embedding dim 1)
print("Loading cycle4_immune embeddings...", flush=True)
emb = np.load(CYCLE4 / "layer_gene_embeddings.npy")

# Load gene list from h5ad (same gene order as embedding)
import anndata as ad
adata = ad.read_h5ad(IMMUNE_H5AD)
vocab = list(adata.var_names)  # 4941 genes
gene2idx = {g: i for i, g in enumerate(vocab)}
print(f"  shape: {emb.shape}, genes: {len(vocab)}", flush=True)

N_LAYERS = emb.shape[0]

# ─── Gene panels ─────────────────────────────────────────────────────────────
# B-cell markers (confirmed in vocab in prior iters)
BCELL_MARKERS = ["MS4A1", "CD19", "CD79A", "BLK", "PRDM1"]
# T-cell TFs and markers
TCELL_TFS = ["FOXP3", "GATA3", "TBX21", "RORC", "RUNX3"]
TCELL_MARKERS = ["CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "IL7R", "CCR7", "SELL", "CD28"]
# Myeloid TFs and markers
MYELOID_TFS = ["SPI1", "CEBPA", "CEBPB", "IRF8", "KLF4"]
MYELOID_MARKERS = ["CD14", "LYZ", "CSF1R", "ITGAM", "FCGR3A", "S100A8", "S100A9", "CD68"]

def filter_vocab(genes):
    found = [g for g in genes if g in gene2idx]
    missing = [g for g in genes if g not in gene2idx]
    return found, missing

bcell_ok, bcell_miss = filter_vocab(BCELL_MARKERS)
tcell_tf_ok, tcell_tf_miss = filter_vocab(TCELL_TFS)
tcell_mk_ok, tcell_mk_miss = filter_vocab(TCELL_MARKERS)
myeloid_tf_ok, myeloid_tf_miss = filter_vocab(MYELOID_TFS)
myeloid_mk_ok, myeloid_mk_miss = filter_vocab(MYELOID_MARKERS)

print(f"\nVocab check:")
print(f"  B-cell markers: {bcell_ok} (missing: {bcell_miss})")
print(f"  T-cell TFs: {tcell_tf_ok} (missing: {tcell_tf_miss})")
print(f"  T-cell markers: {tcell_mk_ok} (missing: {tcell_mk_miss})")
print(f"  Myeloid TFs: {myeloid_tf_ok} (missing: {myeloid_tf_miss})")
print(f"  Myeloid markers: {myeloid_mk_ok} (missing: {myeloid_mk_miss})")

def compute_centroid_ranks(layer_emb, centroid_genes, query_genes):
    """Rank of each query gene among all genes nearest to centroid of centroid_genes."""
    n_genes = layer_emb.shape[0]
    centroid_idx = [gene2idx[g] for g in centroid_genes]
    centroid_vec = layer_emb[centroid_idx].mean(axis=0)
    dists = np.linalg.norm(layer_emb - centroid_vec, axis=1)
    rank_order = np.argsort(dists)
    gene_rank = {vocab[idx]: rank+1 for rank, idx in enumerate(rank_order)}
    return {g: gene_rank.get(g, None) for g in query_genes}


# ═══════════════════════════════════════════════════════════════════════════════
# H01: Multi-lineage attractor screen
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("H01: Multi-lineage attractor screen")
print("="*60)

h01_results = []
for layer in range(N_LAYERS):
    layer_emb = emb[layer]  # [n_genes, 512]

    # T-cell attractor: TF ranks near T-cell marker centroid
    if len(tcell_mk_ok) >= 2 and len(tcell_tf_ok) >= 1:
        tcell_tf_ranks = compute_centroid_ranks(layer_emb, tcell_mk_ok, tcell_tf_ok)
        tcell_mean_rank = np.mean([r for r in tcell_tf_ranks.values() if r is not None])
    else:
        tcell_tf_ranks = {}
        tcell_mean_rank = float('nan')

    # Myeloid attractor: TF ranks near myeloid marker centroid
    if len(myeloid_mk_ok) >= 2 and len(myeloid_tf_ok) >= 1:
        myeloid_tf_ranks = compute_centroid_ranks(layer_emb, myeloid_mk_ok, myeloid_tf_ok)
        myeloid_mean_rank = np.mean([r for r in myeloid_tf_ranks.values() if r is not None])
    else:
        myeloid_tf_ranks = {}
        myeloid_mean_rank = float('nan')

    # B-cell reference (GC-TFs without BCL6)
    gc_tfs = ["BATF", "BACH2", "PAX5"]
    gc_tf_ok = [g for g in gc_tfs if g in gene2idx]
    bcell_gc_ranks = compute_centroid_ranks(layer_emb, bcell_ok, gc_tf_ok)
    bcell_gc_mean = np.mean([r for r in bcell_gc_ranks.values() if r is not None])

    row = {
        "layer": layer,
        "tcell_tf_ranks": tcell_tf_ranks,
        "tcell_mean_rank": float(tcell_mean_rank),
        "myeloid_tf_ranks": myeloid_tf_ranks,
        "myeloid_mean_rank": float(myeloid_mean_rank),
        "bcell_gc_ranks": bcell_gc_ranks,
        "bcell_gc_mean_rank": float(bcell_gc_mean),
    }
    h01_results.append(row)
    print(f"  L{layer:02d}: T-cell TF mean rank={tcell_mean_rank:.0f}  Myeloid TF mean rank={myeloid_mean_rank:.0f}  B-cell GC mean rank={bcell_gc_mean:.0f}")

# Spearman rank correlation (decreasing trend = convergence)
layers = list(range(N_LAYERS))
tcell_ranks = [r["tcell_mean_rank"] for r in h01_results]
myeloid_ranks = [r["myeloid_mean_rank"] for r in h01_results]
bcell_ranks = [r["bcell_gc_mean_rank"] for r in h01_results]

rho_t, p_t = spearmanr(layers, tcell_ranks)
rho_m, p_m = spearmanr(layers, myeloid_ranks)
rho_b, p_b = spearmanr(layers, bcell_ranks)
print(f"\nSpearman convergence (negative = converging):")
print(f"  T-cell TFs:  rho={rho_t:.3f}, p={p_t:.6f}")
print(f"  Myeloid TFs: rho={rho_m:.3f}, p={p_m:.6f}")
print(f"  B-cell GC (ref): rho={rho_b:.3f}, p={p_b:.6f}")

# Find attractor onset (first layer where mean rank < 500)
def find_onset(ranks, threshold=500):
    for i, r in enumerate(ranks):
        if r < threshold:
            return i
    return None

tcell_onset = find_onset(tcell_ranks)
myeloid_onset = find_onset(myeloid_ranks)
bcell_onset = find_onset(bcell_ranks)

print(f"\nAttractor onset (mean rank < 500):")
print(f"  T-cell: L{tcell_onset}")
print(f"  Myeloid: L{myeloid_onset}")
print(f"  B-cell (ref): L{bcell_onset}")

# Save H01
h01_out = {
    "layers": h01_results,
    "spearman_tcell": {"rho": float(rho_t), "p": float(p_t)},
    "spearman_myeloid": {"rho": float(rho_m), "p": float(p_m)},
    "spearman_bcell_ref": {"rho": float(rho_b), "p": float(p_b)},
    "attractor_onset_tcell": tcell_onset,
    "attractor_onset_myeloid": myeloid_onset,
    "attractor_onset_bcell_ref": bcell_onset,
    "tcell_tfs_in_vocab": tcell_tf_ok,
    "tcell_markers_in_vocab": tcell_mk_ok,
    "myeloid_tfs_in_vocab": myeloid_tf_ok,
    "myeloid_markers_in_vocab": myeloid_mk_ok,
}
with open(ITER_DIR / "h01_multilineage_attractor.json", "w") as f:
    json.dump(h01_out, f, indent=2)
print(f"\n  Saved h01_multilineage_attractor.json")


# ═══════════════════════════════════════════════════════════════════════════════
# H02: BCL6 + PAX5 neighborhood characterization
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("H02: BCL6 + PAX5 neighborhood characterization")
print("="*60)

K_NEIGHBORS = 20
LAYERS_OF_INTEREST = [0, 3, 6, 11]

h02_results = {}
for gene in ["BCL6", "PAX5"]:
    if gene not in gene2idx:
        print(f"  {gene}: NOT IN VOCAB")
        continue
    gene_idx = gene2idx[gene]
    h02_results[gene] = {}
    for layer in LAYERS_OF_INTEREST:
        layer_emb = emb[layer]
        gene_vec = layer_emb[gene_idx]
        dists = np.linalg.norm(layer_emb - gene_vec, axis=1)
        dists[gene_idx] = np.inf  # exclude self
        nearest_idx = np.argsort(dists)[:K_NEIGHBORS]
        neighbors = [vocab[i] for i in nearest_idx]
        neighbor_dists = [float(dists[i]) for i in nearest_idx]

        # Count B-cell/T-cell/myeloid content in neighbors
        n_bcell = sum(1 for g in neighbors if g in set(bcell_ok + ["CD19", "MS4A1", "CD79A", "BLK"]))
        n_tcell = sum(1 for g in neighbors if g in set(tcell_mk_ok + tcell_tf_ok))
        n_myeloid = sum(1 for g in neighbors if g in set(myeloid_mk_ok + myeloid_tf_ok))
        n_gc = sum(1 for g in neighbors if g in {"BATF", "BACH2", "PAX5", "BCL6", "SPIB"})

        h02_results[gene][f"L{layer}"] = {
            "top20_neighbors": neighbors,
            "top20_dists": neighbor_dists,
            "n_bcell_in_top20": n_bcell,
            "n_tcell_in_top20": n_tcell,
            "n_myeloid_in_top20": n_myeloid,
            "n_gc_tf_in_top20": n_gc,
        }
        print(f"  {gene} L{layer}: top5={neighbors[:5]}, B-cell={n_bcell}, T-cell={n_tcell}, GC-TF={n_gc}")

# Track BCL6 rank near B-cell centroid across all layers
bcell_centroid_idx = [gene2idx[g] for g in bcell_ok]
bcl6_idx = gene2idx.get("BCL6")
pax5_idx = gene2idx.get("PAX5")

bcl6_bcell_ranks = []
pax5_bcell_ranks = []
for layer in range(N_LAYERS):
    layer_emb = emb[layer]
    centroid = layer_emb[bcell_centroid_idx].mean(axis=0)
    dists = np.linalg.norm(layer_emb - centroid, axis=1)
    rank_order = np.argsort(dists)
    gene_rank = {vocab[idx]: rank+1 for rank, idx in enumerate(rank_order)}
    bcl6_bcell_ranks.append(gene_rank.get("BCL6", None))
    pax5_bcell_ranks.append(gene_rank.get("PAX5", None))

print(f"\n  BCL6 rank near B-cell centroid (all layers): {bcl6_bcell_ranks}")
print(f"  PAX5 rank near B-cell centroid (all layers): {pax5_bcell_ranks}")

h02_results["bcl6_bcell_centroid_ranks_by_layer"] = bcl6_bcell_ranks
h02_results["pax5_bcell_centroid_ranks_by_layer"] = pax5_bcell_ranks

with open(ITER_DIR / "h02_bcl6_pax5_neighborhoods.json", "w") as f:
    json.dump(h02_results, f, indent=2)
print(f"\n  Saved h02_bcl6_pax5_neighborhoods.json")


# ═══════════════════════════════════════════════════════════════════════════════
# H03: TwoNN intrinsic dimensionality by layer
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("H03: TwoNN intrinsic dimensionality")
print("="*60)

def twonn_id(X):
    """TwoNN intrinsic dimensionality estimator (Facco et al., 2017).
    For each point, compute ratio r = d2/d1 (2nd NN / 1st NN distance).
    ID = 1 / (mean(log(r)))
    """
    if X.shape[0] < 4:
        return float('nan')
    dists = cdist(X, X)
    np.fill_diagonal(dists, np.inf)
    sorted_dists = np.sort(dists, axis=1)
    d1 = sorted_dists[:, 0]
    d2 = sorted_dists[:, 1]
    # Avoid division by zero
    valid = (d1 > 0) & (d2 > 0)
    if valid.sum() < 3:
        return float('nan')
    r = d2[valid] / d1[valid]
    log_r = np.log(r)
    if log_r.mean() <= 0:
        return float('nan')
    return float(1.0 / log_r.mean())

# Gene subsets for ID estimation
# Need in-vocab indices
GENE_SUBSETS = {
    "full_set": None,  # will use all 195 genes from prior vocab (top intersect)
    "bcell_markers": bcell_ok,
    "tcell_markers": tcell_mk_ok,
    "myeloid_markers": myeloid_mk_ok,
    "gc_tfs": [g for g in ["BATF", "BACH2", "PAX5", "SPIB"] if g in gene2idx],
    "tcell_tfs": tcell_tf_ok,
    "myeloid_tfs": myeloid_tf_ok,
}

# Build "common" gene set (use same 195-gene set as prior iterations if available)
# Fall back to all in-vocab immune markers
all_immune_genes = list(set(
    bcell_ok + tcell_mk_ok + tcell_tf_ok + myeloid_mk_ok + myeloid_tf_ok +
    [g for g in ["BATF", "BACH2", "PAX5", "BCL6", "SPIB", "IRF4", "JCHAIN", "SDC1", "PRDM1"] if g in gene2idx]
))
GENE_SUBSETS["immune_panel"] = all_immune_genes
print(f"  Immune panel size: {len(all_immune_genes)}")

h03_results = []
for layer in range(N_LAYERS):
    layer_emb = emb[layer]
    row = {"layer": layer}
    for subset_name, gene_list in GENE_SUBSETS.items():
        if gene_list is None:
            # full set - use all 4941 genes (too slow with cdist, skip)
            row[f"id_{subset_name}"] = None
            continue
        if len(gene_list) < 4:
            row[f"id_{subset_name}"] = None
            continue
        idxs = [gene2idx[g] for g in gene_list if g in gene2idx]
        subset_emb = layer_emb[idxs]
        id_est = twonn_id(subset_emb)
        row[f"id_{subset_name}"] = id_est
    h03_results.append(row)
    # Print summary
    bc_id = row.get("id_bcell_markers")
    tc_id = row.get("id_tcell_markers")
    my_id = row.get("id_myeloid_markers")
    ip_id = row.get("id_immune_panel")
    def fv(v): return f"{v:.2f}" if (v is not None and not (isinstance(v, float) and np.isnan(v))) else "N/A"
    print(f"  L{layer:02d}: B-cell ID={fv(bc_id)}  T-cell ID={fv(tc_id)}  Myeloid ID={fv(my_id)}  Panel ID={fv(ip_id)}")

# Summary stats for H03
print("\n  Layer-wise ID summary:")
print(f"  {'Layer':>5}  {'B-cell':>8}  {'T-cell':>8}  {'Myeloid':>8}  {'GC-TFs':>8}  {'Panel':>8}")
for row in h03_results:
    L = row['layer']
    bc = row.get("id_bcell_markers")
    tc = row.get("id_tcell_markers")
    my = row.get("id_myeloid_markers")
    gc = row.get("id_gc_tfs")
    ip = row.get("id_immune_panel")
    def fmt(v): return f"{v:.2f}" if v is not None and not (isinstance(v, float) and np.isnan(v)) else "  N/A"
    print(f"  {L:>5}  {fmt(bc):>8}  {fmt(tc):>8}  {fmt(my):>8}  {fmt(gc):>8}  {fmt(ip):>8}")

# Check if B-cell ID drops at L3
ids_bcell = [r.get("id_bcell_markers") for r in h03_results]
ids_panel = [r.get("id_immune_panel") for r in h03_results]

# Compute Spearman correlation of each subset's ID with layer depth
for name, ids in [("B-cell", ids_bcell), ("Panel", ids_panel)]:
    valid = [(l, v) for l, v in enumerate(ids) if v is not None and not np.isnan(v)]
    if len(valid) >= 4:
        ls, vs = zip(*valid)
        rho, p = spearmanr(ls, vs)
        print(f"\n  {name} ID vs layer: Spearman rho={rho:.3f}, p={p:.4f}")

with open(ITER_DIR / "h03_twonn_intrinsic_dim.json", "w") as f:
    json.dump(h03_results, f, indent=2)
print(f"\n  Saved h03_twonn_intrinsic_dim.json")

# ═══════════════════════════════════════════════════════════════════════════════
# Summary JSON
# ═══════════════════════════════════════════════════════════════════════════════
summary = {
    "iteration": "iter_0041",
    "h01_multilineage": {
        "tcell_spearman_rho": float(rho_t),
        "tcell_spearman_p": float(p_t),
        "myeloid_spearman_rho": float(rho_m),
        "myeloid_spearman_p": float(p_m),
        "bcell_ref_spearman_rho": float(rho_b),
        "bcell_ref_spearman_p": float(p_b),
        "tcell_onset_layer": tcell_onset,
        "myeloid_onset_layer": myeloid_onset,
        "bcell_onset_layer": bcell_onset,
        "tcell_tfs_in_vocab": tcell_tf_ok,
        "myeloid_tfs_in_vocab": myeloid_tf_ok,
    },
    "h02_neighborhoods": {
        "bcl6_in_vocab": "BCL6" in gene2idx,
        "pax5_in_vocab": "PAX5" in gene2idx,
        "bcl6_bcell_ranks": bcl6_bcell_ranks,
        "pax5_bcell_ranks": pax5_bcell_ranks,
    },
    "h03_twonn": {
        "bcell_id_by_layer": ids_bcell,
        "panel_id_by_layer": ids_panel,
    }
}
with open(ITER_DIR / "iter0041_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"\nSaved iter0041_summary.json")
print("\n=== iter_0041 experiments complete ===")
