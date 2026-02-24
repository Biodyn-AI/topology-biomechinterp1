"""
iter_0042 Multi-Hypothesis Screen

H01 (intrinsic_dimensionality / refinement): TwoNN change-point analysis at L3
    Test whether the B-cell ID reduction has an inflection/change-point at L3
    (co-occurring with GC-TF attractor onset). Fit piecewise linear model split
    at L3 vs L5 vs other breakpoints. Compare ΔSLOPE B-cell vs T-cell/myeloid.
    Null: shuffle layer assignment 1000x, get slope-change distribution.

H02 (manifold_distance / new_family): Metabolic isolation generalization
    BCL6 neighbors are dominated by metabolic stress genes (NAMPT, STAT3, PFKFB3, GLUL).
    Test if STAT3, HIF1A, MYC, VEGFA, and HIF1B also cluster near this metabolic neighborhood.
    Method: at each layer, find k=20 nearest neighbors of each candidate gene.
    Count overlap with BCL6-neighborhood set (NAMPT, GLUL, PFKFB3, ACSL1, NIBAN1, FNDC3B, VMP1).
    Compare overlap to random-gene baseline (100 random genes from vocab).

H03 (manifold_distance / new_method): BATF/BACH2 convergence trajectory toward B-cell context
    For GC-TFs BATF and BACH2: at each layer (L0-L11), compute:
    (a) their rank nearest to B-cell marker centroid (MS4A1, CD19, CD79A, BLK, PRDM1)
    (b) their pairwise distance to PAX5 (already B-cell pre-wired from L0)
    (c) their pairwise distance to each other (co-localization test)
    Test if BATF and BACH2 converge toward PAX5 / B-cell centroid from L3 onward.
    Null: 5 random TFs in vocab, track their centroid ranks across layers.
"""

import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist, pdist
import anndata as ad

ROOT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work")
CYCLE4 = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main"
IMMUNE_H5AD = ROOT / "single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad"
ITER_DIR = ROOT / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0042"
ITER_DIR.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(42)

print("Loading cycle4_immune embeddings...", flush=True)
emb = np.load(CYCLE4 / "layer_gene_embeddings.npy")

adata = ad.read_h5ad(IMMUNE_H5AD)
vocab = list(adata.var_names)
gene2idx = {g: i for i, g in enumerate(vocab)}
print(f"  shape: {emb.shape}, genes: {len(vocab)}", flush=True)

N_LAYERS = emb.shape[0]  # 12

# ─── Gene panels ─────────────────────────────────────────────────────────────
BCELL_MARKERS = ["MS4A1", "CD19", "CD79A", "BLK", "PRDM1"]
TCELL_MARKERS = ["CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "IL7R", "CCR7", "SELL", "CD28"]
MYELOID_MARKERS = ["CD14", "LYZ", "ITGAM", "FCGR3A", "S100A8", "S100A9"]

# Filter to vocab
bcell_idx = [gene2idx[g] for g in BCELL_MARKERS if g in gene2idx]
tcell_idx = [gene2idx[g] for g in TCELL_MARKERS if g in gene2idx]
myeloid_idx = [gene2idx[g] for g in MYELOID_MARKERS if g in gene2idx]

print(f"  B-cell markers in vocab: {len(bcell_idx)}")
print(f"  T-cell markers in vocab: {len(tcell_idx)}")
print(f"  Myeloid markers in vocab: {len(myeloid_idx)}")

# ─── H01: TwoNN change-point analysis ────────────────────────────────────────
print("\n=== H01: TwoNN change-point analysis ===", flush=True)

def twonn_id(X):
    """TwoNN intrinsic dimensionality estimator (Facco et al. 2017)."""
    n = X.shape[0]
    if n < 3:
        return float('nan')
    dists = cdist(X, X, metric='euclidean')
    np.fill_diagonal(dists, np.inf)
    sorted_dists = np.sort(dists, axis=1)
    r1 = sorted_dists[:, 0]
    r2 = sorted_dists[:, 1]
    mu = r2 / (r1 + 1e-12)
    # Remove mu <= 1 (can happen due to ties)
    mu = mu[mu > 1]
    if len(mu) < 3:
        return float('nan')
    # Empirical CDF
    mu_sorted = np.sort(mu)
    F = np.arange(1, len(mu_sorted) + 1) / len(mu_sorted)
    # Fit: -log(1 - F) = d * log(mu) => d = slope via OLS through origin
    log_mu = np.log(mu_sorted)
    log_1mF = -np.log(1 - F + 1e-12)
    # OLS: d = sum(log_mu * log_1mF) / sum(log_mu^2)
    d = np.sum(log_mu * log_1mF) / np.sum(log_mu**2)
    return float(d)

# Compute ID per layer for each lineage
bcell_ids = []
tcell_ids = []
myeloid_ids = []
for layer in range(N_LAYERS):
    X_b = emb[layer, bcell_idx, :]
    X_t = emb[layer, tcell_idx, :]
    X_m = emb[layer, myeloid_idx, :]
    bcell_ids.append(twonn_id(X_b))
    tcell_ids.append(twonn_id(X_t))
    myeloid_ids.append(twonn_id(X_m))

print(f"  B-cell IDs:   {[round(x,2) for x in bcell_ids]}", flush=True)
print(f"  T-cell IDs:   {[round(x,2) for x in tcell_ids]}", flush=True)
print(f"  Myeloid IDs:  {[round(x,2) for x in myeloid_ids]}", flush=True)

# Change-point analysis: for each lineage, test all break-points bp in {1..10}
# Fit y = a*layer + b for layers < bp, y = c*layer + d for layers >= bp
# Return breakpoint with min total RSS
def piecewise_linear_fit(ids, layers=None):
    ids = np.array(ids, dtype=float)
    n = len(ids)
    if layers is None:
        layers = np.arange(n, dtype=float)
    best_bp = None
    best_rss = np.inf
    best_slopes = None
    for bp in range(2, n - 2):
        x1, y1 = layers[:bp], ids[:bp]
        x2, y2 = layers[bp:], ids[bp:]
        # OLS for each segment
        def ols_slope_intercept(x, y):
            if len(x) < 2:
                return 0, float(np.mean(y))
            A = np.vstack([x, np.ones(len(x))]).T
            result = np.linalg.lstsq(A, y, rcond=None)
            m, b = result[0]
            return m, b
        m1, b1 = ols_slope_intercept(x1, y1)
        m2, b2 = ols_slope_intercept(x2, y2)
        r1 = y1 - (m1 * x1 + b1)
        r2 = y2 - (m2 * x2 + b2)
        rss = np.sum(r1**2) + np.sum(r2**2)
        if rss < best_rss:
            best_rss = rss
            best_bp = bp
            best_slopes = (m1, m2)
    return best_bp, best_slopes, best_rss

layers_arr = np.arange(N_LAYERS, dtype=float)
bp_b, slopes_b, rss_b = piecewise_linear_fit(bcell_ids, layers_arr)
bp_t, slopes_t, rss_t = piecewise_linear_fit(tcell_ids, layers_arr)
bp_m, slopes_m, rss_m = piecewise_linear_fit(myeloid_ids, layers_arr)

print(f"  B-cell: best breakpoint L{bp_b}, slopes {slopes_b[0]:.3f} -> {slopes_b[1]:.3f}, slope_change={slopes_b[1]-slopes_b[0]:.3f}", flush=True)
print(f"  T-cell: best breakpoint L{bp_t}, slopes {slopes_t[0]:.3f} -> {slopes_t[1]:.3f}, slope_change={slopes_t[1]-slopes_t[0]:.3f}", flush=True)
print(f"  Myeloid: best breakpoint L{bp_m}, slopes {slopes_m[0]:.3f} -> {slopes_m[1]:.3f}, slope_change={slopes_m[1]-slopes_m[0]:.3f}", flush=True)

# Null: shuffle layer indices 500x for B-cell, get slope-change distribution
null_slope_changes = []
for _ in range(500):
    ids_shuffled = rng.permutation(bcell_ids)
    _, sl, _ = piecewise_linear_fit(ids_shuffled, layers_arr)
    null_slope_changes.append(sl[1] - sl[0])
null_slope_changes = np.array(null_slope_changes)
bcell_slope_change = slopes_b[1] - slopes_b[0]
pval_slope_change = float(np.mean(null_slope_changes <= bcell_slope_change))
print(f"  B-cell slope change {bcell_slope_change:.3f} vs null mean {null_slope_changes.mean():.3f} +/- {null_slope_changes.std():.3f}, p={pval_slope_change:.4f}", flush=True)

h01_results = {
    "bcell_ids_by_layer": bcell_ids,
    "tcell_ids_by_layer": tcell_ids,
    "myeloid_ids_by_layer": myeloid_ids,
    "bcell_best_breakpoint": bp_b,
    "tcell_best_breakpoint": bp_t,
    "myeloid_best_breakpoint": bp_m,
    "bcell_slope_pre": float(slopes_b[0]),
    "bcell_slope_post": float(slopes_b[1]),
    "bcell_slope_change": float(bcell_slope_change),
    "tcell_slope_pre": float(slopes_t[0]),
    "tcell_slope_post": float(slopes_t[1]),
    "tcell_slope_change": float(slopes_t[1] - slopes_t[0]),
    "myeloid_slope_pre": float(slopes_m[0]),
    "myeloid_slope_post": float(slopes_m[1]),
    "myeloid_slope_change": float(slopes_m[1] - slopes_m[0]),
    "null_slope_change_mean": float(null_slope_changes.mean()),
    "null_slope_change_std": float(null_slope_changes.std()),
    "null_pval_bcell_slope_change": pval_slope_change,
    "n_null_permutations": 500
}

with open(ITER_DIR / "h01_twonn_changepoint.json", "w") as f:
    json.dump(h01_results, f, indent=2)
print("  H01 saved.", flush=True)


# ─── H02: Metabolic isolation generalization ─────────────────────────────────
print("\n=== H02: Metabolic isolation generalization ===", flush=True)

# BCL6 metabolic neighborhood (stable across layers in iter_0041)
BCL6_METAB_CLUSTER = ["NAMPT", "GLUL", "PFKFB3", "ACSL1", "NIBAN1", "FNDC3B", "VMP1", "STAT3", "CEBPD", "TRIB1"]
metab_idx_set = set(gene2idx[g] for g in BCL6_METAB_CLUSTER if g in gene2idx)
metab_in_vocab = [g for g in BCL6_METAB_CLUSTER if g in gene2idx]
print(f"  Metabolic cluster genes in vocab: {metab_in_vocab}", flush=True)

# Candidate metabolic TFs/oncogenes to test
CANDIDATE_GENES = ["STAT3", "MYC", "HIF1A", "VEGFA", "SLC2A1", "LDHA", "MCM2", "PCNA", "BCL6"]
# Also add some B-cell markers as negative controls
NEG_CONTROLS = ["MS4A1", "CD79A", "PAX5", "BLK"]
ALL_TEST = CANDIDATE_GENES + NEG_CONTROLS

# For each test gene, at each layer, find k=20 NN. Count overlap with metab cluster.
def get_knn_overlap(emb_layer, gene_idx, reference_set, k=20):
    """Get k nearest neighbors of gene_idx and count overlap with reference_set."""
    n_genes = emb_layer.shape[0]
    query = emb_layer[gene_idx:gene_idx+1]
    dists = cdist(query, emb_layer, metric='euclidean')[0]
    dists[gene_idx] = np.inf  # exclude self
    nn_indices = np.argsort(dists)[:k]
    overlap = sum(1 for idx in nn_indices if idx in reference_set)
    return int(overlap), list(nn_indices.tolist())

results_by_gene = {}
for gene in ALL_TEST:
    if gene not in gene2idx:
        print(f"  {gene}: not in vocab", flush=True)
        continue
    gidx = gene2idx[gene]
    layer_overlaps = []
    for layer in range(N_LAYERS):
        overlap, nn_list = get_knn_overlap(emb[layer], gidx, metab_idx_set, k=20)
        layer_overlaps.append(overlap)
    results_by_gene[gene] = layer_overlaps
    print(f"  {gene}: overlaps by layer = {layer_overlaps}", flush=True)

# Random gene baseline: 100 random genes
random_gene_idxs = rng.choice(len(vocab), size=100, replace=False)
random_overlaps_by_layer = []
for layer in range(N_LAYERS):
    layer_overlaps_rand = []
    for gidx in random_gene_idxs:
        # skip if gene is in metab cluster itself
        if gidx in metab_idx_set:
            continue
        overlap, _ = get_knn_overlap(emb[layer], int(gidx), metab_idx_set, k=20)
        layer_overlaps_rand.append(overlap)
    random_overlaps_by_layer.append({
        "mean": float(np.mean(layer_overlaps_rand)),
        "std": float(np.std(layer_overlaps_rand)),
        "p95": float(np.percentile(layer_overlaps_rand, 95))
    })

print(f"  Random baseline overlap (L6): mean={random_overlaps_by_layer[6]['mean']:.2f} std={random_overlaps_by_layer[6]['std']:.2f} p95={random_overlaps_by_layer[6]['p95']:.1f}", flush=True)

h02_results = {
    "metabolic_cluster_genes": BCL6_METAB_CLUSTER,
    "metabolic_cluster_in_vocab": metab_in_vocab,
    "candidate_gene_overlaps": results_by_gene,
    "random_baseline_by_layer": random_overlaps_by_layer
}

with open(ITER_DIR / "h02_metabolic_isolation.json", "w") as f:
    json.dump(h02_results, f, indent=2)
print("  H02 saved.", flush=True)


# ─── H03: BATF/BACH2 convergence trajectory ──────────────────────────────────
print("\n=== H03: BATF/BACH2 convergence trajectory ===", flush=True)

GC_TFS = ["BATF", "BACH2"]
pax5_idx = gene2idx.get("PAX5")

# Compute B-cell centroid at each layer
bcell_centroids = np.array([emb[l, bcell_idx, :].mean(axis=0) for l in range(N_LAYERS)])  # [12, 512]

# For each GC-TF: rank near B-cell centroid, distance to PAX5, distance to each other
gc_results = {}
for gname in GC_TFS:
    if gname not in gene2idx:
        print(f"  {gname}: not in vocab", flush=True)
        continue
    gidx = gene2idx[gname]

    ranks_bcell_centroid = []
    dists_to_pax5 = []

    for layer in range(N_LAYERS):
        # rank near B-cell centroid
        centroid = bcell_centroids[layer]
        all_embs = emb[layer]  # [n_genes, 512]
        dists_to_centroid = np.linalg.norm(all_embs - centroid, axis=1)
        rank = int(np.sum(dists_to_centroid < dists_to_centroid[gidx]))  # rank = number of genes closer
        ranks_bcell_centroid.append(rank)

        # distance to PAX5
        if pax5_idx is not None:
            d_pax5 = float(np.linalg.norm(emb[layer, gidx] - emb[layer, pax5_idx]))
            dists_to_pax5.append(d_pax5)

    gc_results[gname] = {
        "ranks_near_bcell_centroid": ranks_bcell_centroid,
        "dists_to_pax5": dists_to_pax5
    }
    print(f"  {gname} ranks near B-cell centroid: {ranks_bcell_centroid}", flush=True)
    print(f"  {gname} dists to PAX5: {[round(d,2) for d in dists_to_pax5]}", flush=True)

# Pairwise distance between BATF and BACH2 across layers
batf_idx = gene2idx.get("BATF")
bach2_idx = gene2idx.get("BACH2")
if batf_idx is not None and bach2_idx is not None:
    batf_bach2_dists = [float(np.linalg.norm(emb[l, batf_idx] - emb[l, bach2_idx])) for l in range(N_LAYERS)]
    print(f"  BATF-BACH2 pairwise dists: {[round(d,2) for d in batf_bach2_dists]}", flush=True)
    gc_results["BATF_BACH2_pairwise_dists"] = batf_bach2_dists

# Null: 10 random TFs from vocab — track their mean rank near B-cell centroid
# Use genes that are TF-like (just pick random non-B-cell genes as null)
all_indices = list(range(len(vocab)))
bc_set = set(bcell_idx)
non_bcell_indices = [i for i in all_indices if i not in bc_set]
null_tfs_idx = rng.choice(non_bcell_indices, size=10, replace=False)
null_ranks_by_layer = []
for layer in range(N_LAYERS):
    centroid = bcell_centroids[layer]
    all_embs = emb[layer]
    dists_to_centroid = np.linalg.norm(all_embs - centroid, axis=1)
    null_ranks = []
    for gidx in null_tfs_idx:
        rank = int(np.sum(dists_to_centroid < dists_to_centroid[gidx]))
        null_ranks.append(rank)
    null_ranks_by_layer.append({"mean": float(np.mean(null_ranks)), "std": float(np.std(null_ranks))})
gc_results["null_random_gene_ranks_by_layer"] = null_ranks_by_layer
print(f"  Null random gene mean ranks (L3): {null_ranks_by_layer[3]}", flush=True)

# Spearman correlation: rank vs layer for BATF and BACH2
if "BATF" in gc_results and "BACH2" in gc_results:
    rho_batf, p_batf = spearmanr(np.arange(N_LAYERS), gc_results["BATF"]["ranks_near_bcell_centroid"])
    rho_bach2, p_bach2 = spearmanr(np.arange(N_LAYERS), gc_results["BACH2"]["ranks_near_bcell_centroid"])
    gc_results["spearman_batf_rank_vs_layer"] = {"rho": float(rho_batf), "p": float(p_batf)}
    gc_results["spearman_bach2_rank_vs_layer"] = {"rho": float(rho_bach2), "p": float(p_bach2)}
    print(f"  Spearman BATF rank vs layer: rho={rho_batf:.3f} p={p_batf:.4f}", flush=True)
    print(f"  Spearman BACH2 rank vs layer: rho={rho_bach2:.3f} p={p_bach2:.4f}", flush=True)

with open(ITER_DIR / "h03_gc_tf_convergence.json", "w") as f:
    json.dump(gc_results, f, indent=2)
print("  H03 saved.", flush=True)

print("\n=== All hypotheses complete ===", flush=True)
