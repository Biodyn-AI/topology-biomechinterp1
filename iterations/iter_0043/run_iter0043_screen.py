"""
iter_0043 Multi-Hypothesis Screen

H01 (manifold_distance / new_method): GC repression circuit — PRDM1/IRF4 anti-convergence test
    Full pairwise distance matrix for {PAX5, BATF, BACH2, BCL6, PRDM1, IRF4} across L0-L11.
    Test if PRDM1/IRF4 (BCL6 repression targets) show DIVERGENCE from B-cell centroid
    while BATF/BACH2 converge. If so, the model encodes regulatory directionality geometrically.

H02 (manifold_distance / new_method): BCL6 trajectory from B-cell centroid
    Track BCL6 distance to B-cell centroid (MS4A1, CD19, CD79A, BLK, PAX5) across layers.
    Check if BCL6 diverges outward (as it is metabolically isolated) or stays constant.
    Compare to null: 20 random non-B-cell genes.

H03 (intrinsic_dimensionality / new_method): T-cell subtype stratification at L7 ID breakpoint
    The T-cell ID expansion breaks at L7. Test if this is driven by effector vs regulatory vs
    exhausted T-cell subtypes (using CD8A/GZMB for effector, FOXP3/IL2RA for Treg, LAG3/PDCD1
    for exhaustion). Compute per-subtype TwoNN ID at each layer; see which subtype drives the
    L7 expansion.
"""

import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
import anndata as ad

ROOT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work")
CYCLE4 = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main"
IMMUNE_H5AD = ROOT / "single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad"
ITER_DIR = ROOT / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0043"
ITER_DIR.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(42)

print("Loading cycle4_immune embeddings...", flush=True)
emb = np.load(CYCLE4 / "layer_gene_embeddings.npy")
adata = ad.read_h5ad(IMMUNE_H5AD)
vocab = list(adata.var_names)
gene2idx = {g: i for i, g in enumerate(vocab)}
print(f"  shape: {emb.shape}, genes: {len(vocab)}", flush=True)

N_LAYERS = emb.shape[0]  # 12

# =====================================================================
# H01: GC repression circuit — pairwise distance matrix + anti-convergence
# =====================================================================
print("\n=== H01: GC repression circuit ===", flush=True)

GC_CIRCUIT = ["PAX5", "BATF", "BACH2", "BCL6", "PRDM1", "IRF4"]
B_CELL_CENTROID_GENES = ["MS4A1", "CD19", "CD79A", "BLK", "PAX5"]

in_vocab_circuit = {g: gene2idx[g] for g in GC_CIRCUIT if g in gene2idx}
b_centroid_genes_in = [g for g in B_CELL_CENTROID_GENES if g in gene2idx]
print(f"  Circuit in-vocab: {list(in_vocab_circuit.keys())}")
print(f"  B-cell centroid genes: {b_centroid_genes_in}")

h1_results = {}

# Pairwise distance matrix at L0, L3, L6, L11
for layer in [0, 3, 6, 11]:
    E = emb[layer]  # [N_genes, 512]
    circuit_embs = {g: E[idx] for g, idx in in_vocab_circuit.items()}
    genes = list(circuit_embs.keys())
    vecs = np.stack([circuit_embs[g] for g in genes])
    dists = cdist(vecs, vecs, metric='euclidean')
    h1_results[f"L{layer}_pairwise"] = {
        "genes": genes,
        "distances": dists.tolist()
    }
    print(f"  L{layer} pairwise matrix ({len(genes)}x{len(genes)}):")
    for i, g1 in enumerate(genes):
        for j, g2 in enumerate(genes):
            if j > i:
                print(f"    {g1}-{g2}: {dists[i,j]:.3f}")

# B-cell centroid
b_centroid_by_layer = []
for layer in range(N_LAYERS):
    E = emb[layer]
    c = np.mean(np.stack([E[gene2idx[g]] for g in b_centroid_genes_in]), axis=0)
    b_centroid_by_layer.append(c)

# Rank and distance trajectories for each gene
trajectory = {}
for gene, idx in in_vocab_circuit.items():
    ranks = []
    dists_to_centroid = []
    for layer in range(N_LAYERS):
        E = emb[layer]
        gene_vec = E[idx]
        centroid = b_centroid_by_layer[layer]
        d = np.linalg.norm(gene_vec - centroid)
        dists_to_centroid.append(float(d))
        # rank among all genes
        all_dists = np.linalg.norm(E - centroid, axis=1)
        rank = int(np.sum(all_dists < d))
        ranks.append(rank)
    trajectory[gene] = {
        "dists_to_b_centroid": dists_to_centroid,
        "ranks_to_b_centroid": ranks
    }
    rho, pval = spearmanr(np.arange(N_LAYERS), dists_to_centroid)
    print(f"  {gene}: dist L0={dists_to_centroid[0]:.3f} L3={dists_to_centroid[3]:.3f} L6={dists_to_centroid[6]:.3f} L11={dists_to_centroid[11]:.3f} | rho={rho:.3f} p={pval:.4f}")

# Expected: BATF, BACH2 → converge (rho < 0); PRDM1, IRF4 → diverge (rho > 0)
h1_results["trajectories"] = trajectory

# Null: 20 random genes
all_gene_indices = list(range(len(vocab)))
random_genes = rng.choice(all_gene_indices, size=20, replace=False)
null_rhos = []
for idx in random_genes:
    d_series = [float(np.linalg.norm(emb[l][idx] - b_centroid_by_layer[l])) for l in range(N_LAYERS)]
    rho, _ = spearmanr(np.arange(N_LAYERS), d_series)
    null_rhos.append(float(rho))
h1_results["null_rho_mean"] = float(np.mean(null_rhos))
h1_results["null_rho_std"] = float(np.std(null_rhos))
print(f"  Null rho (n=20): mean={np.mean(null_rhos):.3f} ± {np.std(null_rhos):.3f}")

# Summary: compare activators (BATF, BACH2) vs repressors (PRDM1, IRF4) by rho
gene_rhos = {}
for gene in in_vocab_circuit:
    d_series = trajectory[gene]["dists_to_b_centroid"]
    rho, pval = spearmanr(np.arange(N_LAYERS), d_series)
    gene_rhos[gene] = {"rho": float(rho), "pval": float(pval)}

h1_results["gene_rhos"] = gene_rhos

# Save
h1_path = ITER_DIR / "h01_gc_repression_circuit.json"
with open(h1_path, "w") as f:
    json.dump(h1_results, f, indent=2)
print(f"  Saved: {h1_path}")

# =====================================================================
# H02: BCL6 trajectory from B-cell centroid
# =====================================================================
print("\n=== H02: BCL6 divergence from B-cell centroid ===", flush=True)

h2_results = {}

# Check BCL6 vocab membership
BCL6_GENES = ["BCL6", "BATF", "BACH2", "PAX5"]  # reference
bcl6_in = {g: gene2idx[g] for g in BCL6_GENES if g in gene2idx}
print(f"  In vocab: {list(bcl6_in.keys())}")

# B-cell centroid with PAX5 included
B_CENTROID_FULL = ["MS4A1", "CD19", "CD79A", "BLK", "PAX5"]
b_full_in = [g for g in B_CENTROID_FULL if g in gene2idx]
b_centroid_full_by_layer = []
for layer in range(N_LAYERS):
    E = emb[layer]
    c = np.mean(np.stack([E[gene2idx[g]] for g in b_full_in]), axis=0)
    b_centroid_full_by_layer.append(c)

# BCL6 trajectory
if "BCL6" in gene2idx:
    bcl6_idx = gene2idx["BCL6"]
    bcl6_dists = [float(np.linalg.norm(emb[l][bcl6_idx] - b_centroid_full_by_layer[l])) for l in range(N_LAYERS)]
    bcl6_ranks = []
    for layer in range(N_LAYERS):
        E = emb[layer]
        all_dists = np.linalg.norm(E - b_centroid_full_by_layer[layer], axis=1)
        bcl6_ranks.append(int(np.sum(all_dists < bcl6_dists[layer])))
    rho_bcl6, pval_bcl6 = spearmanr(np.arange(N_LAYERS), bcl6_dists)
    print(f"  BCL6 dists: {[f'{d:.2f}' for d in bcl6_dists]}")
    print(f"  BCL6 ranks: {bcl6_ranks}")
    print(f"  BCL6 rho={rho_bcl6:.3f} p={pval_bcl6:.4f}")
    h2_results["BCL6"] = {
        "dists_to_b_centroid": bcl6_dists,
        "ranks_to_b_centroid": bcl6_ranks,
        "rho": float(rho_bcl6),
        "pval": float(pval_bcl6)
    }
else:
    print("  BCL6 not in vocab!")
    h2_results["BCL6"] = "not_in_vocab"

# Null: 20 random non-B-cell genes
null_rhos_h2 = []
null_genes_sampled = rng.choice(len(vocab), size=50, replace=False)
non_b_count = 0
for idx in null_genes_sampled:
    if non_b_count >= 20:
        break
    g = vocab[idx]
    if g not in B_CENTROID_FULL:
        d_s = [float(np.linalg.norm(emb[l][idx] - b_centroid_full_by_layer[l])) for l in range(N_LAYERS)]
        rho, _ = spearmanr(np.arange(N_LAYERS), d_s)
        null_rhos_h2.append(float(rho))
        non_b_count += 1

h2_results["null_rho_mean"] = float(np.mean(null_rhos_h2))
h2_results["null_rho_std"] = float(np.std(null_rhos_h2))
print(f"  Null rho (n=20): mean={np.mean(null_rhos_h2):.3f} ± {np.std(null_rhos_h2):.3f}")

h2_path = ITER_DIR / "h02_bcl6_divergence.json"
with open(h2_path, "w") as f:
    json.dump(h2_results, f, indent=2)
print(f"  Saved: {h2_path}")

# =====================================================================
# H03: T-cell subtype TwoNN ID stratification at L7
# =====================================================================
print("\n=== H03: T-cell subtype ID stratification ===", flush=True)

# T-cell subtypes by gene markers
# Effector: CD8A, GZMB, PRF1
# Treg: FOXP3, IL2RA, CTLA4
# Exhaustion: LAG3, PDCD1, HAVCR2
# Helper: CD4, CXCR5, ICOS
T_SUBTYPES = {
    "effector_cd8": ["CD8A", "GZMB", "PRF1", "NKG7", "GNLY"],
    "treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "exhaustion": ["LAG3", "PDCD1", "HAVCR2", "TIGIT"],
    "helper_tfh": ["CD4", "CXCR5", "ICOS", "BCL6", "SH2D1A"],
    "general_t": ["CD3E", "CD3D", "TRAC", "TRBC1"]
}

h3_results = {}

def twonn_id(X):
    """TwoNN intrinsic dimensionality (Facco et al. 2017)."""
    from scipy.spatial.distance import pdist, squareform
    if X.shape[0] < 3:
        return float("nan")
    D = squareform(pdist(X, metric='euclidean'))
    np.fill_diagonal(D, np.inf)
    sorted_D = np.sort(D, axis=1)
    r1 = sorted_D[:, 0]
    r2 = sorted_D[:, 1]
    valid = (r1 > 0) & (r2 > 0)
    if valid.sum() < 2:
        return float("nan")
    mu = r2[valid] / r1[valid]
    mu_sorted = np.sort(mu)
    n = len(mu_sorted)
    F = np.arange(1, n + 1) / n
    # Linear regression of -log(1-F) ~ d*log(mu)
    log_mu = np.log(mu_sorted)
    neg_log_surv = -np.log(1.0 - F + 1e-10)
    valid2 = log_mu > 0
    if valid2.sum() < 2:
        return float("nan")
    d = np.polyfit(log_mu[valid2], neg_log_surv[valid2], 1)[0]
    return float(d)

for subtype, genes in T_SUBTYPES.items():
    in_v = [g for g in genes if g in gene2idx]
    print(f"  {subtype}: {in_v} ({len(in_v)}/{len(genes)} in vocab)")
    if len(in_v) < 2:
        h3_results[subtype] = {"status": "insufficient_genes", "in_vocab": in_v}
        continue
    ids_per_layer = []
    for layer in range(N_LAYERS):
        E = emb[layer]
        X = np.stack([E[gene2idx[g]] for g in in_v])
        d = twonn_id(X)
        ids_per_layer.append(d)
    # Check if there's an inflection at L7
    pre_l7 = np.mean(ids_per_layer[:7])
    post_l7 = np.mean(ids_per_layer[7:])
    slope_change = post_l7 - pre_l7
    print(f"    IDs: {[f'{d:.1f}' for d in ids_per_layer]}")
    print(f"    Pre-L7 mean: {pre_l7:.2f}, Post-L7 mean: {post_l7:.2f}, change: {slope_change:.2f}")
    h3_results[subtype] = {
        "in_vocab": in_v,
        "ids_per_layer": ids_per_layer,
        "pre_l7_mean": float(pre_l7),
        "post_l7_mean": float(post_l7),
        "l7_slope_change": float(slope_change)
    }

h3_path = ITER_DIR / "h03_tcell_subtype_id.json"
with open(h3_path, "w") as f:
    json.dump(h3_results, f, indent=2)
print(f"  Saved: {h3_path}")

print("\n=== iter_0043 screen complete ===", flush=True)
