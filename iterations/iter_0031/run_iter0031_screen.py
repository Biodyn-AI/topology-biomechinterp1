"""
iter_0031 Multi-Hypothesis Screen — 195 in-vocabulary genes (OOV-corrected)

H01 (graph_topology / new_method): kNN community detection + spectral gap on 195 in-vocab genes
    k=10 kNN graph per layer. Normalized spectral gap (lambda2 / lambda_max).
    Spearman rho(layer, spectral_gap). Compare to shuffled null.

H02 (manifold_distance / refinement-rescue): STRING score → embedding distance on 195 in-vocab genes
    For gene pairs with STRING score >= 0.4, compute L2 embedding distance per layer.
    Spearman rho(STRING_score, L2_dist) per layer. AUROC: high vs low STRING score pairs.
    Slope of STRING→distance correlation across layers (first clean test on 195 in-vocab only).

H03 (intrinsic_dimensionality / new_method): Participation ratio + L2 norm distribution per layer
    For 195 in-vocab genes at each layer: mean-center, compute SVD, participation ratio (PR).
    Also compute L2 norm distribution statistics (mean, std, skewness) per layer.
    Test if PR changes monotonically with layer (Spearman rho).
"""

import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu, skew
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import laplacian
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0031"
ITER_DIR.mkdir(parents=True, exist_ok=True)

CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")
ITER15_DIR = PROJECT / "iterations" / "iter_0015"
STRING_CACHE = ITER15_DIR / "string_ppi_score04_cache.json"

rng = np.random.default_rng(42)

# ─── Load embeddings ──────────────────────────────────────────────────────────
print("Loading embeddings...", flush=True)
emb = np.load(CYCLE1 / "layer_gene_embeddings.npy")   # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb.shape
print(f"  Shape: {emb.shape}", flush=True)

# ─── Load vocab + named genes ─────────────────────────────────────────────────
with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f if line.strip()]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}

import csv
EDGES_FILE = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
                  "/subproject_38_geometric_residual_stream_interpretability"
                  "/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv")
named_genes_set = set()
with open(EDGES_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        named_genes_set.add(row["source"])
        named_genes_set.add(row["target"])
named_genes = sorted(named_genes_set)
print(f"  Total named genes: {len(named_genes)}", flush=True)

# ─── Filter to in-vocabulary (non-zero L0 embedding) ─────────────────────────
OOV_GENES = {"FOS","HLA-A","HLA-DPB1","JUNB","KLF6","LDHA","LGALS1",
             "NCAM1","NCOA3","NR4A3","PAX5","PTGS2","TBXAS1","TNF"}

inv_genes = [g for g in named_genes if g in gene_to_emb_idx and g not in OOV_GENES]
inv_indices = [gene_to_emb_idx[g] for g in inv_genes]
print(f"  In-vocab named genes: {len(inv_genes)}", flush=True)  # expect 195

# Embeddings for in-vocab genes: [12, 195, 512]
emb_inv = emb[:, inv_indices, :]
print(f"  In-vocab emb shape: {emb_inv.shape}", flush=True)

# ─── H01: kNN topology + spectral gap per layer ───────────────────────────────
print("\n=== H01: kNN topology + spectral gap ===", flush=True)

def spectral_gap(X, k=10):
    """Build kNN graph, compute normalized Laplacian eigenvalues, return lambda2/lambda_max."""
    n = X.shape[0]
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(X)
    dists, idxs = nbrs.kneighbors(X)
    # Build symmetric adjacency
    rows, cols = [], []
    for i in range(n):
        for j in idxs[i, 1:]:  # skip self
            rows.append(i); cols.append(j)
            rows.append(j); cols.append(i)
    data = np.ones(len(rows))
    A = csr_matrix((data, (rows, cols)), shape=(n, n))
    A.data = np.ones(A.nnz)  # binarize
    # Normalized Laplacian eigenvalues (sparse approx)
    from scipy.sparse.linalg import eigsh
    L = laplacian(A, normed=True)
    # get smallest non-trivial + largest eigenvalues
    try:
        # Get bottom 4 eigenvalues
        vals_low = eigsh(L, k=4, which='SM', return_eigenvectors=False)
        vals_low = np.sort(vals_low)
        # Get top eigenvalue
        vals_high = eigsh(L, k=2, which='LM', return_eigenvectors=False)
        lambda2 = vals_low[1] if vals_low[1] > 1e-10 else vals_low[2]
        lambda_max = vals_high[-1]
        return lambda2, lambda_max, lambda2 / (lambda_max + 1e-12)
    except Exception as e:
        return np.nan, np.nan, np.nan

h01_results = []
for layer in range(N_LAYERS):
    X = emb_inv[layer]
    l2, lmax, ratio = spectral_gap(X, k=10)
    h01_results.append({"layer": layer, "lambda2": float(l2) if not np.isnan(l2) else None,
                         "lambda_max": float(lmax) if not np.isnan(lmax) else None,
                         "spectral_gap_ratio": float(ratio) if not np.isnan(ratio) else None})
    print(f"  L{layer:02d}: lambda2={l2:.4f}, lmax={lmax:.4f}, ratio={ratio:.4f}", flush=True)

# Spearman rho(layer, spectral_gap_ratio)
layers = np.arange(N_LAYERS)
ratios = np.array([r["spectral_gap_ratio"] for r in h01_results])
valid = ~np.isnan(ratios)
rho_h01, p_h01 = spearmanr(layers[valid], ratios[valid])
print(f"  Spearman rho(layer, spectral_gap_ratio) = {rho_h01:.4f}, p = {p_h01:.3e}", flush=True)

# Null: shuffle gene labels (random permutation of embeddings)
null_ratios = []
X_flat = emb_inv[11].copy()
rng.shuffle(X_flat)  # permute rows
for layer in range(N_LAYERS):
    X_null = rng.standard_normal(X_flat.shape).astype(np.float32)
    _, _, null_r = spectral_gap(X_null, k=10)
    null_ratios.append(null_r)
null_ratios = np.array(null_ratios)
print(f"  Null spectral_gap_ratio mean={np.nanmean(null_ratios):.4f}, std={np.nanstd(null_ratios):.4f}", flush=True)
print(f"  Real spectral_gap_ratio mean={np.nanmean(ratios[valid]):.4f}, std={np.nanstd(ratios[valid]):.4f}", flush=True)

h01_artifact = {
    "per_layer": h01_results,
    "spearman_rho_layer_vs_ratio": float(rho_h01),
    "spearman_p": float(p_h01),
    "null_mean": float(np.nanmean(null_ratios)),
    "null_std": float(np.nanstd(null_ratios)),
    "real_mean": float(np.nanmean(ratios[valid])),
    "real_std": float(np.nanstd(ratios[valid])),
    "n_inv_genes": len(inv_genes)
}
with open(ITER_DIR / "h01_spectral_gap_195.json", "w") as f:
    json.dump(h01_artifact, f, indent=2)
print("  Saved h01_spectral_gap_195.json", flush=True)

# ─── H02: STRING score → embedding distance on 195 in-vocab genes ─────────────
print("\n=== H02: STRING score → embedding distance (195 in-vocab) ===", flush=True)

with open(STRING_CACHE) as f:
    string_data = json.load(f)

# Build mapping: (geneA, geneB) -> score
string_pairs = {}
for key, score in string_data.items():
    parts = key.split("||")
    if len(parts) == 2:
        a, b = parts[0].strip(), parts[1].strip()
        string_pairs[(a, b)] = score
        string_pairs[(b, a)] = score

inv_set = set(inv_genes)

# Collect valid pairs
valid_pairs = []
for (a, b), score in string_pairs.items():
    if a in inv_set and b in inv_set and a < b:
        ia = inv_genes.index(a)
        ib = inv_genes.index(b)
        valid_pairs.append((ia, ib, score))

print(f"  STRING pairs among 195 in-vocab genes: {len(valid_pairs)}", flush=True)

if len(valid_pairs) > 0:
    # Per layer: Spearman rho(STRING_score, L2_dist) and AUROC high vs low pairs
    scores = np.array([p[2] for p in valid_pairs])
    median_score = np.median(scores)
    high_mask = scores >= 0.7
    low_mask = scores < 0.5
    print(f"  Score range: [{scores.min():.3f}, {scores.max():.3f}], median={median_score:.3f}", flush=True)
    print(f"  High (>=0.7): {high_mask.sum()}, Low (<0.5): {low_mask.sum()}", flush=True)

    h02_per_layer = []
    for layer in range(N_LAYERS):
        X = emb_inv[layer]
        dists = []
        for ia, ib, _ in valid_pairs:
            d = np.linalg.norm(X[ia] - X[ib])
            dists.append(d)
        dists = np.array(dists)

        # Spearman: negative expected (higher STRING score → closer in embedding space)
        rho, p = spearmanr(scores, dists)

        # AUROC: high STRING score pairs vs low (test: high pairs have smaller dist)
        if high_mask.sum() > 5 and low_mask.sum() > 5:
            stat, mwp = mannwhitneyu(dists[high_mask], dists[low_mask], alternative='less')
            auroc = stat / (high_mask.sum() * low_mask.sum())
        else:
            auroc, mwp = np.nan, np.nan

        h02_per_layer.append({
            "layer": layer,
            "spearman_rho": float(rho),
            "spearman_p": float(p),
            "auroc_high_vs_low": float(auroc) if not np.isnan(auroc) else None,
            "auroc_p": float(mwp) if not np.isnan(mwp) else None,
            "mean_dist_high": float(dists[high_mask].mean()) if high_mask.sum() > 0 else None,
            "mean_dist_low": float(dists[low_mask].mean()) if low_mask.sum() > 0 else None,
        })
        print(f"  L{layer:02d}: rho={rho:.4f} p={p:.3e} AUROC={auroc:.4f}", flush=True)

    rho_vals = np.array([r["spearman_rho"] for r in h02_per_layer])
    auroc_vals = np.array([r["auroc_high_vs_low"] for r in h02_per_layer if r["auroc_high_vs_low"] is not None])

    # Spearman rho across layers for the STRING→dist correlation
    rho_trend, p_trend = spearmanr(np.arange(N_LAYERS), rho_vals)

    h02_artifact = {
        "n_pairs": len(valid_pairs),
        "n_high": int(high_mask.sum()),
        "n_low": int(low_mask.sum()),
        "per_layer": h02_per_layer,
        "mean_rho": float(np.mean(rho_vals)),
        "layer_trend_rho": float(rho_trend),
        "layer_trend_p": float(p_trend),
        "mean_auroc": float(np.nanmean(auroc_vals)) if len(auroc_vals) > 0 else None,
        "best_layer": int(np.argmin(rho_vals)),  # most negative rho
        "best_rho": float(np.min(rho_vals)),
    }
    with open(ITER_DIR / "h02_string_dist_195.json", "w") as f:
        json.dump(h02_artifact, f, indent=2)
    print("  Saved h02_string_dist_195.json", flush=True)
else:
    print("  WARNING: No STRING pairs found among 195 in-vocab genes!", flush=True)
    h02_artifact = {"error": "no_pairs"}

# ─── H03: Participation ratio + L2 norm distribution per layer ────────────────
print("\n=== H03: Intrinsic dim (PR) + L2 norm distribution ===", flush=True)

h03_per_layer = []
for layer in range(N_LAYERS):
    X = emb_inv[layer]
    # Center
    Xc = X - X.mean(axis=0)
    # SVD
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    # Participation ratio: (sum s^2)^2 / sum(s^4)
    s2 = s**2
    pr = (s2.sum())**2 / (s2**2).sum()
    # Explained variance by PC1
    var_pc1 = s2[0] / s2.sum()
    # L2 norms
    norms = np.linalg.norm(X, axis=1)
    norm_mean = float(norms.mean())
    norm_std = float(norms.std())
    norm_skew = float(skew(norms))
    norm_cv = norm_std / (norm_mean + 1e-12)
    h03_per_layer.append({
        "layer": layer,
        "participation_ratio": float(pr),
        "var_pc1": float(var_pc1),
        "norm_mean": norm_mean,
        "norm_std": norm_std,
        "norm_skew": norm_skew,
        "norm_cv": norm_cv,
    })
    print(f"  L{layer:02d}: PR={pr:.2f}, var_pc1={var_pc1:.3f}, norm_mean={norm_mean:.3f}, cv={norm_cv:.3f}", flush=True)

pr_vals = np.array([r["participation_ratio"] for r in h03_per_layer])
cv_vals = np.array([r["norm_cv"] for r in h03_per_layer])
rho_pr, p_pr = spearmanr(np.arange(N_LAYERS), pr_vals)
rho_cv, p_cv = spearmanr(np.arange(N_LAYERS), cv_vals)
print(f"  Spearman rho(layer, PR) = {rho_pr:.4f}, p = {p_pr:.3e}", flush=True)
print(f"  Spearman rho(layer, norm_cv) = {rho_cv:.4f}, p = {p_cv:.3e}", flush=True)

h03_artifact = {
    "per_layer": h03_per_layer,
    "rho_layer_vs_pr": float(rho_pr),
    "p_layer_vs_pr": float(p_pr),
    "rho_layer_vs_cv": float(rho_cv),
    "p_layer_vs_cv": float(p_cv),
    "n_genes": len(inv_genes),
    "pr_range": [float(pr_vals.min()), float(pr_vals.max())],
    "cv_range": [float(cv_vals.min()), float(cv_vals.max())],
}
with open(ITER_DIR / "h03_intdim_norm_195.json", "w") as f:
    json.dump(h03_artifact, f, indent=2)
print("  Saved h03_intdim_norm_195.json", flush=True)

# ─── Summary ──────────────────────────────────────────────────────────────────
print("\n=== SUMMARY ===", flush=True)
print(f"H01 (spectral gap): rho(layer,gap)={rho_h01:.4f} p={p_h01:.3e}", flush=True)
print(f"     real_mean={np.nanmean(ratios[valid]):.4f} vs null_mean={np.nanmean(null_ratios):.4f}", flush=True)
if len(valid_pairs) > 0:
    print(f"H02 (STRING->dist): mean_rho={h02_artifact['mean_rho']:.4f}, best_rho={h02_artifact['best_rho']:.4f} at L{h02_artifact['best_layer']}", flush=True)
    print(f"     AUROC_mean={h02_artifact['mean_auroc']:.4f}, layer_trend_rho={h02_artifact['layer_trend_rho']:.4f}", flush=True)
print(f"H03 (PR / norm_cv): rho(layer,PR)={rho_pr:.4f} p={p_pr:.3e}, rho(layer,cv)={rho_cv:.4f} p={p_cv:.3e}", flush=True)
