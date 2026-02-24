"""
iter_0029 Multi-Hypothesis Screen

H01 (graph_topology / refinement): Identify 2 kNN connected components at L11,
    test biological enrichment against immune gene families from iter_0026.
    Direct follow-up of iter_0028 H03 (2 components found at every layer).

H02 (manifold_distance / new_method): Gene Trajectory Clustering
    Each named gene has a 12-layer distance-to-centroid profile (trajectory).
    Cluster these 209 trajectories by shape (k-means, k=3..8), test if
    trajectory clusters map to known biology (immune families, cell-type markers,
    STRING hub degree).

H03 (topology_stability / new_method): Spectral Gap k-Robustness Sweep
    Sweep k=5,10,15,20,25,30; for each k compute spectral gap per layer.
    Test if the Spearman rho(layer, gap) < 0 trend is robust across k values.
    Also compare to random-Gaussian null baseline (same n=209 points).
"""

import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings("ignore")

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0029"
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

# ─── Load gene names ──────────────────────────────────────────────────────────
import csv

# Load vocab
with open(CYCLE1 / "gene_list.txt") as f:
    vocab_genes = [line.strip() for line in f if line.strip()]
gene_to_emb_idx = {g: i for i, g in enumerate(vocab_genes)}
print(f"  Total vocab genes: {len(vocab_genes)}", flush=True)

# Load named genes from edge dataset
EDGE_PATH = CYCLE1 / "cycle1_edge_dataset.tsv"
named_gene_set = set()
with open(EDGE_PATH) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        named_gene_set.add(row['source'])
        named_gene_set.add(row['target'])

named_genes = np.array(sorted(g for g in named_gene_set if g in gene_to_emb_idx))
named_idx = np.array([gene_to_emb_idx[g] for g in named_genes])
emb_named = emb[:, named_idx, :]   # [12, N_NAMED, 512]
N_NAMED = len(named_idx)
print(f"  Named genes: {N_NAMED}", flush=True)

# ─── Immune gene families (from iter_0026) ────────────────────────────────────
IMMUNE_FAMILIES = {
    "AP1":      ["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "ATF3"],
    "RUNX":     ["RUNX1", "RUNX2", "RUNX3"],
    "KLF":      ["KLF2", "KLF4", "KLF6"],
    "BCL2fam":  ["BCL2", "BCL2L1", "BCL2L11", "MCL1", "BAX"],
    "HLA_I":    ["HLA-A", "HLA-B", "HLA-C", "B2M"],
    "HLA_II":   ["HLA-DRA", "HLA-DRB1", "HLA-DMA", "HLA-DMB", "HLA-DPA1",
                  "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74"],
    "CCL":      ["CCL2", "CCL3", "CCL4", "CCL5", "CCL20", "CXCL8", "CXCL10"],
    "TNFSF":    ["TNFSF10", "TNFSF11", "TNFSF13B", "FASLG"],
    "IL2_path": ["IL2RA", "IL2RB", "IL2RG", "JAK1", "JAK3", "STAT5A", "STAT5B"],
}

# Cell-type markers
TCELL = ["PRF1", "LCK", "IFNG", "CD8A", "JUN", "JUNB", "GZMB", "ZAP70"]
BCELL = ["CD79A", "CD79B", "MS4A1", "PAX5", "IGHM"]
FIBRO = ["COL1A1", "VIM", "ACTA2"]
EPITHELIAL = ["EPCAM", "KRT19"]

gene2family = {}
for fam, genes in IMMUNE_FAMILIES.items():
    for g in genes:
        gene2family[g] = fam

gene2celltype = {}
for g in TCELL:       gene2celltype[g] = "T_cell"
for g in BCELL:       gene2celltype[g] = "B_cell"
for g in FIBRO:       gene2celltype[g] = "fibroblast"
for g in EPITHELIAL:  gene2celltype[g] = "epithelial"

# ─── STRING degrees (from iter_0015 cache) ────────────────────────────────────
print("Loading STRING degrees...", flush=True)
string_data = json.load(open(STRING_CACHE))
string_degree = defaultdict(int)
for pair in string_data["pairs"]:
    g1, g2 = pair["g1"], pair["g2"]
    string_degree[g1] += 1
    string_degree[g2] += 1

gene2degree = {g: string_degree.get(g, 0) for g in named_genes}
named_degrees = np.array([gene2degree[g] for g in named_genes])
print(f"  Genes with STRING degree > 0: {(named_degrees > 0).sum()}", flush=True)


# ─── Helper: build kNN graph adjacency ────────────────────────────────────────
def build_knn_adj(X, k):
    """Return symmetric kNN adjacency matrix (scipy csr)."""
    n = len(X)
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(X)
    dist_mat, idx_mat = nbrs.kneighbors(X)
    rows, cols = [], []
    for i in range(n):
        for j in idx_mat[i, 1:]:   # skip self
            rows.append(i); cols.append(j)
            rows.append(j); cols.append(i)
    data = np.ones(len(rows))
    A = csr_matrix((data, (rows, cols)), shape=(n, n))
    A.data[:] = 1.0
    return A


def laplacian_spectrum(A, n_eigs=30):
    """Normalized Laplacian eigenvalues (ascending)."""
    n = A.shape[0]
    degrees = np.array(A.sum(axis=1)).flatten()
    d_inv_sqrt = np.where(degrees > 0, 1.0 / np.sqrt(degrees), 0.0)
    # D^{-1/2} A D^{-1/2}
    from scipy.sparse import diags
    D_inv_sqrt = diags(d_inv_sqrt)
    L_sym = csr_matrix(np.eye(n)) - D_inv_sqrt @ A @ D_inv_sqrt
    n_eigs = min(n_eigs, n - 2)
    try:
        vals = eigsh(L_sym, k=n_eigs, which="SM", return_eigenvectors=False, tol=1e-6)
        vals = np.sort(np.real(vals))
    except Exception:
        vals = np.sort(np.linalg.eigvalsh(L_sym.toarray()))[:n_eigs]
    return vals


# ══════════════════════════════════════════════════════════════════════════════
# H01: Connected components at L11, biological enrichment
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== H01: Connected Components Biological Enrichment ===", flush=True)
K_DEFAULT = 10
L_DEEP = 11   # last layer

results_h01 = {}

for layer in [11, 10, 9, 8]:
    X = emb_named[layer]
    A = build_knn_adj(X, K_DEFAULT)
    n_comps, labels = connected_components(A, directed=False)
    comp_sizes = Counter(labels)
    results_h01[f"L{layer}"] = {
        "n_components": int(n_comps),
        "component_sizes": {int(k): int(v) for k, v in sorted(comp_sizes.items())},
    }
    print(f"  L{layer}: {n_comps} components, sizes = {sorted(comp_sizes.values(), reverse=True)}", flush=True)

# Deep dive at L11
layer = 11
X = emb_named[layer]
A = build_knn_adj(X, K_DEFAULT)
n_comps, labels = connected_components(A, directed=False)
print(f"\n  L11 deep dive: {n_comps} components", flush=True)

comp_gene_lists = defaultdict(list)
for gene_i, comp_id in enumerate(labels):
    comp_gene_lists[comp_id].append(named_genes[gene_i])

# For each component, count family memberships and cell-type memberships
enrichment = {}
for comp_id, genes in sorted(comp_gene_lists.items(), key=lambda x: -len(x[1])):
    fam_counts = Counter(gene2family.get(g) for g in genes if g in gene2family)
    ct_counts  = Counter(gene2celltype.get(g) for g in genes if g in gene2celltype)
    enrichment[f"comp_{int(comp_id)}"] = {
        "size": len(genes),
        "genes": sorted(genes),
        "immune_family_counts": {k: v for k, v in fam_counts.items() if k is not None},
        "celltype_counts": {k: v for k, v in ct_counts.items() if k is not None},
        "mean_string_degree": float(np.mean([gene2degree.get(g, 0) for g in genes])),
    }
    print(f"    Comp {comp_id} (n={len(genes)}): fams={dict(fam_counts)}", flush=True)

results_h01["L11_component_enrichment"] = enrichment

# Fisher exact test: is one component enriched for any family?
from scipy.stats import fisher_exact
enrichment_tests = []
if n_comps == 2:
    comp0_genes = set(comp_gene_lists[0])
    comp1_genes = set(comp_gene_lists[1])
    for fam, fam_genes in IMMUNE_FAMILIES.items():
        fam_set = set(fam_genes)
        a = len(fam_set & comp0_genes)   # fam AND comp0
        b = len(comp0_genes) - a          # not-fam AND comp0
        c = len(fam_set & comp1_genes)   # fam AND comp1
        d = len(comp1_genes) - c          # not-fam AND comp1
        if a + c == 0:
            continue
        odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        enrichment_tests.append({
            "family": fam,
            "comp0_count": a, "comp0_total": len(comp0_genes),
            "comp1_count": c, "comp1_total": len(comp1_genes),
            "odds_ratio": float(odds), "p_value": float(p),
        })
    enrichment_tests.sort(key=lambda x: x["p_value"])
    results_h01["fisher_enrichment_tests"] = enrichment_tests
    print("\n  Fisher enrichment (top 5):", flush=True)
    for t in enrichment_tests[:5]:
        print(f"    {t['family']}: comp0={t['comp0_count']}/{t['comp0_total']}, "
              f"comp1={t['comp1_count']}/{t['comp1_total']}, "
              f"OR={t['odds_ratio']:.2f}, p={t['p_value']:.4f}", flush=True)

# Also test Mann-Whitney for STRING degree between components
if n_comps == 2:
    deg0 = [gene2degree.get(g, 0) for g in comp_gene_lists[0]]
    deg1 = [gene2degree.get(g, 0) for g in comp_gene_lists[1]]
    stat, p_deg = mannwhitneyu(deg0, deg1, alternative="two-sided")
    auroc = stat / (len(deg0) * len(deg1))
    results_h01["string_degree_MW"] = {
        "mean_deg_comp0": float(np.mean(deg0)),
        "mean_deg_comp1": float(np.mean(deg1)),
        "auroc": float(auroc),
        "p_value": float(p_deg),
    }
    print(f"\n  STRING degree MW test: comp0_mean={np.mean(deg0):.1f}, "
          f"comp1_mean={np.mean(deg1):.1f}, AUROC={auroc:.3f}, p={p_deg:.4f}", flush=True)

out_h01 = ITER_DIR / "h01_component_enrichment.json"
with open(out_h01, "w") as f:
    json.dump(results_h01, f, indent=2)
print(f"\n  Saved: {out_h01}", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
# H02: Gene Trajectory Clustering (distance-to-centroid profiles)
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== H02: Gene Trajectory Clustering ===", flush=True)

# Compute per-gene, per-layer distance to centroid
traj = np.zeros((N_NAMED, N_LAYERS))
for layer_i in range(N_LAYERS):
    X = emb_named[layer_i]
    centroid = X.mean(axis=0)
    traj[:, layer_i] = np.linalg.norm(X - centroid, axis=1)

print(f"  Trajectory matrix shape: {traj.shape}  (genes x layers)", flush=True)

# Normalize trajectories: standardize per gene (z-score across layers)
traj_norm = (traj - traj.mean(axis=1, keepdims=True)) / (traj.std(axis=1, keepdims=True) + 1e-9)

# k-means clustering for k=2..6
best_k = None
best_sil = -999
sil_scores = {}
from sklearn.metrics import silhouette_score

for k in range(2, 7):
    km = KMeans(n_clusters=k, random_state=42, n_init=20)
    lbls = km.fit_predict(traj_norm)
    if len(set(lbls)) > 1:
        sil = silhouette_score(traj_norm, lbls)
    else:
        sil = -1
    sil_scores[k] = float(sil)
    if sil > best_sil:
        best_sil = sil
        best_k = k
    print(f"  k={k}: silhouette={sil:.4f}", flush=True)

print(f"  Best k={best_k} (sil={best_sil:.4f})", flush=True)

km_best = KMeans(n_clusters=best_k, random_state=42, n_init=20)
cluster_labels = km_best.fit_predict(traj_norm)

# Characterize clusters
results_h02 = {
    "n_genes": int(N_NAMED),
    "n_layers": int(N_LAYERS),
    "silhouette_by_k": sil_scores,
    "best_k": int(best_k),
    "best_silhouette": float(best_sil),
    "clusters": [],
}

for c in range(best_k):
    mask = cluster_labels == c
    genes_in_c = named_genes[mask]
    fam_counts = Counter(gene2family.get(g) for g in genes_in_c if g in gene2family)
    ct_counts  = Counter(gene2celltype.get(g) for g in genes_in_c if g in gene2celltype)
    mean_deg   = float(np.mean([gene2degree.get(g, 0) for g in genes_in_c]))
    mean_traj  = traj[mask].mean(axis=0).tolist()

    results_h02["clusters"].append({
        "cluster_id": int(c),
        "size": int(mask.sum()),
        "genes": sorted(genes_in_c.tolist()),
        "immune_family_counts": {k: v for k, v in fam_counts.items() if k is not None},
        "celltype_counts": {k: v for k, v in ct_counts.items() if k is not None},
        "mean_string_degree": mean_deg,
        "mean_trajectory": mean_traj,
        "mean_traj_slope": float(np.polyfit(np.arange(N_LAYERS), traj[mask].mean(axis=0), 1)[0]),
    })
    print(f"  Cluster {c} (n={mask.sum()}): deg={mean_deg:.1f}, "
          f"fams={dict(fam_counts)}, cts={dict(ct_counts)}", flush=True)

# Kruskal-Wallis test: do clusters differ in STRING degree?
from scipy.stats import kruskal
cluster_groups = [[gene2degree.get(g, 0) for g in c_info["genes"]]
                   for c_info in results_h02["clusters"]]
if best_k > 1:
    kw_stat, kw_p = kruskal(*cluster_groups)
    results_h02["kruskal_wallis_degree"] = {
        "statistic": float(kw_stat), "p_value": float(kw_p)
    }
    print(f"  Kruskal-Wallis (degree across clusters): H={kw_stat:.3f}, p={kw_p:.4f}", flush=True)

# Chi-square: are immune family memberships non-uniform across clusters?
from scipy.stats import chi2_contingency
fam_by_cluster = np.zeros((best_k, len(IMMUNE_FAMILIES)), dtype=int)
for c_i, c_info in enumerate(results_h02["clusters"]):
    for f_i, fam_name in enumerate(IMMUNE_FAMILIES.keys()):
        fam_by_cluster[c_i, f_i] = c_info["immune_family_counts"].get(fam_name, 0)

row_sums = fam_by_cluster.sum(axis=1)
col_sums = fam_by_cluster.sum(axis=0)
if fam_by_cluster.sum() > 0 and (col_sums > 0).sum() > 1 and (row_sums > 0).sum() > 1:
    chi2, p_chi2, dof, expected = chi2_contingency(fam_by_cluster[:, col_sums > 0])
    results_h02["chi2_family_by_cluster"] = {
        "chi2": float(chi2), "p_value": float(p_chi2), "dof": int(dof)
    }
    print(f"  Chi2 (immune families across clusters): chi2={chi2:.3f}, p={p_chi2:.4f}", flush=True)

out_h02 = ITER_DIR / "h02_trajectory_clusters.json"
with open(out_h02, "w") as f:
    json.dump(results_h02, f, indent=2)
print(f"\n  Saved: {out_h02}", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
# H03: Spectral Gap k-Robustness Sweep
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== H03: Spectral Gap k-Robustness ===", flush=True)

K_VALUES = [5, 10, 15, 20, 25, 30]
N_NULL = 50   # random Gaussian nulls for comparison

results_h03 = {
    "k_values": K_VALUES,
    "n_null_runs": N_NULL,
    "layers": list(range(N_LAYERS)),
    "by_k": {},
}

for k in K_VALUES:
    print(f"  k={k}...", flush=True, end="")
    gaps_real = []
    for layer_i in range(N_LAYERS):
        X = emb_named[layer_i]
        A = build_knn_adj(X, k)
        n_comps, _ = connected_components(A, directed=False)
        vals = laplacian_spectrum(A, n_eigs=min(20, N_NAMED - 2))
        # spectral gap = smallest nonzero eigenvalue (Fiedler)
        nonzero = vals[vals > 1e-8]
        gap = float(nonzero[0]) if len(nonzero) > 0 else 0.0
        gaps_real.append(gap)

    rho, p_rho = spearmanr(np.arange(N_LAYERS), gaps_real)

    # Null: random Gaussian data, same shape
    null_rhos = []
    for _ in range(N_NULL):
        X_null = rng.standard_normal((N_NAMED, N_DIM))
        gaps_null_run = []
        for layer_i in range(N_LAYERS):
            A_n = build_knn_adj(X_null, k)
            vals_n = laplacian_spectrum(A_n, n_eigs=min(20, N_NAMED - 2))
            nz = vals_n[vals_n > 1e-8]
            gaps_null_run.append(float(nz[0]) if len(nz) > 0 else 0.0)
        rho_n, _ = spearmanr(np.arange(N_LAYERS), gaps_null_run)
        null_rhos.append(rho_n)
    null_p = float(np.mean(np.array(null_rhos) <= rho))

    results_h03["by_k"][str(k)] = {
        "gaps": [float(g) for g in gaps_real],
        "spearman_rho": float(rho),
        "spearman_p": float(p_rho),
        "null_rho_mean": float(np.mean(null_rhos)),
        "null_rho_std": float(np.std(null_rhos)),
        "null_p": null_p,
    }
    print(f" rho={rho:.4f} (p={p_rho:.2e}), null_p={null_p:.3f}", flush=True)

# Summary: is the trend consistent across k values?
rhos = [results_h03["by_k"][str(k)]["spearman_rho"] for k in K_VALUES]
null_ps = [results_h03["by_k"][str(k)]["null_p"] for k in K_VALUES]
results_h03["summary"] = {
    "rho_mean": float(np.mean(rhos)),
    "rho_min": float(np.min(rhos)),
    "rho_max": float(np.max(rhos)),
    "all_negative": bool(all(r < 0 for r in rhos)),
    "all_null_p_lt_01": bool(all(p < 0.1 for p in null_ps)),
}
print(f"\n  k-sweep rho summary: mean={np.mean(rhos):.4f}, "
      f"range=[{np.min(rhos):.4f}, {np.max(rhos):.4f}], "
      f"all_negative={results_h03['summary']['all_negative']}", flush=True)

out_h03 = ITER_DIR / "h03_spectral_gap_krobustness.json"
with open(out_h03, "w") as f:
    json.dump(results_h03, f, indent=2)
print(f"\n  Saved: {out_h03}", flush=True)

print("\nAll experiments complete.", flush=True)
