"""
iter_0035 Multi-Hypothesis Screen

H01 (module_structure / new_method): TRRUST + GO Jaccard co-annotation vs embedding distance
    On 195 in-vocab genes, build TRRUST co-regulation pairs and GO BP co-annotation pairs.
    Test whether Jaccard similarity (shared GO terms or TRRUST co-targets) predicts
    embedding proximity (lower L2 distance) at L11.
    Mann-Whitney AUROC: high-Jaccard pairs vs random pairs.

H02 (manifold_distance / new_method): B-cell conditional nearest-neighbor profile
    At L11: for each B-cell marker, find its k=10 nearest neighbors in embedding space.
    Compute fraction of neighbors that are also B-cell markers (precision@10).
    Compare to null: random sets of same size, bootstrap (n=1000).
    Also do same at L2 to check for layer emergence.

H03 (intrinsic_dimensionality / new_method): B-cell vs non-B-cell centroid separation trajectory
    Per layer: compute centroid of B-cell markers (n=~14), centroid of non-B-cell genes.
    Measure centroid-centroid L2 distance per layer.
    Test monotonic increase (Spearman rho across layers).
    Normalize by mean pairwise distance to control for shrinkage.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
import warnings
warnings.filterwarnings("ignore")

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

def jdump(obj, path, **kw):
    with open(path, "w") as f:
        json.dump(obj, f, cls=NpEncoder, indent=2, **kw)

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0035"
ITER_DIR.mkdir(parents=True, exist_ok=True)

CYCLE1 = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
              "/subproject_38_geometric_residual_stream_interpretability"
              "/implementation/outputs/cycle1_main")

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

# Filter to in-vocab, drop OOV
invocab_named = sorted([g for g in named_genes_set if g in gene_to_emb_idx])
L0_NORMS = np.linalg.norm(emb[0, [gene_to_emb_idx[g] for g in invocab_named], :], axis=1)
invocab_named = [g for g, n in zip(invocab_named, L0_NORMS) if n > 1e-8]
print(f"  In-vocab OOV-filtered: {len(invocab_named)} genes", flush=True)

# Gene index arrays
gene_idx = np.array([gene_to_emb_idx[g] for g in invocab_named])
N = len(invocab_named)
gene_to_local = {g: i for i, g in enumerate(invocab_named)}

# B-cell markers (from prior iterations)
BCELL_MARKERS = ["MS4A1", "CD19", "CD79A", "CD79B", "PAX5", "BLK", "BANK1",
                 "FCRL1", "CD22", "FCER2", "IGHM", "IGHG1", "IGHA1", "CR2"]
bcell_in_vocab = [g for g in BCELL_MARKERS if g in gene_to_local]
print(f"  B-cell markers in vocab: {len(bcell_in_vocab)} {bcell_in_vocab}", flush=True)

# Embeddings per layer for in-vocab genes: [12, N, 512]
E = emb[:, gene_idx, :]  # [12, N, 512]

# ──────────────────────────────────────────────────────────────────────────────
# H01: TRRUST + GO Jaccard co-annotation vs embedding distance
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H01: TRRUST + GO Jaccard vs embedding distance ===", flush=True)

# TRRUST: use edge dataset — activation/repression pairs as co-regulation signal
# Build TRRUST co-target pairs: genes that share a TF regulator
trrust_pairs_coregulated = set()
tf_targets = {}  # TF -> set of targets
with open(EDGES_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        src, tgt = row["source"], row["target"]
        if src not in gene_to_local or tgt not in gene_to_local:
            continue
        # Assume src=TF, tgt=target gene (based on iter_0027 context)
        tf_targets.setdefault(src, set()).add(tgt)

# Build co-regulated pairs (share at least one TF regulator)
print(f"  TFs with in-vocab targets: {len(tf_targets)}", flush=True)
coregulated = set()
for tf, targets in tf_targets.items():
    targets_list = sorted(targets)
    for i in range(len(targets_list)):
        for j in range(i+1, len(targets_list)):
            coregulated.add((targets_list[i], targets_list[j]))
print(f"  Co-regulated pairs: {len(coregulated)}", flush=True)

# GO BP co-annotation from mygene (use curated gene sets from prior work)
# Build manually from known biology: use immune function gene sets
GO_SETS = {
    "bcell_activation": ["MS4A1", "CD19", "CD79A", "CD79B", "PAX5", "BLK", "BANK1",
                         "FCRL1", "CD22", "FCER2", "IGHM", "IGHG1", "IGHA1", "CR2"],
    "tcell_activation": ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "LCK",
                         "ZAP70", "LAT", "GATA3", "TBX21", "RORC", "FOXP3"],
    "myeloid": ["CD14", "CSF1R", "CD33", "ITGAM", "SPI1", "IRF8", "CEBPA"],
    "interferon": ["STAT1", "STAT2", "IRF3", "IRF7", "MX1", "ISG15", "IFIT1", "IFIT2"],
    "nfkb": ["NFKB1", "NFKB2", "RELA", "RELB", "REL", "IKBKA", "IKBKB"],
    "apoptosis": ["BCL2", "BCL2L1", "BAX", "BAD", "CASP3", "CASP8", "TP53", "MCL1"],
    "cell_cycle": ["CCND1", "CCNE1", "CDK2", "CDK4", "CDK6", "RB1", "E2F1"],
    "tgfb": ["TGFB1", "SMAD2", "SMAD3", "SMAD4", "SMAD7"],
    "wnt": ["CTNNB1", "APC", "AXIN1", "TCF4", "LEF1"],
    "jak_stat": ["JAK1", "JAK2", "STAT3", "STAT5A", "STAT5B", "SOCS1", "SOCS3"],
}

# Build GO co-annotation pairs: share at least one GO set
go_coannot = set()
for gs_name, gs_genes in GO_SETS.items():
    gs_in_vocab = [g for g in gs_genes if g in gene_to_local]
    for i in range(len(gs_in_vocab)):
        for j in range(i+1, len(gs_in_vocab)):
            a, b = sorted([gs_in_vocab[i], gs_in_vocab[j]])
            go_coannot.add((a, b))

print(f"  GO co-annotated pairs: {len(go_coannot)}", flush=True)

# All possible pairs from in-vocab genes
all_pairs = []
for i in range(N):
    for j in range(i+1, N):
        all_pairs.append((invocab_named[i], invocab_named[j]))

# For efficiency, sample random pairs for null
rng2 = np.random.default_rng(42)
random_pair_idx = rng2.choice(len(all_pairs), size=min(5000, len(all_pairs)), replace=False)
random_pairs_sample = [all_pairs[k] for k in random_pair_idx]

# Compute L2 distances for co-regulated pairs at L11
layer = 11
E_L11 = E[layer]  # [N, 512]

def pair_dist(pairs, E_local, g2l):
    dists = []
    for a, b in pairs:
        ia, ib = g2l.get(a), g2l.get(b)
        if ia is None or ib is None:
            continue
        dists.append(np.linalg.norm(E_local[ia] - E_local[ib]))
    return np.array(dists)

d_coreg = pair_dist(list(coregulated), E_L11, gene_to_local)
d_go = pair_dist(list(go_coannot), E_L11, gene_to_local)
d_random = pair_dist(random_pairs_sample, E_L11, gene_to_local)

print(f"  Co-reg pairs with distances: {len(d_coreg)}, mean={d_coreg.mean():.3f}")
print(f"  GO pairs with distances: {len(d_go)}, mean={d_go.mean():.3f}")
print(f"  Random pairs: {len(d_random)}, mean={d_random.mean():.3f}")

# AUROC: fraction of annotated pairs that are closer than random (one-sided)
# lower distance = better co-localization → test: annotated < random
def auroc_lower(focal, background):
    """AUROC for focal distances being LOWER than background."""
    if len(focal) == 0 or len(background) == 0:
        return 0.5, 1.0
    stat, p = mannwhitneyu(focal, background, alternative="less")
    auroc = stat / (len(focal) * len(background))
    return float(auroc), float(p)

auroc_coreg, p_coreg = auroc_lower(d_coreg, d_random)
auroc_go, p_go = auroc_lower(d_go, d_random)

print(f"  Co-reg AUROC(less)={auroc_coreg:.3f} p={p_coreg:.4f}")
print(f"  GO AUROC(less)={auroc_go:.3f} p={p_go:.4f}")

# Also compute per layer for co-reg to test emergence
coreg_aurocs = []
go_aurocs = []
for lyr in range(N_LAYERS):
    E_lyr = E[lyr]
    d_c = pair_dist(list(coregulated), E_lyr, gene_to_local)
    d_g = pair_dist(list(go_coannot), E_lyr, gene_to_local)
    d_r = pair_dist(random_pairs_sample, E_lyr, gene_to_local)
    a_c, _ = auroc_lower(d_c, d_r)
    a_g, _ = auroc_lower(d_g, d_r)
    coreg_aurocs.append(a_c)
    go_aurocs.append(a_g)

rho_coreg, p_rho_coreg = spearmanr(range(N_LAYERS), coreg_aurocs)
rho_go, p_rho_go = spearmanr(range(N_LAYERS), go_aurocs)
print(f"  Co-reg Spearman rho(layer, AUROC)={rho_coreg:.3f} p={p_rho_coreg:.4f}")
print(f"  GO Spearman rho(layer, AUROC)={rho_go:.3f} p={p_rho_go:.4f}")

h01_result = {
    "hypothesis": "trrust_go_jaccard_vs_embedding_distance",
    "layer_tested": 11,
    "n_coregulated_pairs": int(len(d_coreg)),
    "n_go_pairs": int(len(d_go)),
    "n_random_pairs": int(len(d_random)),
    "coregulated_mean_dist_L11": float(d_coreg.mean()) if len(d_coreg) > 0 else None,
    "go_mean_dist_L11": float(d_go.mean()) if len(d_go) > 0 else None,
    "random_mean_dist_L11": float(d_random.mean()),
    "coregulated_auroc_L11": auroc_coreg,
    "coregulated_p_L11": p_coreg,
    "go_auroc_L11": auroc_go,
    "go_p_L11": p_go,
    "coreg_auroc_per_layer": coreg_aurocs,
    "go_auroc_per_layer": go_aurocs,
    "spearman_rho_coreg_layer_auroc": float(rho_coreg),
    "spearman_p_coreg": float(p_rho_coreg),
    "spearman_rho_go_layer_auroc": float(rho_go),
    "spearman_p_go": float(p_rho_go),
}
jdump(h01_result, ITER_DIR / "h01_trrust_go_jaccard.json")
print("  H01 saved.", flush=True)

# ──────────────────────────────────────────────────────────────────────────────
# H02: B-cell conditional nearest-neighbor profile
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H02: B-cell conditional nearest-neighbor profile ===", flush=True)

bcell_local_idx = np.array([gene_to_local[g] for g in bcell_in_vocab])
n_bcell = len(bcell_local_idx)

def precision_at_k_bcell(E_lyr, bcell_idx, k=10):
    """For each B-cell gene, find k nearest neighbors, measure fraction that are B-cell."""
    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=k+1, metric="euclidean", algorithm="ball_tree")
    nbrs.fit(E_lyr)
    # Query only B-cell genes
    E_bcell = E_lyr[bcell_idx]
    _, indices = nbrs.kneighbors(E_bcell)
    # Exclude self (first neighbor)
    bcell_set = set(bcell_idx.tolist())
    precisions = []
    for row in indices:
        neighbors = row[1:]  # exclude self
        frac = sum(1 for nb in neighbors if nb in bcell_set) / k
        precisions.append(frac)
    return np.array(precisions)

# Bootstrap null: random gene sets of same size
def bootstrap_null_precision(E_lyr, n_bcell_local, k=10, n_boot=500, rng_seed=42):
    from sklearn.neighbors import NearestNeighbors
    rng_b = np.random.default_rng(rng_seed)
    nbrs = NearestNeighbors(n_neighbors=k+1, metric="euclidean", algorithm="ball_tree")
    nbrs.fit(E_lyr)
    N_local = E_lyr.shape[0]
    null_means = []
    for _ in range(n_boot):
        rand_idx = rng_b.choice(N_local, size=n_bcell_local, replace=False)
        rand_set = set(rand_idx.tolist())
        E_rand = E_lyr[rand_idx]
        _, indices = nbrs.kneighbors(E_rand)
        precs = []
        for row in indices:
            neighbors = row[1:]
            frac = sum(1 for nb in neighbors if nb in rand_set) / k
            precs.append(frac)
        null_means.append(np.mean(precs))
    return np.array(null_means)

K = 10
h02_per_layer = []
for lyr in [2, 5, 8, 11]:  # key layers
    E_lyr = E[lyr]
    prec = precision_at_k_bcell(E_lyr, bcell_local_idx, k=K)
    obs_mean = float(prec.mean())
    null_dist = bootstrap_null_precision(E_lyr, n_bcell, k=K, n_boot=500, rng_seed=42+lyr)
    null_mean = float(null_dist.mean())
    null_std = float(null_dist.std())
    z = (obs_mean - null_mean) / (null_std + 1e-10)
    emp_p = float(np.mean(null_dist >= obs_mean))
    print(f"  L{lyr:02d}: obs_precision@10={obs_mean:.3f}, null_mean={null_mean:.3f}, z={z:.2f}, emp_p={emp_p:.3f}")
    h02_per_layer.append({
        "layer": lyr,
        "obs_precision_at_10": obs_mean,
        "null_mean": null_mean,
        "null_std": null_std,
        "z_score": z,
        "empirical_p": emp_p,
        "per_gene_precisions": prec.tolist()
    })

h02_result = {
    "hypothesis": "bcell_knn_precision_at_10",
    "k": K,
    "n_bcell_markers": n_bcell,
    "bcell_markers": bcell_in_vocab,
    "n_bootstrap": 500,
    "results_per_layer": h02_per_layer,
    "interpretation": "Higher precision@k for B-cell markers vs null = geometric clustering"
}
jdump(h02_result, ITER_DIR / "h02_bcell_knn_precision.json")
print("  H02 saved.", flush=True)

# ──────────────────────────────────────────────────────────────────────────────
# H03: B-cell vs non-B-cell centroid separation trajectory
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H03: B-cell centroid separation trajectory ===", flush=True)

non_bcell_local_idx = np.array([i for i in range(N) if i not in set(bcell_local_idx.tolist())])

centroid_dists = []
normalized_dists = []
for lyr in range(N_LAYERS):
    E_lyr = E[lyr]
    bcell_centroid = E_lyr[bcell_local_idx].mean(axis=0)
    nonbcell_centroid = E_lyr[non_bcell_local_idx].mean(axis=0)
    dist = float(np.linalg.norm(bcell_centroid - nonbcell_centroid))
    # Normalize by mean pairwise distance (sample 2000 random pairs)
    n_pairs = 2000
    idx_a = rng.integers(0, N, size=n_pairs)
    idx_b = rng.integers(0, N, size=n_pairs)
    mean_pairwise = float(np.mean([np.linalg.norm(E_lyr[a] - E_lyr[b])
                                    for a, b in zip(idx_a, idx_b)]))
    norm_dist = dist / (mean_pairwise + 1e-10)
    centroid_dists.append(dist)
    normalized_dists.append(norm_dist)
    print(f"  L{lyr:02d}: centroid_dist={dist:.3f}, mean_pairwise={mean_pairwise:.3f}, norm={norm_dist:.4f}")

rho_raw, p_rho_raw = spearmanr(range(N_LAYERS), centroid_dists)
rho_norm, p_rho_norm = spearmanr(range(N_LAYERS), normalized_dists)
print(f"  Spearman rho(layer, raw_dist)={rho_raw:.3f} p={p_rho_raw:.4f}")
print(f"  Spearman rho(layer, norm_dist)={rho_norm:.3f} p={p_rho_norm:.4f}")

# Bootstrap null: random sets of same size as B-cell markers
n_boot_h03 = 1000
null_rhos = []
for _ in range(n_boot_h03):
    rand_idx = rng.choice(N, size=n_bcell, replace=False)
    rand_dists = []
    for lyr in range(N_LAYERS):
        E_lyr = E[lyr]
        c1 = E_lyr[rand_idx].mean(axis=0)
        c2_idx = np.array([i for i in range(N) if i not in set(rand_idx.tolist())])
        c2 = E_lyr[c2_idx].mean(axis=0)
        rand_dists.append(np.linalg.norm(c1 - c2))
    rho_r, _ = spearmanr(range(N_LAYERS), rand_dists)
    null_rhos.append(rho_r)

null_rhos = np.array(null_rhos)
z_rho = (rho_raw - null_rhos.mean()) / (null_rhos.std() + 1e-10)
emp_p_rho = float(np.mean(null_rhos >= rho_raw))
print(f"  Bootstrap null rho: mean={null_rhos.mean():.3f}, std={null_rhos.std():.3f}")
print(f"  Observed rho={rho_raw:.3f}, z={z_rho:.2f}, emp_p={emp_p_rho:.3f}")

h03_result = {
    "hypothesis": "bcell_centroid_separation_trajectory",
    "n_bcell": n_bcell,
    "n_non_bcell": int(len(non_bcell_local_idx)),
    "centroid_dists_per_layer": centroid_dists,
    "normalized_dists_per_layer": normalized_dists,
    "spearman_rho_raw": float(rho_raw),
    "spearman_p_raw": float(p_rho_raw),
    "spearman_rho_norm": float(rho_norm),
    "spearman_p_norm": float(p_rho_norm),
    "bootstrap_null_rho_mean": float(null_rhos.mean()),
    "bootstrap_null_rho_std": float(null_rhos.std()),
    "bootstrap_z_score": float(z_rho),
    "bootstrap_empirical_p": emp_p_rho,
    "n_bootstrap": n_boot_h03,
    "interpretation": "If rho > 0 and significant vs null: B-cell centroid progressively separates"
}
jdump(h03_result, ITER_DIR / "h03_bcell_centroid_trajectory.json")
print("  H03 saved.", flush=True)

print("\n=== All experiments complete ===")
print(f"Outputs in: {ITER_DIR}")
