"""
iter_0036 Multi-Hypothesis Screen

H01 (manifold_distance / new_method): Multi-cell-type kNN precision@10 panel
    Expand B-cell markers (5 → as many in-vocab as possible).
    Test B-cell, T-cell, Myeloid, NK cell types at L2 and L11.
    Demonstrates specificity: does kNN clustering occur for ALL cell types
    or selectively for B-cells? Also tests whether z-scores are highest for B-cell.

H02 (null_sensitivity / new_method): Gene-name permutation null
    Shuffle gene↔embedding assignment (random permutation of which embedding
    row belongs to which gene name), then recompute B-cell precision@10.
    Under null: B-cell genes get random embeddings → precision should drop to ~null.
    Demonstrates result is not an artifact of tokenizer/position effects.

H03 (module_structure / new_method): B-cell neighborhood functional characterization
    At L2 (highest z-score layer): for all B-cell markers, collect ALL k=20 nearest
    neighbors. Determine what fraction are: other B-cell markers, TFs, cytokines,
    immune genes vs random expectation. Produces a ranked neighbor list.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.neighbors import NearestNeighbors
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
ITER_DIR = PROJECT / "iterations" / "iter_0036"
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

gene_idx = np.array([gene_to_emb_idx[g] for g in invocab_named])
N = len(invocab_named)
gene_to_local = {g: i for i, g in enumerate(invocab_named)}

# Embeddings per layer for in-vocab genes: [12, N, 512]
E = emb[:, gene_idx, :]  # [12, N, 512]

# ─── Cell-type marker sets (expanded) ─────────────────────────────────────────
# B-cell markers — expanded from literature
BCELL_ALL = [
    "MS4A1", "CD19", "CD79A", "CD79B", "PAX5", "BLK", "BANK1",
    "FCRL1", "CD22", "FCER2", "IGHM", "IGHG1", "IGHA1", "CR2",
    "BLNK", "CD27", "IRF4", "PRDM1", "XBP1", "EBF1",
    "VPREB1", "VPREB3", "IGHD", "IGKC", "IGLC1",
]

# T-cell markers
TCELL_ALL = [
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "LCK",
    "ZAP70", "LAT", "ITK", "GATA3", "TBX21", "RORC", "FOXP3",
    "TRAC", "TRBC1", "TRBC2", "CD5", "CD28", "CD69",
    "IL2RA", "IL7R", "CCR7", "SELL",
]

# Myeloid markers
MYELOID_ALL = [
    "CD14", "CSF1R", "CD33", "ITGAM", "SPI1", "IRF8", "CEBPA",
    "FCGR3A", "FCGR1A", "LYZ", "S100A8", "S100A9", "VCAN",
    "CD68", "MSR1", "MRC1", "VSIG4", "C1QA", "C1QB", "C1QC",
]

# NK cell markers
NK_ALL = [
    "NCAM1", "KLRB1", "KLRD1", "KLRK1", "NKG7", "GNLY",
    "PRF1", "GZMB", "GZMK", "FCGR3A", "NCR1", "NCR3",
    "KIR2DL1", "KIR2DL3", "KIR3DL1",
]

CELL_TYPES = {
    "B-cell": BCELL_ALL,
    "T-cell": TCELL_ALL,
    "Myeloid": MYELOID_ALL,
    "NK": NK_ALL,
}

ct_in_vocab = {}
for ct, markers in CELL_TYPES.items():
    in_v = [g for g in markers if g in gene_to_local]
    ct_in_vocab[ct] = in_v
    print(f"  {ct}: {len(in_v)} in-vocab: {in_v}", flush=True)


# ─── Shared utilities ──────────────────────────────────────────────────────────
def precision_at_k(E_lyr, focal_idx, all_idx_set, k=10):
    """Precision@k: for each focal gene, fraction of k-NN that are also focal."""
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", algorithm="ball_tree")
    nbrs.fit(E_lyr)
    E_focal = E_lyr[focal_idx]
    _, indices = nbrs.kneighbors(E_focal)
    precs = []
    for row in indices:
        neighbors = row[1:]  # exclude self
        frac = sum(1 for nb in neighbors if nb in all_idx_set) / k
        precs.append(frac)
    return np.array(precs)


def bootstrap_null_precision(E_lyr, n_focal, k=10, n_boot=500, seed=42):
    """Bootstrap precision@k for random gene sets of size n_focal."""
    rng_b = np.random.default_rng(seed)
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", algorithm="ball_tree")
    nbrs.fit(E_lyr)
    N_local = E_lyr.shape[0]
    null_means = []
    for _ in range(n_boot):
        rand_idx = rng_b.choice(N_local, size=n_focal, replace=False)
        rand_set = set(rand_idx.tolist())
        E_rand = E_lyr[rand_idx]
        _, indices = nbrs.kneighbors(E_rand)
        precs = [sum(1 for nb in row[1:] if nb in rand_set) / k for row in indices]
        null_means.append(np.mean(precs))
    return np.array(null_means)


# ──────────────────────────────────────────────────────────────────────────────
# H01: Multi-cell-type kNN precision@10 panel
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H01: Multi-cell-type kNN precision@10 panel ===", flush=True)

K = 10
N_BOOT = 500
LAYERS_TO_TEST = [2, 5, 8, 11]

h01_results = {}
for ct, in_v in ct_in_vocab.items():
    if len(in_v) < 3:
        print(f"  {ct}: only {len(in_v)} in-vocab, skipping", flush=True)
        continue
    focal_idx = np.array([gene_to_local[g] for g in in_v])
    focal_set = set(focal_idx.tolist())
    ct_results = []
    for lyr in LAYERS_TO_TEST:
        E_lyr = E[lyr]
        prec = precision_at_k(E_lyr, focal_idx, focal_set, k=K)
        obs_mean = float(prec.mean())
        null_dist = bootstrap_null_precision(E_lyr, len(focal_idx), k=K,
                                             n_boot=N_BOOT, seed=42 + lyr)
        null_mean = float(null_dist.mean())
        null_std = float(null_dist.std())
        z = (obs_mean - null_mean) / (null_std + 1e-10)
        emp_p = float(np.mean(null_dist >= obs_mean))
        print(f"  {ct} L{lyr:02d}: obs={obs_mean:.3f}, null={null_mean:.3f}, z={z:.2f}, p={emp_p:.3f}")
        ct_results.append({
            "layer": lyr,
            "n_markers": len(in_v),
            "markers": in_v,
            "obs_precision_at_10": obs_mean,
            "null_mean": null_mean,
            "null_std": null_std,
            "z_score": float(z),
            "empirical_p": emp_p,
            "per_gene_precisions": prec.tolist(),
        })
    h01_results[ct] = ct_results

h01_out = {
    "hypothesis": "multi_celltype_knn_precision_panel",
    "k": K,
    "n_bootstrap": N_BOOT,
    "layers_tested": LAYERS_TO_TEST,
    "results_by_celltype": h01_results,
}
jdump(h01_out, ITER_DIR / "h01_multi_celltype_knn.json")
print("  H01 saved.", flush=True)

# ──────────────────────────────────────────────────────────────────────────────
# H02: Gene-name permutation null
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H02: Gene-name permutation null ===", flush=True)

# B-cell markers that are in-vocab
bcell_in_vocab = ct_in_vocab["B-cell"]
bcell_local_idx = np.array([gene_to_local[g] for g in bcell_in_vocab])
n_bcell = len(bcell_local_idx)

# Real precision@10 at L2 (highest z-score layer from H02 iter_0035)
TARGET_LAYER = 2
E_L2 = E[TARGET_LAYER]

# Real result
real_prec = precision_at_k(E_L2, bcell_local_idx, set(bcell_local_idx.tolist()), k=K)
real_mean = float(real_prec.mean())
print(f"  Real B-cell precision@10 at L{TARGET_LAYER}: {real_mean:.4f}", flush=True)

# Bootstrap null (same as before)
null_dist_real = bootstrap_null_precision(E_L2, n_bcell, k=K, n_boot=N_BOOT, seed=42)
z_real = (real_mean - null_dist_real.mean()) / (null_dist_real.std() + 1e-10)
print(f"  Real z-score: {z_real:.3f}", flush=True)

# Permutation null: shuffle gene-to-embedding assignment
# Generate N_PERM permutations of the N gene rows; recompute B-cell precision
N_PERM = 200
perm_precisions = []
perm_z_scores = []

rng_perm = np.random.default_rng(123)
for perm_i in range(N_PERM):
    perm_order = rng_perm.permutation(N)
    E_perm = E_L2[perm_order]  # shuffle rows (embeddings) while gene names stay fixed
    # B-cell gene names are still at local indices bcell_local_idx,
    # but their embeddings are now random rows
    perm_prec = precision_at_k(E_perm, bcell_local_idx, set(bcell_local_idx.tolist()), k=K)
    perm_mean = float(perm_prec.mean())
    # Bootstrap null for this permuted dataset
    null_perm = bootstrap_null_precision(E_perm, n_bcell, k=K, n_boot=100, seed=perm_i)
    z_perm = (perm_mean - null_perm.mean()) / (null_perm.std() + 1e-10)
    perm_precisions.append(perm_mean)
    perm_z_scores.append(float(z_perm))

perm_precisions = np.array(perm_precisions)
perm_z_scores_arr = np.array(perm_z_scores)

emp_p_perm = float(np.mean(perm_precisions >= real_mean))
perm_mean_prec = float(perm_precisions.mean())
perm_std_prec = float(perm_precisions.std())
z_against_perm = (real_mean - perm_mean_prec) / (perm_std_prec + 1e-10)

print(f"  Permutation null: mean_prec={perm_mean_prec:.4f}, std={perm_std_prec:.4f}", flush=True)
print(f"  Real vs permutation: z={z_against_perm:.2f}, emp_p={emp_p_perm:.4f}", flush=True)
print(f"  Permutation z-scores: mean={perm_z_scores_arr.mean():.2f}, std={perm_z_scores_arr.std():.2f}", flush=True)

h02_out = {
    "hypothesis": "bcell_precision_permutation_null",
    "layer": TARGET_LAYER,
    "k": K,
    "n_bcell_markers": n_bcell,
    "bcell_markers": bcell_in_vocab,
    "real_precision_at_10": real_mean,
    "real_z_score_vs_bootstrap_null": float(z_real),
    "n_permutations": N_PERM,
    "perm_null_mean_precision": perm_mean_prec,
    "perm_null_std_precision": perm_std_prec,
    "z_real_vs_perm_null": float(z_against_perm),
    "empirical_p_vs_perm_null": emp_p_perm,
    "perm_z_scores_mean": float(perm_z_scores_arr.mean()),
    "perm_z_scores_std": float(perm_z_scores_arr.std()),
    "interpretation": (
        "If z_real_vs_perm_null >> 0 and emp_p << 0.05: "
        "B-cell precision cannot be explained by random embedding geometry"
    ),
}
jdump(h02_out, ITER_DIR / "h02_bcell_perm_null.json")
print("  H02 saved.", flush=True)

# ──────────────────────────────────────────────────────────────────────────────
# H03: B-cell neighborhood functional characterization
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H03: B-cell neighborhood functional characterization ===", flush=True)

K_NEIGHBOR = 20

# At L2: collect k=20 NN for each B-cell marker
nbrs_l2 = NearestNeighbors(n_neighbors=K_NEIGHBOR + 1, metric="euclidean", algorithm="ball_tree")
nbrs_l2.fit(E_L2)
E_bcell = E_L2[bcell_local_idx]
_, indices_l2 = nbrs_l2.kneighbors(E_bcell)

# Build neighbor lists (exclude self)
bcell_set_local = set(bcell_local_idx.tolist())
neighbor_gene_counts = {}  # gene_name -> count of times appearing as neighbor
per_gene_neighbors = {}
all_neighbor_idx = []

for i, (gene_name, row) in enumerate(zip(bcell_in_vocab, indices_l2)):
    neighbors_idx = row[1:]  # exclude self
    neighbor_names = [invocab_named[nb] for nb in neighbors_idx]
    per_gene_neighbors[gene_name] = neighbor_names
    for nb in neighbors_idx:
        all_neighbor_idx.append(nb)
        nb_name = invocab_named[nb]
        neighbor_gene_counts[nb_name] = neighbor_gene_counts.get(nb_name, 0) + 1

# Rank neighbors by frequency
ranked_neighbors = sorted(neighbor_gene_counts.items(), key=lambda x: -x[1])
print(f"  Unique neighbor genes: {len(ranked_neighbors)}", flush=True)
print(f"  Top 30 B-cell neighbors at L{TARGET_LAYER}:", flush=True)
for g, cnt in ranked_neighbors[:30]:
    is_bcell = "(B-cell)" if g in set(bcell_in_vocab) else ""
    print(f"    {g}: {cnt} {is_bcell}", flush=True)

# Compute fraction of neighbor slots that are B-cell markers
total_slots = n_bcell * K_NEIGHBOR
bcell_slots = sum(1 for nb in all_neighbor_idx if nb in bcell_set_local)
frac_bcell_in_nb = bcell_slots / total_slots
expected_frac = n_bcell / N  # random expectation
enrichment = frac_bcell_in_nb / (expected_frac + 1e-10)
print(f"  B-cell fraction in neighbors: {frac_bcell_in_nb:.4f} vs expected {expected_frac:.4f}")
print(f"  Enrichment: {enrichment:.2f}x")

# Categorize neighbor genes using curated sets
IMMUNE_GENES = set(BCELL_ALL + TCELL_ALL + MYELOID_ALL + NK_ALL)
TF_GENES = set(["PAX5", "EBF1", "IRF4", "PRDM1", "XBP1", "BLIMP1", "BCL6",
                "IKZF1", "IKZF3", "NFKB1", "NFKB2", "RELA", "STAT3", "STAT5A",
                "STAT5B", "SPI1", "CEBPA", "GATA3", "TBX21", "RORC", "FOXP3",
                "IRF8", "IRF3", "IRF7"])

# Count neighbor categories
unique_neighbors = set(nb for nb, _ in ranked_neighbors)
n_bcell_in_nb = sum(1 for g in unique_neighbors if g in set(bcell_in_vocab))
n_immune_in_nb = sum(1 for g in unique_neighbors if g in IMMUNE_GENES)
n_tf_in_nb = sum(1 for g in unique_neighbors if g in TF_GENES)

print(f"\n  Neighbor gene categories (unique genes):", flush=True)
print(f"    Total unique: {len(unique_neighbors)}", flush=True)
print(f"    B-cell markers: {n_bcell_in_nb} ({n_bcell_in_nb/len(unique_neighbors):.1%})", flush=True)
print(f"    Immune genes (any type): {n_immune_in_nb} ({n_immune_in_nb/len(unique_neighbors):.1%})", flush=True)
print(f"    Known TFs: {n_tf_in_nb} ({n_tf_in_nb/len(unique_neighbors):.1%})", flush=True)

# Also test at L11 for comparison
E_L11 = E[11]
nbrs_l11 = NearestNeighbors(n_neighbors=K_NEIGHBOR + 1, metric="euclidean", algorithm="ball_tree")
nbrs_l11.fit(E_L11)
_, indices_l11 = nbrs_l11.kneighbors(E_L11[bcell_local_idx])
nb_counts_l11 = {}
all_nb_l11 = []
for row in indices_l11:
    for nb in row[1:]:
        all_nb_l11.append(nb)
        nb_name = invocab_named[nb]
        nb_counts_l11[nb_name] = nb_counts_l11.get(nb_name, 0) + 1

ranked_l11 = sorted(nb_counts_l11.items(), key=lambda x: -x[1])
bcell_slots_l11 = sum(1 for nb in all_nb_l11 if nb in bcell_set_local)
frac_l11 = bcell_slots_l11 / total_slots
enrich_l11 = frac_l11 / (expected_frac + 1e-10)
print(f"\n  L11 comparison: B-cell fraction={frac_l11:.4f}, enrichment={enrich_l11:.2f}x")
print(f"  Top 10 L11 neighbors: {[g for g, _ in ranked_l11[:10]]}")

h03_out = {
    "hypothesis": "bcell_neighborhood_characterization",
    "layer_primary": TARGET_LAYER,
    "k_neighbors": K_NEIGHBOR,
    "n_bcell_markers": n_bcell,
    "bcell_markers": bcell_in_vocab,
    "total_neighbor_slots": total_slots,
    "bcell_slots_in_neighbors": bcell_slots,
    "frac_bcell_in_neighbors": float(frac_bcell_in_nb),
    "expected_frac_random": float(expected_frac),
    "enrichment_over_random": float(enrichment),
    "n_unique_neighbor_genes": len(unique_neighbors),
    "n_bcell_in_unique_neighbors": n_bcell_in_nb,
    "n_immune_in_unique_neighbors": n_immune_in_nb,
    "n_tf_in_unique_neighbors": n_tf_in_nb,
    "top30_neighbors_by_frequency": ranked_neighbors[:30],
    "per_gene_neighbors": per_gene_neighbors,
    "l11_comparison": {
        "bcell_slots": bcell_slots_l11,
        "frac_bcell": float(frac_l11),
        "enrichment": float(enrich_l11),
        "top10_neighbors": [g for g, _ in ranked_l11[:10]],
    },
}
jdump(h03_out, ITER_DIR / "h03_bcell_neighborhood_characterization.json")
print("  H03 saved.", flush=True)

print("\n=== All experiments complete ===")
print(f"Outputs in: {ITER_DIR}")
