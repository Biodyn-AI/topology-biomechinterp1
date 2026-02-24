"""
iter_0037 Multi-Hypothesis Screen

H01 (manifold_distance / new_method): Extended cell-type panel expansion
    Test NK cells, dendritic cells, plasma cells (+ B-cell, T-cell, Myeloid for reference).
    Plasma cells expected near B-cells (differentiation trajectory).
    Also compute centroid-centroid distance matrix across all cell types at L2.

H02 (manifold_distance / new_method): Held-out B-cell marker generalization
    The original 7-marker panel: MS4A1, CD19, CD79A, CD79B, BLK, BANK1, PRDM1.
    Find additional B-cell/plasma TFs/markers in vocab that were NOT in original panel:
    PAX5, IRF4, EBF1, XBP1, BCL6, FOXO1, IKZF1, IRF8, MEF2C, SPIB, BACH2.
    Test if held-out markers cluster near the 7-gene reference centroid at L2.
    Expected: strong positive for B-cell identity TFs (PAX5, EBF1), moderate for plasma (XBP1, IRF4).

H03 (null_sensitivity / refinement): T-cell permutation null + STRING TF neighborhood scoring
    Part A: Run gene-name permutation null for T-cell markers (same method as H02 iter_0036).
    Expected: T-cell z should NOT be significant, confirming the permutation null works differently.
    Part B: For top-50 B-cell neighbors at L2, compute STRING PPI scoring:
    fraction with STRING score >= 400 to at least one B-cell marker.
    Compare against random gene set STRING connectivity.
"""

import numpy as np
import json
import csv
import urllib.request
import urllib.parse
from pathlib import Path
from scipy.stats import mannwhitneyu
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
ITER_DIR = PROJECT / "iterations" / "iter_0037"
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

# Filter to in-vocab, drop zero-norm (OOV)
invocab_named = sorted([g for g in named_genes_set if g in gene_to_emb_idx])
L0_NORMS = np.linalg.norm(emb[0, [gene_to_emb_idx[g] for g in invocab_named], :], axis=1)
invocab_named = [g for g, n in zip(invocab_named, L0_NORMS) if n > 1e-8]
print(f"  In-vocab OOV-filtered: {len(invocab_named)} genes", flush=True)

gene_idx = np.array([gene_to_emb_idx[g] for g in invocab_named])
N = len(invocab_named)
gene_to_local = {g: i for i, g in enumerate(invocab_named)}

# Embeddings per layer for in-vocab genes: [12, N, 512]
E = emb[:, gene_idx, :]  # [12, N, 512]
invocab_set = set(invocab_named)

# ─── Cell-type marker sets ────────────────────────────────────────────────────
# ORIGINAL 7-gene B-cell panel (used in iter_0035/0036)
BCELL_ORIGINAL_7 = ["MS4A1", "CD19", "CD79A", "CD79B", "BLK", "BANK1", "PRDM1"]

# Full B-cell/plasma marker sets (for reference)
BCELL_ALL = [
    "MS4A1", "CD19", "CD79A", "CD79B", "PAX5", "BLK", "BANK1",
    "FCRL1", "CD22", "FCER2", "IGHM", "IGHG1", "IGHA1", "CR2",
    "BLNK", "CD27", "IRF4", "PRDM1", "XBP1", "EBF1",
    "VPREB1", "VPREB3", "IGHD", "IGKC", "IGLC1",
]
PLASMA_MARKERS = ["PRDM1", "XBP1", "IRF4", "CD38", "SDC1", "MZB1", "JCHAIN", "IGHA1", "IGHG1"]
TCELL_ALL = [
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B",
    "TRAC", "TRBC1", "TRBC2", "LCK", "ZAP70", "LAT",
    "IL7R", "SELL", "TCF7", "CCR7", "GZMB", "GZMK",
    "PRF1", "NKG7", "CD44", "FOXP3", "IL2RA",
]
MYELOID_ALL = [
    "CD14", "CD68", "CSF1R", "FCGR3A", "S100A8", "S100A9",
    "LYZ", "MARCO", "CD163", "MRC1",
]
NK_ALL = [
    "NCAM1", "KLRD1", "KLRC1", "KLRB1", "KLRF1",
    "NKG7", "GNLY", "PRF1", "GZMB", "FCGR3A",
    "TYROBP", "FCER1G", "NCR1", "NCR3", "CD56",
]
DC_ALL = [
    "FCER1A", "CLEC9A", "XCR1", "CLEC10A", "CD1C",
    "SIRPA", "ITGAX", "HLA-DRA", "HLA-DRB1", "CD83",
    "LAMP3", "IDO1", "CLEC4C", "IL3RA", "LILRA4",
    "PTPRC", "CCR7", "CD86", "CD80",
]

# ─── Helper functions ──────────────────────────────────────────────────────────
def filter_in_vocab(gene_list):
    return [g for g in gene_list if g in invocab_set]

def precision_at_k(E_layer, marker_local_idx, marker_set_local, k=10):
    """Precision@k: fraction of kNN that are in marker set (excluding self)."""
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", algorithm="ball_tree")
    nbrs.fit(E_layer)
    E_markers = E_layer[marker_local_idx]
    _, indices = nbrs.kneighbors(E_markers)
    precisions = []
    for row in indices:
        neighbors = row[1:]  # exclude self
        p = sum(1 for nb in neighbors if nb in marker_set_local) / k
        precisions.append(p)
    return np.array(precisions)

def bootstrap_null_precision(E_layer, n_markers, k=10, n_boot=500, seed=42):
    """Bootstrap null: random sets of size n_markers, compute mean precision@k."""
    rng_b = np.random.default_rng(seed)
    N_total = E_layer.shape[0]
    null_means = []
    for _ in range(n_boot):
        idx = rng_b.choice(N_total, size=n_markers, replace=False)
        idx_set = set(idx.tolist())
        prec = precision_at_k(E_layer, idx, idx_set, k=k)
        null_means.append(prec.mean())
    return np.array(null_means)

# ─── Common settings ───────────────────────────────────────────────────────────
K = 10
N_BOOT = 500
LAYERS_TO_TEST = [0, 2, 5, 8, 11]

# ──────────────────────────────────────────────────────────────────────────────
# H01: Extended cell-type panel expansion + centroid distance matrix
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H01: Extended cell-type panel expansion ===", flush=True)

cell_types_h01 = {
    "B-cell": BCELL_ORIGINAL_7,
    "T-cell": TCELL_ALL,
    "Myeloid": MYELOID_ALL,
    "NK": NK_ALL,
    "DC": DC_ALL,
    "Plasma": PLASMA_MARKERS,
}

ct_in_vocab = {}
for ct, markers in cell_types_h01.items():
    in_v = filter_in_vocab(markers)
    ct_in_vocab[ct] = in_v
    print(f"  {ct}: {len(in_v)} in-vocab from {len(markers)} total: {in_v}", flush=True)

h01_results = {}
for ct, markers in ct_in_vocab.items():
    if len(markers) < 2:
        print(f"  Skipping {ct} (only {len(markers)} in-vocab markers)", flush=True)
        h01_results[ct] = {"skipped": True, "n_markers": len(markers), "reason": "too_few"}
        continue
    local_idx = np.array([gene_to_local[g] for g in markers])
    local_set = set(local_idx.tolist())
    ct_result = {"n_markers": len(markers), "markers": markers, "layers": {}}
    for layer in LAYERS_TO_TEST:
        E_layer = E[layer]
        prec = precision_at_k(E_layer, local_idx, local_set, k=K)
        obs_mean = float(prec.mean())
        null_dist = bootstrap_null_precision(E_layer, len(markers), k=K, n_boot=N_BOOT, seed=layer*100)
        z = (obs_mean - null_dist.mean()) / (null_dist.std() + 1e-10)
        p_emp = float(np.mean(null_dist >= obs_mean))
        ct_result["layers"][f"L{layer}"] = {
            "obs_precision": obs_mean,
            "null_mean": float(null_dist.mean()),
            "null_std": float(null_dist.std()),
            "z_score": float(z),
            "emp_p": p_emp,
        }
        print(f"  {ct} L{layer}: prec={obs_mean:.4f} z={z:.2f} p={p_emp:.4f}", flush=True)
    h01_results[ct] = ct_result

# Centroid-centroid distance matrix at L2
print("\n  Computing centroid-centroid distance matrix at L2...", flush=True)
E_L2 = E[2]
centroid_matrix = {}
ct_names = [ct for ct, r in h01_results.items() if not r.get("skipped")]
centroids = {}
for ct in ct_names:
    local_idx = np.array([gene_to_local[g] for g in ct_in_vocab[ct]])
    centroids[ct] = E_L2[local_idx].mean(axis=0)

dist_matrix = {}
for ct1 in ct_names:
    dist_matrix[ct1] = {}
    for ct2 in ct_names:
        d = float(np.linalg.norm(centroids[ct1] - centroids[ct2]))
        dist_matrix[ct1][ct2] = d
        if ct1 < ct2:
            print(f"  dist({ct1},{ct2}) = {d:.4f}", flush=True)

h01_out = {
    "hypothesis": "extended_celltype_panel_knn_precision",
    "k": K, "n_bootstrap": N_BOOT, "layers_tested": LAYERS_TO_TEST,
    "results_by_celltype": h01_results,
    "centroid_distance_matrix_L2": dist_matrix,
}
jdump(h01_out, ITER_DIR / "h01_extended_panel_knn.json")
print("  H01 saved.", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
# H02: Held-out B-cell marker generalization
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H02: Held-out B-cell marker generalization ===", flush=True)

# Original 7-gene panel (reference)
BCELL_REF_PANEL = BCELL_ORIGINAL_7
bcell_ref_in_vocab = filter_in_vocab(BCELL_REF_PANEL)
print(f"  Reference panel in-vocab: {len(bcell_ref_in_vocab)} genes: {bcell_ref_in_vocab}", flush=True)

# Held-out B-cell TFs/markers (NOT in original 7-gene panel)
BCELL_HELD_OUT = [
    # Core B-cell identity TFs
    "PAX5",   # master B-cell lineage TF
    "EBF1",   # early B-cell factor
    "SPIB",   # ETS-family TF (was top neighbor in H03 iter_0036)
    "BACH2",  # germinal center B-cell TF (was top neighbor in H03 iter_0036)
    "BATF",   # AP-1 TF (was top neighbor in H03 iter_0036)
    # Plasma cell / late B-cell
    "XBP1",   # plasma cell differentiation TF
    "IRF4",   # plasma cell TF
    "BCL6",   # germinal center TF (represses plasma differentiation)
    # Additional B-cell surface markers
    "CD22",   # B-cell inhibitory receptor
    "FCER2",  # CD23, low-affinity IgE receptor on B-cells
    "BLNK",   # B-cell linker protein
    "CR2",    # CD21, complement receptor on B-cells
    # B-cell signaling
    "IKZF1",  # Ikaros, lymphocyte development TF
    "IKZF3",  # Aiolos, B-cell maturation TF
    "MEF2C",  # B-cell survival/proliferation TF
    "IRF8",   # B-cell commitment TF
    "FOXO1",  # B-cell survival/homing TF
]

held_out_in_vocab = filter_in_vocab(BCELL_HELD_OUT)
print(f"  Held-out in-vocab: {len(held_out_in_vocab)} genes: {held_out_in_vocab}", flush=True)

# Reference centroid at L2
ref_local_idx = np.array([gene_to_local[g] for g in bcell_ref_in_vocab])
E_L2 = E[2]
ref_centroid = E_L2[ref_local_idx].mean(axis=0)

# For each held-out gene: distance to B-cell centroid vs null distribution
# Null: random gene distances to same centroid
all_gene_dists = np.linalg.norm(E_L2 - ref_centroid[np.newaxis, :], axis=1)

# Overall distribution
null_mean = float(all_gene_dists.mean())
null_std = float(all_gene_dists.std())
print(f"  Null distribution (all {N} genes): mean={null_mean:.4f}, std={null_std:.4f}", flush=True)

# Reference genes (lower bound on expected distance)
ref_dists = all_gene_dists[ref_local_idx]
print(f"  Reference genes ({len(ref_local_idx)}): mean dist={ref_dists.mean():.4f}, std={ref_dists.std():.4f}", flush=True)

# Held-out gene distances and z-scores
h02_gene_results = []
for gene in held_out_in_vocab:
    local_i = gene_to_local[gene]
    dist = float(all_gene_dists[local_i])
    z = (dist - null_mean) / (null_std + 1e-10)
    # Negative z = closer to centroid than average (positive signal for clustering)
    pctile = float(np.mean(all_gene_dists > dist))  # fraction farther = percentile of closeness
    h02_gene_results.append({
        "gene": gene,
        "dist_to_bcell_centroid": dist,
        "z_vs_null": float(z),
        "pctile_closeness": pctile,  # high = very close
    })
    print(f"  {gene}: dist={dist:.4f} z={z:.2f} pctile_closeness={pctile:.3f}", flush=True)

# Sort by distance (ascending = closer = more B-cell-like)
h02_gene_results.sort(key=lambda x: x["dist_to_bcell_centroid"])

# Summary: how many held-out genes are in top quartile of closeness?
n_close = sum(1 for r in h02_gene_results if r["pctile_closeness"] >= 0.75)
n_medium = sum(1 for r in h02_gene_results if 0.5 <= r["pctile_closeness"] < 0.75)
print(f"\n  Summary: {n_close}/{len(h02_gene_results)} held-out in top-quartile closeness to B-cell centroid", flush=True)
print(f"  {n_medium}/{len(h02_gene_results)} in 50-75th percentile", flush=True)

# Also test across layers
h02_layer_results = {}
for layer in LAYERS_TO_TEST:
    E_layer = E[layer]
    centroid_l = E_layer[ref_local_idx].mean(axis=0)
    all_dists_l = np.linalg.norm(E_layer - centroid_l[np.newaxis, :], axis=1)
    null_m = all_dists_l.mean()
    null_s = all_dists_l.std()
    held_dists = [float(all_dists_l[gene_to_local[g]]) for g in held_out_in_vocab]
    z_scores = [(d - null_m) / (null_s + 1e-10) for d in held_dists]
    mean_z = float(np.mean(z_scores))
    pctiles = [float(np.mean(all_dists_l > all_dists_l[gene_to_local[g]])) for g in held_out_in_vocab]
    frac_top_quartile = float(np.mean([p >= 0.75 for p in pctiles]))
    h02_layer_results[f"L{layer}"] = {
        "mean_z_held_out": mean_z,
        "mean_pctile_closeness": float(np.mean(pctiles)),
        "frac_top_quartile": frac_top_quartile,
    }
    print(f"  L{layer}: held-out mean_z={mean_z:.3f}, frac_top_quartile={frac_top_quartile:.2f}", flush=True)

h02_out = {
    "hypothesis": "held_out_bcell_marker_generalization",
    "reference_panel": bcell_ref_in_vocab,
    "held_out_genes": held_out_in_vocab,
    "n_held_out_in_vocab": len(held_out_in_vocab),
    "primary_layer": 2,
    "null_mean_dist": null_mean,
    "null_std_dist": null_std,
    "ref_mean_dist": float(ref_dists.mean()),
    "ref_std_dist": float(ref_dists.std()),
    "gene_results_L2": h02_gene_results,
    "n_held_out_top_quartile_L2": n_close,
    "layer_summary": h02_layer_results,
}
jdump(h02_out, ITER_DIR / "h02_held_out_bcell_generalization.json")
print("  H02 saved.", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
# H03: T-cell permutation null + B-cell STRING TF neighborhood scoring
# ──────────────────────────────────────────────────────────────────────────────
print("\n=== H03: T-cell permutation null + STRING TF neighborhood scoring ===", flush=True)

# Part A: T-cell permutation null
print("\n  Part A: T-cell permutation null", flush=True)
tcell_in_vocab = filter_in_vocab(TCELL_ALL)
print(f"  T-cell in-vocab: {len(tcell_in_vocab)} markers: {tcell_in_vocab}", flush=True)

tcell_local_idx = np.array([gene_to_local[g] for g in tcell_in_vocab])
n_tcell = len(tcell_local_idx)
tcell_local_set = set(tcell_local_idx.tolist())

E_L2 = E[2]

# Real T-cell precision@10 at L2
real_tcell_prec = precision_at_k(E_L2, tcell_local_idx, tcell_local_set, k=K)
real_tcell_mean = float(real_tcell_prec.mean())
null_tcell = bootstrap_null_precision(E_L2, n_tcell, k=K, n_boot=N_BOOT, seed=42)
z_tcell_real = (real_tcell_mean - null_tcell.mean()) / (null_tcell.std() + 1e-10)
print(f"  Real T-cell precision@10 L2: {real_tcell_mean:.4f}, z={z_tcell_real:.2f}", flush=True)

# Permutation null for T-cell
N_PERM = 200
perm_prec_tcell = []
perm_z_tcell = []
rng_perm = np.random.default_rng(456)
for perm_i in range(N_PERM):
    perm_order = rng_perm.permutation(N)
    E_perm = E_L2[perm_order]
    p = precision_at_k(E_perm, tcell_local_idx, tcell_local_set, k=K)
    pm = float(p.mean())
    null_p = bootstrap_null_precision(E_perm, n_tcell, k=K, n_boot=100, seed=perm_i)
    z_p = (pm - null_p.mean()) / (null_p.std() + 1e-10)
    perm_prec_tcell.append(pm)
    perm_z_tcell.append(float(z_p))

perm_prec_tcell = np.array(perm_prec_tcell)
perm_z_tcell_arr = np.array(perm_z_tcell)
z_tcell_vs_perm = (real_tcell_mean - perm_prec_tcell.mean()) / (perm_prec_tcell.std() + 1e-10)
emp_p_tcell = float(np.mean(perm_prec_tcell >= real_tcell_mean))

print(f"  T-cell perm null: perm mean_prec={perm_prec_tcell.mean():.4f}", flush=True)
print(f"  T-cell vs perm null: z={z_tcell_vs_perm:.2f}, emp_p={emp_p_tcell:.4f}", flush=True)
print(f"  T-cell perm z-scores: mean={perm_z_tcell_arr.mean():.2f}, std={perm_z_tcell_arr.std():.2f}", flush=True)

# Part B: STRING TF neighborhood scoring for B-cell neighbors
print("\n  Part B: STRING TF neighborhood scoring", flush=True)
# Load STRING interactions from local file (if available)
STRING_FILE = CYCLE1 / "string_interactions.tsv"
string_edges = {}
string_scores = {}
if STRING_FILE.exists():
    print(f"  Loading STRING from {STRING_FILE}...", flush=True)
    with open(STRING_FILE) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            g1, g2 = row.get("source", row.get("gene1", "")), row.get("target", row.get("gene2", ""))
            score = float(row.get("score", row.get("combined_score", 0)))
            if g1 and g2 and score >= 400:
                string_edges.setdefault(g1, set()).add(g2)
                string_edges.setdefault(g2, set()).add(g1)
                string_scores[(g1, g2)] = score
                string_scores[(g2, g1)] = score
    print(f"  STRING loaded: {len(string_edges)} genes with edges", flush=True)
else:
    print(f"  STRING file not found at {STRING_FILE}, using fallback known B-cell TF-marker edges", flush=True)

# B-cell marker neighbors at L2 (from iter_0036 H03 — we recompute here)
bcell_ref_local_idx = np.array([gene_to_local[g] for g in bcell_ref_in_vocab])
bcell_set_local = set(bcell_ref_local_idx.tolist())

K_NEIGHBOR = 50  # top-50 neighbors
nbrs_l2 = NearestNeighbors(n_neighbors=K_NEIGHBOR + 1, metric="euclidean", algorithm="ball_tree")
nbrs_l2.fit(E_L2)
_, indices_l2 = nbrs_l2.kneighbors(E_L2[bcell_ref_local_idx])

# Collect unique neighbors (excluding the B-cell markers themselves)
neighbor_freq = {}
for row in indices_l2:
    for nb in row[1:]:
        if nb not in bcell_set_local:
            nm = invocab_named[nb]
            neighbor_freq[nm] = neighbor_freq.get(nm, 0) + 1

# Top-50 unique neighbors by frequency
top_neighbors = sorted(neighbor_freq.items(), key=lambda x: -x[1])[:50]
top_nb_genes = [g for g, _ in top_neighbors]
print(f"  Top-10 B-cell neighbors: {top_nb_genes[:10]}", flush=True)

# Known B-cell TF-gene regulatory edges (from TRRUST/literature)
# as a fallback for STRING scoring
KNOWN_BCELL_TF_TARGETS = {
    # TF -> targets known from TRRUST/literature
    "PAX5": ["CD19", "CD79A", "MS4A1", "IGHM", "EBF1"],
    "EBF1": ["CD79A", "CD79B", "MS4A1", "BLNK"],
    "PRDM1": ["XBP1", "IRF4"],
    "IRF4": ["PRDM1", "XBP1"],
    "XBP1": ["IGHA1", "IGHG1", "MZB1"],
    "BCL6": ["PRDM1", "IRF4", "XBP1"],
    "SPIB": ["CD79A", "CD79B", "BLNK"],
    "BACH2": ["PRDM1"],
    "IKZF1": ["CD79B", "MS4A1"],
}
# Build set of known B-cell gene pairs
known_bcell_pairs = set()
for tf, targets in KNOWN_BCELL_TF_TARGETS.items():
    for t in targets:
        known_bcell_pairs.add((tf, t))
        known_bcell_pairs.add((t, tf))

# Compute B-cell connectivity score for each top neighbor
# Score: number of B-cell marker genes it has known edges with
bcell_marker_names = set(bcell_ref_in_vocab)
nb_scores = []
for gene in top_nb_genes:
    if string_edges:
        # Use STRING
        string_connected = bcell_marker_names & string_edges.get(gene, set())
        n_connected = len(string_connected)
        connection_type = "STRING"
    else:
        # Use known TF-target fallback
        string_connected = set()
        for bm in bcell_marker_names:
            if (gene, bm) in known_bcell_pairs or (bm, gene) in known_bcell_pairs:
                string_connected.add(bm)
        n_connected = len(string_connected)
        connection_type = "known_TF_target"

    # Also check if gene is a known B-cell TF
    is_bcell_tf = gene in set(BCELL_HELD_OUT + list(KNOWN_BCELL_TF_TARGETS.keys()))
    nb_scores.append({
        "gene": gene,
        "freq_as_neighbor": neighbor_freq[gene],
        "n_bcell_markers_connected": n_connected,
        "connected_bcell_markers": list(string_connected),
        "connection_type": connection_type,
        "is_known_bcell_tf": is_bcell_tf,
    })

# Summary
n_with_any_connection = sum(1 for s in nb_scores if s["n_bcell_markers_connected"] > 0)
n_known_tf = sum(1 for s in nb_scores if s["is_known_bcell_tf"])
print(f"  Top-50 neighbors with any B-cell connection: {n_with_any_connection}/50", flush=True)
print(f"  Top-50 neighbors that are known B-cell TFs: {n_known_tf}/50", flush=True)

# Random baseline: same stats for random 50 genes
rng_baseline = np.random.default_rng(99)
random_baseline_connected = []
for trial in range(200):
    rand_genes = [invocab_named[i] for i in rng_baseline.choice(N, size=50, replace=False)]
    n_conn = 0
    for gene in rand_genes:
        for bm in bcell_marker_names:
            if (gene, bm) in known_bcell_pairs or (bm, gene) in known_bcell_pairs:
                n_conn += 1
                break
    random_baseline_connected.append(n_conn)

rand_mean = float(np.mean(random_baseline_connected))
rand_std = float(np.std(random_baseline_connected))
z_nb_vs_rand = (n_with_any_connection - rand_mean) / (rand_std + 1e-10)
print(f"  Random baseline (200 trials): mean={rand_mean:.2f}, std={rand_std:.2f}", flush=True)
print(f"  z(top-50 neighbors vs random): {z_nb_vs_rand:.2f}", flush=True)

h03_out = {
    "hypothesis": "tcell_perm_null_plus_string_tf_scoring",
    "part_A": {
        "tcell_markers_in_vocab": tcell_in_vocab,
        "n_tcell": n_tcell,
        "real_tcell_precision_L2": real_tcell_mean,
        "z_tcell_vs_bootstrap_null": float(z_tcell_real),
        "n_permutations": N_PERM,
        "perm_null_mean": float(perm_prec_tcell.mean()),
        "perm_null_std": float(perm_prec_tcell.std()),
        "z_tcell_vs_perm_null": float(z_tcell_vs_perm),
        "empirical_p_tcell_vs_perm_null": emp_p_tcell,
        "perm_z_mean": float(perm_z_tcell_arr.mean()),
        "perm_z_std": float(perm_z_tcell_arr.std()),
        "interpretation": "If z_vs_perm_null NOT significant: T-cell has no robust clustering, confirming B-cell specificity"
    },
    "part_B": {
        "k_neighbors": K_NEIGHBOR,
        "bcell_markers_used": bcell_ref_in_vocab,
        "n_top_neighbors_analyzed": len(top_nb_genes),
        "top_neighbors": top_nb_genes,
        "n_with_bcell_connection": n_with_any_connection,
        "n_known_bcell_tf": n_known_tf,
        "random_baseline_mean": rand_mean,
        "random_baseline_std": rand_std,
        "z_vs_random": float(z_nb_vs_rand),
        "neighbor_scores": nb_scores,
        "connection_type_used": "STRING" if string_edges else "known_TF_target",
    },
}
jdump(h03_out, ITER_DIR / "h03_tcell_perm_string_scoring.json")
print("  H03 saved.", flush=True)

print("\n=== All experiments complete ===")
print(f"Outputs in: {ITER_DIR}")
