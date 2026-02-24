"""
iter_0038 Multi-Hypothesis Screen

H01 (manifold_distance / new_method): GC-TF vs plasma-TF centroid separation across layers
    Ask: does scGPT encode B→plasma differentiation as monotonically growing geometric separation
    between germinal-center TF centroid (BATF, SPIB, BACH2) and plasma TF centroid (IRF4, PRDM1)?
    Also compare GC-TF centroid vs B-cell identity centroid distance across layers.

H02 (manifold_distance / new_method): B-cell centroid directional drift
    Track B-cell centroid vector across L0→L11. Test if late-layer displacement from L0 aligns
    with the direction FROM B-cell centroid TOWARD plasma cell centroid.
    Metric: cosine similarity between drift vector and B→plasma direction vector, per layer.

H03 (null_sensitivity + manifold_distance): Master-TF proximity + NK/myeloid specificity screen
    Part A: PAX5/EBF1/BCL6 proximity to B-cell centroid at L2 — master B-cell TFs not in original panel.
    Part B: NK cell markers precision@10 at L2 (expected negative, validates B-cell specificity).
    Part C: Myeloid markers precision@10 at L2 (expected negative).
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr

# ─── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work")
CYCLE1 = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main"
ITER_DIR = ROOT / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0038"
ITER_DIR.mkdir(parents=True, exist_ok=True)

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

print(f"  Gene vocabulary size for analysis: {N}", flush=True)

# ─── Marker sets ─────────────────────────────────────────────────────────────
# B-cell identity markers (reference panel from iter_0036/0037)
BCELL_IDENTITY = ["MS4A1", "CD19", "CD79A", "BLK", "PRDM1"]
bcell_invocab = [g for g in BCELL_IDENTITY if g in invocab_set]
print(f"  B-cell identity in-vocab: {bcell_invocab} (n={len(bcell_invocab)})")

# GC-TF germinal center TFs (strong proximity signal from iter_0037)
GC_TFS = ["BATF", "SPIB", "BACH2"]
gc_invocab = [g for g in GC_TFS if g in invocab_set]
print(f"  GC-TFs in-vocab: {gc_invocab} (n={len(gc_invocab)})")

# Plasma-TF differentiation TFs (low proximity = different geometry)
PLASMA_TFS = ["IRF4", "PRDM1"]  # PRDM1 also in B-cell identity panel
plasma_invocab = [g for g in PLASMA_TFS if g in invocab_set]
print(f"  Plasma-TFs in-vocab: {plasma_invocab} (n={len(plasma_invocab)})")

# Master B-cell TFs for H03 Part A
MASTER_BCELL_TFS = ["PAX5", "EBF1", "BCL6"]
master_invocab = [g for g in MASTER_BCELL_TFS if g in invocab_set]
print(f"  Master B-cell TFs in-vocab: {master_invocab} (n={len(master_invocab)})")

# NK cell markers for H03 Part B
NK_MARKERS = ["KLRB1", "NKG7", "GNLY", "KLRD1", "NCR1", "KLRC1", "XCL1", "GZMB", "PRF1", "CD244"]
nk_invocab = [g for g in NK_MARKERS if g in invocab_set]
print(f"  NK markers in-vocab: {nk_invocab} (n={len(nk_invocab)})")

# Myeloid markers for H03 Part C
MYELOID_MARKERS = ["CD14", "LYZ", "CST3", "FCGR3A", "MS4A7", "CTSS", "SERPINA1", "AIF1", "TYROBP"]
myeloid_invocab = [g for g in MYELOID_MARKERS if g in invocab_set]
print(f"  Myeloid markers in-vocab: {myeloid_invocab} (n={len(myeloid_invocab)})")

# ─── Helper functions ─────────────────────────────────────────────────────────
def get_centroid(layer, gene_list):
    """Compute centroid of gene list at given layer."""
    idxs = [gene_to_local[g] for g in gene_list if g in gene_to_local]
    if not idxs:
        return None
    return E[layer, idxs, :].mean(axis=0)

def l2_dist(a, b):
    return float(np.linalg.norm(a - b))

def precision_at_k(query_genes, all_genes, layer, k=10):
    """Fraction of k-NN of each query gene that are also query genes."""
    if len(query_genes) < 2:
        return None
    query_set = set(query_genes)
    query_local = [gene_to_local[g] for g in query_genes if g in gene_to_local]
    if not query_local:
        return None
    hits = 0
    total = 0
    for qi in query_local:
        q_emb = E[layer, qi, :]
        dists = np.linalg.norm(E[layer] - q_emb, axis=1)
        dists[qi] = np.inf  # exclude self
        nn_idxs = np.argsort(dists)[:k]
        nn_genes = set(all_genes[i] for i in nn_idxs)
        hits += len(nn_genes & query_set)
        total += k
    return hits / total if total > 0 else 0.0

def bootstrap_precision_null(n_query, all_genes, layer, k=10, n_boot=500, seed=42):
    """Bootstrap null distribution for precision@k with random gene sets of same size."""
    rng_boot = np.random.default_rng(seed)
    n_all = len(all_genes)
    precisions = []
    for _ in range(n_boot):
        rand_idxs = rng_boot.choice(n_all, size=n_query, replace=False)
        rand_genes = [all_genes[i] for i in rand_idxs]
        p = precision_at_k(rand_genes, all_genes, layer, k)
        if p is not None:
            precisions.append(p)
    return np.array(precisions)

# ─── H01: GC-TF vs Plasma-TF centroid separation across layers ───────────────
print("\n=== H01: GC-TF vs Plasma-TF centroid separation ===", flush=True)

h01_results = {
    "hypothesis": "H01",
    "description": "GC-TF vs plasma-TF centroid distance across layers L0-L11",
    "gc_tfs_invocab": gc_invocab,
    "plasma_tfs_invocab": plasma_invocab,
    "bcell_identity_invocab": bcell_invocab,
    "layers": []
}

# Pure plasma TFs (excluding PRDM1 which is in B-cell panel)
plasma_pure = [g for g in PLASMA_TFS if g in invocab_set and g != "PRDM1"]
print(f"  Pure plasma TFs (excl. PRDM1): {plasma_pure}")

# Also check IRF4 alone as plasma marker
irf4_in_vocab = "IRF4" in invocab_set
print(f"  IRF4 in vocab: {irf4_in_vocab}")

for layer in range(N_LAYERS):
    gc_centroid = get_centroid(layer, gc_invocab)
    bcell_centroid = get_centroid(layer, bcell_invocab)

    # Per-gene distances to B-cell centroid
    gc_dists_to_bcell = []
    for g in gc_invocab:
        if g in gene_to_local:
            g_emb = E[layer, gene_to_local[g], :]
            gc_dists_to_bcell.append(float(np.linalg.norm(g_emb - bcell_centroid)))

    plasma_dists_to_bcell = []
    for g in plasma_invocab:
        if g in gene_to_local:
            g_emb = E[layer, gene_to_local[g], :]
            plasma_dists_to_bcell.append(float(np.linalg.norm(g_emb - bcell_centroid)))

    # Also compute GC centroid vs plasma centroid distance
    if "IRF4" in gene_to_local and gc_centroid is not None:
        irf4_emb = E[layer, gene_to_local["IRF4"], :]
        gc_to_irf4_dist = float(np.linalg.norm(gc_centroid - irf4_emb))
    else:
        gc_to_irf4_dist = None

    # Null: all 195 genes distances to B-cell centroid
    all_dists_to_bcell = np.linalg.norm(E[layer] - bcell_centroid, axis=1)
    null_mean = float(all_dists_to_bcell.mean())
    null_std = float(all_dists_to_bcell.std())

    gc_mean_dist = float(np.mean(gc_dists_to_bcell)) if gc_dists_to_bcell else None
    plasma_mean_dist = float(np.mean(plasma_dists_to_bcell)) if plasma_dists_to_bcell else None

    gc_z = (gc_mean_dist - null_mean) / null_std if gc_mean_dist is not None else None
    plasma_z = (plasma_mean_dist - null_mean) / null_std if plasma_mean_dist is not None else None

    # GC centroid to plasma centroid distance
    plasma_centroid = get_centroid(layer, plasma_invocab)
    gc_to_plasma_centroid = l2_dist(gc_centroid, plasma_centroid) if (gc_centroid is not None and plasma_centroid is not None) else None
    gc_to_bcell_centroid = l2_dist(gc_centroid, bcell_centroid) if (gc_centroid is not None and bcell_centroid is not None) else None

    layer_data = {
        "layer": layer,
        "gc_mean_dist_to_bcell": gc_mean_dist,
        "gc_z_vs_null": gc_z,
        "plasma_mean_dist_to_bcell": plasma_mean_dist,
        "plasma_z_vs_null": plasma_z,
        "gc_to_bcell_centroid_dist": gc_to_bcell_centroid,
        "gc_to_plasma_centroid_dist": gc_to_plasma_centroid,
        "null_mean": null_mean,
        "null_std": null_std,
        "gc_per_gene": {g: d for g, d in zip(gc_invocab, gc_dists_to_bcell)},
        "plasma_per_gene": {g: d for g, d in zip(plasma_invocab, plasma_dists_to_bcell)},
    }
    h01_results["layers"].append(layer_data)

    print(f"  L{layer:02d}: GC z={gc_z:.2f}" + (f" ({gc_mean_dist:.2f})" if gc_mean_dist else "") +
          f" | Plasma z={plasma_z:.2f}" + (f" ({plasma_mean_dist:.2f})" if plasma_mean_dist else "") +
          f" | null_mean={null_mean:.2f}", flush=True)

# Summarize: is there monotonic trend in separation?
gc_plasma_sep = [row["gc_to_plasma_centroid_dist"] for row in h01_results["layers"] if row["gc_to_plasma_centroid_dist"] is not None]
layers_for_sep = [row["layer"] for row in h01_results["layers"] if row["gc_to_plasma_centroid_dist"] is not None]
if len(gc_plasma_sep) >= 3:
    rho_sep, p_sep = spearmanr(layers_for_sep, gc_plasma_sep)
    h01_results["gc_plasma_separation_spearman_rho"] = float(rho_sep)
    h01_results["gc_plasma_separation_spearman_p"] = float(p_sep)
    print(f"  GC-plasma centroid separation trend: rho={rho_sep:.3f}, p={p_sep:.4f}", flush=True)

gc_bcell_dists = [row["gc_to_bcell_centroid_dist"] for row in h01_results["layers"] if row["gc_to_bcell_centroid_dist"] is not None]
if len(gc_bcell_dists) >= 3:
    rho_gb, p_gb = spearmanr(layers_for_sep, gc_bcell_dists)
    h01_results["gc_bcell_centroid_trend_rho"] = float(rho_gb)
    h01_results["gc_bcell_centroid_trend_p"] = float(p_gb)
    print(f"  GC-bcell centroid distance trend: rho={rho_gb:.3f}, p={p_gb:.4f}", flush=True)

h01_path = ITER_DIR / "h01_gc_plasma_separation.json"
with open(h01_path, "w") as f:
    json.dump(h01_results, f, indent=2)
print(f"  Saved: {h01_path}", flush=True)

# ─── H02: B-cell centroid directional drift ──────────────────────────────────
print("\n=== H02: B-cell centroid directional drift ===", flush=True)

h02_results = {
    "hypothesis": "H02",
    "description": "B-cell centroid drift alignment with B→plasma direction vector",
    "bcell_identity_invocab": bcell_invocab,
    "plasma_tfs_invocab": plasma_invocab,
    "layers": []
}

# B-cell centroid at L0 (anchor)
bcell_L0 = get_centroid(0, bcell_invocab)
bcell_Lmax = get_centroid(N_LAYERS - 1, bcell_invocab)

# B→plasma direction: vector from B-cell centroid to plasma centroid at L11
plasma_L11 = get_centroid(N_LAYERS - 1, plasma_invocab)
bc_to_plasma_vec = plasma_L11 - bcell_Lmax
bc_to_plasma_unit = bc_to_plasma_vec / (np.linalg.norm(bc_to_plasma_vec) + 1e-10)

# Also use GC-TF as "B-cell identity" anchor
gc_L11 = get_centroid(N_LAYERS - 1, gc_invocab)
gc_to_plasma_vec = plasma_L11 - gc_L11
gc_to_plasma_unit = gc_to_plasma_vec / (np.linalg.norm(gc_to_plasma_vec) + 1e-10)

print(f"  B-cell L0→L11 centroid distance: {l2_dist(bcell_L0, bcell_Lmax):.3f}", flush=True)
print(f"  B-cell centroid to plasma centroid at L11: {l2_dist(bcell_Lmax, plasma_L11):.3f}", flush=True)

# Also compute the "B→GC TF" direction as control
gc_L11 = get_centroid(N_LAYERS - 1, gc_invocab)
bc_to_gc_vec = gc_L11 - bcell_Lmax
bc_to_gc_unit = bc_to_gc_vec / (np.linalg.norm(bc_to_gc_vec) + 1e-10)

for layer in range(N_LAYERS):
    bc_layer = get_centroid(layer, bcell_invocab)
    gc_layer = get_centroid(layer, gc_invocab)
    plasma_layer = get_centroid(layer, plasma_invocab)

    # Drift vector from L0
    drift_vec = bc_layer - bcell_L0
    drift_norm = np.linalg.norm(drift_vec)

    # Alignment of drift with B→plasma direction
    cosine_to_plasma = float(np.dot(drift_vec / (drift_norm + 1e-10), bc_to_plasma_unit)) if drift_norm > 1e-10 else 0.0
    # Alignment of drift with B→GC direction
    cosine_to_gc = float(np.dot(drift_vec / (drift_norm + 1e-10), bc_to_gc_unit)) if drift_norm > 1e-10 else 0.0

    # Distance between GC and plasma centroids at this layer
    gc_plasma_dist_layer = l2_dist(gc_layer, plasma_layer)

    # Angle between GC and Plasma from B-cell centroid
    bc_to_gc_here = gc_layer - bc_layer
    bc_to_plasma_here = plasma_layer - bc_layer
    cos_angle_gc_plasma = float(np.dot(bc_to_gc_here / (np.linalg.norm(bc_to_gc_here) + 1e-10),
                                        bc_to_plasma_here / (np.linalg.norm(bc_to_plasma_here) + 1e-10)))
    angle_gc_plasma_deg = float(np.degrees(np.arccos(np.clip(cos_angle_gc_plasma, -1, 1))))

    layer_data = {
        "layer": layer,
        "drift_magnitude": float(drift_norm),
        "cosine_to_plasma_direction": cosine_to_plasma,
        "cosine_to_gc_direction": cosine_to_gc,
        "gc_plasma_dist": gc_plasma_dist_layer,
        "angle_gc_plasma_from_bcell_deg": angle_gc_plasma_deg,
        "bc_centroid_norm": float(np.linalg.norm(bc_layer)),
        "gc_centroid_norm": float(np.linalg.norm(gc_layer)),
        "plasma_centroid_norm": float(np.linalg.norm(plasma_layer)),
    }
    h02_results["layers"].append(layer_data)

    print(f"  L{layer:02d}: drift_mag={drift_norm:.3f}, cos_to_plasma={cosine_to_plasma:.3f}, "
          f"cos_to_gc={cosine_to_gc:.3f}, GC-plasma_dist={gc_plasma_dist_layer:.3f}, "
          f"angle(GC,plasma|BC)={angle_gc_plasma_deg:.1f}°", flush=True)

# Trend analysis
drift_mags = [row["drift_magnitude"] for row in h02_results["layers"]]
cos_to_plasma = [row["cosine_to_plasma_direction"] for row in h02_results["layers"]]
cos_to_gc = [row["cosine_to_gc_direction"] for row in h02_results["layers"]]
gc_plasma_dists = [row["gc_plasma_dist"] for row in h02_results["layers"]]
layers_all = list(range(N_LAYERS))

rho_drift, p_drift = spearmanr(layers_all, drift_mags)
rho_cos_plasma, p_cos_plasma = spearmanr(layers_all, cos_to_plasma)
rho_cos_gc, p_cos_gc = spearmanr(layers_all, cos_to_gc)
rho_gcdist, p_gcdist = spearmanr(layers_all, gc_plasma_dists)

h02_results.update({
    "drift_magnitude_spearman": {"rho": float(rho_drift), "p": float(p_drift)},
    "cosine_to_plasma_trend": {"rho": float(rho_cos_plasma), "p": float(p_cos_plasma)},
    "cosine_to_gc_trend": {"rho": float(rho_cos_gc), "p": float(p_cos_gc)},
    "gc_plasma_dist_trend": {"rho": float(rho_gcdist), "p": float(p_gcdist)},
    "interpretation": "If cosine_to_plasma trend is positive, drift aligns with B→plasma direction across layers."
})

print(f"\n  Drift magnitude trend: rho={rho_drift:.3f}, p={p_drift:.4f}", flush=True)
print(f"  Cosine to plasma trend: rho={rho_cos_plasma:.3f}, p={p_cos_plasma:.4f}", flush=True)
print(f"  Cosine to GC trend: rho={rho_cos_gc:.3f}, p={p_cos_gc:.4f}", flush=True)
print(f"  GC-plasma separation trend: rho={rho_gcdist:.3f}, p={p_gcdist:.4f}", flush=True)

h02_path = ITER_DIR / "h02_directional_drift.json"
with open(h02_path, "w") as f:
    json.dump(h02_results, f, indent=2)
print(f"  Saved: {h02_path}", flush=True)

# ─── H03: Master TF proximity + NK/myeloid specificity screen ────────────────
print("\n=== H03: Master TF proximity + NK/Myeloid specificity screen ===", flush=True)

TARGET_LAYER = 2  # L2 has best B-cell signal
K = 10
N_BOOT = 500

bcell_centroid_L2 = get_centroid(TARGET_LAYER, bcell_invocab)
all_dists_L2 = np.linalg.norm(E[TARGET_LAYER] - bcell_centroid_L2, axis=1)
null_mean_L2 = float(all_dists_L2.mean())
null_std_L2 = float(all_dists_L2.std())
null_pctile_L2 = np.percentile(all_dists_L2, np.arange(0, 101, 1))

h03_results = {
    "hypothesis": "H03",
    "description": "Master B-cell TF proximity + NK/myeloid precision@10 at L2",
    "target_layer": TARGET_LAYER,
    "k": K,
    "bcell_identity_panel": bcell_invocab,
    "null_mean_dist": null_mean_L2,
    "null_std_dist": null_std_L2,
}

# Part A: Master B-cell TF proximity
print(f"\n  Part A: Master B-cell TFs at L{TARGET_LAYER}", flush=True)
master_proximity = {}
for g in MASTER_BCELL_TFS:
    if g in gene_to_local:
        g_emb = E[TARGET_LAYER, gene_to_local[g], :]
        dist = float(np.linalg.norm(g_emb - bcell_centroid_L2))
        z = (dist - null_mean_L2) / null_std_L2
        pctile_rank = float(np.mean(all_dists_L2 > dist))  # fraction farther = closeness pctile
        master_proximity[g] = {
            "dist_to_bcell_centroid": dist,
            "z_vs_null": float(z),
            "closeness_percentile": float(pctile_rank),
            "in_vocab": True
        }
        print(f"    {g}: dist={dist:.2f}, z={z:.2f}, closeness_pctile={pctile_rank:.3f}", flush=True)
    else:
        master_proximity[g] = {"in_vocab": False}
        print(f"    {g}: NOT in vocab", flush=True)

# Also check GC-TFs for reference
print(f"\n  Part A reference: GC-TFs at L{TARGET_LAYER}", flush=True)
gc_proximity = {}
for g in gc_invocab:
    if g in gene_to_local:
        g_emb = E[TARGET_LAYER, gene_to_local[g], :]
        dist = float(np.linalg.norm(g_emb - bcell_centroid_L2))
        z = (dist - null_mean_L2) / null_std_L2
        pctile_rank = float(np.mean(all_dists_L2 > dist))
        gc_proximity[g] = {
            "dist_to_bcell_centroid": dist,
            "z_vs_null": float(z),
            "closeness_percentile": float(pctile_rank),
        }
        print(f"    {g}: dist={dist:.2f}, z={z:.2f}, closeness_pctile={pctile_rank:.3f}", flush=True)

h03_results["master_bcell_tfs"] = master_proximity
h03_results["gc_tfs_reference"] = gc_proximity

# Part B: NK cell precision@10
print(f"\n  Part B: NK cell precision@k={K} at L{TARGET_LAYER}", flush=True)
nk_prec = None
nk_null_dist = None
nk_z = None
if len(nk_invocab) >= 2:
    nk_prec = precision_at_k(nk_invocab, invocab_named, TARGET_LAYER, K)
    nk_null = bootstrap_precision_null(len(nk_invocab), invocab_named, TARGET_LAYER, K, N_BOOT)
    nk_null_mean = float(nk_null.mean())
    nk_null_std = float(nk_null.std())
    nk_z = (nk_prec - nk_null_mean) / (nk_null_std + 1e-10)
    nk_emp_p = float(np.mean(nk_null >= nk_prec))
    print(f"    NK prec@{K}={nk_prec:.4f}, null_mean={nk_null_mean:.4f}, z={nk_z:.2f}, emp_p={nk_emp_p:.3f}", flush=True)
    h03_results["nk_precision_at_k"] = {
        "n_invocab": len(nk_invocab),
        "genes": nk_invocab,
        "precision_at_k": float(nk_prec),
        "null_mean": nk_null_mean,
        "null_std": nk_null_std,
        "z_score": float(nk_z),
        "empirical_p": nk_emp_p,
        "decision": "positive" if nk_z > 2.0 else ("negative" if nk_emp_p > 0.8 else "inconclusive")
    }
else:
    print(f"    NK: only {len(nk_invocab)} in-vocab (need >=2), skipping", flush=True)
    h03_results["nk_precision_at_k"] = {"n_invocab": len(nk_invocab), "status": "skipped"}

# Part C: Myeloid precision@10
print(f"\n  Part C: Myeloid cell precision@k={K} at L{TARGET_LAYER}", flush=True)
myeloid_prec = None
if len(myeloid_invocab) >= 2:
    myeloid_prec = precision_at_k(myeloid_invocab, invocab_named, TARGET_LAYER, K)
    myeloid_null = bootstrap_precision_null(len(myeloid_invocab), invocab_named, TARGET_LAYER, K, N_BOOT)
    myeloid_null_mean = float(myeloid_null.mean())
    myeloid_null_std = float(myeloid_null.std())
    myeloid_z = (myeloid_prec - myeloid_null_mean) / (myeloid_null_std + 1e-10)
    myeloid_emp_p = float(np.mean(myeloid_null >= myeloid_prec))
    print(f"    Myeloid prec@{K}={myeloid_prec:.4f}, null_mean={myeloid_null_mean:.4f}, z={myeloid_z:.2f}, emp_p={myeloid_emp_p:.3f}", flush=True)
    h03_results["myeloid_precision_at_k"] = {
        "n_invocab": len(myeloid_invocab),
        "genes": myeloid_invocab,
        "precision_at_k": float(myeloid_prec),
        "null_mean": myeloid_null_mean,
        "null_std": myeloid_null_std,
        "z_score": float(myeloid_z),
        "empirical_p": float(myeloid_emp_p),
        "decision": "positive" if myeloid_z > 2.0 else ("negative" if myeloid_emp_p > 0.8 else "inconclusive")
    }
else:
    print(f"    Myeloid: only {len(myeloid_invocab)} in-vocab (need >=2), skipping", flush=True)
    h03_results["myeloid_precision_at_k"] = {"n_invocab": len(myeloid_invocab), "status": "skipped"}

# B-cell precision at L2 for reference
bcell_prec_L2 = precision_at_k(bcell_invocab, invocab_named, TARGET_LAYER, K)
bcell_null = bootstrap_precision_null(len(bcell_invocab), invocab_named, TARGET_LAYER, K, N_BOOT)
bcell_null_mean = float(bcell_null.mean())
bcell_null_std = float(bcell_null.std())
bcell_z_L2 = (bcell_prec_L2 - bcell_null_mean) / (bcell_null_std + 1e-10)
print(f"\n  Reference: B-cell prec@{K} at L{TARGET_LAYER}={bcell_prec_L2:.4f}, z={bcell_z_L2:.2f}", flush=True)
h03_results["bcell_reference"] = {
    "precision_at_k": float(bcell_prec_L2),
    "z_score": float(bcell_z_L2),
    "null_mean": bcell_null_mean
}

h03_path = ITER_DIR / "h03_master_tf_nk_myeloid.json"
with open(h03_path, "w") as f:
    json.dump(h03_results, f, indent=2)
print(f"  Saved: {h03_path}", flush=True)

print("\n=== All experiments complete ===", flush=True)
print(f"Artifacts in: {ITER_DIR}", flush=True)
