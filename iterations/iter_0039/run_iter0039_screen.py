"""
iter_0039 Multi-Hypothesis Screen

H01 (manifold_distance / new_method): Drift target identification
    Ask: The B-cell centroid drifts 26.4 units from L0→L11 but not toward plasma direction.
    What IS it drifting toward? Find the top-10 named genes nearest to the "drift endpoint"
    (L0_centroid + displacement_vector) in L11 embedding space.

H02 (manifold_distance / new_method): Extended plasma panel + BCL6/PAX5 GC-TF screen
    Using cycle4_immune embeddings where BCL6 and PAX5 are available:
    B-cell anchor = MS4A1, CD79A, BLK (no PRDM1 overlap artifact)
    Extended GC panel = BATF, BACH2, BCL6, PAX5
    Plasma panel = JCHAIN, SDC1 (cleaner plasma markers, no PRDM1)
    Does the GC/plasma divergence pattern persist with the cleaner anchor and extended GC set?

H03 (null_sensitivity / refinement): Leave-one-out (LOO) ablation on B-cell precision@10
    Using cycle1 main embeddings at L2 and L11:
    For each 4-gene anchor panel (LOO over {MS4A1, CD19, CD79A, BLK, PRDM1}):
    Compute B-cell centroid precision@10 for GC-TF markers (BATF, SPIB, BACH2).
    Report: which gene, if removed, most degrades the GC-TF proximity signal?
    Compare against 200 random 4-gene null sets.
"""

import numpy as np
import json
import csv
from pathlib import Path
from scipy.stats import spearmanr

# ─── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work")
CYCLE1 = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main"
CYCLE4 = ROOT / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main"
IMMUNE_H5AD = ROOT / "single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad"
ITER_DIR = ROOT / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0039"
ITER_DIR.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(42)

# ─── Load cycle1 embeddings ───────────────────────────────────────────────────
print("Loading cycle1 embeddings...", flush=True)
emb1 = np.load(CYCLE1 / "layer_gene_embeddings.npy")   # [12, 4803, 512]
N_LAYERS, N_GENES_TOTAL, N_DIM = emb1.shape
print(f"  cycle1 shape: {emb1.shape}", flush=True)

with open(CYCLE1 / "gene_list.txt") as f:
    vocab1 = [line.strip() for line in f if line.strip()]
gene_to_idx1 = {g: i for i, g in enumerate(vocab1)}

# In-vocab named genes for cycle1
EDGES_FILE1 = CYCLE1 / "cycle1_edge_dataset.tsv"
named_genes1 = set()
with open(EDGES_FILE1) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        named_genes1.add(row["source"])
        named_genes1.add(row["target"])

invocab1 = sorted([g for g in named_genes1 if g in gene_to_idx1])
norms0_1 = np.linalg.norm(emb1[0, [gene_to_idx1[g] for g in invocab1], :], axis=1)
invocab1 = [g for g, n in zip(invocab1, norms0_1) if n > 1e-8]
print(f"  cycle1 in-vocab OOV-filtered: {len(invocab1)} genes", flush=True)

gene_idx1 = np.array([gene_to_idx1[g] for g in invocab1])
E1 = emb1[:, gene_idx1, :]   # [12, N1, 512]
gene_to_local1 = {g: i for i, g in enumerate(invocab1)}
invocab1_set = set(invocab1)

# ─── Load cycle4_immune embeddings ───────────────────────────────────────────
print("Loading cycle4_immune embeddings...", flush=True)
import anndata as ad
adata4 = ad.read_h5ad(IMMUNE_H5AD)
gene_list4 = list(adata4.var_names)
gene_to_idx4 = {g: i for i, g in enumerate(gene_list4)}

emb4 = np.load(CYCLE4 / "layer_gene_embeddings.npy")   # [12, 4941, 512]
print(f"  cycle4 shape: {emb4.shape}", flush=True)

# In-vocab named genes for cycle4
EDGES_FILE4 = CYCLE4 / "cycle1_edge_dataset.tsv"
named_genes4 = set()
with open(EDGES_FILE4) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        named_genes4.add(row["source"])
        named_genes4.add(row["target"])

invocab4 = sorted([g for g in named_genes4 if g in gene_to_idx4])
norms0_4 = np.array([np.linalg.norm(emb4[0, gene_to_idx4[g], :]) for g in invocab4])
invocab4 = [g for g, n in zip(invocab4, norms0_4) if n > 1e-8]
print(f"  cycle4 in-vocab OOV-filtered: {len(invocab4)} genes", flush=True)

gene_idx4 = np.array([gene_to_idx4[g] for g in invocab4])
E4 = emb4[:, gene_idx4, :]   # [12, N4, 512]
gene_to_local4 = {g: i for i, g in enumerate(invocab4)}
invocab4_set = set(invocab4)

print("", flush=True)

# ─── Helper functions ─────────────────────────────────────────────────────────
def centroid(E_layer, gene_list, gene_to_local):
    """Compute centroid of gene embeddings at a given layer."""
    idxs = [gene_to_local[g] for g in gene_list if g in gene_to_local]
    if not idxs:
        return None, []
    found = [g for g in gene_list if g in gene_to_local]
    return E_layer[idxs].mean(axis=0), found

def z_score_dist(E_layer, point, target_genes, gene_to_local, n_perm=500):
    """
    Compute mean L2 distance from point to target_genes,
    z-scored vs null (random n_perm samples of same size from all genes).
    """
    n = E_layer.shape[0]
    idxs = [gene_to_local[g] for g in target_genes if g in gene_to_local]
    if not idxs:
        return None, None, None

    target_dists = np.linalg.norm(E_layer[idxs] - point[None, :], axis=1)
    mean_dist = target_dists.mean()

    # Null distribution: random sets of same size
    k = len(idxs)
    null_means = []
    for _ in range(n_perm):
        rand_idxs = rng.choice(n, size=k, replace=False)
        null_mean = np.linalg.norm(E_layer[rand_idxs] - point[None, :], axis=1).mean()
        null_means.append(null_mean)

    null_mean = np.mean(null_means)
    null_std = np.std(null_means)
    z = (mean_dist - null_mean) / (null_std + 1e-12)
    return mean_dist, z, null_mean

def precision_at_k(E_layer, anchor_genes, target_genes, gene_to_local, k=10):
    """
    Compute precision@k: fraction of k nearest neighbors of anchor centroid
    that are in target_genes set.
    """
    # Compute anchor centroid
    anchor_idxs = [gene_to_local[g] for g in anchor_genes if g in gene_to_local]
    if not anchor_idxs:
        return None
    c = E_layer[anchor_idxs].mean(axis=0)

    # Distance from all genes to centroid
    dists = np.linalg.norm(E_layer - c[None, :], axis=1)

    # Exclude anchor genes themselves from ranking
    for idx in anchor_idxs:
        dists[idx] = np.inf

    nn_idxs = np.argsort(dists)[:k]
    target_set = set(gene_to_local[g] for g in target_genes if g in gene_to_local)
    hits = sum(1 for idx in nn_idxs if idx in target_set)
    return hits / k


# ══════════════════════════════════════════════════════════════════════════════
# H01: Drift Target Identification
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70, flush=True)
print("H01: Drift Target Identification", flush=True)
print("=" * 70, flush=True)

# B-cell identity panel (cycle1)
BCELL_IDENTITY = ["MS4A1", "CD19", "CD79A", "BLK", "PRDM1"]
bcell_invocab1 = [g for g in BCELL_IDENTITY if g in invocab1_set]
print(f"B-cell identity genes available: {bcell_invocab1}", flush=True)

# Compute B-cell centroid at L0 and L11
bc_L0, _ = centroid(E1[0], bcell_invocab1, gene_to_local1)
bc_L11, _ = centroid(E1[11], bcell_invocab1, gene_to_local1)

# Displacement vector L0→L11 (in L0 embedding space)
# Note: L0 and L11 are in the same ambient 512-dim space, so vector arithmetic is valid
displacement = bc_L11 - bc_L0
drift_magnitude = np.linalg.norm(displacement)
print(f"Drift magnitude L0→L11: {drift_magnitude:.4f}", flush=True)

# "Drift endpoint": where would we end up if we started from L0 centroid and added drift?
# This is just bc_L11 itself, but the question is: which L11 gene embeddings are near bc_L11?
# We find top-10 genes in L11 space nearest to bc_L11 (B-cell centroid at L11)

# Get all 195 in-vocab genes' L11 embeddings
all_L11 = E1[11]   # [N1, 512]

# Distance from bc_L11 to all in-vocab genes
dists_to_endpoint = np.linalg.norm(all_L11 - bc_L11[None, :], axis=1)

# Exclude B-cell identity genes from top-k
exclude_idxs = set(gene_to_local1[g] for g in bcell_invocab1)
masked_dists = dists_to_endpoint.copy()
for idx in exclude_idxs:
    masked_dists[idx] = np.inf

# Top-20 nearest genes (excluding anchor)
top20_idxs = np.argsort(masked_dists)[:20]
top20_genes = [(invocab1[i], float(masked_dists[i])) for i in top20_idxs]

print("\nTop-20 genes nearest to B-cell centroid at L11 (drift endpoint):", flush=True)
for rank, (gene, dist) in enumerate(top20_genes, 1):
    print(f"  {rank:2d}. {gene:12s}  dist={dist:.4f}", flush=True)

# Also compute: angle between displacement direction and direction to top genes
displacement_unit = displacement / (np.linalg.norm(displacement) + 1e-12)

# L0 centroid of each top gene as seen from L0
print("\nAngle between drift vector and L0→top_gene direction:", flush=True)
drift_alignment_scores = []
for gene, dist in top20_genes[:10]:
    gene_L0 = E1[0][gene_to_local1[gene]]
    direction_to_gene = gene_L0 - bc_L0
    norm = np.linalg.norm(direction_to_gene)
    if norm > 1e-8:
        cos_sim = np.dot(displacement_unit, direction_to_gene / norm)
    else:
        cos_sim = 0.0
    drift_alignment_scores.append((gene, dist, float(cos_sim)))
    print(f"  {gene:12s}  L11_dist={dist:.4f}  cos_to_drift={cos_sim:.4f}", flush=True)

# Save H01 results
h01_results = {
    "hypothesis": "H01",
    "description": "B-cell centroid drift target: nearest genes to bc_L11 in cycle1",
    "bcell_identity_invocab": bcell_invocab1,
    "drift_magnitude_L0_to_L11": float(drift_magnitude),
    "bc_L0_norm": float(np.linalg.norm(bc_L0)),
    "bc_L11_norm": float(np.linalg.norm(bc_L11)),
    "top20_nearest_genes_at_L11": [
        {"rank": rank+1, "gene": gene, "dist_to_bcell_centroid_L11": float(dist)}
        for rank, (gene, dist) in enumerate(top20_genes)
    ],
    "drift_alignment_top10": [
        {"gene": g, "L11_dist": d, "cos_align_to_drift": c}
        for g, d, c in drift_alignment_scores
    ]
}
with open(ITER_DIR / "h01_drift_target.json", "w") as f:
    json.dump(h01_results, f, indent=2)
print("\nSaved: h01_drift_target.json", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
# H02: Extended Plasma Panel + BCL6/PAX5 GC-TF Screen (cycle4)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70, flush=True)
print("H02: Extended Plasma Panel + BCL6/PAX5 GC-TF Screen (cycle4)", flush=True)
print("=" * 70, flush=True)

# B-cell anchor (plasma-exclusive, no PRDM1 overlap)
BCELL_ANCHOR_C4 = ["MS4A1", "CD79A", "BLK"]
anchor_avail_c4 = [g for g in BCELL_ANCHOR_C4 if g in invocab4_set]
print(f"B-cell anchor (plasma-exclusive): {anchor_avail_c4}", flush=True)

# Extended GC-TF panel
EXT_GC_TFS = ["BATF", "BACH2", "BCL6", "PAX5"]
gc_avail_c4 = [g for g in EXT_GC_TFS if g in invocab4_set]
print(f"Extended GC-TF panel: {gc_avail_c4}", flush=True)

# Plasma panel (cleaner — no PRDM1)
PLASMA_CLEAN = ["JCHAIN", "SDC1"]
plasma_avail_c4 = [g for g in PLASMA_CLEAN if g in invocab4_set]
print(f"Plasma panel (clean): {plasma_avail_c4}", flush=True)

PROBE_LAYERS = [0, 2, 5, 8, 11]
N4 = E4.shape[1]

h02_layers = []
for L in PROBE_LAYERS:
    E4_L = E4[L]   # [N4, 512]

    # Compute B-cell centroid
    bc_c4, bc_found = centroid(E4_L, anchor_avail_c4, gene_to_local4)

    # GC-TF: mean dist to anchor + z-score
    gc_mean, gc_z, gc_null_mean = z_score_dist(E4_L, bc_c4, gc_avail_c4, gene_to_local4, n_perm=200)

    # Plasma: mean dist to anchor + z-score
    pl_mean, pl_z, pl_null_mean = z_score_dist(E4_L, bc_c4, plasma_avail_c4, gene_to_local4, n_perm=200)

    # Per-gene distances
    gc_per_gene = {}
    for g in gc_avail_c4:
        if g in gene_to_local4:
            d = float(np.linalg.norm(E4_L[gene_to_local4[g]] - bc_c4))
            gc_per_gene[g] = d

    pl_per_gene = {}
    for g in plasma_avail_c4:
        if g in gene_to_local4:
            d = float(np.linalg.norm(E4_L[gene_to_local4[g]] - bc_c4))
            pl_per_gene[g] = d

    h02_layers.append({
        "layer": L,
        "gc_mean_dist": float(gc_mean) if gc_mean is not None else None,
        "gc_z_vs_null": float(gc_z) if gc_z is not None else None,
        "plasma_mean_dist": float(pl_mean) if pl_mean is not None else None,
        "plasma_z_vs_null": float(pl_z) if pl_z is not None else None,
        "null_mean": float(gc_null_mean) if gc_null_mean is not None else None,
        "gc_per_gene": gc_per_gene,
        "plasma_per_gene": pl_per_gene
    })
    print(f"  L{L:2d}: GC z={gc_z:.3f}, Plasma z={pl_z:.3f}, null_mean={gc_null_mean:.3f}", flush=True)

# Compute trends
gc_zs = [r["gc_z_vs_null"] for r in h02_layers if r["gc_z_vs_null"] is not None]
pl_zs = [r["plasma_z_vs_null"] for r in h02_layers if r["plasma_z_vs_null"] is not None]
layers_for_trend = [r["layer"] for r in h02_layers if r["gc_z_vs_null"] is not None]

gc_rho, gc_p = spearmanr(layers_for_trend, gc_zs)
pl_rho, pl_p = spearmanr(layers_for_trend, pl_zs)
print(f"\nGC-TF z-trend: rho={gc_rho:.3f}, p={gc_p:.4f}", flush=True)
print(f"Plasma z-trend: rho={pl_rho:.3f}, p={pl_p:.4f}", flush=True)

# Check if divergence (plasma goes positive, GC stays negative)
plasma_final_z = pl_zs[-1] if pl_zs else None
gc_final_z = gc_zs[-1] if gc_zs else None
plasma_initial_z = pl_zs[0] if pl_zs else None
print(f"\nPlasma z: L0={plasma_initial_z:.3f} → L11={plasma_final_z:.3f}", flush=True)
print(f"GC-TF z:  L0={gc_zs[0]:.3f} → L11={gc_final_z:.3f}", flush=True)
divergence_confirmed = (plasma_final_z is not None and gc_final_z is not None and
                        plasma_final_z > 0 and gc_final_z < 0)
print(f"Divergence confirmed (plasma>0, GC<0 at L11): {divergence_confirmed}", flush=True)

h02_results = {
    "hypothesis": "H02",
    "description": "Extended plasma panel + BCL6/PAX5 GC-TF proximity (cycle4_immune, plasma-exclusive anchor)",
    "anchor_genes": anchor_avail_c4,
    "gc_tf_genes": gc_avail_c4,
    "plasma_genes": plasma_avail_c4,
    "layers": h02_layers,
    "gc_z_trend_rho": float(gc_rho),
    "gc_z_trend_p": float(gc_p),
    "plasma_z_trend_rho": float(pl_rho),
    "plasma_z_trend_p": float(pl_p),
    "plasma_z_L0": float(plasma_initial_z) if plasma_initial_z is not None else None,
    "plasma_z_L11": float(plasma_final_z) if plasma_final_z is not None else None,
    "gc_z_L0": float(gc_zs[0]) if gc_zs else None,
    "gc_z_L11": float(gc_final_z) if gc_final_z is not None else None,
    "divergence_confirmed_at_L11": bool(divergence_confirmed)
}
with open(ITER_DIR / "h02_extended_plasma_trajectory.json", "w") as f:
    json.dump(h02_results, f, indent=2)
print("\nSaved: h02_extended_plasma_trajectory.json", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
# H03: LOO Ablation on B-cell Precision@10 (cycle1)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70, flush=True)
print("H03: LOO Ablation on GC-TF Proximity (cycle1)", flush=True)
print("=" * 70, flush=True)

GC_TFS_C1 = ["BATF", "SPIB", "BACH2"]
FULL_BCELL_PANEL = ["MS4A1", "CD19", "CD79A", "BLK", "PRDM1"]
full_avail = [g for g in FULL_BCELL_PANEL if g in invocab1_set]
gc_avail_c1 = [g for g in GC_TFS_C1 if g in invocab1_set]
print(f"Full B-cell panel: {full_avail}", flush=True)
print(f"GC-TFs (targets): {gc_avail_c1}", flush=True)

N_NULL_PERM = 200
LOO_LAYERS = [2, 11]

# Full panel baseline first
loo_results_all = {}
for L in LOO_LAYERS:
    E1_L = E1[L]
    # Full panel precision@10 for GC-TFs
    full_p10 = precision_at_k(E1_L, full_avail, gc_avail_c1, gene_to_local1, k=10)

    # Full panel z-score (mean dist to bc centroid vs null)
    bc_full, _ = centroid(E1_L, full_avail, gene_to_local1)
    full_gc_mean, full_gc_z, full_gc_null = z_score_dist(E1_L, bc_full, gc_avail_c1, gene_to_local1, n_perm=N_NULL_PERM)

    print(f"\nL{L:2d} Full panel: precision@10={full_p10:.3f}, GC_z={full_gc_z:.3f}", flush=True)

    loo_results_all[L] = {
        "full_panel": {
            "genes": full_avail,
            "precision_at_10": float(full_p10) if full_p10 is not None else None,
            "gc_z": float(full_gc_z) if full_gc_z is not None else None
        },
        "loo_panels": []
    }

    # LOO: remove each gene one at a time
    for removed in full_avail:
        loo_panel = [g for g in full_avail if g != removed]
        if len(loo_panel) < 2:
            continue

        loo_p10 = precision_at_k(E1_L, loo_panel, gc_avail_c1, gene_to_local1, k=10)
        bc_loo, _ = centroid(E1_L, loo_panel, gene_to_local1)
        loo_gc_mean, loo_gc_z, loo_null = z_score_dist(E1_L, bc_loo, gc_avail_c1, gene_to_local1, n_perm=N_NULL_PERM)

        delta_p10 = (full_p10 or 0) - (loo_p10 or 0)
        delta_z = (full_gc_z or 0) - (loo_gc_z or 0)

        print(f"  LOO-{removed:8s}: p@10={loo_p10:.3f} (Δ={delta_p10:+.3f}), gc_z={loo_gc_z:.3f} (Δ={delta_z:+.3f})", flush=True)

        loo_results_all[L]["loo_panels"].append({
            "removed_gene": removed,
            "remaining_panel": loo_panel,
            "precision_at_10": float(loo_p10) if loo_p10 is not None else None,
            "gc_z": float(loo_gc_z) if loo_gc_z is not None else None,
            "delta_precision": float(delta_p10),
            "delta_gc_z": float(delta_z)
        })

    # Also: null comparison — random 4-gene panels
    print(f"\n  Null distribution (200 random 4-gene panels at L{L}):", flush=True)
    null_p10s = []
    null_zs = []
    n_inv = len(invocab1)
    for _ in range(N_NULL_PERM):
        rand_panel = list(rng.choice(n_inv, size=4, replace=False))
        rand_genes = [invocab1[i] for i in rand_panel]
        rand_bc, _ = centroid(E1_L, rand_genes, gene_to_local1)
        if rand_bc is None:
            continue
        rand_p10 = precision_at_k(E1_L, rand_genes, gc_avail_c1, gene_to_local1, k=10)
        # z-score vs second-level null would be expensive; just record raw p10
        null_p10s.append(float(rand_p10) if rand_p10 is not None else 0.0)

    null_p10_mean = np.mean(null_p10s)
    null_p10_std = np.std(null_p10s)
    null_p10_pctile = np.mean([p <= (full_p10 or 0) for p in null_p10s])
    print(f"  Null p@10: mean={null_p10_mean:.3f}, std={null_p10_std:.3f}, full_pctile={null_p10_pctile:.3f}", flush=True)

    loo_results_all[L]["null_p10_mean"] = float(null_p10_mean)
    loo_results_all[L]["null_p10_std"] = float(null_p10_std)
    loo_results_all[L]["full_panel_pctile_vs_null"] = float(null_p10_pctile)

h03_results = {
    "hypothesis": "H03",
    "description": "LOO ablation: which B-cell anchor gene most contributes to GC-TF proximity signal?",
    "gc_tf_targets": gc_avail_c1,
    "full_bcell_panel": full_avail,
    "layers": LOO_LAYERS,
    "results_by_layer": loo_results_all
}
with open(ITER_DIR / "h03_loo_ablation.json", "w") as f:
    json.dump(h03_results, f, indent=2)
print("\nSaved: h03_loo_ablation.json", flush=True)

print("\n" + "=" * 70, flush=True)
print("iter_0039 screen complete.", flush=True)
print("=" * 70, flush=True)
