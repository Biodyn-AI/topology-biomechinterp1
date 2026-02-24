"""
iter_0044 Multi-Hypothesis Screen

H01 (intrinsic_dimensionality / new_method): Cross-lineage ID compression law
    Extend effector CD8 compression finding (iter_0043: delta_ID = -45.8) to:
    - NK effector (NCAM1/FCGR3A/KLRB1/KLRD1)
    - Plasma cell (IGHG1/IGKC/JCHAIN/MZB1)
    - Neutrophil effector (S100A8/S100A9/CXCR1/FCGR3B)
    - Pan-myeloid reference (LYZ/CST3/AIF1/FCN1)
    TwoNN ID at each layer; test if effector subtypes compress at L7 while reference expands.

H02 (manifold_distance / new_method): T-cell activation circuit attractor test
    CD28/LAT/LCK/ZAP70 — classical TCR signaling circuit.
    Test if these 4 genes converge to T-cell centroid (CD3E/CD3D/TRAC/TRBC1/CD247) at deep layers.
    Analog of the B-cell GC attractor result. Measure precision@k=10 and mean rank at L0 vs L11.

H03 (manifold_distance / refinement): BCL6 rank-metabolic correlation
    Use existing iter_0042 data: BCL6 nearest-neighbor ranks at L0-L11.
    Correlate BCL6 rank-to-B-cell-centroid with its rank-to-metabolic cluster (NAMPT, GLUL, etc).
    Test: do the two tracks diverge (negative correlation) or co-move?
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
ITER_DIR = ROOT / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0044"
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
# TwoNN intrinsic dimensionality estimator (Facco et al. 2017)
# =====================================================================
def twonn_id(X):
    """Estimate intrinsic dimensionality of point cloud X using TwoNN."""
    n = X.shape[0]
    if n < 3:
        return np.nan
    dists = cdist(X, X, metric='euclidean')
    np.fill_diagonal(dists, np.inf)
    sorted_dists = np.sort(dists, axis=1)
    r1 = sorted_dists[:, 0]  # nearest neighbor
    r2 = sorted_dists[:, 1]  # 2nd nearest neighbor
    # avoid division by zero
    valid = (r1 > 0) & (r2 > 0) & np.isfinite(r1) & np.isfinite(r2)
    mu = r2[valid] / r1[valid]
    mu = mu[mu > 1]  # should always be >= 1 by definition
    if len(mu) < 2:
        return np.nan
    # MLE: ID = -n / sum(ln(mu))
    id_est = len(mu) / np.sum(np.log(mu))
    return float(id_est)


# =====================================================================
# H01: Cross-lineage ID compression law
# =====================================================================
print("\n=== H01: Cross-lineage ID compression law ===", flush=True)

gene_sets = {
    "nk_effector":    ["NCAM1", "FCGR3A", "KLRB1", "KLRD1", "GNLY"],
    "plasma_cell":    ["IGHG1", "IGKC", "JCHAIN", "MZB1", "DERL3"],
    "neutrophil":     ["S100A8", "S100A9", "CXCR1", "FCGR3B", "MPO"],
    "pan_myeloid":    ["LYZ", "CST3", "AIF1", "FCN1", "CD14"],
    "effector_cd8":   ["CD8A", "GZMB", "PRF1", "NKG7", "GNLY"],   # reference from iter_0043
    "general_t":      ["CD3E", "CD3D", "TRAC", "TRBC1"],            # reference pan-T
}

h1_results = {}
for name, genes in gene_sets.items():
    in_vocab = [g for g in genes if g in gene2idx]
    print(f"  {name}: {in_vocab} (n={len(in_vocab)})", flush=True)
    if len(in_vocab) < 3:
        h1_results[name] = {"in_vocab": in_vocab, "status": "insufficient_genes"}
        continue
    ids = []
    for layer in range(N_LAYERS):
        X = emb[layer, [gene2idx[g] for g in in_vocab], :]
        ids.append(twonn_id(X))
    ids = np.array(ids)
    # Measure compression at L7 breakpoint (layers 7-11 vs 0-6)
    pre_mean = float(np.nanmean(ids[:7]))
    post_mean = float(np.nanmean(ids[7:]))
    slope_change = post_mean - pre_mean
    h1_results[name] = {
        "in_vocab": in_vocab,
        "ids_per_layer": [float(x) if np.isfinite(x) else None for x in ids],
        "pre_l7_mean": pre_mean,
        "post_l7_mean": post_mean,
        "l7_slope_change": slope_change,
    }
    print(f"    ID pre-L7: {pre_mean:.1f}, post-L7: {post_mean:.1f}, delta: {slope_change:.1f}", flush=True)

out_path_h1 = ITER_DIR / "h01_crosslineage_id.json"
with open(out_path_h1, "w") as f:
    json.dump(h1_results, f, indent=2)
print(f"  Saved: {out_path_h1}", flush=True)


# =====================================================================
# H02: T-cell activation circuit attractor test
# =====================================================================
print("\n=== H02: T-cell activation circuit attractor test ===", flush=True)

TCR_CIRCUIT = ["CD28", "LAT", "LCK", "ZAP70", "CD247"]
T_CENTROID_GENES = ["CD3E", "CD3D", "TRAC", "TRBC1", "CD247"]

circuit_in_vocab = [g for g in TCR_CIRCUIT if g in gene2idx]
t_centroid_in = [g for g in T_CENTROID_GENES if g in gene2idx]
print(f"  TCR circuit in-vocab: {circuit_in_vocab}", flush=True)
print(f"  T-cell centroid genes: {t_centroid_in}", flush=True)

h2_results = {
    "circuit_genes": circuit_in_vocab,
    "centroid_genes": t_centroid_in,
    "layers": {},
}

if len(circuit_in_vocab) >= 2 and len(t_centroid_in) >= 2:
    circuit_ids = [gene2idx[g] for g in circuit_in_vocab]
    centroid_ids = [gene2idx[g] for g in t_centroid_in]
    all_gene_ids = list(range(len(vocab)))

    for layer in range(N_LAYERS):
        layer_emb = emb[layer]  # [n_genes, 512]
        centroid_vec = layer_emb[centroid_ids].mean(axis=0)  # [512]

        # Distance from each circuit gene to T-cell centroid
        dists_to_centroid = np.linalg.norm(layer_emb - centroid_vec, axis=1)  # [n_genes]

        # Rank each circuit gene among all genes (lower rank = closer)
        ranks = {}
        mean_ranks = []
        for g in circuit_in_vocab:
            gid = gene2idx[g]
            d = dists_to_centroid[gid]
            rank = int(np.sum(dists_to_centroid < d)) + 1
            ranks[g] = rank
            mean_ranks.append(rank)

        # Precision@k=20: how many circuit genes are in top-20 nearest to centroid?
        k = 20
        top_k = set(np.argsort(dists_to_centroid)[:k])
        circuit_in_top_k = sum(1 for g in circuit_in_vocab if gene2idx[g] in top_k)
        precision_at_k = circuit_in_top_k / len(circuit_in_vocab)

        h2_results["layers"][str(layer)] = {
            "ranks": ranks,
            "mean_rank": float(np.mean(mean_ranks)),
            "precision_at_20": float(precision_at_k),
            "n_circuit_in_top20": circuit_in_top_k,
        }
        print(f"  L{layer:02d}: mean_rank={np.mean(mean_ranks):.0f}, prec@20={precision_at_k:.2f}, "
              f"top20 count={circuit_in_top_k}/{len(circuit_in_vocab)}", flush=True)

    # Summary: L0 vs L11 comparison
    l0 = h2_results["layers"]["0"]
    l11 = h2_results["layers"]["11"]
    rank_change = l11["mean_rank"] - l0["mean_rank"]
    h2_results["rank_change_l0_to_l11"] = rank_change
    h2_results["prec_change_l0_to_l11"] = l11["precision_at_20"] - l0["precision_at_20"]
    print(f"\n  Summary: rank change L0->L11: {rank_change:.0f}, prec@20 change: {h2_results['prec_change_l0_to_l11']:.2f}", flush=True)

    # Null: random gene sets of same size, 200 permutations
    n_null = 200
    null_rank_changes = []
    for _ in range(n_null):
        null_genes = rng.choice(all_gene_ids, size=len(circuit_in_vocab), replace=False)
        for lyr in [0, 11]:
            layer_emb = emb[lyr]
            centroid_vec = layer_emb[centroid_ids].mean(axis=0)
            dists = np.linalg.norm(layer_emb - centroid_vec, axis=1)
            null_ranks = [int(np.sum(dists < dists[g])) + 1 for g in null_genes]
            if lyr == 0:
                null_r0 = np.mean(null_ranks)
            else:
                null_r11 = np.mean(null_ranks)
        null_rank_changes.append(null_r11 - null_r0)

    null_arr = np.array(null_rank_changes)
    pval = float(np.mean(null_arr <= rank_change))
    h2_results["null_rank_change_mean"] = float(null_arr.mean())
    h2_results["null_rank_change_std"] = float(null_arr.std())
    h2_results["null_pval_rank_change"] = pval
    print(f"  Null rank change mean={null_arr.mean():.0f} ± {null_arr.std():.0f}, p-val={pval:.3f}", flush=True)
else:
    print("  BLOCKED: insufficient genes in vocab", flush=True)
    h2_results["status"] = "blocked"

out_path_h2 = ITER_DIR / "h02_tcr_circuit_attractor.json"
with open(out_path_h2, "w") as f:
    json.dump(h2_results, f, indent=2)
print(f"  Saved: {out_path_h2}", flush=True)


# =====================================================================
# H03: BCL6 rank-metabolic correlation (refinement from iter_0042)
# =====================================================================
print("\n=== H03: BCL6 rank-metabolic correlation ===", flush=True)

BCL6_GENE = "BCL6"
B_CELL_GENES = ["MS4A1", "CD19", "CD79A", "BLK", "PAX5"]
METABOLIC_GENES = ["NAMPT", "GLUL", "PFKFB3", "ACSL1", "NIBAN1", "FNDC3B", "VMP1", "HMGCS1", "FASN", "FDFT1"]

h3_results = {
    "bcl6_in_vocab": BCL6_GENE in gene2idx,
    "status": "blocked" if BCL6_GENE not in gene2idx else "computed",
}

if BCL6_GENE in gene2idx:
    bcl6_id = gene2idx[BCL6_GENE]
    b_in_vocab = [g for g in B_CELL_GENES if g in gene2idx]
    met_in_vocab = [g for g in METABOLIC_GENES if g in gene2idx]
    print(f"  B-cell centroid genes: {b_in_vocab}", flush=True)
    print(f"  Metabolic genes: {met_in_vocab}", flush=True)

    bcell_ranks = []
    metabolic_ranks = []
    bcell_dists = []
    metabolic_dists = []

    for layer in range(N_LAYERS):
        layer_emb = emb[layer]
        all_dists = np.linalg.norm(layer_emb - layer_emb[bcl6_id], axis=1)

        # Rank to B-cell centroid
        if b_in_vocab:
            b_centroid = layer_emb[[gene2idx[g] for g in b_in_vocab]].mean(axis=0)
            d_bcell = np.linalg.norm(layer_emb - b_centroid, axis=1)
            bcell_rank = int(np.sum(d_bcell < d_bcell[bcl6_id])) + 1
            bcell_ranks.append(bcell_rank)
            bcell_dists.append(float(d_bcell[bcl6_id]))

        # Rank to metabolic cluster centroid
        if met_in_vocab:
            met_centroid = layer_emb[[gene2idx[g] for g in met_in_vocab]].mean(axis=0)
            d_met = np.linalg.norm(layer_emb - met_centroid, axis=1)
            met_rank = int(np.sum(d_met < d_met[bcl6_id])) + 1
            metabolic_ranks.append(met_rank)
            metabolic_dists.append(float(d_met[bcl6_id]))

    bcell_ranks = np.array(bcell_ranks)
    metabolic_ranks = np.array(metabolic_ranks)

    # Spearman correlation between tracks (do they co-move or diverge?)
    if len(bcell_ranks) == len(metabolic_ranks) == N_LAYERS:
        rho, pval = spearmanr(bcell_ranks, metabolic_ranks)
        h3_results["bcell_ranks_by_layer"] = bcell_ranks.tolist()
        h3_results["metabolic_ranks_by_layer"] = metabolic_ranks.tolist()
        h3_results["bcell_dists_by_layer"] = bcell_dists
        h3_results["metabolic_dists_by_layer"] = metabolic_dists
        h3_results["spearman_rho_rank_tracks"] = float(rho)
        h3_results["spearman_pval"] = float(pval)
        h3_results["interpretation"] = (
            "co-move" if rho > 0.3 else ("diverge" if rho < -0.3 else "uncorrelated")
        )
        print(f"  BCL6 rank tracks: B-cell={bcell_ranks.tolist()}", flush=True)
        print(f"  BCL6 rank tracks: Metabolic={metabolic_ranks.tolist()}", flush=True)
        print(f"  Spearman rho={rho:.3f}, p={pval:.4f}, interpretation={h3_results['interpretation']}", flush=True)

out_path_h3 = ITER_DIR / "h03_bcl6_rank_metabolic.json"
with open(out_path_h3, "w") as f:
    json.dump(h3_results, f, indent=2)
print(f"  Saved: {out_path_h3}", flush=True)

print("\nAll done.", flush=True)
