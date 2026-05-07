#!/usr/bin/env python3
"""E8 — Dynamic vs static test for "regulatory dynamics" claim.

Reviewer 1 (#5) and Reviewer 3 (#2): the paper's claim that the model has
"internalized aspects of regulatory dynamics" was based on the across-layer
geometric trajectory of BATF/BACH2 toward PAX5. The reviewers correctly
note that without a comparison to actual cellular dynamics, this could be
a network statement (the GC regulatory circuit is encoded as a structure)
rather than a dynamics claim.

This test compares:
  (A) The across-layer trajectory of each GC-TF: distance to B-cell centroid
      as a function of scGPT layer (from existing layer_gene_embeddings.npy).
  (B) The across-pseudotime trajectory of the same TFs: pseudotime-binned
      mean expression of the TF in B-cell h5ad (subproject_29).

If across-layer convergence is correlated with across-pseudotime expression
rise, the dynamics interpretation has support. If not, the manuscript must
reframe the claim.

Concretely: for each of {BATF, BACH2, PAX5, BCL6, PRDM1, IRF4}, we report:
  - Spearman ρ between layer index and distance-to-PAX5 (from existing
    layer-averaged embeddings).
  - Spearman ρ between pseudotime bin and mean expression (from B-cell h5ad).
  - Whether the two trajectories agree in direction (both negative ⇒
    convergence in space + expression rise; supports dynamics).
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

EMB = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy"
H5AD_IMMUNE = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
BCELL_H5AD = REPO / "subproject_29_pseudotime_directionality_audit/implementation/adata_B_cell.h5ad"

GC_TFS = ["BATF", "BACH2", "PAX5", "BCL6", "PRDM1", "IRF4"]
B_CELL_MARKERS = ["CD19", "CD79A", "MS4A1", "BLK", "VPREB3", "FCRL1", "PAX5"]


def _h5ad_index(p):
    with h5py.File(p, "r") as f:
        raw = f["var"]["_index"][:]
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]


def _spearman(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    if len(x) < 3:
        return float("nan")
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    rxc = rx - rx.mean()
    ryc = ry - ry.mean()
    denom = np.sqrt((rxc * rxc).sum() * (ryc * ryc).sum()) + 1e-12
    return float((rxc * ryc).sum() / denom)


def _bcell_pseudotime_expression(p, gene_list, n_bins=10):
    """Return: (bins, n_genes) array of mean expression per pseudotime bin."""
    with h5py.File(p, "r") as f:
        # Pseudotime
        pt_obj = f["obs"]["dpt_pseudotime"]
        if isinstance(pt_obj, h5py.Group):
            cats = pt_obj["categories"][:]
            codes = pt_obj["codes"][:]
            pt = np.array([float(cats[c]) if cats[c] != b"NA" else np.nan for c in codes])
        else:
            pt = pt_obj[:]
        # Gene index (Ensembl IDs)
        var_idx = f["var"]["_index"][:]
        var_idx = [g.decode() if isinstance(g, bytes) else str(g) for g in var_idx]
        # feature_name (HGNC)
        fn = f["var"]["feature_name"]
        if isinstance(fn, h5py.Group):
            cats = fn["categories"][:]
            codes = fn["codes"][:]
            feature_name = [cats[c].decode() if isinstance(cats[c], bytes) else str(cats[c]) for c in codes]
        else:
            feature_name = [g.decode() for g in fn[:]]
        name_to_var = {g: i for i, g in enumerate(feature_name)}

        # X
        X = f["X"]
        if isinstance(X, h5py.Group):
            from scipy.sparse import csr_matrix
            data = X["data"][:]
            indices = X["indices"][:]
            indptr = X["indptr"][:]
            n_obs = len(indptr) - 1
            n_var = int(np.max(indices) + 1) if len(indices) else 0
            mat = csr_matrix((data, indices, indptr), shape=(n_obs, n_var))
            arr_full = mat
        else:
            arr_full = X

        # Filter cells with valid pseudotime
        valid = ~np.isnan(pt)
        if not valid.any():
            print("  no valid pseudotime values", flush=True)
            return None, None, None
        pt_v = pt[valid]
        # Gene indices
        gene_var_idx = [name_to_var.get(g, -1) for g in gene_list]

        # Bin pseudotime
        bin_edges = np.quantile(pt_v, np.linspace(0, 1, n_bins + 1))
        bin_assign = np.clip(np.digitize(pt_v, bin_edges[1:-1]), 0, n_bins - 1)

        # Compute mean expression per bin per gene
        valid_idx = np.where(valid)[0]
        bin_means = np.zeros((n_bins, len(gene_list)))
        for bi in range(n_bins):
            cells_in_bin = valid_idx[bin_assign == bi]
            if len(cells_in_bin) == 0:
                continue
            for gi, var_idx_g in enumerate(gene_var_idx):
                if var_idx_g < 0:
                    continue
                if isinstance(arr_full, h5py.Dataset):
                    vals = arr_full[cells_in_bin, var_idx_g]
                else:
                    vals = arr_full[cells_in_bin, var_idx_g].toarray().flatten()
                bin_means[bi, gi] = float(vals.mean())

        return bin_means, bin_assign, pt_v


def main():
    print("[E8] Loading…", flush=True)
    gene_index = _h5ad_index(H5AD_IMMUNE)
    n2i = {g: i for i, g in enumerate(gene_index)}

    emb = np.load(EMB, mmap_mode="r")
    n_layers = emb.shape[0]
    print(f"  emb: {emb.shape}", flush=True)

    # Across-layer trajectory: B-cell centroid distance per TF
    # B-cell centroid = mean of B-cell markers
    marker_idx = [n2i[g] for g in B_CELL_MARKERS if g in n2i]
    print(f"  B-cell markers in vocab: {len(marker_idx)}", flush=True)
    layer_dist = np.zeros((len(GC_TFS), n_layers))
    layer_rank = np.zeros((len(GC_TFS), n_layers))
    for li in range(n_layers):
        emb_l = np.asarray(emb[li])
        centroid = emb_l[marker_idx].mean(axis=0)
        # Distance to centroid for every gene
        d_all = np.linalg.norm(emb_l - centroid, axis=1)
        for ti, g in enumerate(GC_TFS):
            if g not in n2i:
                layer_dist[ti, li] = float("nan")
                layer_rank[ti, li] = float("nan")
            else:
                gi = n2i[g]
                layer_dist[ti, li] = d_all[gi]
                layer_rank[ti, li] = (d_all < d_all[gi]).sum()  # rank 0 = closest

    # Across-pseudotime trajectory: mean expression
    print("[E8] Loading B-cell pseudotime + expression…", flush=True)
    bin_means, bin_assign, pt = _bcell_pseudotime_expression(BCELL_H5AD, GC_TFS, n_bins=10)
    if bin_means is None:
        print("  pseudotime data unavailable", flush=True)
        return

    print("[E8] Per-TF trajectory comparison:", flush=True)
    rows = []
    for ti, g in enumerate(GC_TFS):
        layers = np.arange(n_layers)
        bins = np.arange(bin_means.shape[0])
        rho_layer_dist = _spearman(layers, layer_dist[ti])
        rho_layer_rank = _spearman(layers, layer_rank[ti])
        rho_pt_expr = _spearman(bins, bin_means[:, ti])
        # Convergence in space ⇔ rho_layer_dist < 0
        # Expression rises ⇔ rho_pt_expr > 0
        # Dynamics support: both signs as expected (negative * positive = -)
        # sign of correlation between layer-trajectory and pt-trajectory:
        # use cross-correlation between -layer_dist (so larger=closer) and bin_means
        # Use Spearman.
        # We compute Spearman(–layer_dist[ti], interp(bin_means)) — but they have different lengths.
        # Use rank-based: sort layers by -layer_dist, sort pseudotime bins by bin_means,
        # then look at directional agreement. Simplest: report individually + interpret.
        rows.append(
            dict(
                tf=g,
                in_vocab=g in n2i,
                rho_layer_index_vs_dist=rho_layer_dist,
                rho_layer_index_vs_rank=rho_layer_rank,
                rho_pseudotime_bin_vs_expr=rho_pt_expr,
                supports_dynamics=(rho_layer_dist < 0 and rho_pt_expr > 0)
                or (rho_layer_dist > 0 and rho_pt_expr < 0),
                bin_means=bin_means[:, ti].tolist(),
                layer_dist=layer_dist[ti].tolist(),
            )
        )
        print(f"  {g}: rho_layer_dist={rho_layer_dist:.3f}; rho_pt_expr={rho_pt_expr:.3f}; supports_dynamics={(rho_layer_dist < 0 and rho_pt_expr > 0)}", flush=True)

    df = pd.DataFrame([{k: v for k, v in r.items() if k not in ("bin_means", "layer_dist")} for r in rows])
    df.to_csv(OUT / "results.csv", index=False)
    with open(OUT / "results.json", "w") as f:
        json.dump(rows, f, indent=2)
    print("\n[E8] Summary:")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
