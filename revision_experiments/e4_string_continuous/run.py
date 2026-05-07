#!/usr/bin/env python3
"""E4 — Continuous PPI regression (replaces n=5 quintile headline).

For all in-vocabulary STRING pairs (≥0.4, N≈3,092), at every layer:
  - Compute pair-level cosine similarity in SV2 and in the SV2-7 subspace.
  - Compute pair-level co-expression Pearson correlation across cells.
  - Regress cosine similarity on STRING combined score, controlling for
    co-expression, mean expression of the two genes, and the sum of node
    degrees (STRING degree).
  - Report β (standardised coefficient), R², and partial Pearson r with 95%
    bootstrap CI.

The brittle ρ=1.000-on-five-quintiles statistic moves to a sanity check.
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

EMB_IMMUNE = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy"
H5AD_IMMUNE = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
STRING_CACHE = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015/string_ppi_score04_cache.json"


def _h5ad_genes_and_means(p: Path) -> tuple[list[str], np.ndarray]:
    with h5py.File(p, "r") as f:
        raw = f["var"]["_index"][:]
        means = f["var"]["mean"][:] if "mean" in f["var"] else np.array([])
    gene_index = [g.decode() if isinstance(g, bytes) else str(g) for g in raw]
    return gene_index, np.asarray(means)


def _coexpression_for_pairs(h5_path: Path, pairs_idx: np.ndarray, max_cells: int = 8000, rng_seed: int = 42) -> np.ndarray:
    """Return Pearson correlation of expression across cells for each pair.

    Uses a random subsample of cells to keep memory bounded.
    """
    rng = np.random.default_rng(rng_seed)
    with h5py.File(h5_path, "r") as f:
        # Read X (the count matrix)
        X = f["X"]
        if isinstance(X, h5py.Group):
            # CSR
            data = X["data"][:]
            indices = X["indices"][:]
            indptr = X["indptr"][:]
            n_obs = len(indptr) - 1
            from scipy.sparse import csr_matrix
            mat = csr_matrix((data, indices, indptr), shape=(n_obs, len(pairs_idx[0]) if False else 4803))
            # Subsample cells
            if n_obs > max_cells:
                cell_idx = np.sort(rng.choice(n_obs, size=max_cells, replace=False))
                mat = mat[cell_idx]
            mat_dense = mat.toarray().astype(np.float32)
        else:
            mat_dense = X[:].astype(np.float32)
            if mat_dense.shape[0] > max_cells:
                cell_idx = np.sort(rng.choice(mat_dense.shape[0], size=max_cells, replace=False))
                mat_dense = mat_dense[cell_idx]
    # Standardise columns (genes)
    mu = mat_dense.mean(axis=0, keepdims=True)
    sd = mat_dense.std(axis=0, keepdims=True) + 1e-9
    Z = (mat_dense - mu) / sd
    # Pearson correlation = inner product / n_cells
    corrs = np.zeros(len(pairs_idx[0]))
    n = Z.shape[0]
    for k, (i, j) in enumerate(zip(pairs_idx[0], pairs_idx[1])):
        corrs[k] = float((Z[:, i] * Z[:, j]).sum() / n)
    return corrs


def _bootstrap_ci(x: np.ndarray, y: np.ndarray, fn, n_boot: int = 500, alpha: float = 0.05, rng_seed: int = 42) -> tuple[float, float, float]:
    rng = np.random.default_rng(rng_seed)
    n = len(x)
    boot = np.zeros(n_boot)
    for k in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boot[k] = fn(x[idx], y[idx])
    obs = fn(x, y)
    lo, hi = np.quantile(boot, [alpha / 2, 1 - alpha / 2])
    return float(obs), float(lo), float(hi)


def _pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xm = x - x.mean()
    ym = y - y.mean()
    denom = float(np.sqrt((xm * xm).sum() * (ym * ym).sum()) + 1e-12)
    return float((xm * ym).sum() / denom)


def _ols_residualise(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Return residuals of OLS y ~ X (intercept added)."""
    X1 = np.hstack([np.ones((len(X), 1)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return y - X1 @ beta


def main() -> None:
    print("[E4] Loading…")
    obj = json.loads(STRING_CACHE.read_text())
    pairs = obj["pairs"]
    gene_index, gene_means = _h5ad_genes_and_means(H5AD_IMMUNE)
    name_to_i = {g: i for i, g in enumerate(gene_index)}

    # Filter pairs to in-vocab and resolve indices
    rows = []
    pair_ids = []
    for p in pairs:
        i = name_to_i.get(p["g1"])
        j = name_to_i.get(p["g2"])
        if i is None or j is None or i == j:
            continue
        rows.append((i, j, p["score"]))
        pair_ids.append((p["g1"], p["g2"]))
    pi = np.array([r[0] for r in rows])
    pj = np.array([r[1] for r in rows])
    score = np.array([r[2] for r in rows])
    print(f"  N pairs in-vocab: {len(rows)}; score range [{score.min():.2f}, {score.max():.2f}]")

    # Compute STRING-implied node degree at the >=0.4 cutoff
    from collections import Counter
    deg = Counter()
    for r in rows:
        deg[r[0]] += 1
        deg[r[1]] += 1
    pair_degsum = np.array([deg[i] + deg[j] for i, j in zip(pi, pj)])
    pair_meanexpr = (gene_means[pi] + gene_means[pj]) if len(gene_means) else np.zeros(len(rows))

    # Co-expression for these pairs
    print("[E4] Computing co-expression for STRING pairs…")
    coexp = _coexpression_for_pairs(H5AD_IMMUNE, (pi, pj), max_cells=4000)
    print(f"  coexp range: [{coexp.min():.3f}, {coexp.max():.3f}]; mean={coexp.mean():.3f}")

    # Load embeddings
    print("[E4] Loading scGPT immune embeddings…")
    emb = np.load(EMB_IMMUNE, mmap_mode="r")  # (12, 4803, 512)
    n_layers = emb.shape[0]
    print(f"  emb: {emb.shape}")

    layer_results = []
    for li in range(n_layers):
        X = emb[li]
        # SVD; project to SV2 and SV2-7
        Xc = X - X.mean(axis=0, keepdims=True)
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        US = U * S
        # Pair cosine similarity in SV2 (1D) and SV2-7 (6D)
        sv2 = US[:, 1]  # 1D
        sv27 = US[:, 1:7]  # 6D

        # 1-D cosine on SV2 = sign(sv2_i * sv2_j) * |sv2_i*sv2_j| / (|sv2_i||sv2_j|) = sign agreement
        # Use signed product for SV2 (since 1-D cosine is just sign)
        cos_sv2 = (sv2[pi] * sv2[pj]) / (np.abs(sv2[pi]) * np.abs(sv2[pj]) + 1e-12)
        # SV2-7 cosine
        norm = np.linalg.norm(sv27, axis=1) + 1e-12
        cos_sv27 = (sv27[pi] * sv27[pj]).sum(axis=1) / (norm[pi] * norm[pj])

        # OLS regression of cos_sv27 on STRING score, controlling for coexp, degsum, meanexpr
        Xc_design = np.column_stack([
            score,
            coexp,
            pair_degsum.astype(float),
            pair_meanexpr,
        ])
        # Standardise predictors
        Xc_std = (Xc_design - Xc_design.mean(axis=0)) / (Xc_design.std(axis=0) + 1e-9)
        # OLS y = X b + e
        X1 = np.hstack([np.ones((len(Xc_std), 1)), Xc_std])
        beta, *_ = np.linalg.lstsq(X1, cos_sv27, rcond=None)
        # R² for full model
        yhat = X1 @ beta
        ss_res = float(((cos_sv27 - yhat) ** 2).sum())
        ss_tot = float(((cos_sv27 - cos_sv27.mean()) ** 2).sum() + 1e-12)
        r2_full = 1 - ss_res / ss_tot
        # Partial Pearson of (cos_sv27 residualised on coexp+deg+meanexpr) vs (score residualised on same)
        controls = Xc_std[:, 1:]
        cos_resid = _ols_residualise(cos_sv27, controls)
        score_resid = _ols_residualise(score, controls)
        partial_r = _pearson(score_resid, cos_resid)
        # Bootstrap CI for partial r
        obs, lo, hi = _bootstrap_ci(score_resid, cos_resid, _pearson, n_boot=100, rng_seed=42 + li)

        # SV2 alone (no controls): just Pearson(cos_sv2, score)
        r_sv2 = _pearson(cos_sv2, score)
        # And residualised SV2
        sv2_resid = _ols_residualise(cos_sv2, controls)
        r_sv2_resid = _pearson(score_resid, sv2_resid)

        layer_results.append(
            dict(
                layer=li,
                n_pairs=len(rows),
                pearson_score_cos_sv2_raw=r_sv2,
                pearson_score_cos_sv2_residual=r_sv2_resid,
                pearson_score_cos_sv27_raw=_pearson(score, cos_sv27),
                partial_pearson_score_cos_sv27=partial_r,
                partial_pearson_lo95=lo,
                partial_pearson_hi95=hi,
                r2_full_sv27_model=r2_full,
                beta_score_standardised=float(beta[1]),
                beta_coexp_standardised=float(beta[2]),
                beta_degsum_standardised=float(beta[3]),
                beta_meanexpr_standardised=float(beta[4]),
            )
        )
        print(f"  L{li:02d}: r(score, cos_sv2)={r_sv2:.3f}; partial r(score, cos_sv27 | coexp,deg,meanexpr)={partial_r:.3f} [{lo:.3f}, {hi:.3f}]; β_score={float(beta[1]):.3f}, β_coexp={float(beta[2]):.3f}")

    df = pd.DataFrame(layer_results)
    df.to_csv(OUT / "per_layer.csv", index=False)

    summary = dict(
        n_pairs=len(rows),
        partial_pearson_score_cos_sv27_layer_mean=float(df.partial_pearson_score_cos_sv27.mean()),
        partial_pearson_score_cos_sv27_layer_max=float(df.partial_pearson_score_cos_sv27.max()),
        partial_pearson_score_cos_sv27_layer_min=float(df.partial_pearson_score_cos_sv27.min()),
        beta_score_layer_mean=float(df.beta_score_standardised.mean()),
        beta_coexp_layer_mean=float(df.beta_coexp_standardised.mean()),
        all_layers_partial_pos_lo95=int((df.partial_pearson_lo95 > 0).sum()),
    )
    with open(OUT / "results.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n[E4] Summary:")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
