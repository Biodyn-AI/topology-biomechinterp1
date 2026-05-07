#!/usr/bin/env python3
"""E13 — Drug-target candidate ranking from PPI subspace.

For ~30 well-characterised drug targets in the curated 209-gene set, rank
the in-vocabulary genes by cosine similarity in SV2-SV4 (the PPI subspace
identified by the paper). Compare precision@5 and precision@10 against
known STRING ≥ 0.7 interactors:
  (a) PPI subspace ranking (the paper's method).
  (b) Co-expression Pearson ranking.
  (c) Random ranking.
  (d) Full 512-dim cosine similarity (no projection).

Drug targets curated from canonical immune-relevant clinical/oncology pathways.
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
H5AD = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
STRING_CACHE = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015/string_ppi_score04_cache.json"

# Canonical immune-relevant drug targets in the curated 209-gene set.
DRUG_TARGETS = [
    "JAK1", "JAK2", "STAT3", "STAT1", "BCL2", "BCL2L1",
    "EGFR", "MTOR", "MYC", "TP53",
    "NFKB1", "RELA", "BRD4", "BTK",
    "PIK3CA", "PIK3CB", "AKT1",
    "TNFRSF17",  # BCMA
    "CXCR4", "CCR5", "ICAM1", "ITGAL",
    "IRF4", "IRF8", "PRDM1",
    "MEF2C",
    "STAT6", "IL10",
]


def _h5ad_index(p):
    with h5py.File(p, "r") as f:
        raw = f["var"]["_index"][:]
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]


def _load_expression(p, max_cells=4000, rng_seed=42):
    rng = np.random.default_rng(rng_seed)
    with h5py.File(p, "r") as f:
        X = f["X"]
        if isinstance(X, h5py.Group):
            from scipy.sparse import csr_matrix
            data = X["data"][:]
            indices = X["indices"][:]
            indptr = X["indptr"][:]
            n_obs = len(indptr) - 1
            n_var = int(np.max(indices) + 1) if len(indices) else 0
            mat = csr_matrix((data, indices, indptr), shape=(n_obs, n_var))
            if n_obs > max_cells:
                idx = np.sort(rng.choice(n_obs, size=max_cells, replace=False))
                mat = mat[idx]
            arr = mat.toarray().astype(np.float32)
        else:
            arr = X[:].astype(np.float32)
            if arr.shape[0] > max_cells:
                idx = np.sort(rng.choice(arr.shape[0], size=max_cells, replace=False))
                arr = arr[idx]
    return arr


def precision_at_k(ranked_genes, true_set, k):
    top_k = ranked_genes[:k]
    return sum(1 for g in top_k if g in true_set) / k


def main():
    print("[E13] Loading…", flush=True)
    gene_index = _h5ad_index(H5AD)
    n2i = {g: i for i, g in enumerate(gene_index)}

    obj = json.loads(STRING_CACHE.read_text())
    string_pairs = [(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= 0.7]
    # Build per-gene true-interactor set
    interactors = {}
    for g1, g2 in string_pairs:
        interactors.setdefault(g1, set()).add(g2)
        interactors.setdefault(g2, set()).add(g1)
    print(f"  STRING ≥ 0.7 pairs: {len(string_pairs)}", flush=True)

    # Filter drug targets to in-vocab + in-STRING
    targets_in_vocab = [g for g in DRUG_TARGETS if g in n2i and g in interactors]
    print(f"  Drug targets in vocab + STRING: {len(targets_in_vocab)}", flush=True)

    # Load scGPT L11 (PPI subspace is at L11 per paper)
    emb = np.load(EMB, mmap_mode="r")
    X = np.asarray(emb[11])
    Xc = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    US = U * S
    sv24 = US[:, 1:4]  # PPI subspace per paper
    sv24_norm = sv24 / (np.linalg.norm(sv24, axis=1, keepdims=True) + 1e-12)
    full_norm = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)

    # Co-expression
    expr = _load_expression(H5AD, max_cells=3000)
    expr_z = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-9)
    n_cells = expr_z.shape[0]

    rng = np.random.default_rng(42)
    rows = []
    for tg in targets_in_vocab:
        ti = n2i[tg]
        true_int = interactors[tg]

        # Rank by PPI subspace (SV2-4)
        sims_sv24 = sv24_norm @ sv24_norm[ti]
        sims_sv24[ti] = -np.inf
        order_sv24 = np.argsort(-sims_sv24)
        ranked_sv24 = [gene_index[i] for i in order_sv24]

        # Rank by full 512-dim cosine
        sims_full = full_norm @ full_norm[ti]
        sims_full[ti] = -np.inf
        order_full = np.argsort(-sims_full)
        ranked_full = [gene_index[i] for i in order_full]

        # Rank by co-expression
        coexp_to_tg = (expr_z * expr_z[:, ti:ti+1]).sum(axis=0) / n_cells
        coexp_to_tg[ti] = -np.inf
        order_coexp = np.argsort(-coexp_to_tg)
        ranked_coexp = [gene_index[i] for i in order_coexp]

        # Random
        rand_perm = rng.permutation(len(gene_index))
        rand_perm = rand_perm[rand_perm != ti]
        ranked_rand = [gene_index[i] for i in rand_perm]

        for method, ranked in [
            ("ppi_sv2-4", ranked_sv24),
            ("scgpt_full512", ranked_full),
            ("coexpression", ranked_coexp),
            ("random", ranked_rand),
        ]:
            rows.append(
                dict(
                    target=tg,
                    method=method,
                    n_known_interactors=len(true_int),
                    p_at_5=precision_at_k(ranked, true_int, 5),
                    p_at_10=precision_at_k(ranked, true_int, 10),
                    p_at_25=precision_at_k(ranked, true_int, 25),
                )
            )

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_target.csv", index=False)

    # Aggregate by method
    agg = df.groupby("method").agg(
        n_targets=("target", "count"),
        mean_p_at_5=("p_at_5", "mean"),
        mean_p_at_10=("p_at_10", "mean"),
        mean_p_at_25=("p_at_25", "mean"),
        median_n_interactors=("n_known_interactors", "median"),
    ).reset_index()
    agg.to_csv(OUT / "results.csv", index=False)
    print("\n[E13] Aggregated results:")
    print(agg.to_string(index=False))

    summary = dict(
        n_drug_targets=len(targets_in_vocab),
        per_method=agg.to_dict(orient="records"),
    )
    with open(OUT / "results.json", "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
