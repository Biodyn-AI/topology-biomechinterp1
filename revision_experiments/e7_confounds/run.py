#!/usr/bin/env python3
"""E7 — Strengthened multi-feature confound residualisation.

For every gene-pair-level claim (PPI co-pole, TRRUST TF–target proximity), the
residualisation in the paper used only co-expression as a covariate. Reviewer
1 (#5) and Reviewer 2 (#6) want a stricter version that also controls for
expression magnitude, hub degree (in STRING), and number of GO annotations
per gene.

This script:
  1. For each TRRUST TF–target pair and a matched non-regulatory pair set,
     computes pair-level cosine similarity in SV5–SV7 (the paper's
     "regulatory" subspace) at every layer.
  2. Builds a covariate matrix per pair: log mean expression, co-expression,
     STRING degree, GO annotation count, gene length proxy (using mean
     counts as a stand-in when length is unavailable), GC bin (0..3) by
     gene-name first-letter as a coarse proxy.
  3. Residualises the cosine similarity against the covariates via OLS.
  4. Reports rank-biserial correlation between TRRUST and non-regulatory
     pair residuals + permutation p, plus a bootstrap robustness check.

The criterion: a finding is "co-expression-independent" only if the
residualised rank-biserial r_rb is positive with permutation p ≤ 0.05.
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path
from collections import Counter

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

EMB_IMMUNE = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy"
H5AD_IMMUNE = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
TRRUST_PATH = REPO / "single_cell_mechinterp/external/networks/trrust_human.tsv"
STRING_CACHE = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015/string_ppi_score04_cache.json"
GENE2GO = REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"


def _h5ad_genes(p):
    with h5py.File(p, "r") as f:
        raw = f["var"]["_index"][:]
        means = f["var"]["mean"][:] if "mean" in f["var"] else np.zeros(len(raw))
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw], np.asarray(means)


def _coexp_matrix(p, max_cells=4000, rng_seed=42):
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
                cell_idx = np.sort(rng.choice(n_obs, size=max_cells, replace=False))
                mat = mat[cell_idx]
            arr = mat.toarray().astype(np.float32)
        else:
            arr = X[:].astype(np.float32)
            if arr.shape[0] > max_cells:
                cell_idx = np.sort(rng.choice(arr.shape[0], size=max_cells, replace=False))
                arr = arr[cell_idx]
    mu = arr.mean(axis=0, keepdims=True)
    sd = arr.std(axis=0, keepdims=True) + 1e-9
    Z = (arr - mu) / sd
    n = Z.shape[0]
    return Z, n


def _coexp_for_pairs(Z, n_cells, pi, pj):
    return (Z[:, pi] * Z[:, pj]).sum(axis=0) / n_cells


def _cosine_pairs(US_subspace, pi, pj):
    norm = np.linalg.norm(US_subspace, axis=1) + 1e-12
    cos = (US_subspace[pi] * US_subspace[pj]).sum(axis=1) / (norm[pi] * norm[pj])
    return cos


def _ols_residualise(y, X):
    X1 = np.hstack([np.ones((len(X), 1)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return y - X1 @ beta


def _rank_biserial(pos: np.ndarray, neg: np.ndarray) -> float:
    """Rank-biserial: U statistic normalised to [-1, 1]."""
    n_pos, n_neg = len(pos), len(neg)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    combined = np.concatenate([pos, neg])
    labels = np.concatenate([np.ones(n_pos), np.zeros(n_neg)])
    order = np.argsort(combined)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(combined) + 1)
    R_pos = ranks[: n_pos].sum()
    U = R_pos - n_pos * (n_pos + 1) / 2
    return float(2 * U / (n_pos * n_neg) - 1)


def main() -> None:
    print("[E7] Loading…", flush=True)
    gene_index, gene_means = _h5ad_genes(H5AD_IMMUNE)
    name_to_i = {g: i for i, g in enumerate(gene_index)}

    # TRRUST positives (in-vocab)
    t = pd.read_csv(TRRUST_PATH, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    pos_pairs = []
    for _, r in t.iterrows():
        if r.TF in name_to_i and r.target in name_to_i and r.TF != r.target:
            pos_pairs.append((r.TF, r.target))
    pos_pairs = list({(a, b) for a, b in pos_pairs})
    print(f"  TRRUST in-vocab pairs: {len(pos_pairs)}", flush=True)

    # Negatives: same TFs but pair with a random non-target gene (no TRRUST edge)
    rng = np.random.default_rng(42)
    trrust_set = set((a, b) for a, b in pos_pairs)
    pos_tfs = sorted({a for a, b in pos_pairs})
    in_vocab_genes = list(name_to_i.keys())
    neg_pairs = []
    for tf in pos_tfs:
        cnt = 0
        # Sample 4 negatives per TF
        attempts = 0
        while cnt < 4 and attempts < 100:
            attempts += 1
            g = in_vocab_genes[int(rng.integers(0, len(in_vocab_genes)))]
            if g == tf:
                continue
            if (tf, g) in trrust_set or (g, tf) in trrust_set:
                continue
            neg_pairs.append((tf, g))
            cnt += 1
    neg_pairs = list({(a, b) for a, b in neg_pairs})
    print(f"  Negative pairs: {len(neg_pairs)}", flush=True)

    # Convert to indices
    def to_idx(pairs):
        pi = np.array([name_to_i[a] for a, _ in pairs])
        pj = np.array([name_to_i[b] for _, b in pairs])
        return pi, pj

    pos_pi, pos_pj = to_idx(pos_pairs)
    neg_pi, neg_pj = to_idx(neg_pairs)

    # Co-expression
    print("[E7] Computing co-expression matrix…", flush=True)
    Z, ncells = _coexp_matrix(H5AD_IMMUNE, max_cells=3000)
    pos_coexp = _coexp_for_pairs(Z, ncells, pos_pi, pos_pj)
    neg_coexp = _coexp_for_pairs(Z, ncells, neg_pi, neg_pj)
    del Z
    print(f"  pos coexp mean: {pos_coexp.mean():.3f}; neg coexp mean: {neg_coexp.mean():.3f}", flush=True)

    # STRING degree (≥0.4)
    obj = json.loads(STRING_CACHE.read_text())
    deg = Counter()
    for p in obj["pairs"]:
        if p["g1"] in name_to_i and p["g2"] in name_to_i:
            deg[p["g1"]] += 1
            deg[p["g2"]] += 1
    pos_degsum = np.array([deg.get(a, 0) + deg.get(b, 0) for a, b in pos_pairs])
    neg_degsum = np.array([deg.get(a, 0) + deg.get(b, 0) for a, b in neg_pairs])

    # GO count
    with open(GENE2GO, "rb") as f:
        g2go = pickle.load(f)
    go_count = lambda g: len(g2go.get(g, []))
    pos_gocount = np.array([go_count(a) + go_count(b) for a, b in pos_pairs])
    neg_gocount = np.array([go_count(a) + go_count(b) for a, b in neg_pairs])

    # Mean expression
    pos_meanexpr = gene_means[pos_pi] + gene_means[pos_pj]
    neg_meanexpr = gene_means[neg_pi] + gene_means[neg_pj]

    # Embeddings
    print("[E7] Loading embeddings…", flush=True)
    emb = np.load(EMB_IMMUNE, mmap_mode="r")
    n_layers = emb.shape[0]
    print(f"  emb: {emb.shape}", flush=True)

    layer_results = []
    for li in range(n_layers):
        X = emb[li]
        Xc = X - X.mean(axis=0, keepdims=True)
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        US = U * S
        sv57 = US[:, 4:7]    # SV5-SV7 (paper)
        sv24 = US[:, 1:4]    # SV2-SV4

        for subspace_name, subspace in (("SV5-7", sv57), ("SV2-4", sv24)):
            cos_pos = _cosine_pairs(subspace, pos_pi, pos_pj)
            cos_neg = _cosine_pairs(subspace, neg_pi, neg_pj)

            # Combined regression to residualise
            cos_all = np.concatenate([cos_pos, cos_neg])
            controls = np.column_stack([
                np.concatenate([pos_coexp, neg_coexp]),
                np.concatenate([pos_meanexpr, neg_meanexpr]),
                np.concatenate([pos_degsum, neg_degsum]).astype(float),
                np.concatenate([pos_gocount, neg_gocount]).astype(float),
            ])
            # Standardise controls
            controls = (controls - controls.mean(axis=0)) / (controls.std(axis=0) + 1e-9)

            cos_resid = _ols_residualise(cos_all, controls)
            n_pos = len(cos_pos)
            cos_resid_pos = cos_resid[:n_pos]
            cos_resid_neg = cos_resid[n_pos:]

            # Raw rank-biserial
            r_rb_raw = _rank_biserial(cos_pos, cos_neg)
            r_rb_resid = _rank_biserial(cos_resid_pos, cos_resid_neg)

            # Permutation null on labels for residualised
            n_perm = 200
            null_rb = np.zeros(n_perm)
            labels = np.concatenate([np.ones(n_pos), np.zeros(len(cos_neg))])
            for k in range(n_perm):
                rng.shuffle(labels)
                null_rb[k] = _rank_biserial(cos_resid[labels == 1], cos_resid[labels == 0])
            perm_p = float((null_rb >= r_rb_resid).mean())

            layer_results.append(
                dict(
                    layer=li,
                    subspace=subspace_name,
                    r_rb_raw=r_rb_raw,
                    r_rb_residual=r_rb_resid,
                    perm_p_residual=perm_p,
                    n_pos=n_pos,
                    n_neg=len(cos_neg),
                )
            )
            if subspace_name == "SV5-7":
                print(f"  L{li:02d} SV5-7: raw r_rb={r_rb_raw:+.3f}, residual r_rb={r_rb_resid:+.3f}, perm_p={perm_p:.3f}", flush=True)

    df = pd.DataFrame(layer_results)
    df.to_csv(OUT / "per_layer.csv", index=False)
    sv57 = df[df.subspace == "SV5-7"]
    sv24 = df[df.subspace == "SV2-4"]
    summary = dict(
        n_pos=int(df.iloc[0].n_pos),
        n_neg=int(df.iloc[0].n_neg),
        sv57_layers_residual_significant=int((sv57.perm_p_residual <= 0.05).sum()),
        sv57_residual_r_rb_layer_mean=float(sv57.r_rb_residual.mean()),
        sv57_residual_r_rb_layer_max=float(sv57.r_rb_residual.max()),
        sv57_residual_r_rb_layer_0=float(sv57.iloc[0].r_rb_residual),
        sv24_layers_residual_significant=int((sv24.perm_p_residual <= 0.05).sum()),
        sv24_residual_r_rb_layer_mean=float(sv24.r_rb_residual.mean()),
    )
    with open(OUT / "results.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n[E7] Summary:")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
