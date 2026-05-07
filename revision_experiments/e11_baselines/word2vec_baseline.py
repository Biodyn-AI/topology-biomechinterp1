#!/usr/bin/env python3
"""E11 (Word2Vec baseline) — train Word2Vec on cell-tokenised expression.

For each cell, tokenise as a list of gene names ranked by expression magnitude
(top-K most-expressed genes, in order). Train Word2Vec (skip-gram, 256-dim) on
these "documents." Per-gene Word2Vec vectors then serve as a baseline gene
embedding that captures co-expression relationships from the same data scGPT
trained on, but without the transformer architecture.

Compare these Word2Vec gene embeddings against scGPT-L11 on the same headline
metrics (SV1 GO:0005615 enrichment, SV2 PPI co-pole, TF/target AUROC).

This addresses Reviewer 2 (#8): "compare to simpler baselines."
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

H5AD = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
TRRUST_PATH = REPO / "single_cell_mechinterp/external/networks/trrust_human.tsv"
STRING_CACHE = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015/string_ppi_score04_cache.json"
GENE2GO = REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"
GO_EXTRACELLULAR = "GO:0005615"


def _h5ad_genes(p):
    with h5py.File(p, "r") as f:
        raw = f["var"]["_index"][:]
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]


def _load_expression(p, max_cells=8000, rng_seed=42):
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
            return mat  # keep sparse
        return X[:].astype(np.float32)


def fisher_2x2(a, b, c, d):
    from math import lgamma, exp
    n = a + b + c + d
    r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    if r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0:
        return 1.0, 1.0
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    def lpx(x):
        return (lgamma(r1 + 1) - lgamma(x + 1) - lgamma(r1 - x + 1)
                + lgamma(r2 + 1) - lgamma(c1 - x + 1) - lgamma(c2 - r1 + x + 1)
                - lgamma(n + 1) + lgamma(c1 + 1) + lgamma(c2 + 1))
    obs = lpx(a)
    p = 0.0
    for x in range(max(0, c1 - r2), min(c1, r1) + 1):
        if lpx(x) <= obs + 1e-12:
            p += exp(lpx(x))
    return float(odds), float(min(p, 1.0))


def go_pole_enrichment(proj, gene_names, gene2go, go_term, pole_frac=0.27):
    n = len(gene_names)
    K = max(1, int(pole_frac * n))
    order = np.argsort(proj)
    in_top = set(order[-K:].tolist())
    in_bot = set(order[:K].tolist())
    a_top = sum(1 for i in in_top if go_term in gene2go.get(gene_names[i], set()))
    out_with = sum(1 for i in range(n) if i not in in_top and go_term in gene2go.get(gene_names[i], set()))
    or_top, p_top = fisher_2x2(a_top, K - a_top, out_with, (n - K) - out_with)
    a_bot = sum(1 for i in in_bot if go_term in gene2go.get(gene_names[i], set()))
    out_with2 = sum(1 for i in range(n) if i not in in_bot and go_term in gene2go.get(gene_names[i], set()))
    or_bot, p_bot = fisher_2x2(a_bot, K - a_bot, out_with2, (n - K) - out_with2)
    return max(or_top, or_bot), min(p_top, p_bot)


def copole_z(proj, gene_names, pairs, n_perm=100, pole_frac=0.27, rng=None):
    rng = rng or np.random.default_rng(42)
    n2i = {g: i for i, g in enumerate(gene_names)}
    n = len(gene_names)
    K = max(1, int(pole_frac * n))
    order = np.argsort(proj)
    top = set(order[-K:].tolist())
    bot = set(order[:K].tolist())
    co = total = 0
    for g1, g2 in pairs:
        if g1 not in n2i or g2 not in n2i:
            continue
        i, j = n2i[g1], n2i[g2]
        total += 1
        if (i in top and j in top) or (i in bot and j in bot):
            co += 1
    obs = co / total if total else 0.0
    null_rates = np.zeros(n_perm)
    perm_idx = np.arange(n)
    for k in range(n_perm):
        rng.shuffle(perm_idx)
        n2i_p = {gene_names[perm_idx[i]]: i for i in range(n)}
        co_n = total_n = 0
        for g1, g2 in pairs:
            if g1 not in n2i_p or g2 not in n2i_p:
                continue
            i, j = n2i_p[g1], n2i_p[g2]
            total_n += 1
            if (i in top and j in top) or (i in bot and j in bot):
                co_n += 1
        null_rates[k] = co_n / total_n if total_n else 0.0
    return dict(obs=obs, null_mean=float(null_rates.mean()), z=float((obs - null_rates.mean()) / max(null_rates.std(ddof=1), 1e-9)), n_pairs=total)


def tf_target_auroc(proj_subspace, gene_names, tfs, targets_only):
    name_to_i = {g: i for i, g in enumerate(gene_names)}
    tf_idx = [name_to_i[g] for g in tfs if g in name_to_i]
    target_idx = [name_to_i[g] for g in targets_only if g in name_to_i]
    if len(tf_idx) < 5 or len(target_idx) < 5:
        return float("nan"), len(tf_idx), len(target_idx)
    norm_proj = proj_subspace / (np.linalg.norm(proj_subspace, axis=1, keepdims=True) + 1e-12)
    sim = norm_proj @ norm_proj[tf_idx].T
    score = sim.mean(axis=1)
    pos, neg = score[tf_idx], score[target_idx]
    all_s = np.concatenate([pos, neg])
    labels = np.array([1] * len(pos) + [0] * len(neg))
    order = np.argsort(-all_s)
    cn = ar = 0
    for r in labels[order]:
        if r == 1:
            ar += cn
        else:
            cn += 1
    return float(1 - ar / (len(pos) * len(neg))), len(tf_idx), len(target_idx)


def main():
    print("[E11-w2v] Loading…", flush=True)
    gene_names = _h5ad_genes(H5AD)
    n_genes = len(gene_names)

    print("[E11-w2v] Loading expression…", flush=True)
    X_sparse = _load_expression(H5AD, max_cells=6000)
    n_cells = X_sparse.shape[0]
    print(f"  expr shape: {n_cells} x {X_sparse.shape[1]}", flush=True)

    # Tokenise: for each cell, top-K genes by expression
    print("[E11-w2v] Tokenising cells (top-200 genes per cell, ranked)…", flush=True)
    K = 200
    sentences = []
    if hasattr(X_sparse, "toarray"):
        # CSR format - process row by row to save memory
        for i in range(n_cells):
            row = X_sparse.getrow(i).toarray().flatten()
            top = np.argsort(-row)[:K]
            top = [t for t in top if row[t] > 0]
            if len(top) >= 5:
                sentences.append([gene_names[t] for t in top])
    else:
        for i in range(n_cells):
            row = X_sparse[i]
            top = np.argsort(-row)[:K]
            top = [t for t in top if row[t] > 0]
            if len(top) >= 5:
                sentences.append([gene_names[t] for t in top])
    print(f"  N sentences: {len(sentences)}", flush=True)

    # Train Word2Vec
    print("[E11-w2v] Training Word2Vec (256-dim, skip-gram, window=10)…", flush=True)
    try:
        from gensim.models import Word2Vec
    except ImportError:
        print("  gensim not available; falling back to PMI matrix", flush=True)
        return _pmi_baseline(sentences, gene_names, n_genes)

    model = Word2Vec(
        sentences=sentences,
        vector_size=256,
        window=10,
        min_count=2,
        sg=1,           # skip-gram
        workers=4,
        epochs=5,
        seed=42,
    )

    # Build per-gene embedding matrix in vocab order
    w2v_emb = np.zeros((n_genes, 256), dtype=np.float32)
    n_in = 0
    for i, g in enumerate(gene_names):
        if g in model.wv:
            w2v_emb[i] = model.wv[g]
            n_in += 1
    print(f"  Word2Vec coverage: {n_in}/{n_genes} genes", flush=True)

    # Annotations
    t = pd.read_csv(TRRUST_PATH, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    tfs = set(t.TF.astype(str))
    targets_only = set(t.target.astype(str)) - tfs
    obj = json.loads(STRING_CACHE.read_text())
    string_pairs = [(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= 0.7]
    with open(GENE2GO, "rb") as f:
        gene2go = pickle.load(f)

    # SVD on Word2Vec embeddings
    Xc = w2v_emb - w2v_emb.mean(axis=0, keepdims=True)
    print("[E11-w2v] SVD on Word2Vec embeddings…", flush=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    US = U * S
    sv1 = US[:, 0]
    sv2 = US[:, 1]
    sv27 = US[:, 1:7]

    rng = np.random.default_rng(42)
    or_extra, p_extra = go_pole_enrichment(sv1, gene_names, gene2go, GO_EXTRACELLULAR)
    ppi = copole_z(sv2, gene_names, string_pairs, n_perm=100, rng=rng)
    au, n_tf, n_tgt = tf_target_auroc(sv27, gene_names, tfs, targets_only)

    result = dict(
        embedding="word2vec_256_skipgram",
        n_sentences=len(sentences),
        n_in_vocab=n_in,
        or_GO0005615=or_extra,
        p_GO0005615=p_extra,
        ppi_obs_rate=ppi["obs"],
        ppi_null_rate=ppi["null_mean"],
        ppi_z=ppi["z"],
        ppi_n_pairs=ppi["n_pairs"],
        tf_target_auroc=au,
        n_tf=n_tf,
        n_target=n_tgt,
    )
    print("\n[E11-w2v] Result:")
    print(json.dumps(result, indent=2))
    with open(OUT / "word2vec_results.json", "w") as f:
        json.dump(result, f, indent=2)


def _pmi_baseline(sentences, gene_names, n_genes):
    """Fallback: PMI (pointwise mutual information) co-occurrence baseline if
    gensim is unavailable. Gene-gene co-occurrence within sentences, then SVD
    of PMI matrix."""
    print("[E11-w2v PMI fallback] Building co-occurrence…", flush=True)
    name_to_i = {g: i for i, g in enumerate(gene_names)}
    cooc = np.zeros((n_genes, n_genes), dtype=np.float32)
    for s in sentences:
        idx = [name_to_i.get(g, -1) for g in s]
        idx = [i for i in idx if i >= 0]
        for a in idx:
            for b in idx:
                if a != b:
                    cooc[a, b] += 1
    rowsum = cooc.sum(axis=1) + 1e-9
    total = cooc.sum() + 1e-9
    pmi = np.log((cooc * total) / (rowsum[:, None] * rowsum[None, :] + 1e-9) + 1e-9)
    pmi = np.maximum(pmi, 0)  # PPMI
    Xc = pmi - pmi.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return None


if __name__ == "__main__":
    main()
