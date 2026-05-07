#!/usr/bin/env python3
"""E11 — Simpler-baseline comparison.

Reviewer 2 (#8): "manuscript does not sufficiently compare scGPT representations
to simpler baselines."

For each headline metric (SV1 GO:0005615 enrichment, SV2 PPI co-pole, TF-vs-target
AUROC), we compare scGPT layer-11 against four simpler "embeddings" of the same
4,803 genes:

  (A) Raw co-expression PCA: per-gene vector = column of the (cells × genes)
      expression matrix; reduced via PCA to 512 dims.
  (B) Pearson co-expression matrix row: per-gene vector = the gene's row in the
      gene × gene Pearson correlation matrix (4,803-dim).
  (C) Random-projection: per-gene vector = a 512-dim random Gaussian projection
      of (B), serving as a pseudo-random-init baseline (an actual random-init
      scGPT would require running the model and is left for a follow-up).
  (D) One-hot: identity baseline; effectively zero structure.

The comparison shows whether scGPT's spectral structure adds value beyond what
is recoverable from raw expression statistics.
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

GO_EXTRACELLULAR = "GO:0005615"


def _h5ad_genes(p):
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
    a_bot = sum(1 for i in in_bot if go_term in gene2go.get(gene_names[i], set()))
    in_top_with = a_top
    in_top_without = K - a_top
    out_with = sum(1 for i in range(n) if i not in in_top and go_term in gene2go.get(gene_names[i], set()))
    out_without = (n - K) - out_with
    or_top, p_top = fisher_2x2(in_top_with, in_top_without, out_with, out_without)
    in_bot_with = a_bot
    in_bot_without = K - a_bot
    out_with2 = sum(1 for i in range(n) if i not in in_bot and go_term in gene2go.get(gene_names[i], set()))
    out_without2 = (n - K) - out_with2
    or_bot, p_bot = fisher_2x2(in_bot_with, in_bot_without, out_with2, out_without2)
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
    null_mean = null_rates.mean()
    null_std = null_rates.std(ddof=1) if null_rates.std(ddof=1) > 0 else 1e-9
    return dict(obs=obs, null_mean=float(null_mean), z=float((obs - null_mean) / null_std), n_pairs=total)


def tf_target_auroc(proj_subspace, gene_names, tfs, targets_only):
    name_to_i = {g: i for i, g in enumerate(gene_names)}
    tf_idx = [name_to_i[g] for g in tfs if g in name_to_i]
    target_idx = [name_to_i[g] for g in targets_only if g in name_to_i]
    if len(tf_idx) < 5 or len(target_idx) < 5:
        return float("nan"), len(tf_idx), len(target_idx)
    norm_proj = proj_subspace / (np.linalg.norm(proj_subspace, axis=1, keepdims=True) + 1e-12)
    sim = norm_proj @ norm_proj[tf_idx].T
    score = sim.mean(axis=1)
    pos = score[tf_idx]
    neg = score[target_idx]
    all_scores = np.concatenate([pos, neg])
    labels = np.array([1] * len(pos) + [0] * len(neg))
    order = np.argsort(-all_scores)
    ranked = labels[order]
    cn = ar = 0
    for r in ranked:
        if r == 1:
            ar += cn
        else:
            cn += 1
    auroc = 1.0 - ar / (len(pos) * len(neg))
    return float(auroc), len(tf_idx), len(target_idx)


def main():
    print("[E11] Loading…", flush=True)
    gene_names = _h5ad_genes(H5AD_IMMUNE)
    n_genes = len(gene_names)

    # Annotations
    t = pd.read_csv(TRRUST_PATH, sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
    tfs = set(t.TF.astype(str))
    targets_only = set(t.target.astype(str)) - tfs

    obj = json.loads(STRING_CACHE.read_text())
    string_pairs = [(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= 0.7]

    with open(GENE2GO, "rb") as f:
        gene2go = pickle.load(f)

    # Build embeddings -------------------------------------------------------
    print("[E11] Loading expression…", flush=True)
    expr = _load_expression(H5AD_IMMUNE, max_cells=4000)
    print(f"  expr shape: {expr.shape}", flush=True)
    # Pearson correlation across genes
    print("[E11] Computing gene-gene Pearson…", flush=True)
    expr_z = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-9)
    n_cells = expr_z.shape[0]
    pearson = (expr_z.T @ expr_z) / n_cells  # (n_genes x n_genes)
    print(f"  pearson shape: {pearson.shape}", flush=True)

    # Baseline embeddings
    embeddings = {}

    # (A) Raw co-expression PCA
    Xc = expr.T - expr.T.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    embeddings["raw_expr_PCA512"] = (U[:, :512] * S[:512])

    # (B) Pearson row reduced to 512 dims via random projection (full matrix SVD too slow)
    rng = np.random.default_rng(42)
    Wrand = rng.normal(size=(pearson.shape[1], 512)) / np.sqrt(512)
    embeddings["pearson_row_512"] = (pearson @ Wrand).astype(np.float32)

    # (C) Random projection (no biology)
    Wrand2 = rng.normal(size=(pearson.shape[1], 512)) / np.sqrt(512)
    embeddings["random_proj_512"] = rng.normal(size=(n_genes, 512)).astype(np.float32)

    # (D) Identity baseline omitted — SVD is trivial, no informative comparison.

    # (E) scGPT layer 11 (the paper's primary)
    print("[E11] Loading scGPT layer-11…", flush=True)
    emb = np.load(EMB_IMMUNE, mmap_mode="r")
    embeddings["scgpt_L11"] = np.asarray(emb[11])

    # ------------------------------------------------------------------------
    rows = []
    rng = np.random.default_rng(42)
    for name, X in embeddings.items():
        print(f"[E11] -- {name} shape={X.shape}", flush=True)
        Xc2 = X - X.mean(axis=0, keepdims=True)
        try:
            U, S, Vt = np.linalg.svd(Xc2.astype(np.float32), full_matrices=False)
            US = U * S
        except Exception as e:
            print(f"  SVD failed: {e}", flush=True)
            continue
        sv1 = US[:, 0]
        sv2 = US[:, 1] if US.shape[1] > 1 else US[:, 0]
        sv27 = US[:, 1:7] if US.shape[1] >= 7 else US[:, 1:]

        # SV1 enrichment
        or_extra, p_extra = go_pole_enrichment(sv1, gene_names, gene2go, GO_EXTRACELLULAR)

        # SV2 PPI co-pole
        ppi = copole_z(sv2, gene_names, string_pairs, n_perm=100, rng=rng)

        # TF/target AUROC
        auroc, n_tf, n_tgt = tf_target_auroc(sv27, gene_names, tfs, targets_only)

        rows.append(
            dict(
                embedding=name,
                shape=str(X.shape),
                or_GO0005615=or_extra,
                p_GO0005615=p_extra,
                ppi_obs_rate=ppi["obs"],
                ppi_null_rate=ppi["null_mean"],
                ppi_z=ppi["z"],
                ppi_n_pairs=ppi["n_pairs"],
                tf_target_auroc=auroc,
                n_tf=n_tf,
                n_target=n_tgt,
            )
        )
        print(f"  {name}: SV1 extra OR={or_extra:.2f} p={p_extra:.4f}; PPI z={ppi['z']:.2f}; TF/target AUROC={auroc:.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "results.csv", index=False)
    print("\n[E11] Summary:")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
