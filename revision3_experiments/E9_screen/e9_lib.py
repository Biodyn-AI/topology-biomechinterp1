"""Shared loaders for E9 (screen provenance). Read-only access to the screen's inputs.

Two gene->row mappings were used by the screen on cycle1_main (= kidney_A):
  * TRUE      : edge-dataset source_idx/target_idx (iters 5-18), or gene_list.txt read WITH its blank lines
                (iters 19-26, 28). gene_list.txt line i == h5ad var index i, so both are correct.
  * MISALIGNED: gene_list.txt read with blank lines DROPPED and then enumerated (iters 4 [GO-CC part], 27,
                29-40): the k-th named gene (in var order) is assigned matrix row k, i.e. the 209 labels sit
                on rows 0..208. The 14 'OOV' genes of iter_0030 are exactly the labels that land on
                all-zero rows under this mapping.
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import contexts as X          # noqa: E402
import stats3 as S            # noqa: E402

P41 = R3.parent
ITER = P41 / "iterations"
CYC1 = X.S38 / "cycle1_main"
CYC4 = X.S38 / "cycle4_immune_main"


def edge_dataset(run_dir):
    rows = list(csv.DictReader(open(run_dir / "cycle1_edge_dataset.tsv"), delimiter="\t"))
    for r in rows:
        r["label"] = int(r["label"]); r["source_idx"] = int(r["source_idx"]); r["target_idx"] = int(r["target_idx"])
    return rows


def named_genes_true(run_dir=CYC1):
    """{gene: true row} for the edge-dataset genes (209 for cycle1_main, 361 for cycle4_immune_main)."""
    g2i = {}
    for r in edge_dataset(run_dir):
        g2i.setdefault(r["source"], r["source_idx"]); g2i.setdefault(r["target"], r["target_idx"])
    return g2i


def named_genes_misaligned():
    """{gene: row} as produced by `[l.strip() for l in f if l.strip()]` + enumerate on gene_list.txt."""
    vocab = [l.strip() for l in open(CYC1 / "gene_list.txt") if l.strip()]
    g2row = {g: i for i, g in enumerate(vocab)}
    named = set(named_genes_true())
    return {g: g2row[g] for g in sorted(named) if g in g2row}


def nonzero_rows(emb):
    return np.linalg.norm(np.asarray(emb[0], dtype=np.float64), axis=1) > 0


def svd_rows(M):
    """Mean-centred thin SVD of the given rows (float64)."""
    Xc = np.asarray(M, dtype=np.float64)
    Xc = Xc - Xc.mean(axis=0, keepdims=True)
    return np.linalg.svd(Xc, full_matrices=False)


def er_sq(S):
    """Effective rank with p = s^2/sum s^2 (the screen's and the preprint's stated convention)."""
    p = S ** 2 / (S ** 2).sum(); p = p[p > 0]
    return float(np.exp(-(p * np.log(p)).sum()))


def pr(S):
    v = S ** 2
    return float(v.sum() ** 2 / (v ** 2).sum())


def principal_cosines(A, B):
    """Cosines of the principal angles between row-spaces of A (k x d) and B (m x d) (orthonormal rows)."""
    return np.clip(np.linalg.svd(A @ B.T, compute_uv=False), 0, 1)


def trrust_raw_rows():
    out = []
    for line in open(X.NET / "trrust_human.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3: out.append((p[0].strip(), p[1].strip(), p[2].strip()))
    return out
