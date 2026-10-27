"""Shared data layer and estimators for the round-2 revision experiments.

Every estimator here is a faithful re-implementation of the round-1 code in
`revision_experiments/*/run.py`, so that re-run numbers are directly comparable
to the ones currently in the manuscript. The only intended differences are:

  * permutation counts are raised to N >= 10,000 (round 1 used 50-1,000), and
  * empirical p-values use the (b+1)/(N+1) estimator, never a bare 0.

Canonical data layer (all paths verified present):
  embeddings  cycle1_main/layer_gene_embeddings.npy      (12, 4803, 512) float32
  gene axis   tabula_sapiens_processed.h5ad var index    4803 genes
  TRRUST      external/networks/trrust_human.tsv
  STRING      iter_0015/string_ppi_score04_cache.json
  GO          data/perturb/gene2go_all.pkl
"""
from __future__ import annotations

import json
import pickle
from math import exp, lgamma
from pathlib import Path

import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent

EMB_IMMUNE = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/layer_gene_embeddings.npy"
EMB_SEED43 = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_seed43/layer_gene_embeddings.npy"
EMB_SEED44 = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_seed44/layer_gene_embeddings.npy"
EMB_KIDNEY = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0067/kidney_scgpt_run/layer_gene_embeddings.npy"
H5AD = REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"
H5AD_KIDNEY = REPO / "single_cell_mechinterp/outputs/invariant_causal_edges/kidney/processed.h5ad"
TRRUST_PATH = REPO / "single_cell_mechinterp/external/networks/trrust_human.tsv"
STRING_CACHE = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0015/string_ppi_score04_cache.json"
GENE2GO = REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"
ESM2_NPZ = REPO / "subproject_60_uce_esm2_control_gpl_replication/embeddings/uce_esm2_human_gene_embeddings.npz"
GF_V1 = REPO / "subproject_55_geneformer_v1_10m_gpl_replication/embeddings/geneformer_v1_10m_gene_embeddings.npz"
GF_V2 = REPO / "subproject_56_geneformer_v2_104m_gpl_replication/embeddings/geneformer_v2_104m_gene_embeddings.npz"

GO_EXTRACELLULAR = "GO:0005615"
GO_CYTOSOL = "GO:0005829"
POLE_FRAC = 0.27          # round-1 constant: K = 0.27 n  (n=195 -> K=52, as in Methods)
BONFERRONI_ALPHA = 0.05 / 183


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def h5ad_genes(path=H5AD) -> list[str]:
    with h5py.File(path, "r") as f:
        raw = f["var"]["_index"][:]
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]


def load_emb(path=EMB_IMMUNE) -> np.ndarray:
    return np.load(path, mmap_mode="r")


def load_trrust(path=TRRUST_PATH):
    """Return (pairs, tf_set, target_only_set, signed_pairs)."""
    pairs, tfs, targets, signed = [], set(), set(), []
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        tf, tgt, mode = p[0].strip(), p[1].strip(), p[2].strip()
        if not tf or not tgt:
            continue
        pairs.append((tf, tgt))
        signed.append((tf, tgt, mode))
        tfs.add(tf)
        targets.add(tgt)
    return pairs, tfs, targets - tfs, signed


def load_string(path=STRING_CACHE, min_score=0.7):
    """Return list of (g1, g2, score) at or above min_score."""
    obj = json.loads(Path(path).read_text())
    rows = obj["pairs"] if isinstance(obj, dict) and "pairs" in obj else obj
    out = []
    for r in rows:
        if isinstance(r, dict):
            g1, g2, s = r.get("a") or r.get("g1"), r.get("b") or r.get("g2"), r.get("score")
        else:
            g1, g2, s = r[0], r[1], r[2]
        s = float(s)
        if s > 1.0:
            s /= 1000.0
        if s >= min_score:
            out.append((g1, g2, s))
    return out


def load_gene2go(path=GENE2GO):
    with open(path, "rb") as f:
        return pickle.load(f)


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #
def layer_us(emb: np.ndarray, layer: int, gene_idx=None) -> np.ndarray:
    """Mean-centred thin SVD of one layer; returns U*S (genes x components).

    ROUND-1 SEMANTICS, preserved deliberately: the SVD is computed on the FULL
    gene x 512 matrix and the result is *then* restricted to `gene_idx`. Running
    the SVD on the subset instead changes the axes and every downstream number,
    so this ordering matters.
    """
    X = np.asarray(emb[layer], dtype=np.float64)
    X = X - X.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    US = U * S
    return US if gene_idx is None else US[np.asarray(gene_idx)]


def curated_genes(gene_names):
    """The 209-gene curated set: genes named in the cycle1 TRRUST-filtered
    immune regulatory edge list, intersected with the model vocabulary."""
    import pandas as pd
    tsv = (REPO / "subproject_38_geometric_residual_stream_interpretability"
           / "implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv")
    e = pd.read_csv(tsv, sep="\t")
    named = set(e["source"].astype(str)) | set(e["target"].astype(str))
    return sorted([i for i, g in enumerate(gene_names) if g in named])


def nonzero_mask(emb: np.ndarray, layer: int) -> np.ndarray:
    return np.linalg.norm(np.asarray(emb[layer], dtype=np.float64), axis=1) > 0


# --------------------------------------------------------------------------- #
# statistics
# --------------------------------------------------------------------------- #
def emp_p(observed: float, null: np.ndarray, tail: str = "greater") -> tuple[float, int, int]:
    """Empirical permutation p using the (b+1)/(N+1) estimator.

    Returns (p, b, N). Never returns 0 -- the floor is 1/(N+1), which is the
    resolution limit of an N-permutation test and must be reported as such.
    """
    null = np.asarray(null, dtype=np.float64)
    N = null.size
    if tail == "greater":
        b = int(np.sum(null >= observed))
    elif tail == "less":
        b = int(np.sum(null <= observed))
    else:
        b = int(np.sum(np.abs(null) >= abs(observed)))
    return (b + 1) / (N + 1), b, N


def p_floor(n_perm: int) -> float:
    """Smallest p an n_perm permutation test can report."""
    return 1.0 / (n_perm + 1)


def fisher_2x2(a, b, c, d):
    """Two-sided Fisher exact test; analytic, so not resolution-limited."""
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


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Rank-based AUROC, ties handled by average rank."""
    pos, neg = np.asarray(pos), np.asarray(neg)
    allv = np.concatenate([pos, neg])
    order = allv.argsort()
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.arange(1, allv.size + 1)
    # average ranks for ties
    _, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(cnt.size)
    np.add.at(sums, inv, ranks)
    ranks = (sums / cnt)[inv]
    n1, n0 = pos.size, neg.size
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def stouffer(pvals, weights=None):
    """Combine independent one-sided p-values (Stouffer's Z)."""
    from scipy.stats import norm
    p = np.clip(np.asarray(pvals, dtype=np.float64), 1e-300, 1 - 1e-16)
    z = norm.isf(p)
    w = np.ones_like(z) if weights is None else np.asarray(weights, dtype=np.float64)
    Z = float((w * z).sum() / np.sqrt((w ** 2).sum()))
    return Z, float(norm.sf(Z))


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals, dtype=np.float64)
    n = p.size
    order = p.argsort()
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        val = p[i] * n / (rank + 1)
        prev = min(prev, val)
        q[i] = min(prev, 1.0)
    return q


# --------------------------------------------------------------------------- #
# estimators (round-1 semantics preserved)
# --------------------------------------------------------------------------- #
def poles(proj: np.ndarray, pole_frac: float = POLE_FRAC):
    """Top-K / bottom-K index sets by projection. K = int(pole_frac * n)."""
    n = proj.size
    K = max(1, int(pole_frac * n))
    order = np.argsort(proj)
    return set(order[-K:].tolist()), set(order[:K].tolist()), K


def pole_labels(proj: np.ndarray, pole_frac: float = POLE_FRAC) -> np.ndarray:
    """+1 top pole, -1 bottom pole, 0 neither. Vectorised co-pole substrate."""
    n = proj.size
    K = max(1, int(pole_frac * n))
    order = np.argsort(proj)
    lab = np.zeros(n, dtype=np.int8)
    lab[order[-K:]] = 1
    lab[order[:K]] = -1
    return lab


def copole_rate(lab: np.ndarray, pi: np.ndarray, pj: np.ndarray) -> float:
    """Fraction of pairs with both genes in the same (non-empty) pole."""
    li, lj = lab[pi], lab[pj]
    return float(np.mean((li == lj) & (li != 0)))


def copole_null(lab: np.ndarray, pi: np.ndarray, pj: np.ndarray, n_perm: int,
                rng: np.random.Generator, batch: int = 512) -> np.ndarray:
    """Gene-label-shuffle null for the co-pole rate, batched for speed."""
    n = lab.size
    out = np.empty(n_perm, dtype=np.float64)
    done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        # each row is an independent permutation of gene labels
        perm = np.argsort(rng.random((b, n)), axis=1)
        pl = lab[perm]                      # (b, n)
        li = np.take_along_axis(pl, np.broadcast_to(pi, (b, pi.size)), axis=1)
        lj = np.take_along_axis(pl, np.broadcast_to(pj, (b, pj.size)), axis=1)
        out[done:done + b] = np.mean((li == lj) & (li != 0), axis=1)
        done += b
    return out


def tf_target_scores(US: np.ndarray, comps: slice, tf_pos: np.ndarray,
                     exclude_self: bool = True) -> np.ndarray:
    """Mean cosine similarity of every gene to the TF set.

    Round 1 included each TF's self-similarity (=1) in its own mean, which adds
    1/|TF| to every TF score and to no target score -- a systematic upward bias
    in the TF/target AUROC. `exclude_self=True` removes it (leave-one-out);
    `exclude_self=False` reproduces the round-1 estimator exactly.
    """
    E = US[:, comps]
    E = E / (np.linalg.norm(E, axis=1, keepdims=True) + 1e-12)
    sim = E @ E[tf_pos].T                       # (n_genes, n_tf)
    k = tf_pos.size
    if not exclude_self:
        return sim.mean(axis=1)
    tot = sim.sum(axis=1)
    self_sim = np.zeros(E.shape[0])
    self_sim[tf_pos] = sim[tf_pos, np.arange(k)]
    is_tf = np.zeros(E.shape[0], dtype=bool)
    is_tf[tf_pos] = True
    denom = np.where(is_tf, k - 1, k)
    return (tot - self_sim) / np.maximum(denom, 1)


def tf_target_auroc(US, comps, tf_pos, tgt_pos, exclude_self=True):
    s = tf_target_scores(US, comps, tf_pos, exclude_self)
    return auroc(s[tf_pos], s[tgt_pos])
