"""F1 shared helpers: compartment classes, GO cellular-component design, class-pair-preserving
Maslov-Sneppen chain, interval helper with both scalings.

Reuses E2's helpers (../E2_ppi/e2_common.py: context loading, observed-gene SVD coordinates, vectorised
stats3.assortativity, MS / permutation runners, max-T summary) and ../stats3.py. The only new estimator code is
(i) the class-pair-preserving double-edge swap kernel (a copy of stats3._ms_swaps with one extra acceptance
condition) and (ii) the delete-d rescaling of the stats3.vertex_subsample_ci draws (CONVENTION3_ERRATA item 3).
"""
from __future__ import annotations

import json, pickle, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3 / "E2_ppi"))
import e2_common as E          # noqa: E402  (also puts R3 on sys.path)

S = E.S
X = E.X
CKPT = HERE / "_ckpt"
GENE2GO = X.REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"
GO_OBO = X.REPO / "subproject_61_two_axes_revision/data/go-basic.obo"
CLASS_ORDER = ["extracellular", "plasma_membrane", "ER_Golgi", "mitochondrion", "nucleus", "cytosol",
               "other_annotated", "unannotated"]
CLASS_TERMS = [("extracellular", {"GO:0005615"}), ("plasma_membrane", {"GO:0005886"}),
               ("ER_Golgi", {"GO:0005783", "GO:0005788", "GO:0005794"}), ("mitochondrion", {"GO:0005739"}),
               ("nucleus", {"GO:0005634"}), ("cytosol", {"GO:0005829"})]
N_CC = 30
E2_TOP = {"lung": [(3, 1), (6, 2), (8, 3)], "lung_external": [(8, 2), (3, 1), (7, 2)],
          "kidney_B": [(9, 3), (8, 3), (7, 2)]}          # E2 h1_summary.csv top-3 cells by |z| (layer, k)


def col_of(layer, k):
    return layer * E.K_SV + (k - 1)


def load_gene2go():
    with open(GENE2GO, "rb") as f:
        return pickle.load(f)


def go_terms_info():
    """id -> (namespace, name), alt_ids mapped to their primary term."""
    info, cur, alts = {}, None, []
    for line in open(GO_OBO):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {}; alts = []
        elif cur is not None and line.startswith("id: GO:"):
            cur["id"] = line[4:]
        elif cur is not None and line.startswith("name: "):
            cur["name"] = line[6:]
        elif cur is not None and line.startswith("namespace: "):
            cur["ns"] = line[11:]
        elif cur is not None and line.startswith("alt_id: "):
            alts.append(line[8:])
        elif line == "" and cur is not None and "id" in cur:
            info[cur["id"]] = (cur.get("ns", ""), cur.get("name", ""))
            for a in alts: info.setdefault(a, info[cur["id"]])
            cur = None
    return info


def assign_classes(genes, g2go):
    """One class per gene by the pre-specified priority. other_annotated = >= 1 direct GO annotation of any
    namespace but none of the six priority terms; unannotated = no direct annotation in gene2go_all.pkl."""
    cls = np.empty(len(genes), dtype=np.int64)
    for i, g in enumerate(genes):
        ann = g2go.get(g, ())
        ann = set(ann)
        c = None
        for ci, (name, terms) in enumerate(CLASS_TERMS):
            if ann & terms:
                c = ci; break
        if c is None:
            c = CLASS_ORDER.index("other_annotated") if ann else CLASS_ORDER.index("unannotated")
        cls[i] = c
    return cls


def cc_design(genes, g2go, info, n_terms=N_CC):
    """Indicators for the n_terms most frequent direct GO cellular-component annotations among V
    (ties broken by GO id). Returns (term_ids, names, counts, indicator matrix n x n_terms)."""
    from collections import Counter
    cnt = Counter()
    per_gene = []
    for g in genes:
        cc = {t for t in g2go.get(g, ()) if info.get(t, ("", ""))[0] == "cellular_component"}
        per_gene.append(cc); cnt.update(cc)
    top = sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[:n_terms]
    ids = [t for t, _ in top]
    M = np.array([[t in cc for t in ids] for cc in per_gene], dtype=np.float64)
    n_cc_annot = int(sum(len(cc) > 0 for cc in per_gene))
    return ids, [info[t][1] for t in ids], [int(c) for _, c in top], M, n_cc_annot


def residualise(Xc, D):
    """OLS residuals of every column of Xc on [1, D] (QR). Returns (residuals, R^2 per column)."""
    Z = np.column_stack([np.ones(Xc.shape[0]), D])
    Q, R = np.linalg.qr(Z)
    rank = int((np.abs(np.diag(R)) > 1e-10 * np.abs(R).max()).sum())
    assert rank == Z.shape[1], f"design rank {rank} < {Z.shape[1]}"
    Res = Xc - Q @ (Q.T @ Xc)
    tss = ((Xc - Xc.mean(0)) ** 2).sum(0)
    return Res, 1.0 - (Res ** 2).sum(0) / tss


def class_centre(Xc, cls):
    out = Xc.copy()
    for c in np.unique(cls):
        k = cls == c; out[k] -= out[k].mean(0)
    return out


def same_class_frac(cls, a, b):
    return float(np.mean(cls[a] == cls[b]))


def class_pair_counts(cls, a, b, n_cls):
    lo = np.minimum(cls[a], cls[b]); hi = np.maximum(cls[a], cls[b])
    return np.bincount(lo * n_cls + hi, minlength=n_cls * n_cls)


def config_same_class_expectation(cls, a, b):
    """Stub-matching expectation of the same-class edge fraction: sum_c (stub share of class c)^2."""
    st = cls[np.concatenate([a, b])]
    sh = np.bincount(st) / st.size
    return float((sh ** 2).sum())


# ------------------------------------------------------------------ class-pair-preserving MS chain
import numba as _nb                         # noqa: E402


@_nb.njit(cache=True)
def _ms_swaps_cp(a, b, cls, n_nodes, n_swaps, seed, keyset, max_factor):
    """stats3._ms_swaps with one extra acceptance condition: the double-edge swap (u,v),(x,y) -> (u,y),(x,v)
    is accepted only if cls[v] == cls[y] (each edge keeps its own class pair) or cls[u] == cls[x] (the two
    edges exchange class pairs). Either way every node's degree and the multiset of edge class pairs
    {cls_i, cls_j} are preserved exactly; no self-loops or multi-edges."""
    np.random.seed(seed); m = a.size; acc = 0; tried = 0
    while acc < n_swaps:
        tried += 1
        if tried > max_factor * n_swaps + 1000: break
        e1 = np.random.randint(m); e2 = np.random.randint(m)
        if e1 == e2: continue
        u, v, x, y = a[e1], b[e1], a[e2], b[e2]
        if np.random.random() < 0.5: x, y = y, x
        if u == y or x == v: continue
        if not (cls[v] == cls[y] or cls[u] == cls[x]): continue
        k1 = min(u, y) * n_nodes + max(u, y); k2 = min(x, v) * n_nodes + max(x, v)
        if k1 in keyset or k2 in keyset or k1 == k2: continue
        o1 = min(u, v) * n_nodes + max(u, v); o2 = min(x, y) * n_nodes + max(x, y)
        keyset.pop(o1); keyset.pop(o2); keyset[k1] = 1; keyset[k2] = 1
        a[e1], b[e1], a[e2], b[e2] = u, y, x, v
        acc += 1
    return acc, tried


def ms_classpair_chain(pi, pj, cls, n_nodes, n_draws, seed=0, burn_factor=10, thin_factor=1, max_factor=400):
    """Like stats3.maslov_sneppen_chain (burn burn_factor*m, thin thin_factor*m ACCEPTED swaps) but restricted to
    class-pair-preserving swaps. Yields (a, b, acc_rate, completed) where completed is False if the try cap
    was hit before the requested number of accepted swaps."""
    from numba import types
    from numba.typed import Dict
    a = np.minimum(pi, pj).astype(np.int64).copy(); b = np.maximum(pi, pj).astype(np.int64).copy()
    cls = np.asarray(cls, np.int64)
    ks = Dict.empty(key_type=types.int64, value_type=types.int64)
    for i, j in zip(a, b): ks[int(i) * n_nodes + int(j)] = 1
    if len(ks) != a.size: raise ValueError("duplicate edges")
    m = a.size
    acc, tried = _ms_swaps_cp(a, b, cls, n_nodes, burn_factor * m, seed, ks, max_factor)
    burn_ok = acc >= burn_factor * m
    for d in range(n_draws):
        acc, tried = _ms_swaps_cp(a, b, cls, n_nodes, thin_factor * m, seed + 1 + d, ks, max_factor)
        yield a.copy(), b.copy(), acc / max(tried, 1), bool(burn_ok and acc >= thin_factor * m)


def run_chain(chain, cols_list, extra=None, log=None, every=500, tag=""):
    """Evaluate E.assort_cols for each matrix in cols_list on every chain draw; extra(a, b) -> float vector."""
    R = [[] for _ in cols_list]; ex = []; acc = []; ok = []
    for d, item in enumerate(chain):
        a, b, ar = item[0], item[1], item[2]
        ok.append(item[3] if len(item) > 3 else True)
        for i, M in enumerate(cols_list): R[i].append(E.assort_cols(M, a, b))
        if extra is not None: ex.append(extra(a, b))
        acc.append(ar)
        if log and (d + 1) % every == 0: log(f"    {tag} draw {d+1} acc={ar:.3f}")
    return [np.array(r) for r in R], (np.array(ex) if extra is not None else None), np.array(acc), np.array(ok)


# ------------------------------------------------------------------ intervals
def vs_ci_both(stat_fn, n, pi, pj, B, seed, frac=0.5, level=0.95):
    """stats3.vertex_subsample_ci (pre-registered sqrt(m/n) scaling) plus the delete-d scaling sqrt(m/(n-m))
    applied to the SAME draws (CONVENTION3_ERRATA item 3)."""
    rng = np.random.default_rng(seed)
    est, lo, hi, draws = S.vertex_subsample_ci(stat_fn, n, pi, pj, B, rng, frac=frac, level=level)
    m = int(frac * n)
    dev = (draws - np.nanmean(draws)) * np.sqrt(m / (n - m))
    a = (1 - level) / 2
    lo_d, hi_d = est - np.nanquantile(dev, 1 - a), est - np.nanquantile(dev, a)
    return dict(est=float(est), ci_lo_prereg=float(lo), ci_hi_prereg=float(hi), ci_lo_deleted=float(lo_d),
                ci_hi_deleted=float(hi_d), se_prereg=float(np.nanstd(draws, ddof=1) * np.sqrt(m / n)),
                se_deleted=float(np.nanstd(draws, ddof=1) * np.sqrt(m / (n - m))), B=int(B),
                n_nan=int(np.isnan(draws).sum()))


def save(path, **kw):
    np.savez(path, **kw)


def jdump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
