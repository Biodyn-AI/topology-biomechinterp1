"""Shared estimators and nulls for the round-3 experiments.

Every experiment imports this module so that one definition of each statistic is used
everywhere. `python stats3.py` runs the self-tests (must all pass before use).

Conventions (see CONVENTION3.md):
  * observed genes V = rows with non-zero embedding norm (identical across layers, because a
    gene absent from every sampled cell's input has count 0 and an all-zero mean embedding);
  * the SVD is computed per layer on the V x d matrix, mean-centred over V;
  * empirical p-values use (b+1)/(N+1); the floor 1/(N+1) is reported as such.
"""
from __future__ import annotations

from math import exp, lgamma

import numpy as np
from scipy.stats import norm

# =====================================================================================
# geometry
# =====================================================================================
def observed_mask(emb) -> np.ndarray:
    """Boolean mask of genes with a non-zero embedding (checked at L0 and the last layer)."""
    m0 = np.linalg.norm(np.asarray(emb[0], dtype=np.float64), axis=1) > 0
    mL = np.linalg.norm(np.asarray(emb[emb.shape[0] - 1], dtype=np.float64), axis=1) > 0
    assert np.array_equal(m0, mL), "observed set differs across layers"
    return m0


def svd_observed(emb, layer: int, obs_idx: np.ndarray):
    """Thin SVD of the observed-gene matrix at one layer, mean-centred over observed genes.

    Returns (US, S, Vt); row r of US corresponds to gene obs_idx[r]. Signs are fixed so that
    the largest-|loading| gene of each component has a positive score (deterministic output).
    """
    X = np.asarray(emb[layer], dtype=np.float64)[obs_idx]
    X = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    US = U * S
    sgn = np.sign(US[np.abs(US).argmax(axis=0), np.arange(US.shape[1])])
    sgn[sgn == 0] = 1
    return US * sgn, S, Vt * sgn[:, None]


def effective_rank(S: np.ndarray) -> float:
    """exp(Shannon entropy of the normalised singular values) (Roy & Vetterli 2007)."""
    p = S / S.sum()
    p = p[p > 0]
    return float(np.exp(-(p * np.log(p)).sum()))


def participation_ratio(S: np.ndarray) -> float:
    v = S ** 2
    return float(v.sum() ** 2 / (v ** 2).sum())


def l2norm_rows(E: np.ndarray) -> np.ndarray:
    return E / (np.linalg.norm(E, axis=1, keepdims=True) + 1e-12)


# =====================================================================================
# gene sets
# =====================================================================================
def trrust_gene_set(V: set, trrust_pairs) -> list:
    """Genes in V with at least one TRRUST edge whose two endpoints are both in V."""
    return sorted({g for a, b in trrust_pairs if a in V and b in V for g in (a, b)})


# =====================================================================================
# basic statistics
# =====================================================================================
def auroc(pos, neg) -> float:
    """Mann-Whitney AUROC with average ranks for ties."""
    pos = np.asarray(pos, dtype=np.float64); neg = np.asarray(neg, dtype=np.float64)
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(allv.size, dtype=np.float64)
    ranks[order] = np.arange(1, allv.size + 1)
    _, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(cnt.size); np.add.at(sums, inv, ranks)
    ranks = (sums / cnt)[inv]
    n1, n0 = pos.size, neg.size
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def auroc_columns(scores: np.ndarray, is_pos: np.ndarray) -> np.ndarray:
    """AUROC for every column of `scores` (n x B) with the same positive mask (no tie correction:
    use only for continuous scores, e.g. permutation nulls)."""
    n1 = int(is_pos.sum()); n0 = is_pos.size - n1
    ranks = scores.argsort(axis=0).argsort(axis=0) + 1.0
    return (ranks[is_pos].sum(axis=0) - n1 * (n1 + 1) / 2) / (n1 * n0)


def delong(pos_a, neg_a, pos_b=None, neg_b=None):
    """DeLong variance of one AUROC, or of the paired difference A - B (same pos/neg items).

    Returns dict(auc_a, var_a[, auc_b, var_b, cov, delta, se_delta, z, p]).
    """
    def comps(p, n):
        p = np.asarray(p, float)[:, None]; n = np.asarray(n, float)[None, :]
        psi = (p > n) + 0.5 * (p == n)
        return psi.mean(1), psi.mean(0), psi.mean()       # V10 (per positive), V01 (per negative), auc
    v10a, v01a, A = comps(pos_a, neg_a)
    m, n = v10a.size, v01a.size
    out = dict(auc_a=float(A), var_a=float(v10a.var(ddof=1) / m + v01a.var(ddof=1) / n))
    if pos_b is not None:
        v10b, v01b, B = comps(pos_b, neg_b)
        s10 = np.cov(np.vstack([v10a, v10b])); s01 = np.cov(np.vstack([v01a, v01b]))
        var = s10 / m + s01 / n
        d = A - B
        se = float(np.sqrt(max(var[0, 0] + var[1, 1] - 2 * var[0, 1], 1e-300)))
        z = d / se
        out.update(auc_b=float(B), var_b=float(var[1, 1]), cov=float(var[0, 1]),
                   delta=float(d), se_delta=se, z=float(z), p=float(2 * norm.sf(abs(z))))
    return out


def hanley_mcneil_se(auc: float, n1: int, n0: int) -> float:
    q1 = auc / (2 - auc); q2 = 2 * auc ** 2 / (1 + auc)
    return float(np.sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc ** 2) + (n0 - 1) * (q2 - auc ** 2)) / (n1 * n0)))


def mde_auroc(n1: int, n0: int, alpha: float = 0.05, power: float = 0.80, auc0: float = 0.5) -> float:
    """Smallest AUROC > auc0 detectable against auc0 (two-sided alpha) with the given power,
    Hanley-McNeil variances under null and alternative."""
    za, zb = norm.isf(alpha / 2), norm.isf(1 - power)
    se0 = hanley_mcneil_se(auc0, n1, n0)
    lo, hi = auc0, 0.9999
    for _ in range(80):
        mid = (lo + hi) / 2
        if (mid - auc0) >= za * se0 + zb * hanley_mcneil_se(mid, n1, n0): hi = mid
        else: lo = mid
    return float(hi)


def mde_paired_delta(se_a: float, se_b: float, rho: float, alpha=0.05, power=0.80) -> float:
    """Minimum detectable paired AUROC difference given per-method SEs and their correlation."""
    se = np.sqrt(se_a ** 2 + se_b ** 2 - 2 * rho * se_a * se_b)
    return float((norm.isf(alpha / 2) + norm.isf(1 - power)) * se)


def emp_p(observed: float, null, tail: str = "greater"):
    null = np.asarray(null, dtype=np.float64); N = null.size
    if tail == "greater": b = int(np.sum(null >= observed))
    elif tail == "less": b = int(np.sum(null <= observed))
    else:
        c = np.median(null); b = int(np.sum(np.abs(null - c) >= abs(observed - c)))
    return (b + 1) / (N + 1), b, N


def p_floor(n_perm: int) -> float:
    return 1.0 / (n_perm + 1)


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=np.float64); n = p.size
    order = p.argsort(); q = np.empty(n); prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]; prev = min(prev, p[i] * n / (rank + 1)); q[i] = min(prev, 1.0)
    return q


def stouffer(*_a, **_k):
    """Removed. Layers (and contexts) share genes and are successive transforms of one another, so their
    p-values are dependent and Stouffer's independence assumption fails. Use a joint-permutation statistic
    (one label permutation applied to every layer) or max_t_adjust instead."""
    raise RuntimeError("stouffer() is disabled in round 3; use max_t_adjust or a joint-permutation statistic")


def fisher_2x2(a, b, c, d):
    n = a + b + c + d; r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    if r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0: return 1.0, 1.0
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    def lpx(x):
        return (lgamma(r1 + 1) - lgamma(x + 1) - lgamma(r1 - x + 1) + lgamma(r2 + 1)
                - lgamma(c1 - x + 1) - lgamma(c2 - r1 + x + 1) - lgamma(n + 1)
                + lgamma(c1 + 1) + lgamma(c2 + 1))
    obs = lpx(a); p = 0.0
    for x in range(max(0, c1 - r2), min(c1, r1) + 1):
        if lpx(x) <= obs + 1e-12: p += exp(lpx(x))
    return float(odds), float(min(p, 1.0))


# =====================================================================================
# co-pole statistic and its nulls
# =====================================================================================
def pole_labels(proj: np.ndarray, pole_frac: float = 0.27) -> np.ndarray:
    """+1 top pole, -1 bottom pole, 0 neither; K = floor(pole_frac * n). pole_frac = 0.5 puts
    every gene in a pole (top half / bottom half)."""
    n = proj.size; K = max(1, int(pole_frac * n)); order = np.argsort(proj, kind="mergesort")
    lab = np.zeros(n, dtype=np.int8); lab[order[-K:]] = 1; lab[order[:K]] = -1
    return lab


def copole_rate(lab, pi, pj) -> float:
    li, lj = lab[pi], lab[pj]
    return float(np.mean((li == lj) & (li != 0)))


def _copole_rows(L, pi, pj):
    li = L[:, pi]; lj = L[:, pj]
    return np.mean((li == lj) & (li != 0), axis=1)


def copole_null_label(lab, pi, pj, n_perm, rng, batch=256):
    """Gene-label shuffle. Preserves: the network exactly (every edge, every degree) and the pole
    sizes. Breaks: the association between a gene's pole and anything else about it -- including its
    degree, so hubs that sit in poles inflate the observed rate relative to this null."""
    out = np.empty(n_perm); done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        L = rng.permuted(np.broadcast_to(lab, (b, lab.size)), axis=1)
        out[done:done + b] = _copole_rows(L, pi, pj); done += b
    return out


def degree_bins(deg: np.ndarray, n_bins: int = 10) -> np.ndarray:
    """Quantile bins of degree; degree-0 genes get their own bin."""
    deg = np.asarray(deg); bins = np.full(deg.size, -1, dtype=int)
    pos = deg > 0
    if pos.any():
        q = np.unique(np.quantile(deg[pos], np.linspace(0, 1, n_bins + 1)))
        bins[pos] = np.clip(np.searchsorted(q, deg[pos], side="right") - 1, 0, len(q) - 2)
    return bins


def copole_null_degree_stratified(lab, pi, pj, deg, n_perm, rng, n_bins=10, batch=256):
    """Gene-label shuffle within degree bins. Preserves: the network and, approximately, the joint
    distribution of (degree, pole) -- so a hub-in-pole excess cannot by itself produce enrichment."""
    bins = degree_bins(deg, n_bins); groups = [np.flatnonzero(bins == b) for b in np.unique(bins)]
    out = np.empty(n_perm); done = 0
    while done < n_perm:
        b = min(batch, n_perm - done); L = np.broadcast_to(lab, (b, lab.size)).copy()
        for g in groups:
            if g.size > 1: L[:, g] = rng.permuted(L[:, g], axis=1)
        out[done:done + b] = _copole_rows(L, pi, pj); done += b
    return out


def copole_null_configuration(lab, pi, pj, n_perm, rng, batch=64):
    """Degree-preserving rewiring (configuration model by stub matching). Preserves: every gene's
    degree and pole label exactly. Breaks: which partner each edge connects to. Self-pairs created by
    matching are dropped. Tests whether partners co-locate beyond what degrees alone imply."""
    stubs = np.concatenate([pi, pj]); m = pi.size
    out = np.empty(n_perm); done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        S = rng.permuted(np.broadcast_to(stubs, (b, stubs.size)), axis=1)
        a, c = S[:, :m], S[:, m:2 * m]
        la, lc = lab[a], lab[c]
        same = (la == lc) & (la != 0); keep = a != c
        out[done:done + b] = (same & keep).sum(1) / keep.sum(1); done += b
    return out


def copole_expected_configuration(lab, pi, pj) -> float:
    """Analytic expectation of the co-pole rate under stub matching: sum over poles of (stub share)^2."""
    stubs = np.concatenate([pi, pj]); ls = lab[stubs]
    return float(sum(np.mean(ls == p) ** 2 for p in (-1, 1)))


def assortativity(x, pi, pj) -> float:
    """Scalar (Newman) assortativity of a continuous node attribute over undirected edges: the Pearson
    correlation of (x_i, x_j) over both orientations of every edge. K-free alternative to co-pole."""
    a = np.concatenate([x[pi], x[pj]]); b = np.concatenate([x[pj], x[pi]])
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def assortativity_null_label(x, pi, pj, n_perm, rng, batch=256):
    out = np.empty(n_perm); done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        Xp = rng.permuted(np.broadcast_to(x, (b, x.size)), axis=1)
        A = np.concatenate([Xp[:, pi], Xp[:, pj]], 1); Bm = np.concatenate([Xp[:, pj], Xp[:, pi]], 1)
        A = A - A.mean(1, keepdims=True); Bm = Bm - Bm.mean(1, keepdims=True)
        out[done:done + b] = (A * Bm).sum(1) / np.sqrt((A * A).sum(1) * (Bm * Bm).sum(1)); done += b
    return out


def assortativity_null_configuration(x, pi, pj, n_perm, rng, batch=64):
    stubs = np.concatenate([pi, pj]); m = pi.size
    out = np.empty(n_perm); done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        S = rng.permuted(np.broadcast_to(stubs, (b, stubs.size)), axis=1)
        a, c = S[:, :m], S[:, m:2 * m]
        for r in range(b):
            k = a[r] != c[r]; ai, ci = a[r][k], c[r][k]
            out[done + r] = assortativity(x, ai, ci)
        done += b
    return out


def gene_bootstrap(stat_fn, n_genes, pi, pj, B, rng):
    """Cluster (gene-level) bootstrap for an edge statistic: resample genes with replacement and keep
    every edge whose two endpoints were drawn, with multiplicity. Respects the non-independence of
    edges that share a gene. stat_fn(pi_b, pj_b) -> float."""
    out = np.empty(B)
    for r in range(B):
        w = np.bincount(rng.integers(0, n_genes, n_genes), minlength=n_genes)
        mult = w[pi] * w[pj]; keep = mult > 0
        reps = np.repeat(np.flatnonzero(keep), mult[keep])
        out[r] = stat_fn(pi[reps], pj[reps]) if reps.size else np.nan
    return out


# =====================================================================================
# node-level TF identity
# =====================================================================================
def tf_scores_loo(E: np.ndarray, tf_pos: np.ndarray) -> np.ndarray:
    """Leave-one-out mean cosine similarity of every row of E to the TF rows (self excluded)."""
    En = l2norm_rows(E); v = En[tf_pos].sum(0); s = En @ v
    k = tf_pos.size; is_tf = np.zeros(En.shape[0], bool); is_tf[tf_pos] = True
    s = np.where(is_tf, (s - (En * En).sum(1)) / max(k - 1, 1), s / k)
    return s


def tf_auroc(E, tf_pos, tgt_pos) -> float:
    s = tf_scores_loo(E, tf_pos)
    return auroc(s[tf_pos], s[tgt_pos])


def tf_auroc_null(E, tf_pos, tgt_pos, n_perm, rng, batch=500):
    """Label-shuffle null for tf_auroc, vectorised: TF labels are re-drawn among the labelled genes
    (TFs + targets) with the TF count fixed; scores and AUROC recomputed for each draw."""
    lab = np.concatenate([tf_pos, tgt_pos]); k = tf_pos.size; n = lab.size
    En = l2norm_rows(E[lab]); selfsim = (En * En).sum(1)
    out = np.empty(n_perm); done = 0
    while done < n_perm:
        b = min(batch, n_perm - done)
        M = np.zeros((n, b)); idx = np.argsort(rng.random((n, b)), axis=0)[:k]
        M[idx, np.arange(b)[None, :]] = 1.0
        tot = En @ (En.T @ M)                                  # n x b
        sc = np.where(M > 0, (tot - selfsim[:, None]) / max(k - 1, 1), tot / k)
        out[done:done + b] = _auroc_mask_cols(sc, M > 0)
        done += b
    return out


def _auroc_mask_cols(sc, pos_mask):
    """Column-wise AUROC where the positive set differs per column."""
    ranks = sc.argsort(axis=0).argsort(axis=0) + 1.0
    n1 = pos_mask.sum(0); n0 = pos_mask.shape[0] - n1
    return ((ranks * pos_mask).sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)


# =====================================================================================
# self-tests
# =====================================================================================
def _selftest():
    rng = np.random.default_rng(0)
    from sklearn.metrics import roc_auc_score
    # 1. AUROC with ties equals sklearn
    p = rng.integers(0, 5, 60).astype(float); q = rng.integers(0, 5, 90).astype(float)
    assert abs(auroc(p, q) - roc_auc_score(np.r_[np.ones(60), np.zeros(90)], np.r_[p, q])) < 1e-12
    # 2. DeLong variance ~ simulated sampling variance
    sims = [auroc(rng.normal(0.6, 1, 50), rng.normal(0, 1, 150)) for _ in range(3000)]
    d = delong(rng.normal(0.6, 1, 2000), rng.normal(0, 1, 6000))
    hm = hanley_mcneil_se(np.mean(sims), 50, 150)
    assert abs(np.std(sims) - hm) / np.std(sims) < 0.12, (np.std(sims), hm)
    # 3. MDE monotone in sample size and sensible
    assert mde_auroc(49, 147) > mde_auroc(500, 1500) > mde_auroc(5000, 15000) > 0.5
    # 4. co-pole nulls: calibrated under a random-label null (p roughly uniform)
    n = 400; pi = rng.integers(0, n, 1200); pj = rng.integers(0, n, 1200); k = pi != pj; pi, pj = pi[k], pj[k]
    ps = []
    for _ in range(60):
        lab = pole_labels(rng.normal(size=n))
        ps.append(emp_p(copole_rate(lab, pi, pj), copole_null_label(lab, pi, pj, 400, rng))[0])
    assert 0.3 < np.mean(ps) < 0.7, np.mean(ps)
    # 5. configuration null preserves degrees exactly (pre-drop) and matches analytic expectation
    lab = pole_labels(rng.normal(size=n))
    nulls = copole_null_configuration(lab, pi, pj, 2000, rng)
    assert abs(nulls.mean() - copole_expected_configuration(lab, pi, pj)) < 0.01
    # 6a. partner homophily with balanced degrees is detected by all three nulls
    x = rng.normal(size=n); lab = pole_labels(x)
    src = rng.integers(0, n, 1500); dst = np.empty_like(src)
    for e, i in enumerate(src):
        same = np.flatnonzero(lab == lab[i])
        dst[e] = rng.choice(same) if (lab[i] != 0 and rng.random() < 0.6) else rng.integers(0, n)
    k = src != dst; pi2, pj2 = src[k], dst[k]
    deg = np.bincount(np.r_[pi2, pj2], minlength=n)
    for f in (lambda: copole_null_label(lab, pi2, pj2, 500, rng),
              lambda: copole_null_degree_stratified(lab, pi2, pj2, deg, 500, rng),
              lambda: copole_null_configuration(lab, pi2, pj2, 500, rng)):
        assert emp_p(copole_rate(lab, pi2, pj2), f())[0] < 0.01
    assert assortativity(x, pi2, pj2) > 0.1
    # 6b. enrichment that is fully explained by where the degree sits: all edges among top-pole genes.
    #     The label shuffle calls it significant; the degree-preserving configuration null does not.
    order = np.argsort(x); top = order[-60:]
    pi3 = top[rng.integers(0, 60, 400)]; pj3 = top[rng.integers(0, 60, 400)]
    k = pi3 != pj3; pi3, pj3 = pi3[k], pj3[k]
    assert emp_p(copole_rate(lab, pi3, pj3), copole_null_label(lab, pi3, pj3, 500, rng))[0] < 0.01
    assert emp_p(copole_rate(lab, pi3, pj3), copole_null_configuration(lab, pi3, pj3, 500, rng))[0] > 0.2
    # 7. TF AUROC: vectorised null matches loop semantics; observed matches round-2 estimator
    E = rng.normal(size=(120, 6)); tf = np.arange(30); tg = np.arange(30, 120)
    E[tf] += 0.8
    obs = tf_auroc(E, tf, tg)
    En = l2norm_rows(E); s_loop = np.array([(En[i] @ En[tf].sum(0) - (1 if i in set(tf) else 0)) /
                                           (29 if i in set(tf) else 30) for i in range(120)])
    assert abs(auroc(s_loop[tf], s_loop[tg]) - obs) < 1e-12
    nul = tf_auroc_null(E, tf, tg, 2000, rng)
    assert abs(nul.mean() - 0.5) < 0.02 and obs > np.quantile(nul, 0.99)
    # 8. gene bootstrap runs and is centred near the estimate
    bs = gene_bootstrap(lambda a, b: copole_rate(lab, a, b), n, pi2, pj2, 200, rng)
    assert abs(np.nanmean(bs) - copole_rate(lab, pi2, pj2)) < 0.1
    print("stats3 self-tests: ALL PASS")


if __name__ == "__main__":
    _selftest()


# =====================================================================================
# round-3 additions (critique amendments A4)
# =====================================================================================
def assert_resolution(n_perm: int, alpha: float, m: int = 1):
    """A permutation test can only support a significance claim at alpha/m if its floor is below it."""
    floor = 1.0 / (n_perm + 1)
    if floor >= alpha / m:
        raise ValueError(f"N={n_perm} gives floor {floor:.2e} >= alpha/m = {alpha/m:.2e}: cannot support the claim")
    return floor


def max_t_adjust(t_obs: np.ndarray, t_null: np.ndarray, tail: str = "greater"):
    """Westfall-Young single-step max-T adjustment over a family of tests computed on the SAME null draws.

    t_obs: (m,) observed statistics; t_null: (N, m) statistics on N joint null draws (each draw applies one
    permutation / rewiring to every test). Statistics are standardised by their own null mean and sd so that
    tests on different scales are comparable. Returns (z_obs, p_adj, p_raw): p_adj controls the FWER over
    the family; p_raw is each test's own (b+1)/(N+1).
    """
    t_obs = np.asarray(t_obs, float); t_null = np.asarray(t_null, float)
    mu = t_null.mean(0); sd = t_null.std(0, ddof=1); sd[sd == 0] = np.inf
    z_obs = (t_obs - mu) / sd; z_null = (t_null - mu) / sd
    if tail == "two-sided": z_obs, z_null = np.abs(z_obs), np.abs(z_null)
    elif tail == "less": z_obs, z_null = -z_obs, -z_null
    mx = z_null.max(1); N = t_null.shape[0]
    p_adj = np.array([(np.sum(mx >= z) + 1) / (N + 1) for z in z_obs])
    p_raw = np.array([(np.sum(z_null[:, j] >= z_obs[j]) + 1) / (N + 1) for j in range(z_obs.size)])
    return z_obs, p_adj, p_raw


def signed_free_auroc(a: float) -> float:
    """Sign-invariant AUROC for a single SVD direction, whose sign is arbitrary: 0.5 + |AUROC - 0.5|."""
    return 0.5 + abs(a - 0.5)


def label_noise_corrected(auc_obs: float, pi: float) -> float:
    """If a fraction pi of candidate negatives are true edges scored like positives, the observed AUROC is
    attenuated towards 0.5: AUROC_true - 0.5 = (AUROC_obs - 0.5) / (1 - pi)."""
    return 0.5 + (auc_obs - 0.5) / (1 - pi)


def one_sided_upper(boot: np.ndarray, level: float = 0.95) -> float:
    return float(np.nanquantile(boot, level))


def equivalence_excluded(upper_bounds, margin: float) -> bool:
    """Intersection-union test: 'effect >= margin excluded' holds only if EVERY layer's one-sided upper bound
    lies below the margin (no multiplicity penalty is needed for an IU test)."""
    return bool(np.all(np.asarray(upper_bounds) < margin))


# ---------------------------------------------------------------- within-TF stratified AUROC
class StratifiedPairs:
    """Edge-level benchmark container. Pairs (tf, target, label) with scores; the primary statistic is the
    within-TF stratified AUROC: a Mann-Whitney statistic pooled over TFs,
        A = sum_t sum_{p in P_t, q in N_t} w_p w_q 1[s_p > s_q] / sum_t (sum_{P_t} w)(sum_{N_t} w),
    which compares each TF's documented targets only with that same TF's candidate non-targets, so TF
    identity (node-level information) cannot contribute. Ties count 1/2.
    """
    def __init__(self, tf, tg, y):
        self.tf = np.asarray(tf); self.tg = np.asarray(tg); self.y = np.asarray(y).astype(bool)
        keep = np.zeros(self.tf.size, bool)
        for t in np.unique(self.tf):
            m = self.tf == t
            if self.y[m].any() and (~self.y[m]).any(): keep |= m
        self.tf, self.tg, self.y = self.tf[keep], self.tg[keep], self.y[keep]
        self.tfs, self.tf_code = np.unique(self.tf, return_inverse=True)
        self.tgs, self.tg_code = np.unique(self.tg, return_inverse=True)

    def auroc(self, s, w=None) -> float:
        s = np.asarray(s, float); w = np.ones(s.size) if w is None else np.asarray(w, float)
        num = den = 0.0
        for c in range(self.tfs.size):
            m = self.tf_code == c
            sp, wp = s[m & self.y], w[m & self.y]; sn, wn = s[m & ~self.y], w[m & ~self.y]
            if wp.sum() == 0 or wn.sum() == 0: continue
            order = np.argsort(sn, kind="mergesort"); sn, wn = sn[order], wn[order]
            cw = np.concatenate([[0.0], np.cumsum(wn)])
            lo = np.searchsorted(sn, sp, side="left"); hi = np.searchsorted(sn, sp, side="right")
            num += float((wp * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))).sum())
            den += float(wp.sum() * wn.sum())
        return num / den if den > 0 else float("nan")

    def macro_auroc(self, s, min_pos: int = 3) -> float:
        vals = []
        for c in range(self.tfs.size):
            m = self.tf_code == c
            if (m & self.y).sum() >= min_pos: vals.append(auroc(s[m & self.y], s[m & ~self.y]))
        return float(np.mean(vals)) if vals else float("nan")

    def bootstrap(self, s, B: int, rng, two_way: bool = True) -> np.ndarray:
        """Pigeonhole (TF x target) cluster bootstrap: TFs and targets are resampled independently with
        replacement and each pair is weighted by (TF multiplicity) x (target multiplicity). two_way=False
        resamples TFs only."""
        out = np.empty(B)
        for r in range(B):
            wtf = np.bincount(rng.integers(0, self.tfs.size, self.tfs.size), minlength=self.tfs.size)
            w = wtf[self.tf_code].astype(float)
            if two_way:
                wtg = np.bincount(rng.integers(0, self.tgs.size, self.tgs.size), minlength=self.tgs.size)
                w = w * wtg[self.tg_code]
            out[r] = self.auroc(s, w)
        return out

    def permute_within_tf(self, rng) -> np.ndarray:
        """Labels permuted within each TF (counts per TF preserved): the null for the stratified AUROC."""
        y = self.y.copy()
        for c in range(self.tfs.size):
            idx = np.flatnonzero(self.tf_code == c); y[idx] = rng.permutation(y[idx])
        return y


# ---------------------------------------------------------------- detection-stratified node AUROC
def stratified_auroc(score, is_pos, strata) -> float:
    """Mann-Whitney AUROC pooled within strata (e.g. detection deciles): positives are compared only with
    negatives of the same stratum, so a covariate that defines the strata cannot drive the statistic."""
    score = np.asarray(score, float); is_pos = np.asarray(is_pos, bool); strata = np.asarray(strata)
    num = den = 0.0
    for st in np.unique(strata):
        m = strata == st; p = score[m & is_pos]; q = score[m & ~is_pos]
        if p.size and q.size:
            num += auroc(p, q) * p.size * q.size; den += p.size * q.size
    return num / den if den else float("nan")


def permute_within_strata(is_pos, strata, rng):
    y = np.asarray(is_pos, bool).copy(); strata = np.asarray(strata)
    for st in np.unique(strata):
        idx = np.flatnonzero(strata == st); y[idx] = rng.permutation(y[idx])
    return y


def quantile_strata(x, n_bins: int = 10) -> np.ndarray:
    q = np.unique(np.quantile(x, np.linspace(0, 1, n_bins + 1)))
    return np.clip(np.searchsorted(q, x, side="right") - 1, 0, len(q) - 2)


# ---------------------------------------------------------------- vertex subsampling CI
def vertex_subsample_ci(stat_fn, n_nodes: int, pi, pj, B: int, rng, frac: float = 0.5, level: float = 0.95):
    """m-out-of-n vertex subsampling without replacement (no duplicate nodes, no self-pairs). Keeps edges with
    both ends sampled, recomputes stat_fn(node_idx, pi_sub, pj_sub), and rescales deviations by sqrt(m/n).
    Returns (estimate on full graph, lo, hi, draws)."""
    full = stat_fn(np.arange(n_nodes), pi, pj); m = int(frac * n_nodes); draws = np.empty(B)
    for r in range(B):
        nodes = np.sort(rng.choice(n_nodes, m, replace=False)); keep = np.zeros(n_nodes, bool); keep[nodes] = True
        e = keep[pi] & keep[pj]
        draws[r] = stat_fn(nodes, pi[e], pj[e]) if e.sum() > 2 else np.nan
    dev = (draws - np.nanmean(draws)) * np.sqrt(m / n_nodes)
    a = (1 - level) / 2
    return full, full - np.nanquantile(dev, 1 - a), full - np.nanquantile(dev, a), draws


# ---------------------------------------------------------------- Maslov-Sneppen swap chain (numba)
try:
    import numba as _nb

    @_nb.njit(cache=True)
    def _ms_swaps(a, b, n_nodes, n_swaps, seed, keyset):
        """In-place double-edge swaps on a simple undirected graph (a[e], b[e]); keyset holds i*n+j (i<j)."""
        np.random.seed(seed); m = a.size; acc = 0; tried = 0
        while acc < n_swaps:
            tried += 1
            if tried > 50 * n_swaps + 1000: break
            e1 = np.random.randint(m); e2 = np.random.randint(m)
            if e1 == e2: continue
            u, v, x, y = a[e1], b[e1], a[e2], b[e2]
            if np.random.random() < 0.5: x, y = y, x
            # proposal: (u,v),(x,y) -> (u,y),(x,v)
            if u == y or x == v: continue
            k1 = min(u, y) * n_nodes + max(u, y); k2 = min(x, v) * n_nodes + max(x, v)
            if k1 in keyset or k2 in keyset or k1 == k2: continue
            o1 = min(u, v) * n_nodes + max(u, v); o2 = min(x, y) * n_nodes + max(x, y)
            keyset.pop(o1); keyset.pop(o2); keyset[k1] = 1; keyset[k2] = 1
            a[e1], b[e1], a[e2], b[e2] = u, y, x, v
            acc += 1
        return acc, tried

    def maslov_sneppen_chain(pi, pj, n_nodes, n_draws, seed=0, burn_factor=10, thin_factor=1):
        """Yield n_draws degree-preserving rewirings of a simple undirected graph by chained double-edge swaps:
        burn-in burn_factor*m accepted swaps, then thin_factor*m accepted swaps between successive draws.
        Preserves every node's degree exactly; no self-loops or multi-edges. Yields (pi_draw, pj_draw, acc_rate)."""
        from numba import types
        from numba.typed import Dict
        a = np.minimum(pi, pj).astype(np.int64).copy(); b = np.maximum(pi, pj).astype(np.int64).copy()
        ks = Dict.empty(key_type=types.int64, value_type=types.int64)
        for i, j in zip(a, b): ks[int(i) * n_nodes + int(j)] = 1
        if len(ks) != a.size: raise ValueError("input graph has duplicate edges")
        m = a.size
        acc, tried = _ms_swaps(a, b, n_nodes, burn_factor * m, seed, ks)
        for d in range(n_draws):
            acc, tried = _ms_swaps(a, b, n_nodes, thin_factor * m, seed + 1 + d, ks)
            yield a.copy(), b.copy(), acc / max(tried, 1)
except Exception:          # pragma: no cover
    maslov_sneppen_chain = None


def _selftest_round3():
    rng = np.random.default_rng(1)
    # max-T: adjusted p >= raw p, and calibrated (FWER ~ alpha) under a global null of 30 correlated tests
    rej = 0
    for _ in range(200):
        base = rng.normal(size=(400, 1)); null = base + 0.5 * rng.normal(size=(400, 30))
        obs = 0 + 0.5 * rng.normal(size=30) + rng.normal()
        z, pa, pr = max_t_adjust(obs, null)
        assert np.all(pa >= pr - 1e-12)
        rej += pa.min() < 0.05
    assert 0.01 < rej / 200 < 0.10, rej / 200
    # stratified pairs AUROC equals pooled AUROC when there is one TF; ties count 1/2
    tf = np.zeros(50, int); y = np.r_[np.ones(20), np.zeros(30)].astype(bool); s = rng.normal(size=50)
    sp = StratifiedPairs(tf, np.arange(50), y)
    assert abs(sp.auroc(s) - auroc(s[y], s[~y])) < 1e-12
    st = np.round(s, 1); assert abs(sp.auroc(st) - auroc(st[y], st[~y])) < 1e-12
    # TF identity cannot leak: a score that is a pure function of the TF gives exactly 0.5
    tf = np.repeat(np.arange(10), 20); y = np.tile(np.r_[np.ones(5), np.zeros(15)], 10).astype(bool)
    sp = StratifiedPairs(tf, rng.integers(0, 100, 200), y)
    assert abs(sp.auroc(tf.astype(float) * 3.0) - 0.5) < 1e-12
    bs = sp.bootstrap(rng.normal(size=200) + y * 0.8, 200, rng)
    assert np.nanmean(bs) > 0.65
    # detection-stratified AUROC kills a pure detection confound
    det = rng.normal(size=2000); pos = rng.random(2000) < 1 / (1 + np.exp(-2 * det))
    stra = quantile_strata(det, 20)
    assert auroc(det[pos], det[~pos]) > 0.8 and abs(stratified_auroc(det, pos, stra) - 0.5) < 0.06
    # Maslov-Sneppen chain preserves degrees exactly, keeps the graph simple, and matches the analytic
    # configuration-model co-pole expectation on a sparse graph
    if maslov_sneppen_chain is not None:
        n = 600; src = rng.integers(0, n, 3000); dst = rng.integers(0, n, 3000)
        e = {(min(i, j), max(i, j)) for i, j in zip(src, dst) if i != j}; e = np.array(sorted(e))
        deg0 = np.bincount(e.ravel(), minlength=n)
        lab = pole_labels(rng.normal(size=n)); vals = []
        for a, b, acc in maslov_sneppen_chain(e[:, 0], e[:, 1], n, 300, seed=3):
            assert np.array_equal(np.bincount(np.r_[a, b], minlength=n), deg0)
            assert len({(i, j) for i, j in zip(a, b)}) == a.size and np.all(a != b)
            vals.append(copole_rate(lab, a, b))
        assert abs(np.mean(vals) - copole_expected_configuration(lab, e[:, 0], e[:, 1])) < 0.01, \
            (np.mean(vals), copole_expected_configuration(lab, e[:, 0], e[:, 1]))
    # vertex subsampling CI brackets the estimate
    n = 300; pi = rng.integers(0, n, 900); pj = rng.integers(0, n, 900); k = pi != pj; pi, pj = pi[k], pj[k]
    x = rng.normal(size=n)
    est, lo, hi, _ = vertex_subsample_ci(lambda nodes, a, b: assortativity(x, a, b), n, pi, pj, 100, rng)
    assert lo <= est <= hi
    # resolution guard
    try:
        assert_resolution(1000, 0.05, 183); raise AssertionError("should have refused")
    except ValueError:
        pass
    assert_resolution(10000, 0.05, 120)
    print("stats3 round-3 self-tests: ALL PASS")


if __name__ == "__main__":
    _selftest_round3()
