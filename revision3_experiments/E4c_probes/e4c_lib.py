"""E4c shared code: data, CV regimes, pair features, batched L2 logistic solver, fold metric, rewiring.

Sensitivity analysis (CONVENTION3 §5, E4 extras: "supervised probes (C1 random pairs, C2 held-out TFs,
C3 pair-disjoint; rewired-label null)"). Nothing here is an H1-H4 primary analysis.

Estimators taken from stats3 (not re-implemented): StratifiedPairs.auroc (within-TF stratified AUROC,
used for every reported fold metric), auroc (pooled fallback), assert_resolution, max_t_adjust.
A vectorised within-TF AUROC (`fast_within_tf_auroc`) is used ONLY to pick hyperparameters on inner
validation splits; it is asserted equal to stats3.StratifiedPairs.auroc on every call's first column.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import scipy.linalg as _sla

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import stats3 as S  # noqa: E402

CACHE = R3 / "cache"
EMB_LOCAL = Path.home() / "biomechinterp-local/subproject41_r3/emb"
REPO = R3.parents[1]
ESM2_NPZ = REPO / "subproject_60_uce_esm2_control_gpl_replication/embeddings/uce_esm2_human_gene_embeddings.npz"
STATIC_NPZ = REPO / "subproject_53_scgpt_gpl_replication/embeddings/scgpt_gene_embeddings.npz"
SEED = 20261026
N_FOLDS = 5
N_REPEATS = 5
MIN_POS = 3          # within-TF metric: TF needs >= 3 positives and >= 3 non-edges in the test fold
MIN_NEG = 3
ALPHAS = np.array([1e-4, 1e-3, 1e-2, 1e-1])   # M1 L2 grid (mean log-loss + alpha/2 ||w||^2, standardised features)


def conv_hash():
    return (R3 / "FROZEN.txt").read_text().split()[2]


# ===================================================================================== data
class Ctx:
    def __init__(self, name):
        self.name = name
        d = CACHE / name
        self.genes = json.loads((d / "obs_genes.json").read_text())
        self.obs_idx = np.load(d / "obs_idx.npy")
        z = np.load(d / "h3_pairs.npz")
        self.tf = z["tf"].astype(np.int64); self.tg = z["tg"].astype(np.int64); self.y = z["y"].astype(bool)
        self.n = self.y.size
        self.G = len(self.genes)
        nl = np.load(d / "node_labels.npz")
        self.is_tf_lambert = nl["is_tf_lambert"]; self.is_trrust_tf = nl["is_trrust_tf"]
        self.emb = np.load(EMB_LOCAL / f"{name}.npy", mmap_mode="r")
        assert self.emb.shape[1] > self.obs_idx.max()

    def layer(self, l):
        """V x 512 embedding of the observed genes, mean-centred over V (CONVENTION3 §2)."""
        X = np.asarray(self.emb[l], dtype=np.float64)[self.obs_idx]
        assert np.all(np.linalg.norm(X, axis=1) > 0), "never-observed row present"
        X = X - X.mean(0, keepdims=True)
        return X.astype(np.float32)

    def coexpr_gene_stats(self):
        """Per-gene z-scored expression on the extraction cells (expr_extraction.npy), mean and detection."""
        M = np.load(CACHE / self.name / "expr_extraction.npy").astype(np.float64)   # cells x V
        mu = M.mean(0); det = (M > 0).mean(0); sd = M.std(0)
        Z = np.where(sd > 0, (M - mu) / np.where(sd > 0, sd, 1), 0.0)
        return Z, mu, det

    def coexpr_features(self):
        Z, mu, det = self.coexpr_gene_stats()
        corr = np.einsum("ci,ci->i", Z[:, self.tf], Z[:, self.tg]) / Z.shape[0]
        F = np.column_stack([np.abs(corr), corr, mu[self.tf], mu[self.tg], det[self.tf], det[self.tg]])
        return F.astype(np.float32)

    def external_matrix(self, which):
        """Gene x d matrix for ESM2 ('esm2') or scGPT static token embeddings ('static'); missing -> NaN rows."""
        z = np.load(ESM2_NPZ if which == "esm2" else STATIC_NPZ)
        names = [str(s) for s in z["gene_names"]]; pos = {g: i for i, g in enumerate(names)}
        E = z["embeddings"]
        out = np.full((self.G, E.shape[1]), np.nan, dtype=np.float32)
        have = np.array([g in pos for g in self.genes])
        out[have] = E[[pos[g] for g, h in zip(self.genes, have) if h]]
        return out, have


# ===================================================================================== CV regimes
def h3_genes(ctx):
    return np.unique(np.r_[ctx.tf, ctx.tg])


def gene_partition(genes, tf_set, k, rng):
    """Stratified random partition of genes into k groups: TFs (genes that are h3 TFs) and the remaining
    genes are each shuffled and dealt round-robin, so every group holds ~1/k of the TFs."""
    lab = {}
    a = np.array([g for g in genes if g in tf_set]); b = np.array([g for g in genes if g not in tf_set])
    for arr in (a, b):
        arr = rng.permutation(arr); off = int(rng.integers(k))
        for r, g in enumerate(arr): lab[int(g)] = (r + off) % k
    return lab


def make_folds(ctx, regime, repeat):
    """List of (train_idx, test_idx) for one 5-fold repetition. Seeds are deterministic in (regime, repeat).
    C1: stratified random pairs. C2: TFs dealt to folds (balanced by pair count, GroupKFold semantics).
    C3: genes partitioned; train = both ends training genes, test = both ends held-out genes."""
    rid = {"C1": 1, "C2": 2, "C3": 3}[regime]
    rng = np.random.default_rng([SEED, rid, repeat])
    n = ctx.n; folds = []
    if regime == "C1":
        f = np.empty(n, np.int64)
        for cls in (True, False):
            idx = rng.permutation(np.flatnonzero(ctx.y == cls)); off = int(rng.integers(N_FOLDS))
            f[idx] = (np.arange(idx.size) + off) % N_FOLDS
        for k in range(N_FOLDS): folds.append((np.flatnonzero(f != k), np.flatnonzero(f == k)))
    elif regime == "C2":
        tfs, cnt = np.unique(ctx.tf, return_counts=True)
        order = rng.permutation(tfs.size)
        order = order[np.argsort(-cnt[order], kind="stable")]       # largest first, random tie order
        load = np.zeros(N_FOLDS); assign = {}
        for i in order:
            k = int(np.argmin(load + rng.random(N_FOLDS) * 1e-9)); assign[int(tfs[i])] = k; load[k] += cnt[i]
        f = np.array([assign[int(t)] for t in ctx.tf])
        for k in range(N_FOLDS): folds.append((np.flatnonzero(f != k), np.flatnonzero(f == k)))
    elif regime == "C3":
        lab = gene_partition(h3_genes(ctx), set(np.unique(ctx.tf).tolist()), N_FOLDS, rng)
        ft = np.array([lab[int(g)] for g in ctx.tf]); fg = np.array([lab[int(g)] for g in ctx.tg])
        for k in range(N_FOLDS):
            folds.append((np.flatnonzero((ft != k) & (fg != k)), np.flatnonzero((ft == k) & (fg == k))))
    return folds


def inner_split(ctx, regime, train_idx, seed):
    """One regime-matched inner train/validation split of a training fold (hyperparameter selection only)."""
    rng = np.random.default_rng([SEED, 99, seed])
    tf, tg, y = ctx.tf[train_idx], ctx.tg[train_idx], ctx.y[train_idx]
    if regime == "C1":
        v = np.zeros(train_idx.size, bool)
        for cls in (True, False):
            idx = rng.permutation(np.flatnonzero(y == cls)); v[idx[: idx.size // 4]] = True
        return np.flatnonzero(~v), np.flatnonzero(v)
    if regime == "C2":
        tfs = rng.permutation(np.unique(tf)); val = set(tfs[: max(1, tfs.size // 4)].tolist())
        v = np.array([t in val for t in tf])
        return np.flatnonzero(~v), np.flatnonzero(v)
    genes = np.unique(np.r_[tf, tg])
    lab = gene_partition(genes, set(np.unique(tf).tolist()), 4, rng)
    a = np.array([lab[int(g)] for g in tf]); b = np.array([lab[int(g)] for g in tg])
    return np.flatnonzero((a != 0) & (b != 0)), np.flatnonzero((a == 0) & (b == 0))


# ===================================================================================== features
FEATURES = ("concat", "hadamard", "absdiff", "bilinear")


def pca_basis(E, genes, k):
    """PCA of the rows `genes` of E (fitted on those genes only). Returns (mean, k x d components)."""
    X = np.asarray(E[genes], dtype=np.float64); mu = X.mean(0)
    _, _, Vt = np.linalg.svd(X - mu, full_matrices=False)
    return mu, Vt[:k]


def pair_features(kind, E, tf, tg, train_genes=None, extra=None):
    """Pair features from a gene x d matrix E (rows already V-centred / PCA-reduced as required).
    concat [x_i, x_j]; hadamard x_i*x_j; absdiff |x_i-x_j|; bilinear vec(p_i p_j^T) with p = PCA-16 of E
    fitted on train_genes only; tf_only [x_i]; tg_only [x_j]. `extra`: per-gene columns appended for both
    ends (missing-embedding indicators)."""
    if kind == "concat": F = np.hstack([E[tf], E[tg]])
    elif kind == "hadamard": F = E[tf] * E[tg]
    elif kind == "absdiff": F = np.abs(E[tf] - E[tg])
    elif kind == "bilinear":
        mu, C = pca_basis(E, train_genes, 16)
        Pj = ((E.astype(np.float64) - mu) @ C.T).astype(np.float32)
        F = (Pj[tf][:, :, None] * Pj[tg][:, None, :]).reshape(tf.size, 256)
    elif kind in ("tg_hadamard", "tg_bilinear"):          # sensitivity: target-only block + interaction block
        F = np.hstack([E[tg], pair_features(kind[3:], E, tf, tg, train_genes=train_genes)])
    elif kind == "tf_only": F = E[tf]
    elif kind == "tg_only": F = E[tg]
    else: raise ValueError(kind)
    if extra is not None: F = np.hstack([F, extra[tf], extra[tg]])
    return np.ascontiguousarray(F, dtype=np.float32)


def standardise(F, tr):
    mu = F[tr].mean(0); sd = F[tr].std(0); sd[sd == 0] = 1.0
    return ((F - mu) / sd).astype(np.float32), mu, sd


def gene_block(E, idx_train):
    """Gene-level block standardised by the pair-weighted training mean / sd of E[idx_train]
    (identical to standardising the materialised pair columns)."""
    Xb = E[idx_train].astype(np.float64); mu = Xb.mean(0); sd = Xb.std(0); sd[sd == 0] = 1.0
    return ((E - mu) / sd).astype(np.float32)


# ===================================================================================== solver
def _softplus(z):
    return np.logaddexp(0, z)


class Dense:
    """Design matrix held explicitly (n x d float32)."""
    def __init__(self, X):
        self.X = np.ascontiguousarray(X, dtype=np.float32); self.n, self.d = X.shape
    def mv(self, W): return self.X @ W.astype(np.float32)
    def rmv(self, R): return (self.X.T @ R).astype(np.float64)
    def gram(self, dv):
        Xd = self.X * dv[:, None].astype(np.float32)
        return (self.X.T @ Xd).astype(np.float64), Xd.sum(0, dtype=np.float64)
    def rows(self, idx): return Dense(self.X[idx])


class Gather:
    """Design whose rows are concatenations of gene-level blocks: row p = [B_1[idx_1[p]], ..., B_m[idx_m[p]]]
    (e.g. concat = [x_tf, x_tg]; tf_only = [x_tf]). Products are computed at gene level (exact, cheaper)."""
    def __init__(self, blocks, idxs, G):
        import scipy.sparse as sp
        self.B = [np.ascontiguousarray(b, dtype=np.float32) for b in blocks]; self.I = [np.asarray(i) for i in idxs]
        self.G = G; self.n = self.I[0].size; self.dims = [b.shape[1] for b in self.B]; self.d = sum(self.dims)
        self.off = np.r_[0, np.cumsum(self.dims)]
        self.Sm = [sp.csr_matrix((np.ones(self.n, np.float32), (i, np.arange(self.n))), shape=(G, self.n)) for i in self.I]
    def mv(self, W):
        W = W.astype(np.float32); Z = 0
        for b, (B, I) in enumerate(zip(self.B, self.I)): Z = Z + (B @ W[self.off[b]:self.off[b + 1]])[I]
        return Z
    def rmv(self, R):
        R = np.asarray(R, dtype=np.float32)
        return np.vstack([(B.T @ (Sm @ R)).astype(np.float64) for B, Sm in zip(self.B, self.Sm)])
    def gram(self, dv):
        import scipy.sparse as sp
        dv = dv.astype(np.float32); m = len(self.B); H = np.zeros((self.d, self.d)); cs = []
        for a in range(m):
            cs.append((self.B[a].T @ (self.Sm[a] @ dv)).astype(np.float64))
            for c in range(a, m):
                M = (self.Sm[a] @ sp.diags(dv) @ self.Sm[c].T).tocsr()          # G x G
                blk = (self.B[a].T @ (M @ self.B[c])).astype(np.float64)
                H[self.off[a]:self.off[a + 1], self.off[c]:self.off[c + 1]] = blk
                if c != a: H[self.off[c]:self.off[c + 1], self.off[a]:self.off[a + 1]] = blk.T
        return H, np.concatenate(cs)
    def rows(self, idx): return Gather(self.B, [i[idx] for i in self.I], self.G)


# ---- fused elementwise kernels (numba): the solver's cost is dominated by elementwise work on n x K arrays
_KERN = {}


def _kernels():
    if _KERN: return _KERN
    import numba as nb

    @nb.njit(parallel=True, fastmath=False, cache=True)
    def loss_multi(Z, U, Y, etas):
        n, K = Z.shape; E = etas.size; nch = 64
        acc = np.zeros((nch, E, K))
        step = (n + nch - 1) // nch
        for ch in nb.prange(nch):
            lo = ch * step; hi = min(n, lo + step)
            for i in range(lo, hi):
                for k in range(K):
                    z = Z[i, k]; u = U[i, k]; y = Y[i, k]
                    for e in range(E):
                        t = z + etas[e] * u
                        sp = max(t, 0.0) + np.log1p(np.exp(-abs(t)))
                        acc[ch, e, k] += sp - y * t
        return acc.sum(0) / n

    @nb.njit(parallel=True, cache=True)
    def resid(Z, Y):
        n, K = Z.shape; R = np.empty((n, K), np.float32); d = np.empty(n)
        for i in nb.prange(n):
            s = 0.0
            for k in range(K):
                p = 1.0 / (1.0 + np.exp(-Z[i, k])); R[i, k] = (p - Y[i, k]) / n; s += p * (1 - p)
            d[i] = s / K
        return R, d

    @nb.njit(parallel=True, cache=True)
    def axpy_cols(Z, U, beta):
        n, K = Z.shape
        for i in nb.prange(n):
            for k in range(K): Z[i, k] += U[i, k] * beta[k]
    _KERN.update(loss_multi=loss_multi, resid=resid, axpy_cols=axpy_cols)
    return _KERN


def fit_logistic_batched(D, Y, alpha, max_iter=40, tol=1e-7, refresh=5):
    """L2 logistic regression for K label columns at once on a shared design D (Dense, Gather or array).
    Objective per column: mean_i[softplus(z) - y z] + alpha_k/2 ||w_k||^2 (intercept unpenalised).
    Quasi-Newton: Hessian shared across columns, X^T diag(dbar) X / n with dbar the column-averaged
    sigma(1-sigma), refreshed every `refresh` iterations; per-column step length chosen from
    {1.5,1,0.5,0.25,0.1} by exact loss evaluation (monotone decrease). Returns W (d x K), b (K,), n_iter."""
    if isinstance(D, np.ndarray): D = Dense(D)
    n, d = D.n, D.d; K = Y.shape[1]; alpha = np.broadcast_to(np.asarray(alpha, float), (K,)).copy()
    Y = np.ascontiguousarray(Y, dtype=np.float32)
    p0 = np.clip(Y.mean(0), 1e-6, 1 - 1e-6)
    W = np.zeros((d, K), np.float64); b = np.log(p0 / (1 - p0))
    Z = np.broadcast_to(b.astype(np.float32), (n, K)).copy()
    kn = _kernels()
    L = kn["loss_multi"](Z, np.zeros_like(Z), Y, np.zeros(1))[0] + 0.5 * alpha * (W ** 2).sum(0)
    ualphas = np.unique(alpha); active = np.ones(K, bool)
    etas = np.array([1.5, 1.0, 0.5, 0.25, 0.1]); chol = {}
    it = 0
    for it in range(1, max_iter + 1):
        R, dbar = kn["resid"](Z, Y)                                             # n x K, n
        gW = D.rmv(R) + alpha * W; gb = R.sum(0, dtype=np.float64)
        if (it - 1) % refresh == 0:
            H, hx = D.gram(dbar); H /= n; hx /= n; hbb = dbar.sum() / n; chol = {}
            for a in ualphas:
                Ha = np.empty((d + 1, d + 1)); Ha[:d, :d] = H + a * np.eye(d); Ha[:d, d] = hx; Ha[d, :d] = hx
                Ha[d, d] = hbb; Ha[np.diag_indices(d + 1)] += 1e-10
                try: chol[a] = ("c", _sla.cho_factor(Ha, lower=True, check_finite=False))
                except np.linalg.LinAlgError: chol[a] = ("p", np.linalg.pinv(Ha))
        DW = np.zeros_like(W); Db = np.zeros_like(b)
        for a in ualphas:
            cols = np.flatnonzero((alpha == a) & active)
            if cols.size == 0: continue
            g = np.vstack([gW[:, cols], gb[cols][None, :]]); kind, c = chol[a]
            sol = _sla.cho_solve(c, g, check_finite=False) if kind == "c" else c @ g
            DW[:, cols] = -sol[:d]; Db[cols] = -sol[d]
        U = np.ascontiguousarray(D.mv(DW) + Db.astype(np.float32), dtype=np.float32)   # n x K
        Le = kn["loss_multi"](Z, U, Y, etas)                                     # n_eta x K
        Le = Le + 0.5 * alpha * ((W[None] + etas[:, None, None] * DW[None]) ** 2).sum(1)
        best = L.copy(); beta = np.zeros(K)
        for e in range(etas.size):
            imp = Le[e] < best; best[imp] = Le[e][imp]; beta[imp] = etas[e]
        rel = (L - best) / np.maximum(np.abs(L), 1e-12)
        W += beta * DW; b += beta * Db; kn["axpy_cols"](Z, U, beta.astype(np.float32)); L = best
        active = rel > tol
        if not active.any(): break
    return W, b, it


def predict_linear(D, W, b):
    if isinstance(D, np.ndarray): D = Dense(D)
    return D.mv(W) + b.astype(np.float32)


# ===================================================================================== fold metric
def fold_metric(tf, tg, y, s):
    """Reported fold metric (stats3 estimators). Within-TF stratified AUROC over the TFs with >= 3 positives
    and >= 3 candidate non-edges in the fold; if no TF qualifies, pooled AUROC over all fold pairs.
    Returns (auc, metric_type, n_tf_used, n_pos_used, n_neg_used)."""
    tfs, inv, = np.unique(tf, return_inverse=True)
    npos = np.bincount(inv, weights=y, minlength=tfs.size); nneg = np.bincount(inv, weights=~y, minlength=tfs.size)
    q = (npos >= MIN_POS) & (nneg >= MIN_NEG)
    if q.any():
        m = q[inv]
        sp = S.StratifiedPairs(tf[m], tg[m], y[m])
        return sp.auroc(s[m]), "within_tf", int(q.sum()), int(y[m].sum()), int((~y[m]).sum())
    if y.any() and (~y).any():
        return S.auroc(s[y], s[~y]), "pooled", int(tfs.size), int(y.sum()), int((~y).sum())
    return float("nan"), "none", 0, int(y.sum()), int((~y).sum())


def fast_within_tf_auroc(tf, Ys, Sc):
    """Vectorised within-TF AUROC (same rule as fold_metric) for score columns Sc (n x K) with label columns
    Ys (n x K, bool). Used only for inner hyperparameter selection."""
    from scipy.stats import rankdata
    tfs, inv = np.unique(tf, return_inverse=True); K = Sc.shape[1]
    num = np.zeros(K); den = np.zeros(K); pnum = None
    order = np.argsort(inv, kind="stable"); bounds = np.r_[0, np.cumsum(np.bincount(inv, minlength=tfs.size))]
    for c in range(tfs.size):
        idx = order[bounds[c]:bounds[c + 1]]
        yc = Ys[idx]; n1 = yc.sum(0); n0 = idx.size - n1
        ok = (n1 >= MIN_POS) & (n0 >= MIN_NEG)
        if not ok.any(): continue
        r = rankdata(Sc[idx], axis=0)
        num += np.where(ok, (r * yc).sum(0) - n1 * (n1 + 1) / 2, 0); den += np.where(ok, n1 * n0, 0)
    out = np.full(K, np.nan); has = den > 0; out[has] = num[has] / den[has]
    if (~has).any():                                    # pooled fallback
        for k in np.flatnonzero(~has):
            yk = Ys[:, k]
            if yk.any() and (~yk).any(): out[k] = S.auroc(Sc[yk, k], Sc[~yk, k])
    return out


def check_fast(tf, tg, Ys, Sc, fast):
    a = fold_metric(tf, tg, Ys[:, 0], Sc[:, 0].astype(np.float64))[0]
    assert (np.isnan(a) and np.isnan(fast[0])) or abs(a - fast[0]) < 1e-6, (a, fast[0])


# ===================================================================================== rewiring
def _rewire_numba():
    import numba as nb

    @nb.njit(cache=True)
    def chain(tf, tg, y, keys, kidx, G, n_acc, max_try, seed):
        np.random.seed(seed)
        pos = np.flatnonzero(y); npos = pos.size; acc = 0; tried = 0
        while acc < n_acc and tried < max_try:
            tried += 1
            i1 = np.random.randint(npos); i2 = np.random.randint(npos)
            p1 = pos[i1]; p2 = pos[i2]
            a = tf[p1]; b = tg[p1]; c = tf[p2]; d = tg[p2]
            if a == c or b == d: continue
            k1 = a * G + d; k2 = c * G + b
            j1 = np.searchsorted(keys, k1); j2 = np.searchsorted(keys, k2)
            if j1 >= keys.size or keys[j1] != k1 or j2 >= keys.size or keys[j2] != k2: continue
            q1 = kidx[j1]; q2 = kidx[j2]
            if y[q1] or y[q2]: continue
            y[p1] = False; y[p2] = False; y[q1] = True; y[q2] = True
            pos[i1] = q1; pos[i2] = q2; acc += 1
        return acc, tried
    return chain


def rewire_draws(ctx, n_draws, seed, burn=10, thin=1):
    """Degree-preserving rewired label vectors on the h3 candidate universe: MCMC swaps
    (tf1,t1),(tf2,t2) -> (tf1,t2),(tf2,t1), accepted only when both new pairs are candidate pairs of h3 that
    are currently unlabelled; every TF's and every target's positive count is preserved exactly.
    Burn-in burn*n_pos accepted swaps, thin*n_pos accepted swaps between successive draws."""
    chain = _rewire_numba()
    keys = ctx.tf * ctx.G + ctx.tg; assert np.unique(keys).size == keys.size
    o = np.argsort(keys); ks = keys[o]; kidx = o.astype(np.int64)
    y = ctx.y.copy(); npos = int(y.sum()); out = []; stats = []
    acc, tried = chain(ctx.tf, ctx.tg, y, ks, kidx, ctx.G, burn * npos, 2000 * burn * npos, seed)
    stats.append(("burn", acc, tried))
    for d in range(n_draws):
        acc, tried = chain(ctx.tf, ctx.tg, y, ks, kidx, ctx.G, thin * npos, 2000 * thin * npos, seed + 1 + d)
        out.append(y.copy()); stats.append((d, acc, tried))
    return np.array(out), stats


def log(msg, fh=None):
    s = time.strftime("%H:%M:%S ") + msg
    print(s, flush=True)
    if fh is not None:
        with open(fh, "a") as f: f.write(s + "\n")
