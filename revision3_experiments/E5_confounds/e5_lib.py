"""E5 shared code: data assembly, an exact (numerically stabilised) conditional-logit engine, bootstrap and
permutation draws, within-TF residualisation.  Read-only use of contexts.py, stats3.py and cache/.

Conditional logistic regression (stratum = TF)
---------------------------------------------
For a stratum with n candidate pairs of which k are TRRUST edges, the conditional likelihood is
    L = exp(sum_{cases} eta_i) / e_k(exp(eta_1), ..., exp(eta_n)),     eta = X beta,
with e_k the k-th elementary symmetric polynomial (the denominator sums over all k-subsets). This is exactly the
likelihood of statsmodels' ConditionalLogit. statsmodels evaluates e_k by an unscaled memoised Python recursion
of depth n; for the lung / lung_external / kidney_B strata (n up to 1,580, k up to 367) the unscaled e_k
exceeds 1e308 (C(1579,367) ~ e^853) and the recursion exceeds Python's depth limit, so statsmodels cannot fit
these contexts. The engine below evaluates the SAME likelihood with the binomially normalised recursion
    f_i[j] = ((i-j)/i) f_{i-1}[j] + (j/i) w_i f_{i-1}[j-1],      e_k = C(n,k) f_n[k],
(a convex combination, so it neither overflows nor underflows for moderate eta spread; an explicit rescale guard
is added), and carries first and second derivatives with respect to beta forward through the recursion
(forward-mode), giving the exact gradient and exact Hessian. It is validated against statsmodels on every
context/layer where statsmodels runs (immune, kidney_A, and the small-stratum subset of the large contexts):
see e5_validate.py / validation_statsmodels.csv.

Stratum weights (TF bootstrap multiplicity) multiply a stratum's log-likelihood; target bootstrap multiplicity is
implemented by physically replicating rows inside the stratum (a target drawn twice is two candidate pairs),
which is what refitting on the pigeonhole-resampled data means for a conditional likelihood.
"""
from __future__ import annotations
import hashlib, json, pickle, sys, math
from pathlib import Path
import numpy as np
import numba as nb

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import contexts as X   # noqa: E402
import stats3 as S     # noqa: E402

LOCAL_EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
GENE2GO = X.REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl"
CTXS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
CTX_ID = {c: i for i, c in enumerate(CTXS)}
N_LAYERS = 12
SEED = 20261026
B_BOOT = 500
N_PERM = 500
COEXP_EPS = 1e-3          # log(|r| + eps): floor for exactly-zero correlations
COV_NAMES = ["log_abs_coexp", "log_mean_expr_target", "log_nobs_target", "string_deg_target", "go_count_target"]
TF_COV_NAMES = ["log_mean_expr_tf", "log_nobs_tf", "string_deg_tf", "go_count_tf"]   # constant within stratum


def convention_hash():
    txt = (R3 / "FROZEN.txt").read_text().split()
    frozen = txt[txt.index("sha256") + 1]
    actual = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
    return frozen, actual, frozen == actual


def load_emb(ctx):
    p = LOCAL_EMB / f"{ctx}.npy"
    return np.load(p if p.exists() else X.CONTEXTS[ctx]["emb"], mmap_mode="r")


# ----------------------------------------------------------------------------------------------- node covariates
def string_degree_genomewide():
    """Genome-wide STRING v12 physical (combined >= 700) degree per symbol; cached in this directory."""
    p = HERE / "string_physical700_degree_genomewide.json"
    if p.exists():
        return json.loads(p.read_text())
    pairs = X.load_string_genomewide(700)
    deg = {}
    for a, b in pairs:
        deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
    p.write_text(json.dumps(dict(n_pairs=len(pairs), n_genes=len(deg), degree=deg)))
    return json.loads(p.read_text())


def go_counts(genes):
    with open(GENE2GO, "rb") as f:
        g2g = pickle.load(f)
    return np.array([len(g2g.get(g, ())) for g in genes], dtype=float)


def load_context(ctx):
    """Pairs, node/pair covariates (raw), and the extraction-cell expression for one context."""
    d = R3 / "cache" / ctx
    genes = json.loads((d / "obs_genes.json").read_text())
    obs_idx = np.load(d / "obs_idx.npy")
    nobs = np.load(d / "ncells.npy").astype(float)
    P = np.load(d / "h3_pairs.npz")
    tf, tg, y = P["tf"].astype(np.int64), P["tg"].astype(np.int64), P["y"].astype(np.int8)
    E = np.load(d / "expr_extraction.npy").astype(np.float64)          # extraction cells x observed genes
    assert E.shape[1] == len(genes)
    mean_expr = E.mean(0)
    assert np.all(mean_expr > 0), "gene with zero mean expression on the extraction cells"
    Z = (E - E.mean(0)) / E.std(0)
    r = np.einsum("ij,ij->j", Z[:, tf], Z[:, tg]) / E.shape[0]        # Pearson r on extraction cells
    sd = string_degree_genomewide()["degree"]
    deg = np.array([sd.get(g, 0) for g in genes], dtype=float)
    go = go_counts(genes)
    node = dict(log_mean_expr=np.log(mean_expr), log_nobs=np.log(nobs), string_deg=deg, go_count=go)
    cov_raw = np.column_stack([np.log(np.abs(r) + COEXP_EPS), node["log_mean_expr"][tg], node["log_nobs"][tg],
                               deg[tg], go[tg]])
    tf_cov_raw = np.column_stack([node["log_mean_expr"][tf], node["log_nobs"][tf], deg[tf], go[tf]])
    return dict(ctx=ctx, genes=genes, obs_idx=obs_idx, tf=tf, tg=tg, y=y, coexp_r=r, cov_raw=cov_raw,
                tf_cov_raw=tf_cov_raw, node=node, n_cells_extraction=int(E.shape[0]))


def layer_scores(emb, obs_idx, L, tf, tg):
    """Signed cosine of the V-centred full-dimensional embeddings (CONVENTION3 §2)."""
    Em = np.asarray(emb[L], dtype=np.float64)[obs_idx]
    Em = Em - Em.mean(0, keepdims=True)
    En = S.l2norm_rows(Em)
    return np.einsum("ij,ij->i", En[tf], En[tg])


def zscore(x):
    x = np.asarray(x, float)
    if x.ndim == 1:
        return (x - x.mean()) / x.std()
    return (x - x.mean(0)) / x.std(0)


# ----------------------------------------------------------------------------------------------- clogit engine
@nb.njit(cache=True)
def _group_eval(Xg, yg, beta, want_hess):
    n, p = Xg.shape
    k = 0
    for i in range(n):
        k += yg[i]
    eta = Xg @ beta
    m = eta.max()
    s = 0.0
    for i in range(n):
        s += math.exp(eta[i] - m)
    c = m + math.log(s / n)
    f = np.zeros(k + 1); f[0] = 1.0
    g = np.zeros((k + 1, p))
    H = np.zeros((k + 1, p, p))
    logscale = 0.0
    for i in range(1, n + 1):
        w = math.exp(eta[i - 1] - c)
        x = Xg[i - 1]
        jmax = i if i < k else k
        for j in range(jmax, 0, -1):
            a = (i - j) / i; b = j / i * w
            fj1 = f[j - 1]
            if want_hess:
                for u in range(p):
                    gu = g[j - 1, u]; xu = x[u]
                    for v in range(u, p):
                        H[j, u, v] = a * H[j, u, v] + b * (H[j - 1, u, v] + xu * g[j - 1, v] + gu * x[v] + xu * x[v] * fj1)
            for u in range(p):
                g[j, u] = a * g[j, u] + b * (g[j - 1, u] + x[u] * fj1)
            f[j] = a * f[j] + b * fj1
        mx = 0.0
        for j in range(jmax + 1):
            if f[j] > mx: mx = f[j]
        if mx > 1e100 or mx < 1e-100:
            f /= mx; g /= mx
            if want_hess: H /= mx
            logscale += math.log(mx)
    fk = f[k]
    lbin = math.lgamma(n + 1.0) - math.lgamma(k + 1.0) - math.lgamma(n - k + 1.0)
    log_ek = lbin + math.log(fk) + logscale + k * c
    ll = -log_ek
    grad = np.zeros(p)
    for i in range(n):
        if yg[i]:
            ll += eta[i]
            for u in range(p):
                grad[u] += Xg[i, u]
    mu = g[k] / fk
    grad -= mu
    hess = np.zeros((p, p))
    if want_hess:
        for u in range(p):
            for v in range(u, p):
                hv = -(H[k, u, v] / fk - mu[u] * mu[v])
                hess[u, v] = hv; hess[v, u] = hv
    return ll, grad, hess


@nb.njit(parallel=True, cache=True)
def _total_eval(Xs, ys, gstart, gw, beta, want_hess):
    G = gstart.size - 1; p = Xs.shape[1]
    lls = np.zeros(G); grads = np.zeros((G, p)); hesss = np.zeros((G, p, p))
    for gi in nb.prange(G):
        a = gstart[gi]; b = gstart[gi + 1]
        ll, gr, he = _group_eval(Xs[a:b], ys[a:b], beta, want_hess)
        lls[gi] = ll * gw[gi]; grads[gi] = gr * gw[gi]; hesss[gi] = he * gw[gi]
    return lls.sum(), grads.sum(0), hesss.sum(0)


class CLogitData:
    """Rows sorted by stratum; strata without both labels dropped. gw = stratum weights."""
    def __init__(self, X, y, groups, gw_by_group=None):
        X = np.asarray(X, float)
        if X.ndim == 1: X = X[:, None]
        y = np.asarray(y).astype(np.int64); groups = np.asarray(groups)
        order = np.argsort(groups, kind="mergesort")
        X, y, groups = X[order], y[order], groups[order]
        ug, start = np.unique(groups, return_index=True)
        end = np.r_[start[1:], groups.size]
        keep_rows = np.zeros(groups.size, bool); kept_g = []; bounds = []
        for gid, a, b in zip(ug, start, end):
            s = y[a:b].sum()
            if 0 < s < (b - a):
                keep_rows[a:b] = True; kept_g.append(gid)
        self.X = np.ascontiguousarray(X[keep_rows]); self.y = y[keep_rows]; g2 = groups[keep_rows]
        _, st = np.unique(g2, return_index=True)
        self.gstart = np.r_[st, g2.size].astype(np.int64)
        self.groups = np.array(kept_g)
        if gw_by_group is None: self.gw = np.ones(self.groups.size)
        else: self.gw = np.asarray([gw_by_group[g] for g in self.groups], float)
        self.p = self.X.shape[1]

    def eval(self, beta, hess=True):
        return _total_eval(self.X, self.y, self.gstart, self.gw, np.asarray(beta, float), hess)

    def ll0(self):
        return self.eval(np.zeros(self.p), hess=False)[0]


def clogit_fit(D: CLogitData, beta0=None, tol=1e-9, maxit=60):
    """Newton-Raphson with step halving on the exact conditional log-likelihood."""
    beta = np.zeros(D.p) if beta0 is None else np.asarray(beta0, float).copy()
    ll, gr, he = D.eval(beta)
    conv = False; it = 0
    for it in range(1, maxit + 1):
        try:
            step = np.linalg.solve(-he, gr)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(-he, gr, rcond=None)[0]
        t = 1.0
        for _ in range(40):
            nb_ = beta + t * step
            ll_n, gr_n, he_n = D.eval(nb_)
            if np.isfinite(ll_n) and ll_n >= ll - 1e-10: break
            t /= 2
        dmax = np.max(np.abs(nb_ - beta)); dll = ll_n - ll
        beta, ll, gr, he = nb_, ll_n, gr_n, he_n
        if dmax < tol or (abs(dll) < 1e-12 and np.max(np.abs(gr)) < 1e-6):
            conv = True; break
    try:
        cov = np.linalg.inv(-he)
    except np.linalg.LinAlgError:
        cov = np.full((D.p, D.p), np.nan)
    return dict(beta=beta, ll=float(ll), cov=cov, se=np.sqrt(np.clip(np.diag(cov), 0, None)),
                converged=conv, iters=it, max_grad=float(np.max(np.abs(gr))))


# ----------------------------------------------------------------------------------------------- resampling
def tf_target_codes(tf, tg):
    tfs, tfc = np.unique(tf, return_inverse=True)
    tgs, tgc = np.unique(tg, return_inverse=True)
    return tfs, tfc, tgs, tgc


def boot_weights(ctx, b, n_tf, n_tg):
    """Pigeonhole two-way draw b (identical for every layer of the context)."""
    rng = np.random.default_rng([SEED, 5, CTX_ID[ctx], b])
    wtf = np.bincount(rng.integers(0, n_tf, n_tf), minlength=n_tf)
    wtg = np.bincount(rng.integers(0, n_tg, n_tg), minlength=n_tg)
    return wtf, wtg


def perm_labels(ctx, r, sp):
    """Within-TF label permutation r (identical for every layer of the context)."""
    rng = np.random.default_rng([SEED, 55, CTX_ID[ctx], r])
    return sp.permute_within_tf(rng)


def expand(X, y, tfc, wtf, wtg, tgc):
    """Rows replicated by target multiplicity; strata weighted by TF multiplicity."""
    mult = wtg[tgc] * (wtf[tfc] > 0)
    idx = np.repeat(np.arange(y.size), mult)
    return X[idx], y[idx], tfc[idx]


# ----------------------------------------------------------------------------------------------- residualisation
def within_center(v, codes, w=None):
    v = np.asarray(v, float); w = np.ones(codes.size) if w is None else np.asarray(w, float)
    G = codes.max() + 1
    sw = np.bincount(codes, weights=w, minlength=G)
    if v.ndim == 1:
        m = np.bincount(codes, weights=w * v, minlength=G) / np.where(sw > 0, sw, 1)
        return v - m[codes]
    out = np.empty_like(v)
    for c in range(v.shape[1]):
        m = np.bincount(codes, weights=w * v[:, c], minlength=G) / np.where(sw > 0, sw, 1)
        out[:, c] = v[:, c] - m[codes]
    return out


def residualise_within_tf(score, cov, codes, w=None):
    """OLS of the within-TF-centred score on the within-TF-centred covariates (weighted when w is given);
    returns the residual (TF-level covariates are absorbed by the centring)."""
    w = np.ones(score.size) if w is None else np.asarray(w, float)
    sc = within_center(score, codes, w); Xc = within_center(cov, codes, w)
    sw = np.sqrt(w)
    beta = np.linalg.lstsq(Xc * sw[:, None], sc * sw, rcond=None)[0]
    return sc - Xc @ beta, beta


def clogit_fit_fast(D: CLogitData, beta0=None, tol=1e-7, maxit=60):
    """Point estimate only (no SE): exact Hessian at the start point, then chord (fixed-Hessian Newton) steps that
    need only the gradient (about a quarter of the cost of a Hessian evaluation). Falls back to clogit_fit if a chord
    step decreases the log-likelihood or the iteration does not converge. Converges to the same maximiser as
    clogit_fit (checked in e5_check_fast.py); used for bootstrap and permutation refits after layer 2 of the large
    contexts (see SUMMARY)."""
    beta = np.zeros(D.p) if beta0 is None else np.asarray(beta0, float).copy()
    ll, gr, he = D.eval(beta, hess=True)
    try:
        Hinv = np.linalg.inv(-he)
    except np.linalg.LinAlgError:
        return clogit_fit(D, beta0)
    for it in range(1, maxit + 1):
        step = Hinv @ gr
        nb_ = beta + step
        ll_n, gr_n, _ = D.eval(nb_, hess=False)
        if not np.isfinite(ll_n) or ll_n < ll - 1e-10:
            return clogit_fit(D, beta)
        beta, ll, gr = nb_, ll_n, gr_n
        if np.max(np.abs(step)) < tol:
            return dict(beta=beta, ll=float(ll), converged=True, iters=it, max_grad=float(np.max(np.abs(gr))))
    return clogit_fit(D, beta)
