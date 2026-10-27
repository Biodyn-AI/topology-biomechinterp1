"""E3 shared helpers: data loading and vectorised versions of the stats3 node statistics.

The definitions are stats3's (tf_scores_loo, stratified_auroc, permute_within_strata, max_t_adjust,
signed_free_auroc). The column-wise versions below exist only for speed inside permutation nulls and are
checked against the stats3 scalar functions on the observed labels in every run (assert |diff| < 1e-9).
"""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
sys.dont_write_bytecode = True
import stats3 as S          # noqa: E402

CACHE = R3 / "cache"
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
STATIC_SCGPT = R3.parents[1] / "subproject_53_scgpt_gpl_replication/embeddings/scgpt_gene_embeddings.npz"
STATIC_ESM2 = R3.parents[1] / "subproject_60_uce_esm2_control_gpl_replication/embeddings/uce_esm2_human_gene_embeddings.npz"
CONTEXTS = ["lung", "lung_external", "kidney_B", "immune", "kidney_A", "immune_seed43", "immune_seed44",
            "immune_TSP14", "immune_TSP21", "immune_TSP25", "immune_Bcell", "immune_CD4T", "immune_CD8T",
            "kidney_ckpt_whole_human", "kidney_ckpt_kidney", "kidney_ckpt_brain"]
DISCOVERY = {"immune", "kidney_A"}


def conv_hash():
    return hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()


def load(ctx):
    d = CACHE / ctx
    genes = json.loads((d / "obs_genes.json").read_text())
    obs_idx = np.load(d / "obs_idx.npy")
    nl = np.load(d / "node_labels.npz")
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    ee = np.load(d / "expr_extraction.npy", mmap_mode="r"); ex = np.load(d / "expr.npy", mmap_mode="r")
    out = dict(genes=genes, obs_idx=obs_idx, emb=emb, nobs=np.load(d / "ncells.npy").astype(float),
               log_nobs=nl["log_nobs"], det=nl["det_decile"], is_tf=nl["is_tf_lambert"].astype(bool),
               is_trrust_tf=nl["is_trrust_tf"].astype(bool), is_trrust_tgt=nl["is_trrust_target_only"].astype(bool),
               mean_extr=np.asarray(ee, dtype=np.float64).mean(0), mean_expr=np.asarray(ex, dtype=np.float64).mean(0))
    assert len(genes) == obs_idx.size == out["is_tf"].size
    return out


def centred(emb, layer, obs_idx, sub=None):
    """V-centred layer matrix (float64); if sub is given, the rows of sub only, centred over sub."""
    X = np.asarray(emb[layer], dtype=np.float64)[obs_idx]
    if sub is not None: X = X[sub]
    return X - X.mean(0, keepdims=True)


def log_mean(m):
    return np.log(m + 1e-3)


# ------------------------------------------------------------------ column-wise stratified AUROC
class Strata:
    def __init__(self, strata, is_pos):
        strata = np.asarray(strata); is_pos = np.asarray(is_pos, bool)
        self.groups = []; self.w = []
        for st in np.unique(strata):
            idx = np.flatnonzero(strata == st); n1 = int(is_pos[idx].sum()); n0 = idx.size - n1
            if n1 and n0: self.groups.append((idx, n1, n0)); self.w.append(n1 * n0)
        self.W = float(sum(self.w))


def strat_cols_varying(SC, Y, st: Strata):
    """Pooled within-stratum AUROC per column; scores and labels both vary by column (continuous scores,
    no tie correction). Counts per stratum are those of st (within-stratum permutation preserves them)."""
    num = np.zeros(SC.shape[1])
    for idx, n1, n0 in st.groups:
        sub = SC[idx]; ysub = Y[idx]
        order = np.argsort(sub, axis=0, kind="stable")
        ypos = np.take_along_axis(ysub, order, axis=0)
        rsum = (ypos * np.arange(1, idx.size + 1, dtype=np.float64)[:, None]).sum(0)
        num += rsum - n1 * (n1 + 1) / 2
    return num / st.W


def strat_ranks_fixed(score, st: Strata, n):
    """Within-stratum average ranks of a fixed score (tie-correct)."""
    r = np.zeros(n)
    for idx, n1, n0 in st.groups:
        r[idx] = rankdata(score[idx])
    return r


def strat_cols_fixed(ranks, Y, st: Strata):
    """Pooled within-stratum AUROC of a fixed score under label matrices Y (n x B); ties handled via
    average ranks. Exact Mann-Whitney per stratum."""
    num = np.zeros(Y.shape[1])
    for idx, n1, n0 in st.groups:
        num += ranks[idx] @ Y[idx] - n1 * (n1 + 1) / 2
    return num / st.W


def loo_scores_cols(En, Y, k):
    """stats3.tf_scores_loo for every column of the label matrix Y (n x B, 0/1), En row-normalised."""
    tot = En @ (En.T @ Y)
    self_ = (En * En).sum(1)[:, None]
    return np.where(Y > 0, (tot - self_) / max(k - 1, 1), tot / k)


def perm_batch(is_pos, strata, rng, b):
    return np.stack([S.permute_within_strata(is_pos, strata, rng) for _ in range(b)], axis=1)


def null_summary(null, obs=None):
    null = np.asarray(null, float)
    d = dict(N=int(null.size), mean=float(null.mean()), sd=float(null.std(ddof=1)),
             q=dict(zip(["0.5", "2.5", "50", "97.5", "99.5"], map(float, np.quantile(null, [.005, .025, .5, .975, .995])))),
             max=float(null.max()))
    if obs is not None:
        d["b"] = int(np.sum(null >= obs))
    return d


def boot_strat(score, is_pos, strata, B, rng):
    """Gene bootstrap (scores fixed, genes resampled with replacement) of stats3.stratified_auroc and auroc."""
    n = score.size; out = np.empty((B, 2))
    for r in range(B):
        i = rng.integers(0, n, n)
        out[r, 0] = S.stratified_auroc(score[i], is_pos[i], strata[i])
        p = is_pos[i]; out[r, 1] = S.auroc(score[i][p], score[i][~p])
    return out


def bands12():
    return [(i, j) for i in range(12) for j in range(i, 12)]       # 0-based inclusive; 78 bands


# ------------------------------------------------------------------ numba kernel (2x faster than numpy)
import numba as _nb


@_nb.njit(parallel=True, cache=False)
def _strat_nb(SCT, YT, starts, idxcat, n1s):
    B = SCT.shape[0]; out = np.zeros(B)
    for b in _nb.prange(B):
        tot = 0.0
        for g in range(starts.size - 1):
            idx = idxcat[starts[g]:starts[g + 1]]
            v = SCT[b, idx]; o = np.argsort(v)
            r = 0.0
            for p in range(o.size):
                if YT[b, idx[o[p]]] > 0: r += p + 1
            tot += r - n1s[g] * (n1s[g] + 1) / 2
        out[b] = tot
    return out


class FastStrata(Strata):
    """Strata + numba column AUROC. Same numbers as strat_cols_varying (checked in self_check)."""
    def __init__(self, strata, is_pos):
        super().__init__(strata, is_pos)
        self.idxcat = np.concatenate([g[0] for g in self.groups]).astype(np.int64)
        self.starts = np.cumsum([0] + [g[0].size for g in self.groups]).astype(np.int64)
        self.n1s = np.array([g[1] for g in self.groups], dtype=np.float64)

    def cols(self, SC, Y):
        """SC, Y: (n x B). Returns (B,) pooled within-stratum AUROC."""
        SCT = np.ascontiguousarray(SC.T); YT = np.ascontiguousarray(np.broadcast_to(Y, SC.shape).T, dtype=np.float64)
        return _strat_nb(SCT, YT, self.starts, self.idxcat, self.n1s) / self.W


# ------------------------------------------------------------------ vectorised gene bootstrap
def boot_weights(n, B, rng):
    """Multinomial gene-bootstrap weights: row r = bincount of n draws with replacement (identical in law,
    and draw-for-draw, to indexing with rng.integers(0, n, n))."""
    W = np.empty((B, n))
    for r in range(B):
        W[r] = np.bincount(rng.integers(0, n, n), minlength=n)
    return W


def _weighted_mw(score, is_pos, W):
    """Weighted Mann-Whitney numerator and denominator for every bootstrap row of W (B x n); ties count 1/2.
    Equals the AUROC of the resampled data (duplicates as repeated observations, average ranks)."""
    o = np.argsort(score, kind="mergesort"); s = score[o]; p = is_pos[o]; Wo = W[:, o]
    _, start = np.unique(s, return_index=True)
    Wn = np.where(~p[None, :], Wo, 0.0); Wp = np.where(p[None, :], Wo, 0.0)
    gneg = np.add.reduceat(Wn, start, axis=1); gpos = np.add.reduceat(Wp, start, axis=1)
    below = np.cumsum(gneg, axis=1) - gneg
    num = (gpos * (below + 0.5 * gneg)).sum(1)
    den = gpos.sum(1) * gneg.sum(1)
    return num, den


def boot_strat_fast(score, is_pos, strata, W):
    """(B x 2): column 0 = stratified AUROC (stats3.stratified_auroc semantics: strata lacking a class in the
    resample are skipped), column 1 = unstratified AUROC, for bootstrap weight rows W."""
    score = np.asarray(score, float); is_pos = np.asarray(is_pos, bool); strata = np.asarray(strata)
    num = np.zeros(W.shape[0]); den = np.zeros(W.shape[0])
    for st in np.unique(strata):
        m = np.flatnonzero(strata == st)
        if is_pos[m].all() or (~is_pos[m]).all(): continue
        a, b = _weighted_mw(score[m], is_pos[m], W[:, m]); num += a; den += b
    a, b = _weighted_mw(score, is_pos, W)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.column_stack([num / den, a / b])
