"""Shared helpers for E10 (cross-context variance decomposition and diagnosis; SENSITIVITY analysis).

Read-only use of contexts.py, stats3.py and cache/. Point estimates always come from stats3 functions
(tf_scores_loo, stratified_auroc, auroc, StratifiedPairs.auroc, svd_observed, maslov_sneppen_chain).
The vectorised *weighted* AUROCs below exist only to evaluate many bootstrap weight vectors at once; at
unit weights they reproduce the stats3 values exactly (checked by selftest()).
"""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import contexts as X   # noqa: E402
import stats3 as S     # noqa: E402

LOCAL_EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
CTXS = list(X.CONTEXTS.keys())          # all 16
N_LAYERS = 12
SEED = 20261026
BAND = slice(1, 7)                       # SV2..SV7 of the observed-gene SVD (0-based columns 1..6)
LAMBERT = X.REPO / "subproject_61_two_axes_revision/data/humantfs_v1.01.txt"


def convention_hash():
    txt = (R3 / "FROZEN.txt").read_text().split()
    frozen = txt[txt.index("sha256") + 1]
    actual = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
    return dict(frozen=frozen, actual=actual, match=frozen == actual)


def load_emb(ctx):
    p = LOCAL_EMB / f"{ctx}.npy"
    return np.load(p if p.exists() else X.CONTEXTS[ctx]["emb"], mmap_mode="r")


def load_cache(ctx):
    d = R3 / "cache" / ctx
    genes = json.loads((d / "obs_genes.json").read_text())
    obs_idx = np.load(d / "obs_idx.npy")
    nl = dict(np.load(d / "node_labels.npz"))
    nobs = np.load(d / "ncells.npy").astype(float)
    return dict(genes=genes, obs_idx=obs_idx, nl=nl, nobs=nobs, d=d)


def centred(emb, layer, obs_idx):
    """V x 512 layer matrix, mean-centred over observed genes V (CONVENTION3 §2)."""
    Xm = np.asarray(emb[layer], dtype=np.float64)[obs_idx]
    return Xm - Xm.mean(0, keepdims=True)


def unit(E):
    return S.l2norm_rows(E)


def null_summary(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    q = np.quantile(x, [0.025, 0.05, 0.5, 0.95, 0.975])
    return dict(N=int(x.size), mean=float(x.mean()), sd=float(x.std(ddof=1)), q025=float(q[0]),
                q05=float(q[1]), q50=float(q[2]), q95=float(q[3]), q975=float(q[4]))


# ------------------------------------------------------------------------------------------------
# vectorised weighted AUROCs (bootstrap replicates as columns)
# ------------------------------------------------------------------------------------------------
def wauroc_fixed(s, is_pos, W):
    """Weighted Mann-Whitney AUROC for FIXED scores s (n,) and weight matrix W (n, B); ties count 1/2.
    Returns (B,) AUROCs. Pair weight = w_p * w_q."""
    s = np.asarray(s, float); is_pos = np.asarray(is_pos, bool)
    sp, sn = s[is_pos], s[~is_pos]; Wp, Wn = W[is_pos], W[~is_pos]
    o = np.argsort(sn, kind="mergesort"); sn = sn[o]; Wn = Wn[o]
    cw = np.vstack([np.zeros((1, W.shape[1])), np.cumsum(Wn, 0)])          # (nneg+1, B)
    lo = np.searchsorted(sn, sp, "left"); hi = np.searchsorted(sn, sp, "right")
    num = (Wp * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))).sum(0)
    den = Wp.sum(0) * Wn.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def wauroc_var(Sc, is_pos, W):
    """Weighted AUROC when scores differ per column: Sc (n, B), W (n, B). Ties ignored (continuous scores)."""
    is_pos = np.asarray(is_pos, bool)
    o = np.argsort(Sc, axis=0, kind="mergesort")
    Ws = np.take_along_axis(W, o, 0); P = is_pos[o]
    wn = np.where(P, 0.0, Ws); wp = np.where(P, Ws, 0.0)
    cneg_before = np.cumsum(wn, 0) - wn
    num = (wp * cneg_before).sum(0); den = wp.sum(0) * wn.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def wstrat_auroc_var(Sc, is_pos, strata, W):
    """Stratified weighted AUROC (pooled within strata, weights = W_pos * W_neg), columns = replicates."""
    num = np.zeros(W.shape[1]); den = np.zeros(W.shape[1])
    for st in np.unique(strata):
        m = strata == st
        p = is_pos[m]
        if p.sum() == 0 or (~p).sum() == 0: continue
        a = wauroc_var(Sc[m], p, W[m])
        d = W[m][p].sum(0) * W[m][~p].sum(0)
        a = np.nan_to_num(a); num += a * d; den += d
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def loo_scores_weighted(En, tf_pos, W):
    """Bootstrap version of stats3.tf_scores_loo: every gene's mean cosine to the (weighted) TF set, leaving
    the GENE out: all w_i copies of a resampled TF are removed from its own reference set (removing only one
    copy would leave duplicate self-similarities of 1 in the score and bias every bootstrap replicate upwards).
    At unit weights this is exactly stats3.tf_scores_loo. En (n, d) unit rows; W (n, B) multiplicities."""
    C = En @ En[tf_pos].T                                   # (n, k)
    Wt = W[tf_pos]                                          # (k, B)
    tot = C @ Wt                                            # (n, B)
    ktot = Wt.sum(0, keepdims=True)                         # (1, B)
    is_tf = np.zeros(En.shape[0], bool); is_tf[tf_pos] = True
    selfc = (En * En).sum(1)[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        Wself = np.where(is_tf[:, None], W, 0.0)
        out = np.where(is_tf[:, None], (tot - Wself * selfc) / np.maximum(ktot - Wself, 1), tot / np.maximum(ktot, 1))
    return out


def w_strat_pairs(sp: "S.StratifiedPairs", s, Wtf, Wtg):
    """Two-way (TF x target) bootstrap of the within-TF stratified AUROC for fixed scores s.
    Wtf (n_tf_codes, B), Wtg (n_tg_codes, B): multiplicities. Pair weight = wtf(t) * wtg(g), as in
    stats3.StratifiedPairs.bootstrap(two_way=True)."""
    s = np.asarray(s, float); B = Wtf.shape[1]
    num = np.zeros(B); den = np.zeros(B)
    for c in range(sp.tfs.size):
        m = sp.tf_code == c; mp = m & sp.y; mn = m & ~sp.y
        if not mp.any() or not mn.any(): continue
        spv, snv = s[mp], s[mn]
        M = (spv[:, None] > snv[None, :]) + 0.5 * (spv[:, None] == snv[None, :])     # (P, N)
        wp = Wtg[sp.tg_code[mp]]; wn = Wtg[sp.tg_code[mn]]                          # (P,B), (N,B)
        t2 = Wtf[c] ** 2
        num += t2 * np.einsum("pb,pn,nb->b", wp, M, wn)
        den += t2 * wp.sum(0) * wn.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def h3_obj(ctx_cache):
    z = np.load(ctx_cache["d"] / "h3_pairs.npz")
    return S.StratifiedPairs(z["tf"], z["tg"], z["y"])


def rewired_nonpairs(edges, n_nodes, seed):
    """One Maslov-Sneppen degree-preserving rewiring (burn-in 10m accepted swaps, then m more) of the
    V-induced STRING graph; rewired pairs that coincide with a real edge are dropped (a NON-pair set)."""
    a = np.minimum(edges[:, 0], edges[:, 1]); b = np.maximum(edges[:, 0], edges[:, 1])
    key = np.unique(a.astype(np.int64) * n_nodes + b)
    a, b = key // n_nodes, key % n_nodes
    gen = S.maslov_sneppen_chain(a, b, n_nodes, 1, seed=seed, burn_factor=10, thin_factor=1)
    ra, rb, acc = next(gen)
    rk = ra * n_nodes + rb
    keep = ~np.isin(rk, key)
    return np.c_[a, b], np.c_[ra[keep], rb[keep]], float(acc)


def lambert_all():
    import pandas as pd
    lam = pd.read_csv(LAMBERT, sep="\t")
    lam = lam[lam["Is TF?"].astype(str).str.strip() == "Yes"]
    return lam[["Ensembl ID", "HGNC symbol"]].astype(str)


def selftest():
    rng = np.random.default_rng(0)
    # fixed-score weighted AUROC == stats3.auroc at unit weights (with ties) and with integer weights
    s = np.round(rng.normal(size=300), 1); y = rng.random(300) < 0.3
    assert abs(wauroc_fixed(s, y, np.ones((300, 1)))[0] - S.auroc(s[y], s[~y])) < 1e-12
    w = rng.poisson(1, 300)
    rep = np.repeat(np.arange(300), w)
    assert abs(wauroc_fixed(s, y, w[:, None].astype(float))[0] - S.auroc(s[rep][y[rep]], s[rep][~y[rep]])) < 1e-12
    # LOO weighted scores == stats3.tf_scores_loo at unit weights; stratified version == stats3.stratified_auroc
    E = unit(rng.normal(size=(400, 16))); tf = np.sort(rng.choice(400, 40, replace=False))
    sc = loo_scores_weighted(E, tf, np.ones((400, 1)))[:, 0]
    assert np.allclose(sc, S.tf_scores_loo(E, tf))
    w = rng.poisson(1, 400).astype(float); scw = loo_scores_weighted(E, tf, w[:, None])[:, 0]
    for i in (tf[0], tf[1], 5, 7):
        others = [j for j in tf if j != i]
        ref = sum(w[j] * E[i] @ E[j] for j in others) / max(sum(w[j] for j in others), 1)
        assert abs(scw[i] - ref) < 1e-10
    pos = np.zeros(400, bool); pos[tf] = True; st = rng.integers(0, 10, 400)
    a1 = wstrat_auroc_var(sc[:, None], pos, st, np.ones((400, 1)))[0]
    assert abs(a1 - S.stratified_auroc(sc, pos, st)) < 1e-12, (a1, S.stratified_auroc(sc, pos, st))
    # two-way pairs bootstrap == StratifiedPairs.auroc(s, w) for one replicate
    tfc = np.repeat(np.arange(8), 30); tg = rng.integers(0, 60, 240); yy = rng.random(240) < 0.25
    sp = S.StratifiedPairs(tfc, tg, yy); s2 = rng.normal(size=sp.tf.size)
    wtf = rng.poisson(1, sp.tfs.size).astype(float); wtg = rng.poisson(1, sp.tgs.size).astype(float)
    ref = sp.auroc(s2, wtf[sp.tf_code] * wtg[sp.tg_code])
    got = w_strat_pairs(sp, s2, wtf[:, None], wtg[:, None])[0]
    assert abs(ref - got) < 1e-12, (ref, got)
    assert abs(w_strat_pairs(sp, s2, np.ones((sp.tfs.size, 1)), np.ones((sp.tgs.size, 1)))[0] - sp.auroc(s2)) < 1e-12
    print("e10 selftest: ALL PASS")


if __name__ == "__main__":
    selftest()
