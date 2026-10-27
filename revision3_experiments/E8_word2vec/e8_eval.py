"""E8 step 3 - evaluate the Word2Vec instrument on the three tasks and compare it, pair for pair, with the
scGPT residual stream (sensitivity analysis; CONVENTION3 §5 'E8 Word2Vec'; no pass/fail).

Gene set: V (cache/<ctx>/obs_genes.json) intersected with every evaluated model's vocabulary (asserted equal
to V for the models produced by e8_train.py). Each representation (Word2Vec model, scGPT layer) is
mean-centred over that gene set; SVD = stats3.svd_observed on the centred matrix.

Tasks (per representation):
  H1a  scalar assortativity (stats3.assortativity) of SV1..SV10 over STRING-physical>=700 edges within V;
       Maslov-Sneppen degree-preserving null (stats3.maslov_sneppen_chain, burn-in 10m, thinning m),
       N = 1,000 joint draws shared by every representation; Westfall-Young max-T (stats3.max_t_adjust,
       tail 'greater' = assortative) over the 10 SVs of each Word2Vec model, and over SV1-10 x 12 layers
       (120) for scGPT (sensitivity at N = 1,000; the primary H1 uses N = 10,000).
  H1b  basis-free pair-cosine AUROC: cosine of the centred full-dimensional vectors on STRING edges vs one
       Maslov-Sneppen rewired non-pair set (original edges removed); interval and paired delta by vertex
       subsampling (stats3.vertex_subsample_ci, m = n/2, B = 2,000).
  H2   detection-stratified LOO TF AUROC (stats3.tf_scores_loo + stats3.stratified_auroc; Lambert TF vs every
       other gene in V; strata = det_decile); within-stratum label permutation p (stats3.permute_within_strata,
       N = 2,000, joint across representations); paired delta by gene bootstrap (B = 2,000, scores fixed).
  H3   within-TF stratified AUROC of signed cosine on h3_pairs (stats3.StratifiedPairs); two-way TF x target
       cluster bootstrap (B = 2,000) for the interval and the paired delta (identical weights for every
       representation: each bootstrap call uses the same seed).
Word2Vec results: per seed, and mean / SD over the 5 seeds. Paired delta = mean over seeds of the Word2Vec
statistic minus the scGPT statistic at layer l, computed inside each bootstrap replicate.

Writes eval/<ctx>.json and eval/<ctx>_w2v.csv, eval/<ctx>_paired.csv, eval/<ctx>_h1_scgpt.csv (checkpoint per
context). Usage: python e8_eval.py <ctx> [--models DIR]
"""
from __future__ import annotations
import argparse, hashlib, json, re, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import stats3 as S  # noqa: E402

LOCAL_EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
N_MS, N_PERM, B = 1000, 2000, 2000
N_SV, N_LAYERS = 10, 12
SEED = 20261026
ALPHA = 0.05


def convention_hash():
    txt = (R3 / "FROZEN.txt").read_text().split()
    frozen = txt[txt.index("sha256") + 1]
    actual = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
    return frozen, actual


def null_summary(x):
    x = np.asarray(x, float)
    q = np.quantile(x, [0.025, 0.5, 0.95, 0.975, 0.99])
    return dict(N=int(x.size), mean=float(x.mean()), sd=float(x.std(ddof=1)), q025=float(q[0]), q50=float(q[1]),
                q95=float(q[2]), q975=float(q[3]), q99=float(q[4]))


def centred(M):
    M = np.asarray(M, np.float64)
    return M - M.mean(0, keepdims=True)


def pair_cos(E, a, b):
    En = S.l2norm_rows(E)
    return (En[a] * En[b]).sum(1)


def load_models(ctx, mdir):
    reps = {}
    for p in sorted(mdir.glob(f"{ctx}__*.npy")):
        m = re.match(rf"{ctx}__(\w+?)__w(\d+)__e(\d+)__s(\d+)$", p.stem)
        if not m: continue
        corpus, w, e, s = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
        reps[(corpus, w, e, s)] = np.load(p)
    return reps


def group_key(k):
    corpus, w, e, s = k
    return f"{corpus}_w{w}_e{e}"


class IndexedStratifiedPairs(S.StratifiedPairs):
    """stats3.StratifiedPairs with the per-TF boolean masks replaced by precomputed index arrays (speed only).
    auroc() is the stats3 formula line for line; bootstrap() and the constructor are inherited unchanged.
    Equality with stats3.StratifiedPairs.auroc is asserted on the observed scores and on bootstrap weights."""
    def __init__(self, tf, tg, y):
        super().__init__(tf, tg, y)
        self._pi = [np.flatnonzero((self.tf_code == c) & self.y) for c in range(self.tfs.size)]
        self._ni = [np.flatnonzero((self.tf_code == c) & ~self.y) for c in range(self.tfs.size)]

    def auroc(self, s, w=None) -> float:
        s = np.asarray(s, float); w = np.ones(s.size) if w is None else np.asarray(w, float)
        num = den = 0.0
        for c in range(self.tfs.size):
            ip, iq = self._pi[c], self._ni[c]
            sp, wp = s[ip], w[ip]; sn, wn = s[iq], w[iq]
            if wp.sum() == 0 or wn.sum() == 0: continue
            order = np.argsort(sn, kind="mergesort"); sn, wn = sn[order], wn[order]
            cw = np.concatenate([[0.0], np.cumsum(wn)])
            lo = np.searchsorted(sn, sp, side="left"); hi = np.searchsorted(sn, sp, side="right")
            num += float((wp * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))).sum())
            den += float(wp.sum() * wn.sum())
        return num / den if den > 0 else float("nan")


def h2_null_scores(En, M):
    """stats3.tf_scores_loo for b label sets at once (columns of the n x b indicator M): identical formula,
    En @ (En.T @ M) replaces the per-set En @ v (speed only; equality asserted in main)."""
    k = M[:, 0].sum(); selfsim = (En * En).sum(1)[:, None]
    tot = En @ (En.T @ M.astype(np.float64))
    return np.where(M, (tot - selfsim) / max(k - 1, 1), tot / k)


def strat_auroc_cols(sc, M, strata):
    """stats3.stratified_auroc for every column (positives = M[:, j]); within-stratum positive counts are equal
    across columns (within-stratum permutation). Continuous scores: no tie correction (equality asserted)."""
    num = np.zeros(sc.shape[1]); den = 0.0
    for st in np.unique(strata):
        m = strata == st; s = sc[m]; P = M[m]
        n1 = int(P[:, 0].sum()); n0 = int(m.sum()) - n1
        if n1 == 0 or n0 == 0: continue
        ranks = s.argsort(0).argsort(0) + 1.0
        a = ((ranks * P).sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)
        num += a * n1 * n0; den += n1 * n0
    return num / den


def replay_ci(values_per_call, n_nodes, pi, pj, seed, nb):
    """stats3.vertex_subsample_ci on pre-computed per-draw values (identical node draws: same seed)."""
    it = iter(values_per_call)
    return S.vertex_subsample_ci(lambda nodes, a, b: next(it), n_nodes, pi, pj, nb, np.random.default_rng(seed))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("ctx"); ap.add_argument("--models", default=str(HERE / "models"))
    ap.add_argument("--out", default=str(HERE / "eval")); ap.add_argument("--n-ms", type=int, default=N_MS)
    ap.add_argument("--n-perm", type=int, default=N_PERM); ap.add_argument("--B", type=int, default=B)
    a = ap.parse_args(); ctx = a.ctx; out = Path(a.out); out.mkdir(exist_ok=True)
    n_ms, n_perm, nb = a.n_ms, a.n_perm, a.B
    t0 = time.time()
    log = lambda *x: print(time.strftime("%H:%M:%S"), ctx, *x, flush=True)
    d = R3 / "cache" / ctx
    genes = json.loads((d / "obs_genes.json").read_text()); obs_idx = np.load(d / "obs_idx.npy")
    nl = np.load(d / "node_labels.npz"); e_st = np.load(d / "edges_string.npy"); h3 = np.load(d / "h3_pairs.npz")
    w2v = load_models(ctx, Path(a.models))
    assert w2v, "no models"
    inv = np.ones(len(genes), bool)
    for W in w2v.values(): inv &= np.isfinite(W).all(1)
    n_missing = int((~inv).sum())
    assert n_missing == 0, f"{n_missing} V genes missing from a Word2Vec vocabulary"   # expected 0 (step-1 check)
    n = len(genes)
    emb = np.load(LOCAL_EMB / f"{ctx}.npy", mmap_mode="r")
    reps = {}
    for k in sorted(w2v): reps[("w2v",) + k] = centred(w2v[k])
    for l in range(N_LAYERS): reps[("scgpt", l)] = centred(np.asarray(emb[l], np.float64)[obs_idx])
    keys = list(reps)
    log("representations", len(keys), "genes", n)

    # ------------------------------------------------------------------ SVD coordinates
    svc = {}
    for k in keys:
        if k[0] == "w2v":
            US, Sv, _ = S.svd_observed(reps[k][None], 0, np.arange(n))
        else:
            US, Sv, _ = S.svd_observed(emb, k[1], obs_idx)
            assert np.allclose(US @ _ , reps[k], atol=1e-6 * np.abs(reps[k]).max())
        svc[k] = US[:, :N_SV]
    Xall = np.concatenate([svc[k] for k in keys], 1)          # n x (R*10)

    # ------------------------------------------------------------------ H1a assortativity + MS null
    pi, pj = e_st[:, 0].astype(np.int64), e_st[:, 1].astype(np.int64)
    r_obs = np.array([S.assortativity(Xall[:, c], pi, pj) for c in range(Xall.shape[1])])
    # Under degree-preserving rewiring the stub mean and stub variance of each attribute are invariant, so
    # r = (mean over edges of x_i x_j - mu^2) / var_stub; verified against stats3.assortativity on the first draw.
    stubs = np.concatenate([pi, pj]); mu = Xall[stubs].mean(0); var = Xall[stubs].var(0)
    r_null = np.empty((n_ms, Xall.shape[1])); acc = []
    for t, (ra, rb, ar) in enumerate(S.maslov_sneppen_chain(pi, pj, n, n_ms, seed=SEED + 1, burn_factor=10, thin_factor=1)):
        r_null[t] = ((Xall[ra] * Xall[rb]).mean(0) - mu ** 2) / var; acc.append(ar)
        if t == 0:
            chk = np.array([S.assortativity(Xall[:, c], ra, rb) for c in range(0, Xall.shape[1], 7)])
            assert np.allclose(chk, r_null[0, ::7], atol=1e-9), "vectorised assortativity mismatch"
    r_chk = ((Xall[pi] * Xall[pj]).mean(0) - mu ** 2) / var
    assert np.allclose(r_chk, r_obs, atol=1e-9)
    S.assert_resolution(n_ms, ALPHA, 1); S.assert_resolution(n_ms, ALPHA, N_SV)
    log("H1a MS null done", round(time.time() - t0))
    rows, h1s = [], []
    col = {k: slice(i * N_SV, (i + 1) * N_SV) for i, k in enumerate(keys)}
    for k in keys:
        if k[0] != "w2v": continue
        z, pa, pr = S.max_t_adjust(r_obs[col[k]], r_null[:, col[k]], tail="greater")
        for j in range(N_SV):
            ns = null_summary(r_null[:, col[k]][:, j])
            rows.append(dict(context=ctx, corpus=k[1], window=k[2], epochs=k[3], seed=k[4], task="H1a_assortativity",
                             sv=j + 1, value=r_obs[col[k]][j], z=z[j], p_raw=pr[j], p_maxT=pa[j], family=N_SV,
                             null_N=n_ms, null_mean=ns["mean"], null_sd=ns["sd"], null_q975=ns["q975"],
                             n_pos=int(pi.size)))
    sc_cols = np.concatenate([np.arange(N_SV) + i * N_SV for i, k in enumerate(keys) if k[0] == "scgpt"])
    z, pa, pr = S.max_t_adjust(r_obs[sc_cols], r_null[:, sc_cols], tail="greater")
    for q, c in enumerate(sc_cols):
        l = keys[c // N_SV][1]
        h1s.append(dict(context=ctx, layer=l, sv=c % N_SV + 1, r=r_obs[c], z=z[q], p_raw=pr[q], p_maxT_120=pa[q],
                        null_N=n_ms, null_mean=r_null[:, c].mean(), null_sd=r_null[:, c].std(ddof=1)))
    pd.DataFrame(h1s).to_csv(out / f"{ctx}_h1_scgpt.csv", index=False)

    # ------------------------------------------------------------------ H1b pair-cosine AUROC vs rewired non-pairs
    ra, rb, _ = next(S.maslov_sneppen_chain(pi, pj, n, 1, seed=SEED + 2, burn_factor=10, thin_factor=1))
    ek = set(map(tuple, np.sort(e_st, 1).tolist()))
    keep = np.array([(min(x, y), max(x, y)) not in ek for x, y in zip(ra, rb)])
    ni, nj = ra[keep], rb[keep]
    np.save(out / f"{ctx}_rewired_nonpairs.npy", np.stack([ni, nj], 1))
    cos_p = {k: pair_cos(reps[k], pi, pj) for k in keys}; cos_n = {k: pair_cos(reps[k], ni, nj) for k in keys}
    groups = sorted({group_key(k[1:]) for k in keys if k[0] == "w2v"})
    gmem = {g: [k for k in keys if k[0] == "w2v" and group_key(k[1:]) == g] for g in groups}
    sc_keys = [k for k in keys if k[0] == "scgpt"]
    nidx = {k: i for i, k in enumerate(keys)}
    store = []

    def h1b_stats(nodes, a_, b_):
        m = np.zeros(n, bool); m[nodes] = True
        ep = m[pi] & m[pj]; en = m[ni] & m[nj]
        v = np.array([S.auroc(cos_p[k][ep], cos_n[k][en]) for k in keys]); store.append(v); return v[0]
    S.vertex_subsample_ci(h1b_stats, n, pi, pj, nb, np.random.default_rng(SEED + 3))
    store = np.array(store)                         # (1 + B) x R ; row 0 = full graph
    log("H1b vertex subsampling done", round(time.time() - t0))

    def lvl_ci(vals):
        est, lo, hi, _ = replay_ci(list(vals), n, pi, pj, SEED + 3, nb); return est, lo, hi
    h1b = {}
    for k in keys:
        est, lo, hi = lvl_ci(store[:, nidx[k]]); h1b[k] = (est, lo, hi)
        if k[0] == "w2v":
            rows.append(dict(context=ctx, corpus=k[1], window=k[2], epochs=k[3], seed=k[4], task="H1b_pair_cos_auroc",
                             value=est, ci_lo=lo, ci_hi=hi, ci_method="vertex_subsample_m=n/2", B=nb,
                             n_pos=int(pi.size), n_neg=int(ni.size)))
    paired = []
    for g in groups:
        gm = np.array([nidx[k] for k in gmem[g]])
        wm = store[:, gm].mean(1)
        est, lo, hi = lvl_ci(wm)
        rows.append(dict(context=ctx, corpus=g, seed="mean", task="H1b_pair_cos_auroc", value=est, ci_lo=lo, ci_hi=hi,
                         sd_seeds=float(store[0, gm].std(ddof=1)) if gm.size > 1 else np.nan, n_seeds=int(gm.size),
                         ci_method="vertex_subsample_m=n/2", B=nb, n_pos=int(pi.size), n_neg=int(ni.size)))
        for k in sc_keys:
            dv = wm - store[:, nidx[k]]
            de, dlo, dhi = lvl_ci(dv)
            paired.append(dict(context=ctx, w2v=g, task="H1b_pair_cos_auroc", layer=k[1], w2v_mean=wm[0],
                               w2v_sd_seeds=float(store[0, gm].std(ddof=1)) if gm.size > 1 else np.nan,
                               w2v_seed_min=float(store[0, gm].min()), w2v_seed_max=float(store[0, gm].max()),
                               scgpt=store[0, nidx[k]], scgpt_ci_lo=h1b[k][1], scgpt_ci_hi=h1b[k][2],
                               delta=de, delta_ci_lo=dlo, delta_ci_hi=dhi, ci_method="vertex_subsample_m=n/2",
                               B=nb, n_genes=n, n_pos=int(pi.size), n_neg=int(ni.size)))

    # ------------------------------------------------------------------ H2 detection-stratified LOO TF AUROC
    is_tf = nl["is_tf_lambert"].astype(bool); strata = nl["det_decile"]
    tf_pos = np.flatnonzero(is_tf)
    sc2 = {k: S.tf_scores_loo(reps[k], tf_pos) for k in keys}
    h2 = {k: S.stratified_auroc(sc2[k], is_tf, strata) for k in keys}
    En = {k: S.l2norm_rows(reps[k]) for k in keys}
    rng = np.random.default_rng(SEED + 4)
    Yp = np.stack([S.permute_within_strata(is_tf, strata, rng) for _ in range(n_perm)], 1)    # n x N
    h2_null = np.empty((n_perm, len(keys)))
    for i, k in enumerate(keys):
        for b0 in range(0, n_perm, 250):
            Mb = Yp[:, b0:b0 + 250]
            h2_null[b0:b0 + Mb.shape[1], i] = strat_auroc_cols(h2_null_scores(En[k], Mb), Mb, strata)
        if i == 0:      # the batched computation equals stats3 on the first permutations
            for t in range(3):
                tp = np.flatnonzero(Yp[:, t])
                ref = S.stratified_auroc(S.tf_scores_loo(reps[k], tp), Yp[:, t], strata)
                assert abs(ref - h2_null[t, i]) < 1e-9, (ref, h2_null[t, i])
            sc0 = h2_null_scores(En[k], is_tf[:, None])[:, 0]
            assert np.allclose(sc0, sc2[k], atol=1e-10) and abs(strat_auroc_cols(sc0[:, None], is_tf[:, None], strata)[0] - h2[k]) < 1e-9
    S.assert_resolution(n_perm, ALPHA, 1)
    log("H2 permutation null done", round(time.time() - t0))
    # gene bootstrap (fixed scores), identical draws for every representation
    rng = np.random.default_rng(SEED + 5); h2_boot = np.empty((nb, len(keys)))
    for r in range(nb):
        ix = rng.integers(0, n, n)
        for i, k in enumerate(keys):
            h2_boot[r, i] = S.stratified_auroc(sc2[k][ix], is_tf[ix], strata[ix])
    log("H2 gene bootstrap done", round(time.time() - t0))
    sc_i = np.array([nidx[k] for k in sc_keys])
    zs, pa_sc, pr_sc = S.max_t_adjust(np.array([h2[k] for k in sc_keys]), h2_null[:, sc_i], tail="greater")
    h2_scgpt_p = {k: (pr_sc[q], pa_sc[q]) for q, k in enumerate(sc_keys)}
    for k in keys:
        if k[0] != "w2v": continue
        i = nidx[k]; p = float(S.emp_p(h2[k], h2_null[:, i], "greater")[0])
        ns = null_summary(h2_null[:, i])
        rows.append(dict(context=ctx, corpus=k[1], window=k[2], epochs=k[3], seed=k[4], task="H2_strat_tf_auroc",
                         value=h2[k], ci_lo=np.quantile(h2_boot[:, i], 0.025), ci_hi=np.quantile(h2_boot[:, i], 0.975),
                         ci_method="gene_bootstrap", B=nb, p_raw=p, null_N=n_perm, null_mean=ns["mean"],
                         null_sd=ns["sd"], null_q975=ns["q975"], n_pos=int(is_tf.sum()), n_neg=int((~is_tf).sum())))

    # ------------------------------------------------------------------ H3 within-TF stratified AUROC
    sp = IndexedStratifiedPairs(h3["tf"], h3["tg"], h3["y"]); sp_ref = S.StratifiedPairs(h3["tf"], h3["tg"], h3["y"])
    cos3 = {k: pair_cos(reps[k], sp.tf, sp.tg) for k in keys}
    h3v = {k: sp.auroc(cos3[k]) for k in keys}
    k0 = keys[0]; assert abs(h3v[k0] - sp_ref.auroc(cos3[k0])) < 1e-12
    b_fast = sp.bootstrap(cos3[k0], 3, np.random.default_rng(7)); b_ref = sp_ref.bootstrap(cos3[k0], 3, np.random.default_rng(7))
    assert np.allclose(b_fast, b_ref, atol=1e-12), (b_fast, b_ref)
    boots = Parallel(n_jobs=3)(delayed(sp.bootstrap)(cos3[k], nb, np.random.default_rng(SEED + 6)) for k in keys)
    h3_boot = np.stack(boots, 1)
    log("H3 two-way bootstrap done", round(time.time() - t0))
    for k in keys:
        if k[0] != "w2v": continue
        i = nidx[k]
        rows.append(dict(context=ctx, corpus=k[1], window=k[2], epochs=k[3], seed=k[4], task="H3_withinTF_auroc",
                         value=h3v[k], ci_lo=np.nanquantile(h3_boot[:, i], 0.025),
                         ci_hi=np.nanquantile(h3_boot[:, i], 0.975), upper95_1sided=S.one_sided_upper(h3_boot[:, i]),
                         ci_method="two_way_cluster_bootstrap", B=nb, n_pos=int(sp.y.sum()), n_neg=int((~sp.y).sum()),
                         n_tfs=int(sp.tfs.size)))

    # ------------------------------------------------------------------ seed means and paired deltas for H2, H3
    for task, obs, boot in (("H2_strat_tf_auroc", h2, h2_boot), ("H3_withinTF_auroc", h3v, h3_boot)):
        meth = "gene_bootstrap" if task.startswith("H2") else "two_way_cluster_bootstrap"
        for g in groups:
            gm = np.array([nidx[k] for k in gmem[g]]); vals = np.array([obs[k] for k in gmem[g]])
            bm = boot[:, gm].mean(1)
            extra = {}
            if task.startswith("H2"):
                nm = h2_null[:, gm].mean(1)
                extra = dict(p_raw=float(S.emp_p(vals.mean(), nm, "greater")[0]), null_N=n_perm,
                             null_mean=float(nm.mean()), null_sd=float(nm.std(ddof=1)))
            rows.append(dict(context=ctx, corpus=g, seed="mean", task=task, value=vals.mean(),
                             sd_seeds=float(vals.std(ddof=1)) if vals.size > 1 else np.nan, n_seeds=int(vals.size),
                             ci_lo=np.nanquantile(bm, 0.025), ci_hi=np.nanquantile(bm, 0.975), ci_method=meth, B=nb,
                             upper95_1sided=S.one_sided_upper(bm) if task.startswith("H3") else np.nan, **extra))
            for k in sc_keys:
                i = nidx[k]; dv = bm - boot[:, i]
                rec = dict(context=ctx, w2v=g, task=task, layer=k[1], w2v_mean=vals.mean(),
                           w2v_sd_seeds=float(vals.std(ddof=1)) if vals.size > 1 else np.nan,
                           w2v_seed_min=float(vals.min()), w2v_seed_max=float(vals.max()), scgpt=obs[k],
                           scgpt_ci_lo=np.nanquantile(boot[:, i], 0.025), scgpt_ci_hi=np.nanquantile(boot[:, i], 0.975),
                           delta=vals.mean() - obs[k], delta_ci_lo=np.nanquantile(dv, 0.025),
                           delta_ci_hi=np.nanquantile(dv, 0.975), ci_method=meth, B=nb, n_genes=n)
                if task.startswith("H2"):
                    # the LOO estimator is not centred at 0.5 under the null (a member's score omits one
                    # same-stratum neighbour); the null-centred delta removes each representation's own null mean
                    off_ = h2_null[:, gm].mean(1).mean() - h2_null[:, i].mean()
                    rec.update(scgpt_p_raw=h2_scgpt_p[k][0], scgpt_p_maxT12=h2_scgpt_p[k][1],
                               w2v_null_mean=float(h2_null[:, gm].mean(1).mean()), scgpt_null_mean=float(h2_null[:, i].mean()),
                               delta_nullcentred=rec["delta"] - off_, delta_nullcentred_ci_lo=rec["delta_ci_lo"] - off_,
                               delta_nullcentred_ci_hi=rec["delta_ci_hi"] - off_)
                else:
                    rec.update(n_pos=int(sp.y.sum()), n_neg=int((~sp.y).sum()), n_tfs=int(sp.tfs.size))
                paired.append(rec)
    # seed-mean H1a summary rows (per SV: mean / SD of r, and the fraction of seeds with max-T p < 0.05)
    df = pd.DataFrame(rows)
    h1a = df[df.task == "H1a_assortativity"].copy(); h1a["grp"] = [group_key((c, int(w), int(e), 0)) for c, w, e in
                                                                   zip(h1a.corpus, h1a.window, h1a.epochs)]
    for (g, sv), sub in h1a.groupby(["grp", "sv"]):
        rows.append(dict(context=ctx, corpus=g, seed="mean", task="H1a_assortativity", sv=int(sv), value=sub.value.mean(),
                         sd_seeds=sub.value.std(ddof=1) if len(sub) > 1 else np.nan, n_seeds=len(sub),
                         z=sub.z.mean(), p_maxT=sub.p_maxT.max(), frac_seeds_maxT_lt_005=float((sub.p_maxT < 0.05).mean()),
                         family=N_SV, null_N=n_ms))
    # scGPT levels for reference (same estimators, same genes / pairs)
    for k in sc_keys:
        i = nidx[k]
        paired.append(dict(context=ctx, w2v="(scgpt level)", task="H1a_assortativity_maxSV", layer=k[1],
                           scgpt=float(r_obs[col[k]].max()), scgpt_sv=int(r_obs[col[k]].argmax() + 1),
                           scgpt_p_maxT120=float(min(h1["p_maxT_120"] for h1 in h1s if h1["layer"] == k[1]))))
    for g in groups:
        sub = [r for r in rows if r.get("task") == "H1a_assortativity" and r.get("seed") != "mean"
               and group_key((r["corpus"], r["window"], r["epochs"], 0)) == g]
        best = {}
        for r in sub: best.setdefault(r["seed"], []).append(r["value"])
        mx = np.array([max(v) for v in best.values()])
        for k in sc_keys:
            paired.append(dict(context=ctx, w2v=g, task="H1a_assortativity_maxSV", layer=k[1], w2v_mean=mx.mean(),
                               w2v_sd_seeds=mx.std(ddof=1) if mx.size > 1 else np.nan, scgpt=float(r_obs[col[k]].max()),
                               delta=mx.mean() - float(r_obs[col[k]].max()), ci_method="none (descriptive; SV index not matched)"))
    pd.DataFrame(rows).to_csv(out / f"{ctx}_w2v.csv", index=False)
    pd.DataFrame(paired).to_csv(out / f"{ctx}_paired.csv", index=False)
    frozen, actual = convention_hash()
    res = dict(context=ctx, convention3_sha256_frozen=frozen, convention3_sha256_actual=actual,
               n_genes=n, n_missing_from_vocab=n_missing, n_string_edges=int(pi.size), n_rewired_nonpairs=int(ni.size),
               n_lambert_tf=int(is_tf.sum()), h3_pos=int(sp.y.sum()), h3_neg=int((~sp.y).sum()), h3_tfs=int(sp.tfs.size),
               models=[group_key(k[1:]) + f"_s{k[4]}" for k in keys if k[0] == "w2v"],
               ms_null=dict(N=n_ms, burn_factor=10, thin_factor=1, seed=SEED + 1, mean_accept_rate=float(np.mean(acc))),
               h2_perm=dict(N=n_perm, seed=SEED + 4), B=nb,
               null_summaries=dict(
                   h1a={f"{group_key(k[1:])}_s{k[4]}" if k[0] == "w2v" else f"scgpt_L{k[1]}":
                        {f"SV{j+1}": null_summary(r_null[:, col[k]][:, j]) for j in range(N_SV)} for k in keys},
                   h2={f"{group_key(k[1:])}_s{k[4]}" if k[0] == "w2v" else f"scgpt_L{k[1]}": null_summary(h2_null[:, nidx[k]])
                       for k in keys}),
               seconds=round(time.time() - t0, 1))
    (out / f"{ctx}.json").write_text(json.dumps(res, indent=1, default=float))
    log("DONE", round(time.time() - t0))


if __name__ == "__main__":
    main()
