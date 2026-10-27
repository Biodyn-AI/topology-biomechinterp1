"""Independent adversarial re-computation of E7 headline numbers from raw inputs.

Does NOT import e7lib. Uses only: the canonical embedding files (contexts.CONTEXTS paths, not the local copies),
cache/<ctx>/ (obs_idx, obs_genes, node_labels, h3_pairs, edges_string), and stats3.maslov_sneppen_chain
(the pre-specified rewiring tool; seed 20261026, burn_factor 10). SVD, cosines, LOO TF score, stratified AUROCs,
vertex subsampling, two-way bootstrap and the joint permutation null are re-implemented here with sklearn /
scipy primitives and different code paths / seeds.

  python e7_verify.py select        # k* from immune_seed43 (all three tasks, full T curve)
  python e7_verify.py lung          # lung point estimates + PPI vertex CI (both scalings) + EDGE 2-way CI
  python e7_verify.py anyk_lungext  # lung_external NODE any-k joint permutation p (own seed)
  python e7_verify.py static        # one static model spot check
Outputs: verify/verify_<stage>.json
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent; E7 = HERE.parent; R3 = E7.parent
sys.path.insert(0, str(R3))
import contexts as X            # paths only
import stats3 as S              # maslov_sneppen_chain only

CACHE = R3 / "cache"; SEED = 20261026
K_T = [1, 2, 3, 5, 7, 10, 16, 32, 64, 128, 256, 512]
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)


def load(ctx):
    d = CACHE / ctx
    genes = json.loads((d / "obs_genes.json").read_text()); oi = np.load(d / "obs_idx.npy")
    emb = np.load(X.CONTEXTS[ctx]["emb"], mmap_mode="r")
    allg = X.h5ad_genes(X.CONTEXTS[ctx]["h5ad"])
    assert [allg[i] for i in oi] == genes, "gene alignment mismatch"
    nz = np.flatnonzero(np.linalg.norm(np.asarray(emb[0], np.float64), axis=1) > 0)
    nz11 = np.flatnonzero(np.linalg.norm(np.asarray(emb[11], np.float64), axis=1) > 0)
    assert np.array_equal(nz, oi) and np.array_equal(nz11, oi), "obs_idx != nonzero rows"
    nl = np.load(d / "node_labels.npz"); z = np.load(d / "h3_pairs.npz")
    E = np.load(d / "edges_string.npy")
    a = np.minimum(E[:, 0], E[:, 1]); b = np.maximum(E[:, 0], E[:, 1]); n = len(genes)
    key = np.unique(a.astype(np.int64) * n + b); a, b = key // n, key % n
    ra, rb, acc = next(S.maslov_sneppen_chain(a, b, n, 1, seed=SEED, burn_factor=10, thin_factor=1))
    deg0 = np.bincount(np.r_[a, b], minlength=n); deg1 = np.bincount(np.r_[ra, rb], minlength=n)
    assert np.array_equal(deg0, deg1)
    keep = ~np.isin(ra * n + rb, key)
    return dict(n=n, emb=emb, oi=oi, is_tf=nl["is_tf_lambert"].astype(bool), det=nl["det_decile"].astype(int),
                tf=z["tf"], tg=z["tg"], y=z["y"].astype(bool), pa=np.r_[a, ra[keep]], pb=np.r_[b, rb[keep]],
                py=np.r_[np.ones(a.size, bool), np.zeros(int(keep.sum()), bool)], n_edges=int(a.size),
                n_nonpairs=int(keep.sum()))


def us(D, L):
    Xm = np.asarray(D["emb"][L], np.float64)[D["oi"]]; Xm = Xm - Xm.mean(0)
    U, s, _ = np.linalg.svd(Xm, full_matrices=False)
    return U * s, Xm


def unit(M):
    return M / np.linalg.norm(M, axis=1, keepdims=True)


def ppi(En, D):
    s = (En[D["pa"]] * En[D["pb"]]).sum(1); return roc_auc_score(D["py"], s), s


def loo_scores(En, ytf):
    tot = En @ En[ytf].sum(0); k = ytf.sum()
    return np.where(ytf, (tot - 1.0) / (k - 1), tot / k)       # rows are unit norm: self-similarity = 1


def strat_auc(s, y, st):
    num = den = 0.0
    for v in np.unique(st):
        m = st == v
        if y[m].any() and (~y[m]).any():
            w = y[m].sum() * (~y[m]).sum(); num += roc_auc_score(y[m], s[m]) * w; den += w
    return num / den


def node(En, D):
    return strat_auc(loo_scores(En, D["is_tf"]), D["is_tf"], D["det"])


def edge(En, D):
    s = (En[D["tf"]] * En[D["tg"]]).sum(1); return strat_auc(s, D["y"], D["tf"]), s


def stage_select():
    D = load("immune_seed43"); out = {t: {k: [] for k in K_T} for t in ("PPI", "NODE", "EDGE")}
    for L in range(12):
        US, _ = us(D, L)
        for k in K_T:
            En = unit(US[:, :k])
            out["PPI"][k].append(ppi(En, D)[0]); out["NODE"][k].append(node(En, D)); out["EDGE"][k].append(edge(En, D)[0])
        log("select layer", L)
    res = {}
    for t in out:
        lm = {k: float(np.mean(v)) for k, v in out[t].items()}
        best = max(lm.values()); ks = min(k for k, v in lm.items() if v == best)
        res[t] = dict(k_star=ks, layer_mean_curve=lm)
    kj = json.loads((E7 / "kstar.json").read_text())
    res["compare_kstar_json"] = {t: dict(mine=res[t]["k_star"], e7=kj["tasks"][t]["k_star"],
                                         max_abs_curve_diff=max(abs(res[t]["layer_mean_curve"][k] -
                                                                    kj["tasks"][t]["layer_mean_curve"][str(k)]) for k in K_T))
                                 for t in ("PPI", "NODE", "EDGE")}
    (HERE / "verify_select.json").write_text(json.dumps(res, indent=1)); log(json.dumps(res["compare_kstar_json"]))


def wstrat_tf(s, y, tfc, W):
    """Within-TF weighted MW AUROC for continuous scores (no ties); W (B, n) pair weights."""
    o = np.lexsort((s, tfc)); y_, t_, W_ = y[o], tfc[o], W[:, o]
    seg = np.r_[0, np.flatnonzero(t_[1:] != t_[:-1]) + 1]
    Wn = W_ * (~y_); Wp = W_ * y_
    C = np.cumsum(Wn, 1); segstart = np.repeat(seg, np.diff(np.r_[seg, y_.size]))
    Cb = np.c_[np.zeros((W.shape[0], 1)), C][:, segstart]       # cum neg weight before segment start
    below = C - Wn - Cb                                          # neg weight strictly below within TF
    num = (Wp * below).sum(1)
    P = np.add.reduceat(Wp, seg, 1); N = np.add.reduceat(Wn, seg, 1)
    return num / (P * N).sum(1)


def stage_lung():
    kj = json.loads((E7 / "kstar.json").read_text()); kp, ke = kj["tasks"]["PPI"]["k_star"], kj["tasks"]["EDGE"]["k_star"]
    D = load("lung"); n = D["n"]; res = dict(n_V=n, n_edges=D["n_edges"], n_nonpairs=D["n_nonpairs"], layers={})
    sp_full, sp_k, se_full, se_k = [], [], [], []
    node_full, node_r5 = [], []
    for L in range(12):
        US, Xm = us(D, L)
        Ef = unit(US); Ec = unit(Xm)
        cmax = float(np.abs((Ef[:300] @ Ef[:300].T) - (Ec[:300] @ Ec[:300].T)).max())
        a_f, s_f = ppi(Ef, D); a_k, s_k = ppi(unit(US[:, :kp]), D)
        e_f, t_f = edge(Ef, D); e_k, t_k = edge(unit(US[:, :ke]), D)
        nf = node(Ef, D); nr5 = node(unit(US[:, 5:]), D)
        sp_full.append(s_f); sp_k.append(s_k); se_full.append(t_f); se_k.append(t_k); node_full.append(nf); node_r5.append(nr5)
        res["layers"][L] = dict(ppi_full=a_f, ppi_k=a_k, edge_full=e_f, edge_k=e_k, node_full=nf, node_R5=nr5, full_vs_centred=cmax)
        log("lung layer", L)
    lay = res["layers"]
    res["layer_mean"] = {k: float(np.mean([lay[L][k] for L in range(12)])) for k in lay[0]}
    res["delta_ppi"] = res["layer_mean"]["ppi_k"] - res["layer_mean"]["ppi_full"]
    res["delta_edge"] = res["layer_mean"]["edge_k"] - res["layer_mean"]["edge_full"]
    res["diff_node_R5"] = res["layer_mean"]["node_R5"] - res["layer_mean"]["node_full"]
    # PPI vertex subsampling of the paired layer-mean difference, own rng
    rng = np.random.default_rng(777); B = 2000; m = n // 2; dr = np.full(B, np.nan)
    for r in range(B):
        keep = np.zeros(n, bool); keep[rng.permutation(n)[:m]] = True
        e = keep[D["pa"]] & keep[D["pb"]]; yy = D["py"][e]
        dr[r] = np.mean([roc_auc_score(yy, sp_k[L][e]) - roc_auc_score(yy, sp_full[L][e]) for L in range(12)])
        if r % 500 == 0: log("ppi draw", r)
    est = res["delta_ppi"]; ci = {}
    for nm, sc in (("prereg_sqrt_m_over_n", np.sqrt(m / n)), ("errata_delete_d", np.sqrt(m / (n - m)))):
        dev = (dr - dr.mean()) * sc; ci[nm] = [est - float(np.quantile(dev, .975)), est - float(np.quantile(dev, .025))]
    res["ppi_ci"] = ci
    # EDGE two-way TF x target bootstrap of the paired layer-mean difference, own rng
    tfu, tfc = np.unique(D["tf"], return_inverse=True); tgu, tgc = np.unique(D["tg"], return_inverse=True)
    rng = np.random.default_rng(778); B = 2000; bd = np.empty(B); bs = 100
    for a0 in range(0, B, bs):
        wt = np.stack([np.bincount(rng.integers(0, tfu.size, tfu.size), minlength=tfu.size) for _ in range(bs)])
        wg = np.stack([np.bincount(rng.integers(0, tgu.size, tgu.size), minlength=tgu.size) for _ in range(bs)])
        W = (wt[:, tfc] * wg[:, tgc]).astype(np.float64)
        acc = np.zeros(bs)
        for L in range(12):
            acc += wstrat_tf(se_k[L], D["y"], tfc, W) - wstrat_tf(se_full[L], D["y"], tfc, W)
        bd[a0:a0 + bs] = acc / 12
    # sanity: unit weights reproduce the point estimate
    one = np.ones((1, D["y"].size))
    assert abs(wstrat_tf(se_full[0], D["y"], tfc, one)[0] - lay[0]["edge_full"]) < 1e-9
    res["edge_ci"] = [float(np.nanquantile(bd, .025)), float(np.nanquantile(bd, .975))]
    (HERE / "verify_lung.json").write_text(json.dumps(res, indent=1, default=float))
    log(json.dumps({k: res[k] for k in ("delta_ppi", "delta_edge", "diff_node_R5", "ppi_ci", "edge_ci")}))


def stage_anyk_lungext():
    D = load("lung_external"); n = D["n"]; N = 500; rng = np.random.default_rng(999)
    y = D["is_tf"]; st = D["det"]; Y = np.empty((n, N), bool)
    for r in range(N):
        yy = y.copy()
        for v in np.unique(st):
            idx = np.flatnonzero(st == v); yy[idx] = y[idx][rng.permutation(idx.size)]
        Y[:, r] = yy
    assert np.all(Y.sum(0) == y.sum())
    groups = [np.flatnonzero(st == v) for v in np.unique(st)]; k = y.sum()
    obs = np.empty((len(K_T), 12)); nul = np.empty((N, len(K_T), 12))
    M = Y.astype(np.float64)
    for L in range(12):
        US, _ = us(D, L)
        for i, kk in enumerate(K_T):
            En = unit(US[:, :kk]); obs[i, L] = node(En, D)
            tot = En @ (En.T @ M); sc = np.where(Y, (tot - 1.0) / (k - 1), tot / k)
            num = np.zeros(N); den = np.zeros(N)
            for g in groups:
                R = rankdata(sc[g], axis=0); P = Y[g]; n1 = P.sum(0); n0 = g.size - n1
                num += (R * P).sum(0) - n1 * (n1 + 1) / 2; den += n1 * n0
            nul[:, i, L] = num / den
        log("anyk layer", L)
    iF = K_T.index(512); red = [i for i in range(len(K_T)) if i != iF]
    g_obs = obs.mean(1)[red] - obs.mean(1)[iF]; g_null = nul.mean(2)[:, red] - nul.mean(2)[:, [iF]]
    T = float(g_obs.max()); Tn = g_null.max(1); p = (np.sum(Tn >= T) + 1) / (N + 1)
    res = dict(observed_max_gap=T, argmax_k=K_T[red[int(np.argmax(g_obs))]], p=float(p), null_q95=float(np.quantile(Tn, .95)))
    (HERE / "verify_anyk_lungext.json").write_text(json.dumps(res, indent=1)); log(res)


def stage_static():
    gs = np.load(R3 / "E6_crossmodel/work/geneset_tierB_lung.npz")
    genes = [str(g) for g in gs["genes"]]
    z = np.load(X.REPO / "subproject_56_geneformer_v2_104m_gpl_replication/embeddings/geneformer_v2_104m_gene_embeddings.npz")
    g = np.array([str(x) for x in z["gene_names"]]); E = np.asarray(z["embeddings"], np.float64)
    nz = np.abs(E).sum(1) > 0; g, E = g[nz], E[nz]; idx = {x: i for i, x in enumerate(g)}
    Xm = E[[idx[x] for x in genes]]; Xm = Xm - Xm.mean(0); U, s, _ = np.linalg.svd(Xm, full_matrices=False); US = U * s
    py = np.r_[np.ones(gs["pi"].size, bool), np.zeros(gs["ni"].size, bool)]
    pa = np.r_[gs["pi"], gs["ni"]]; pb = np.r_[gs["pj"], gs["nj"]]
    res = {}
    for kk in (7, US.shape[1]):
        En = unit(US[:, :kk]); res[f"ppi_k{kk}"] = float(roc_auc_score(py, (En[pa] * En[pb]).sum(1)))
    (HERE / "verify_static.json").write_text(json.dumps(res, indent=1)); log(res)


if __name__ == "__main__":
    dict(select=stage_select, lung=stage_lung, anyk_lungext=stage_anyk_lungext, static=stage_static)[sys.argv[1]]()
