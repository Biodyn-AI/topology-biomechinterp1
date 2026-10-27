"""E7 (H4, rank reduction vs full dimension): shared data layer, representations, metrics, resampling.

Point estimates ALWAYS call stats3 directly:
  PPI   stats3.auroc(cos of STRING-physical>=700 edges, cos of ONE Maslov-Sneppen rewired non-pair set)
  NODE  stats3.stratified_auroc(stats3.tf_scores_loo(E, Lambert TFs), is_tf_lambert, det_decile)   (= H2 statistic)
  EDGE  stats3.StratifiedPairs(h3_pairs).auroc(signed cosine)                                     (= H3 statistic)
Representations: U_k S_k columns of stats3.svd_observed (observed genes V, mean-centred over V); k = 512 is
the full representation (cosines of the full U S rows equal cosines of the V-centred embedding; checked).

Vectorised paths exist ONLY for resampling / permutation loops and are verified against stats3 on live data
every run (verify_* functions; mismatch > 1e-9 aborts):
  node_boot_scores + wstrat_auroc_ties   gene bootstrap of the H2 statistic (weighted, leave-gene-out, ties 1/2);
                                          equals stats3 at unit weights and equals stats3 on the explicitly
                                          row-duplicated gene set for a bootstrap replicate;
  FastStrat (E4a_edges/e4a_lib.py, imported read-only)   two-way TF x target bootstrap with the SAME integer
                                          draws as StratifiedPairs.bootstrap, and the within-TF permutation null;
  ppi_vertex_draws                        node subsets drawn exactly as stats3.vertex_subsample_ci draws them
                                          (same rng calls), statistic re-evaluated for many score vectors.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import contexts as X  # noqa: E402
import stats3 as S    # noqa: E402

WORK = HERE / "work"; WORK.mkdir(exist_ok=True)
CACHE = R3 / "cache"
LOCAL_EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
SEED = 20261026
N_LAYERS = 12
SELECTION_CTX = "immune_seed43"
EVAL_CTX = ["lung", "lung_external", "kidney_B", "immune"]
EXTRA_CTX = ["kidney_A"]                                   # curves only (where the legacy mapping is defined)
K_T = [1, 2, 3, 5, 7, 10, 16, 32, 64, 128, 256, 512]
J_R = [1, 2, 3, 5, 10]
B_BOOT = 2000
N_PERM_ANYK = 500
TASKS = ["PPI", "NODE", "EDGE"]

# family N: (name, first SV, last SV), 1-based inclusive, in the observed basis of the context
BANDS_NOMINAL = [("SV2-4", 2, 4), ("SV5-7", 5, 7), ("SV2-7", 2, 7)]
# E9_screen/legacy_band_mapping.json, kidney_A 'legacy_full_zero' (submitted-manuscript convention):
# legacy SV2-4 -> observed SV1-3, SV5-7 -> SV4-6, SV2-7 -> SV1-6. On immune the screen's bands map to themselves.
BANDS_KIDNEYA_RULE = [("SV2-4", 1, 3), ("SV5-7", 4, 6), ("SV2-7", 1, 6)]


def convention_hash():
    txt = (R3 / "FROZEN.txt").read_text().split()
    frozen = txt[txt.index("sha256") + 1]
    actual = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
    return dict(frozen=frozen, actual=actual, match=frozen == actual)


def file_sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_faststrat():
    """E4a's FastStrat, imported read-only from its file (no bytecode written)."""
    sys.dont_write_bytecode = True
    p = R3 / "E4a_edges" / "e4a_lib.py"
    spec = importlib.util.spec_from_file_location("e4a_lib_ro", p)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod, file_sha(p)


def check_mapping():
    """Read the fixed legacy band mapping and confirm the rule hard-coded above."""
    m = json.loads((R3 / "E9_screen" / "legacy_band_mapping.json").read_text())
    got = {}
    for b in m["bands"]:
        got[(b["context"], b["legacy_convention"], b["legacy_band"])] = b["mapped_observed_band"]
    for (nm, a, z) in BANDS_KIDNEYA_RULE:
        assert got[("kidney_A", "legacy_full_zero", nm)] == f"SV{a}-{z}", (nm, got)
    for nm in ("SV2-4", "SV5-7", "SV2-7"):
        assert got[("immune", "screen_nonzero", nm)] == nm
    return {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in got.items()}


# =====================================================================================
# data
# =====================================================================================
def load_emb(ctx):
    p = LOCAL_EMB / f"{ctx}.npy"
    return np.load(p if p.exists() else X.CONTEXTS[ctx]["emb"], mmap_mode="r")


def rewired_nonpairs(edges, n, seed):
    """ONE Maslov-Sneppen degree-preserving rewiring (stats3.maslov_sneppen_chain, burn-in 10m accepted swaps,
    then m more) of the V-induced STRING graph; rewired pairs equal to a true edge are dropped."""
    a = np.minimum(edges[:, 0], edges[:, 1]).astype(np.int64); b = np.maximum(edges[:, 0], edges[:, 1]).astype(np.int64)
    key = np.unique(a * n + b); a, b = key // n, key % n
    ra, rb, acc = next(S.maslov_sneppen_chain(a, b, n, 1, seed=seed, burn_factor=10, thin_factor=1))
    keep = ~np.isin(ra * n + rb, key)
    return np.c_[a, b], np.c_[ra[keep], rb[keep]], float(acc), int((~keep).sum())


def load_ctx(ctx):
    d = CACHE / ctx
    genes = json.loads((d / "obs_genes.json").read_text()); obs_idx = np.load(d / "obs_idx.npy")
    nl = dict(np.load(d / "node_labels.npz")); n = len(genes)
    is_tf = nl["is_tf_lambert"].astype(bool); strata = nl["det_decile"].astype(int)
    z = np.load(d / "h3_pairs.npz"); sp = S.StratifiedPairs(z["tf"], z["tg"], z["y"])
    edges = np.load(d / "edges_string.npy")
    epos, eneg, acc, n_drop = rewired_nonpairs(edges, n, SEED)
    pa = np.r_[epos[:, 0], eneg[:, 0]]; pb = np.r_[epos[:, 1], eneg[:, 1]]
    py = np.r_[np.ones(len(epos), bool), np.zeros(len(eneg), bool)]
    return dict(ctx=ctx, genes=genes, obs_idx=obs_idx, n=n, is_tf=is_tf, tf=np.flatnonzero(is_tf), strata=strata,
                sp=sp, pa=pa, pb=pb, py=py, ms_acc=acc, ms_true_dropped=n_drop,
                meta=dict(n_V=n, n_tf=int(is_tf.sum()), n_string=int(len(epos)), n_rewired_nonpairs=int(len(eneg)),
                          rewired_equal_true_dropped=n_drop, ms_acceptance=acc, h3_pairs=int(sp.y.size),
                          h3_pos=int(sp.y.sum()), h3_tfs=int(sp.tfs.size)))


def layer_svds(emb, obs_idx):
    out = []
    for L in range(N_LAYERS):
        US, Sv, _ = S.svd_observed(emb, L, obs_idx); out.append((US, Sv))
    return out


def check_full_equals_centred(emb, obs_idx, US, L, n_check=400):
    Xm = np.asarray(emb[L], dtype=np.float64)[obs_idx]; Xm = Xm - Xm.mean(0, keepdims=True)
    A = S.l2norm_rows(Xm[:n_check]); Bm = S.l2norm_rows(US[:n_check])
    return float(np.max(np.abs(A @ A.T - Bm @ Bm.T)))


# =====================================================================================
# representations
# =====================================================================================
def representations(ctx):
    """List of (rep_id, family, param, col_slice). Columns are 0-based into U S (512 columns)."""
    R = [(f"T{k}", "T", k, slice(0, k)) for k in K_T]
    R.append(("T6", "T_aux", 6, slice(0, 6)))                         # same-width partner of 6-wide bands
    R += [(f"R{j}", "R", j, slice(j, 512)) for j in J_R]
    R += [(f"Nnom_{nm}", "N_nominal", f"SV{a}-{z}", slice(a - 1, z)) for nm, a, z in BANDS_NOMINAL]
    R += [(f"Nmap_{nm}", "N_kidneyA_rule", f"SV{a}-{z}", slice(a - 1, z)) for nm, a, z in BANDS_KIDNEYA_RULE]
    return R


def unit_rep(US, cols):
    """Exact row normalisation (verifier fix, 2026-10-26). stats3.l2norm_rows divides by (norm + 1e-12); at k = 1
    that maps x to sign(x) * (1 - 1e-12/|x|), so the cosines of a 1-D representation, which are exactly +-1 and
    should tie, were instead ordered by |SV1| magnitude, shifting k = 1 PPI / EDGE AUROCs by ~0.01. Rows of U_k S_k
    are never exactly zero here (asserted); for k >= 2 the change is below 1e-12 and does not alter any rank."""
    M = np.ascontiguousarray(US[:, cols]); nr = np.linalg.norm(M, axis=1, keepdims=True)
    assert np.all(nr > 0), "zero row in U_k S_k"
    En = M / nr
    if M.shape[1] == 1:
        En = np.sign(M)                                    # exact +-1 (division can leave 1 ulp noise)
    return En


# =====================================================================================
# point metrics (stats3)
# =====================================================================================
def pair_cos(En, a, b, chunk=200000):
    out = np.empty(a.size)
    for i in range(0, a.size, chunk):
        out[i:i + chunk] = np.einsum("ij,ij->i", En[a[i:i + chunk]], En[b[i:i + chunk]])
    return out


def metric_ppi(En, D, return_scores=False):
    s = pair_cos(En, D["pa"], D["pb"]); v = S.auroc(s[D["py"]], s[~D["py"]])
    return (v, s) if return_scores else v


def metric_node(En, D, return_scores=False):
    s = S.tf_scores_loo(En, D["tf"]); v = S.stratified_auroc(s, D["is_tf"], D["strata"])
    return (v, s) if return_scores else v


def metric_edge(En, D, return_scores=False):
    sp = D["sp"]; s = pair_cos(En, sp.tf, sp.tg); v = sp.auroc(s)
    return (v, s) if return_scores else v


METRIC = dict(PPI=metric_ppi, NODE=metric_node, EDGE=metric_edge)


# =====================================================================================
# NODE: weighted gene bootstrap (leave-gene-out), tie-aware weighted stratified AUROC
# =====================================================================================
def wstrat_auroc_ties(Sc, pos, strata, W):
    """Weighted Mann-Whitney AUROC pooled within strata, ties 1/2; columns = replicates.
    Sc (n, b) scores, pos (n,) bool, strata (n,), W (n, b) nonnegative weights; pair weight w_p * w_q."""
    b = Sc.shape[1]; num = np.zeros(b); den = np.zeros(b)
    for st in np.unique(strata):
        m = strata == st; p = pos[m]
        if p.sum() == 0 or (~p).sum() == 0: continue
        Ss = Sc[m]; Ws = W[m]; ns = Ss.shape[0]
        o = np.argsort(Ss, axis=0, kind="mergesort")
        Ss = np.take_along_axis(Ss, o, 0); Ws = np.take_along_axis(Ws, o, 0); P = p[o]
        wn = np.where(P, 0.0, Ws); wp = np.where(P, Ws, 0.0)
        cneg = np.cumsum(wn, 0)
        ar = np.arange(ns)[:, None]
        diff = Ss[1:] != Ss[:-1]
        newg = np.vstack([np.ones((1, b), bool), diff])
        start = np.maximum.accumulate(np.where(newg, ar, 0), axis=0)
        endf = np.vstack([diff, np.ones((1, b), bool)])
        end = np.minimum.accumulate(np.where(endf, ar, ns - 1)[::-1], axis=0)[::-1]
        before = np.take_along_axis(cneg - wn, start, 0)
        through = np.take_along_axis(cneg, end, 0)
        num += (wp * (before + 0.5 * (through - before))).sum(0)
        den += wp.sum(0) * wn.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def node_boot_scores(En, tf, W):
    """Bootstrap version of stats3.tf_scores_loo for gene multiplicities W (n, b): TF centroid rebuilt from the
    drawn TF copies; every gene scored leave-GENE-out (all copies of itself excluded)."""
    C = En @ En[tf].T                                 # (n, ntf)
    Wt = W[tf]                                        # (ntf, b)
    tot = C @ Wt                                      # (n, b)
    K = Wt.sum(0, keepdims=True)
    selfs = (En * En).sum(1)
    out = tot / np.maximum(K, 1)
    out[tf] = (tot[tf] - Wt * selfs[tf, None]) / np.maximum(K - Wt, 1)
    return out


def gene_boot_weights(n, B, seed):
    rng = np.random.default_rng(seed); W = np.empty((n, B), np.uint16)
    for r in range(B): W[:, r] = np.bincount(rng.integers(0, n, n), minlength=n)
    return W


def node_boot(En, D, W, batch=200):
    out = np.empty(W.shape[1])
    for a in range(0, W.shape[1], batch):
        Wb = W[:, a:a + batch].astype(np.float64)
        out[a:a + batch] = wstrat_auroc_ties(node_boot_scores(En, D["tf"], Wb), D["is_tf"], D["strata"], Wb)
    return out


def verify_node_boot(En, D, W):
    """(1) unit weights == stats3; (2) replicate 0 == stats3.stratified_auroc on the explicitly duplicated set,
    with scores computed as the mean cosine to every drawn TF copy of a DIFFERENT gene."""
    n = D["n"]; one = np.ones((n, 1))
    a1 = wstrat_auroc_ties(node_boot_scores(En, D["tf"], one), D["is_tf"], D["strata"], one)[0]
    ref1 = metric_node(En, D)
    assert abs(a1 - ref1) < 1e-9, ("node unit-weight mismatch", a1, ref1)
    c = W[:, 0].astype(int); idx = np.repeat(np.arange(n), c); Ed = En[idx]
    istf_d = D["is_tf"][idx]; tfrows = np.flatnonzero(istf_d)
    Cm = Ed @ Ed[tfrows].T; same = idx[:, None] == idx[tfrows][None, :]
    Cm[same] = 0.0; den = (~same).sum(1); sc = Cm.sum(1) / np.maximum(den, 1)
    ref2 = S.stratified_auroc(sc, istf_d, D["strata"][idx])
    got2 = node_boot(En, D, W[:, :1])[0]
    assert abs(ref2 - got2) < 1e-9, ("node bootstrap replicate mismatch", ref2, got2)
    return dict(unit_absdiff=float(abs(a1 - ref1)), rep0_absdiff=float(abs(ref2 - got2)))


# =====================================================================================
# PPI: vertex subsampling (stats3.vertex_subsample_ci draws), paired across score vectors
# =====================================================================================
def ppi_vertex_draws(D, score_list, B, seed, frac=0.5):
    """draws (len(score_list), B): AUROC of every score vector on the SAME vertex subsets, drawn with exactly the
    rng calls of stats3.vertex_subsample_ci (np.sort(rng.choice(n, m, replace=False)) per draw)."""
    n = D["n"]; m = int(frac * n); rng = np.random.default_rng(seed)
    pa, pb, y = D["pa"], D["pb"], D["py"]
    out = np.full((len(score_list), B), np.nan)
    for r in range(B):
        nodes = np.sort(rng.choice(n, m, replace=False)); keep = np.zeros(n, bool); keep[nodes] = True
        e = keep[pa] & keep[pb]
        if e.sum() <= 2: continue
        ep = e & y; en = e & ~y
        if ep.sum() == 0 or en.sum() == 0: continue
        for v, s in enumerate(score_list):
            out[v, r] = S.auroc(s[ep], s[en])
    return out, m


def vertex_ci_from_draws(full, draws, n, m, level=0.95):
    dev = (draws - np.nanmean(draws)) * np.sqrt(m / n); a = (1 - level) / 2
    return float(full - np.nanquantile(dev, 1 - a)), float(full - np.nanquantile(dev, a))


def verify_ppi_draws(D, s1, s2, seed, B=15):
    """Compare our paired difference draws/CI with stats3.vertex_subsample_ci applied to the difference."""
    pa, pb, y = D["pa"], D["pb"], D["py"]; n = D["n"]

    def stat(nodes, a, b):
        keep = np.zeros(n, bool); keep[nodes] = True; e = keep[pa] & keep[pb]
        return S.auroc(s1[e & y], s1[e & ~y]) - S.auroc(s2[e & y], s2[e & ~y])
    full, lo, hi, dr = S.vertex_subsample_ci(stat, n, pa, pb, B, np.random.default_rng(seed))
    mine, m = ppi_vertex_draws(D, [s1, s2], B, seed)
    d = mine[0] - mine[1]
    lo2, hi2 = vertex_ci_from_draws(full, d, n, m)
    assert np.nanmax(np.abs(d - dr)) < 1e-12 and abs(lo - lo2) < 1e-12 and abs(hi - hi2) < 1e-12, "ppi vertex mismatch"
    return dict(draw_max_absdiff=float(np.nanmax(np.abs(d - dr))), ci_absdiff=float(max(abs(lo - lo2), abs(hi - hi2))))


# =====================================================================================
# permutation nulls for the descriptive 'any k beats full' test (joint across k and layers)
# =====================================================================================
def node_perm_labels(D, N, seed):
    rng = np.random.default_rng(seed); Y = np.empty((D["n"], N), bool)
    for r in range(N): Y[:, r] = S.permute_within_strata(D["is_tf"], D["strata"], rng)
    return Y


def node_perm_null(En, D, Y, batch=250):
    """H2 statistic for every permuted TF label vector (columns of Y; within-detection-stratum permutations,
    so the per-stratum TF counts are fixed). Tie-aware (average ranks)."""
    n = D["n"]; ntf = int(D["is_tf"].sum()); selfs = (En * En).sum(1)
    strata = D["strata"]; out = np.empty(Y.shape[1])
    groups = [np.flatnonzero(strata == st) for st in np.unique(strata)]
    for a in range(0, Y.shape[1], batch):
        M = Y[:, a:a + batch].astype(np.float64)
        tot = En @ (En.T @ M)
        sc = np.where(M > 0, (tot - selfs[:, None]) / max(ntf - 1, 1), tot / ntf)
        num = np.zeros(M.shape[1]); den = np.zeros(M.shape[1])
        for g in groups:
            P = M[g] > 0; n1 = P.sum(0); n0 = g.size - n1
            ok = (n1 > 0) & (n0 > 0)
            if not ok.any(): continue
            R = rankdata(sc[g], axis=0)
            A = ((R * P).sum(0) - n1 * (n1 + 1) / 2) / np.where(ok, n1 * n0, 1)
            num += np.where(ok, A * n1 * n0, 0); den += np.where(ok, n1 * n0, 0)
        out[a:a + batch] = num / den
    return out


def verify_node_perm(En, D, Y):
    got = node_perm_null(En, D, Y[:, :3])
    ref = []
    for r in range(3):
        t = np.flatnonzero(Y[:, r]); s = S.tf_scores_loo(En, t); ref.append(S.stratified_auroc(s, Y[:, r], D["strata"]))
    d = float(np.max(np.abs(got - np.array(ref)))); assert d < 1e-9, ("node perm mismatch", got, ref)
    return d


def ppi_perm_labels(D, N, seed):
    rng = np.random.default_rng(seed); Y = np.empty((D["py"].size, N), bool)
    for r in range(N): Y[:, r] = rng.permutation(D["py"])
    return Y


def ppi_perm_null(score_cols, Y):
    """AUROC (midranks = stats3.auroc tie rule) of every score column for every permuted pair-label column."""
    R = rankdata(score_cols, axis=0)
    n1 = Y[:, 0].sum(); n0 = Y.shape[0] - n1
    return (Y.T.astype(np.float64) @ R - n1 * (n1 + 1) / 2) / (n1 * n0)


def summarise_null(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    q = np.quantile(x, [0.025, 0.05, 0.5, 0.95, 0.975, 0.99])
    return dict(N=int(x.size), mean=float(x.mean()), sd=float(x.std(ddof=1)), q025=float(q[0]), q05=float(q[1]),
                q50=float(q[2]), q95=float(q[3]), q975=float(q[4]), q99=float(q[5]), max=float(x.max()))
