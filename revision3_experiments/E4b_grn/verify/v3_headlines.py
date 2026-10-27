"""Verifier check 3: independent recomputation of headline numbers from raw inputs (no stats3 / e4b code).

AUROC: per-TF sklearn roc_auc_score (ties 1/2), pooled with weights n_pos*n_neg (= within-TF stratified Mann-Whitney).
Bootstrap: own pigeonhole TF x target implementation, B = 500, own seed (Monte-Carlo comparison, not bit-equality).
Permutation: own within-TF label permutation, N = 2,000, own seed.
Direction: own lookup of importance(i->j) vs importance(j->i) from the partial npz files, scipy binomtest.
Writes verify/v3_headlines.json.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.stats import binomtest
from sklearn.metrics import roc_auc_score

V = Path(__file__).resolve().parent; E = V.parent; CACHE = E.parent / "cache"
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"


def pairs(ctx):
    h = np.load(CACHE / ctx / "h3_pairs.npz"); return h["tf"].astype(int), h["tg"].astype(int), h["y"].astype(bool)


class Strat:
    def __init__(self, tf, tg, y):
        self.tf, self.tg, self.y = tf, tg, y
        self.utf, self.ctf = np.unique(tf, return_inverse=True); self.utg, self.ctg = np.unique(tg, return_inverse=True)
        self.groups = [np.flatnonzero(self.ctf == c) for c in range(self.utf.size)]

    def auc(self, s, w=None):
        num = den = 0.0
        for g in self.groups:
            yy = self.y[g]; ww = None if w is None else w[g]
            if ww is not None:
                k = ww > 0; g2 = g[k]; yy = yy[k]; ww = ww[k]
            else:
                g2 = g
            if yy.all() or (~yy).any() is False or yy.sum() == 0: continue
            wp = yy.sum() if ww is None else ww[yy].sum(); wn = (~yy).sum() if ww is None else ww[~yy].sum()
            a = roc_auc_score(yy, s[g2], sample_weight=ww)
            num += a * wp * wn; den += wp * wn
        return num / den

    def boot(self, scores, B, seed):
        rng = np.random.default_rng(seed); out = np.empty((B, len(scores)))
        for r in range(B):
            w = np.bincount(rng.integers(0, self.utf.size, self.utf.size), minlength=self.utf.size)[self.ctf] * \
                np.bincount(rng.integers(0, self.utg.size, self.utg.size), minlength=self.utg.size)[self.ctg]
            w = w.astype(float)
            out[r] = [self.auc(s, w) for s in scores]
        return out


def ens_score(ctx, tier, method, seeds, src, dst):
    mats = []
    for sd in seeds:
        z = np.load(E / "partial" / f"{ctx}_{tier}_{method}_s{sd}.npz"); mats.append(z["imp"].astype(np.float64))
        regs, tg = z["regs"], z["targets"]
    M = np.mean(mats, 0); rp = {int(r): k for k, r in enumerate(regs)}; tp = {int(t): k for k, t in enumerate(tg)}
    per = [np.array([m[tp[int(b)], rp[int(a)]] for a, b in zip(src, dst)]) for m in mats]
    return np.array([M[tp[int(b)], rp[int(a)]] for a, b in zip(src, dst)]), per


def abscor(X, a, b):
    X = np.asarray(X, np.float64)
    out = np.empty(a.size)
    for k, (i, j) in enumerate(zip(a, b)):
        xi, xj = X[:, i], X[:, j]
        out[k] = 0.0 if xi.std() == 0 or xj.std() == 0 else abs(np.corrcoef(xi, xj)[0, 1])
    return out


def cosL(ctx, layer, a, b):
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r"); oi = np.load(CACHE / ctx / "obs_idx.npy")
    M = np.asarray(emb[layer], np.float64)[oi]; M = M - M.mean(0)
    return np.sum(M[a] * M[b], 1) / (np.linalg.norm(M[a], axis=1) * np.linalg.norm(M[b], axis=1))


res = {}
B = 500
# (1) lung GENIE3 T1 ensemble + per-seed SD of GRNBoost2 T1
tf, tg, y = pairs("lung"); st = Strat(tf, tg, y)
s, per = ens_score("lung", "T1", "genie3", range(5), tf, tg)
bt = st.boot([s], B, 999)[:, 0]
res["lung_genie3_T1_ens"] = dict(auroc=st.auc(s), ci=[float(np.quantile(bt, .025)), float(np.quantile(bt, .975))],
                                 upper95=float(np.quantile(bt, .95)), frac_zero=float(np.mean(s == 0)))
_, perb = ens_score("lung", "T1", "grnboost2", range(5), tf, tg)
pa = [st.auc(p) for p in perb]
res["lung_grnboost2_T1_per_seed"] = dict(aurocs=pa, sd=float(np.std(pa, ddof=1)))
print(res, flush=True)
# (2) lung_external |coexpr| T2 and GRNBoost2 T2 zero fraction
tf, tg, y = pairs("lung_external"); st = Strat(tf, tg, y)
c2 = abscor(np.load(CACHE / "lung_external" / "expr.npy"), tf, tg)
sb, _ = ens_score("lung_external", "T2", "grnboost2", [0, 1], tf, tg)
bt = st.boot([c2], B, 998)[:, 0]
res["lung_external_abscoexpr_T2"] = dict(auroc=st.auc(c2), ci=[float(np.quantile(bt, .025)), float(np.quantile(bt, .975))])
res["lung_external_grnboost2_T2_ens_frac_zero"] = float(np.mean(sb == 0))
print(res, flush=True)
# (3) kidney_B scGPT signed cosine L0 and L11
tf, tg, y = pairs("kidney_B"); st = Strat(tf, tg, y)
res["kidney_B_scgpt_L0"] = st.auc(cosL("kidney_B", 0, tf, tg)); res["kidney_B_scgpt_L11"] = st.auc(cosL("kidney_B", 11, tf, tg))
print(res, flush=True)
# (4) immune GENIE3 T2 minus |coexpr| T2 (paired bootstrap, same weights) + |coexpr| T1 permutation p
tf, tg, y = pairs("immune"); st = Strat(tf, tg, y)
g2, _ = ens_score("immune", "T2", "genie3", [0, 1], tf, tg)
c2 = abscor(np.load(CACHE / "immune" / "expr.npy"), tf, tg); c1 = abscor(np.load(CACHE / "immune" / "expr_extraction.npy"), tf, tg)
bt = st.boot([g2, c2], 1000, 997); d = bt[:, 0] - bt[:, 1]
res["immune_genie3_T2_minus_abscoexpr_T2"] = dict(a=st.auc(g2), b=st.auc(c2), delta=st.auc(g2) - st.auc(c2),
                                                  ci=[float(np.quantile(d, .025)), float(np.quantile(d, .975))])
rng = np.random.default_rng(996); a0 = st.auc(c1); nul = np.empty(2000)
for r in range(2000):
    yy = y.copy()
    for g in st.groups: yy[g] = rng.permutation(yy[g])
    st.y = yy; nul[r] = st.auc(c1)
st.y = y
bt1 = st.boot([c1], B, 995)[:, 0]
res["immune_abscoexpr_T1"] = dict(auroc=a0, perm_p=float((np.sum(nul >= a0) + 1) / 2001), null_sd=float(nul.std(ddof=1)),
                                  ci=[float(np.quantile(bt1, .025)), float(np.quantile(bt1, .975))], boot_sd=float(bt1.std(ddof=1)))
print(res, flush=True)
# (5) direction test kidney_B GENIE3 T2 ensemble and lung GRNBoost2 T1 ensemble
for ctx, tier, method, seeds in [("kidney_B", "T2", "genie3", [0, 1]), ("lung", "T1", "grnboost2", range(5))]:
    vj = json.loads((V / "v1_inputs.json").read_text())  # pair counts only
    sys.path.insert(0, str(E)); import e4b_common as C
    dp = C.load_sets(ctx)["dir_pairs"]  # equality with the raw-file rebuild was checked in v1
    f, _ = ens_score(ctx, tier, method, seeds, dp[:, 0], dp[:, 1]); r_, _ = ens_score(ctx, tier, method, seeds, dp[:, 1], dp[:, 0])
    k = int((f > r_).sum()); n = int((f != r_).sum()); bt_ = binomtest(k, n, 0.5)
    ci = bt_.proportion_ci(0.95, "exact")
    res[f"direction_{ctx}_{method}_{tier}"] = dict(n_pairs=int(dp.shape[0]), n_untied=n, k=k, frac=k / n, p=float(bt_.pvalue),
                                                 cp=[float(ci.low), float(ci.high)])
    print(res, flush=True)
(V / "v3_headlines.json").write_text(json.dumps(res, indent=1, default=float))
