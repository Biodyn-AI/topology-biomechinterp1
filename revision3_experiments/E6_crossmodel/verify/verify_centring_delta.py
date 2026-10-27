"""Verifier: paired gene-bootstrap delta (ESM2-15B minus scGPT static / Geneformer V2-104M) for the TRRUST TF/target-only
AUROC under the two centring conventions (Tier A; and Tier B lung incl. scGPT contextual L1/L11).
Own implementation: labelled genes resampled with replacement (same draws for every embedding), leave-gene-out centroid.
Writes verify/centring_delta.csv."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_e6 as V

def normed(X):
    X = np.asarray(X, np.float64); X = X - X.mean(0); return X / np.linalg.norm(X, axis=1, keepdims=True)

def boot(Xn_lab, ytf, B, seed):
    """Xn_lab: unit rows of labelled genes (already centred by the chosen convention); returns (estimate, draws)."""
    tf = np.flatnonzero(ytf); G = Xn_lab @ Xn_lab[tf].T; n = ytf.size; k = tf.size
    s = np.where(ytf, (G.sum(1) - 1) / (k - 1), G.sum(1) / k); est = roc_auc_score(ytf, s)
    rng = np.random.default_rng(seed); out = np.empty(B)
    for r in range(B):
        c = np.bincount(rng.integers(0, n, n), minlength=n); ct = c[tf]; K = ct.sum(); tot = G @ ct
        sc = np.where(ytf, (tot - c * 1.0) / np.maximum(K - c, 1), tot / K)   # leave-gene-out (self copies excluded)
        w = c > 0; yy = np.repeat(ytf[w], c[w]); ss = np.repeat(sc[w], c[w]); out[r] = roc_auc_score(yy, ss)
    return est, out

voc = {k: V.load_model(k) for k in ("esm2_uce", "scgpt_static", "gf_v2_104m")}
allv = [set(V.load_model(k)[0]) for k in V.MODELS]; A = sorted(set.intersection(*allv))
T = V.trrust(); TF = set(T.TF); TGO = set(T.target) - TF
rows = []
def run(setname, genes, embs):
    ttf = np.array([g in TF for g in genes]); lab = ttf | np.array([g in TGO for g in genes])
    res = {}
    for name, X in embs.items():
        for conv in ("centre_set", "centre_labelled"):
            Xn = normed(X)[lab] if conv == "centre_set" else normed(X[lab])
            res[(name, conv)] = boot(Xn, ttf[lab], 1000, 2026)
    for conv in ("centre_set", "centre_labelled"):
        for a in embs:
            for b in embs:
                if a == b or a not in ("esm2_uce", "gf_v2_104m"): continue
                if a == "gf_v2_104m" and b == "esm2_uce": continue
                ea, da = res[(a, conv)]; eb, db = res[(b, conv)]; d = ea - eb; dd = da - db; sh = d - dd.mean()
                rows.append(dict(set=setname, convention=conv, model_a=a, model_b=b, value_a=ea, value_b=eb, delta=d,
                                 ci_lo=np.quantile(dd, .025) + sh, ci_hi=np.quantile(dd, .975) + sh, n_tf=int(ttf[lab].sum()), n_lab=int(lab.sum())))
                print(rows[-1], flush=True)
run("tierA", A, {k: voc[k][1][[voc[k][0][g] for g in A]] for k in voc})
lg, li = V.ctx_obs("lung"); emb = np.load(V.EMB / "lung.npy", mmap_mode="r"); posl = {g: i for i, g in enumerate(lg)}
Bg = sorted(set(lg) & set(A)); r = np.array([posl[g] for g in Bg])
embs = {k: voc[k][1][[voc[k][0][g] for g in Bg]] for k in voc}
for L in (1, 11): embs[f"lung@L{L}"] = np.asarray(emb[L], np.float64)[li][r]
run("tierB_lung", Bg, embs)
pd.DataFrame(rows).to_csv(V.HERE / "centring_delta.csv", index=False)
