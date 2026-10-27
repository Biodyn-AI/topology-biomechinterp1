"""Compartment axis: is it specific to the expression-trained residual stream?
Compares, on one common gene set per context (V ∩ ESM2 ∩ scGPT static vocabulary):
  scGPT residual stream L0, L3, L11 | scGPT static input-token embedding | ESM2-15B per-protein embedding (UCE file)
Two read-outs per representation and term (extracellular space GO:0005615, cytosol GO:0005829):
  (a) unsupervised: sign-free AUROC of each of the first 12 directions of the representation's own centred SVD on the
      common genes; reported for SV1 and for the best of SV1-12 (max over 12, so optimistic by construction)
  (b) supervised: out-of-fold AUROC of an L2 logistic probe on the first 256 principal components (5 outer folds,
      C chosen by inner 5-fold CV), with a 2,000-replicate gene bootstrap CI; detection-stratified AUROC of the same
      out-of-fold scores; paired bootstrap difference vs scGPT L3.
Sensitivity analysis (CONVENTION3 §5); no pass/fail. Writes static_baselines.csv, static_baselines_paired.csv."""
import json, sys, time, warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
from e11_common import load_emb, load_cache, load_gene2go, convention_hash
import stats3 as S
E3P = HERE.parent / "E3_node" / "partial"
CTXS = ["lung", "lung_external", "kidney_B", "immune"]
TERMS = {"GO:0005615": "extracellular_space", "GO:0005829": "cytosol"}
SEED, B = 20261026, 2000


def auc(score, y):
    """AUROC of `score` for boolean labels y (stats3.auroc takes (pos, neg))."""
    score = np.asarray(score, float); return S.auroc(score[y], score[~y])


def signfree(x, y):
    a = auc(x, y); return max(a, 1 - a)


def strat_auc(score, y, strata):
    num = den = 0.0
    for s in np.unique(strata):
        m = strata == s; yp, yn = score[m & y], score[m & ~y]
        if yp.size and yn.size:
            r = np.argsort(np.argsort(np.r_[yp, yn])) + 1.0
            num += r[:yp.size].sum() - yp.size * (yp.size + 1) / 2; den += yp.size * yn.size
    return num / den


def oof_probe(Xf, y, seed):
    Z = PCA(n_components=min(256, Xf.shape[1]), random_state=seed).fit_transform(StandardScaler().fit_transform(Xf))
    out = np.zeros(len(y))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(Z, y):
        m = LogisticRegressionCV(Cs=np.logspace(-4, 1, 8), cv=5, scoring="roc_auc", max_iter=4000).fit(Z[tr], y[tr])
        out[te] = m.decision_function(Z[te])
    return out


def main():
    assert convention_hash()[2]
    g2go = load_gene2go(); rows, prow = [], []
    for ctx in CTXS:
        t0 = time.time()
        genes, obs_idx, nl, _ = load_cache(ctx)
        ze, zs = np.load(E3P / f"static_esm2_{ctx}.npz"), np.load(E3P / f"static_scgpt_{ctx}.npz")
        common = np.intersect1d(ze["idx"], zs["idx"])
        ie = {k: r for r, k in enumerate(ze["idx"])}; isg = {k: r for r, k in enumerate(zs["idx"])}
        emb = load_emb(ctx)
        reps = {f"scGPT_L{L}": np.asarray(emb[L][obs_idx[common]], float) for L in (0, 3, 11)}
        reps["scGPT_static_token"] = zs["emb"][[isg[k] for k in common]].astype(float)
        reps["ESM2_15B"] = ze["emb"][[ie[k] for k in common]].astype(float)
        strata = nl["det_decile"][common]
        rng = np.random.default_rng(SEED)
        boot = [rng.integers(0, common.size, common.size) for _ in range(B)]
        for tid, tname in TERMS.items():
            y = np.array([tid in g2go.get(genes[k], ()) for k in common])
            oof = {}
            for rname, M in reps.items():
                Mc = M - M.mean(0); U, s, _ = np.linalg.svd(Mc, full_matrices=False); C = U[:, :12] * s[:12]
                au = [signfree(C[:, k], y) for k in range(12)]
                sc = oof_probe(M, y, SEED); oof[rname] = sc
                a = auc(sc, y); bs = np.array([auc(sc[b], y[b]) for b in boot])
                rows.append(dict(context=ctx, term=tname, representation=rname, n_genes=common.size, n_pos=int(y.sum()),
                                 sv1_auroc=au[0], best_sv=int(np.argmax(au)) + 1, best_sv_auroc=max(au),
                                 probe_auroc=a, probe_lo=np.quantile(bs, .025), probe_hi=np.quantile(bs, .975),
                                 probe_strat_auroc=strat_auc(sc, y, strata)))
                print(f"[{time.time()-t0:5.0f}s] {ctx} {tname} {rname}: SV1 {au[0]:.3f} best SV{np.argmax(au)+1} {max(au):.3f} "
                      f"probe {a:.3f} strat {rows[-1]['probe_strat_auroc']:.3f}", flush=True)
            for rname in reps:
                if rname == "scGPT_L3": continue
                d = np.array([auc(oof[rname][b], y[b]) - auc(oof["scGPT_L3"][b], y[b]) for b in boot])
                prow.append(dict(context=ctx, term=tname, representation=rname, vs="scGPT_L3",
                                 delta=auc(oof[rname], y) - auc(oof["scGPT_L3"], y),
                                 lo=np.quantile(d, .025), hi=np.quantile(d, .975)))
        pd.DataFrame(rows).to_csv(HERE / "static_baselines.csv", index=False)
        pd.DataFrame(prow).to_csv(HERE / "static_baselines_paired.csv", index=False)
    print("done")


if __name__ == "__main__":
    main()
