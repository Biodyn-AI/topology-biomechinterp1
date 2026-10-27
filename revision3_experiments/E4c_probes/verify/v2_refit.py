"""E4c independent verifier, part 2: re-fit headline M1 probes with sklearn from raw embeddings.

Independent of e4c_lib's features, solver and metric: own V-centring, own PCA-16 (numpy SVD on training-pair
genes), own standardisation, sklearn LogisticRegression (lbfgs; C = 1 / (n_train * alpha), equivalent to the
implementer's mean-log-loss + alpha/2 ||w||^2 with unpenalised intercept), own within-TF Mann-Whitney (scipy rankdata).
Only the fold assignment (e4c_lib.make_folds; must be identical) and the per-fold alpha chosen by the implementer's
inner split (read from the checkpoint) are taken from the implementation.
Units (C3 repeat 0, 5 folds each):
  A lung          scgpt bilinear    L2  (headline lung probe, 0.597 over 5 repeats)
  B lung          scgpt tg_only     L2  (its target-only control)
  C kidney_B      scgpt tg_bilinear L0  (headline kidney probe, 0.589 over 5 repeats)
  D lung_external scgpt concat      L11 (+ direction-test forward-minus-reverse scores vs direction_L11_r0.csv)
usage: python v2_refit.py A|B|C|D    -> verify/v2_<unit>.csv
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression

HERE = Path(__file__).resolve().parent; E4C = HERE.parent; R3 = E4C.parent
sys.path.insert(0, str(E4C))
import e4c_lib as L   # make_folds only
UNITS = {"A": ("lung", "bilinear", 2), "B": ("lung", "tg_only", 2), "C": ("kidney_B", "tg_bilinear", 0),
         "D": ("lung_external", "concat", 11)}
unit = sys.argv[1]; ctx, feat, layer = UNITS[unit]


def within_tf_auc(tf, y, s):
    num = den = 0.0
    for t in np.unique(tf):
        m = tf == t; yy = y[m]; n1 = yy.sum(); n0 = m.sum() - n1
        if n1 < 3 or n0 < 3: continue
        r = rankdata(s[m]); num += r[yy].sum() - n1 * (n1 + 1) / 2; den += n1 * n0
    return num / den


d = R3 / "cache" / ctx
oi = np.load(d / "obs_idx.npy"); z = np.load(d / "h3_pairs.npz")
tf, tg, y = z["tf"].astype(int), z["tg"].astype(int), z["y"].astype(bool)
emb = np.load(Path.home() / f"biomechinterp-local/subproject41_r3/emb/{ctx}.npy", mmap_mode="r")
E = np.asarray(emb[layer], np.float64)[oi]; E -= E.mean(0)
ck = pd.read_csv(E4C / "partial" / "m1" / ctx / f"scgpt_{feat}_L{layer}_C3_r0.csv")
c = L.Ctx(ctx)
assert np.array_equal(c.tf, tf) and np.array_equal(c.tg, tg)
rows, drows = [], []
if unit == "D":
    te_ = np.load(d / "edges_trrust.npy"); tset = {(int(a), int(b)) for a, b in te_}
    lam = np.load(d / "node_labels.npz")["is_tf_lambert"]
    DIR = np.array([(a, b) for a, b in sorted(tset) if lam[a] and lam[b] and (b, a) not in tset])
for k, (tr, te) in enumerate(L.make_folds(c, "C3", 0)):
    trg = np.unique(np.r_[tf[tr], tg[tr]])
    if feat in ("bilinear", "tg_bilinear"):
        mu = E[trg].mean(0); _, _, Vt = np.linalg.svd(E[trg] - mu, full_matrices=False)
        P = (E - mu) @ Vt[:16].T
        B = np.einsum("pa,pb->pab", P[tf], P[tg]).reshape(tf.size, 256)
        Fm = B if feat == "bilinear" else np.hstack([E[tg], B])
    elif feat == "tg_only": Fm = E[tg]
    else: Fm = np.hstack([E[tf], E[tg]])
    m_ = Fm[tr].mean(0); s_ = Fm[tr].std(0); s_[s_ == 0] = 1
    Xtr = (Fm[tr] - m_) / s_; Xte = (Fm[te] - m_) / s_
    alpha = float(ck[(ck.fold == k)].alpha.values[0])
    clf = LogisticRegression(C=1.0 / (tr.size * alpha), max_iter=5000, tol=1e-8).fit(Xtr, y[tr])
    s = clf.decision_function(Xte)
    a = within_tf_auc(tf[te], y[te], s)
    rows.append(dict(context=ctx, feature=feat, layer=layer, fold=k, alpha=alpha, auc_sklearn=a,
                     auc_implementer=float(ck[ck.fold == k].auc.values[0])))
    print(rows[-1], flush=True)
    if unit == "D":
        ok = ~np.isin(DIR[:, 0], trg) & ~np.isin(DIR[:, 1], trg)
        a_, b_ = DIR[ok, 0], DIR[ok, 1]
        mt, st_ = m_[:512], s_[:512]; mg, sg = m_[512:], s_[512:]
        def score(i, j): return ((E[i] - mt) / st_) @ clf.coef_[0][:512] + ((E[j] - mg) / sg) @ clf.coef_[0][512:]
        dd = score(a_, b_) - score(b_, a_)
        for x, w, v in zip(a_, b_, dd): drows.append(dict(fold=k, i=int(x), j=int(w), diff_sklearn=float(v)))
R = pd.DataFrame(rows); R["absdiff"] = (R.auc_sklearn - R.auc_implementer).abs()
R.to_csv(HERE / f"v2_{unit}.csv", index=False)
print(unit, "mean sklearn", R.auc_sklearn.mean(), "mean implementer", R.auc_implementer.mean(), "max absdiff", R.absdiff.max())
if unit == "D":
    Dm = pd.DataFrame(drows); Di = pd.read_csv(E4C / "partial" / "m1" / ctx / "direction_L11_r0.csv")
    J = Dm.merge(Di, on=["fold", "i", "j"])
    J.to_csv(HERE / "v2_D_direction.csv", index=False)
    print("direction pairs", len(Dm), len(Di), "merged", len(J), "sign agreement", float((np.sign(J.diff_sklearn) == np.sign(J["diff"])).mean()),
          "frac fwd sklearn (per pair, r0)", float((J.groupby(['i', 'j']).diff_sklearn.mean() > 0).mean()),
          "implementer", float((J.groupby(['i', 'j'])["diff"].mean() > 0).mean()))
