"""E10 step 6 (task 5): transfer of a linear TF-vs-non-TF read-out across contexts in the shared 512-d space.

Train: immune genes (V_immune), Lambert TF vs non-TF, per layer; features = unit-normalised V-centred embedding
(the cosine geometry used by S2); sklearn L2 logistic regression, class_weight='balanced', C chosen by 5-fold
stratified CV AUROC on immune from {0.01, 0.1, 1, 10} (seeded). Reference: immune 5-fold out-of-fold AUROC.
Apply the immune-trained model to kidney_B, lung, lung_external (and kidney_A as an extra) on
  novel  = target-context genes NOT observed in immune (never seen in training), and
  shared = genes observed in both (their immune embedding was in training; the target embedding was not).
Report plain AUROC and detection-stratified AUROC (target det_decile), gene-bootstrap 95% CI (B = 2,000, fixed
scores), and the in-context reference (5-fold out-of-fold within the target context, same layer, same C grid)
evaluated on the same gene subsets (the in-context model uses the C chosen on immune). Sensitivity: un-centred embeddings (no per-context centring).
Output: transfer.csv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e10_common as K
S = K.S
TARGETS = ["kidney_B", "lung", "lung_external", "kidney_A"]
CS = [0.01, 0.1, 1.0, 10.0]
B = 2000


def feats(emb, L, obs_idx, centre=True):
    Xm = np.asarray(emb[L], dtype=np.float64)[obs_idx]
    if centre: Xm = Xm - Xm.mean(0, keepdims=True)
    return K.unit(Xm)


def oof(Xf, y, C, seed):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed); s = np.zeros(y.size)
    for tr, te in skf.split(Xf, y):
        m = LogisticRegression(C=C, class_weight="balanced", max_iter=3000).fit(Xf[tr], y[tr])
        s[te] = m.decision_function(Xf[te])
    return s


def pick_C(Xf, y, seed):
    sc = {C: roc_auc_score(y, oof(Xf, y, C, seed)) for C in CS}
    return max(sc, key=sc.get), sc


def main():
    rows = []
    ci = K.load_cache("immune"); ei = K.load_emb("immune"); yi = ci["nl"]["is_tf_lambert"].astype(bool)
    imm = set(ci["genes"])
    tg = {t: (K.load_cache(t), K.load_emb(t)) for t in TARGETS}
    for centre in (True, False):
        for L in range(K.N_LAYERS):
            Xi = feats(ei, L, ci["obs_idx"], centre)
            C, grid = pick_C(Xi, yi, K.SEED + L)
            s_oof = oof(Xi, yi, C, K.SEED + 100 + L)
            rows.append(dict(centred=centre, layer=L, train="immune", test="immune", subset="oof_5fold", C=C,
                             n_pos=int(yi.sum()), n_neg=int((~yi).sum()), auroc=S.auroc(s_oof[yi], s_oof[~yi]),
                             auroc_detstrat=S.stratified_auroc(s_oof, yi, ci["nl"]["det_decile"])))
            model = LogisticRegression(C=C, class_weight="balanced", max_iter=3000).fit(Xi, yi)
            for t, (c, e) in tg.items():
                Xt = feats(e, L, c["obs_idx"], centre); y = c["nl"]["is_tf_lambert"].astype(bool); st = c["nl"]["det_decile"]
                s = model.decision_function(Xt)
                s_in = oof(Xt, y, C, K.SEED + 200 + L) if centre else None
                novel = np.array([g not in imm for g in c["genes"]])
                for nm, msk in (("novel", novel), ("shared", ~novel), ("all", np.ones(novel.size, bool))):
                    yy = y[msk]
                    if yy.sum() < 5 or (~yy).sum() < 5: continue
                    rng = np.random.default_rng(K.SEED + 300 + L)
                    m = msk.sum(); W = rng.multinomial(m, np.full(m, 1 / m), size=B).T.astype(float)
                    bt = K.wauroc_fixed(s[msk], yy, W)
                    r = dict(centred=centre, layer=L, train="immune", test=t, subset=nm, C=C, n_pos=int(yy.sum()),
                             n_neg=int((~yy).sum()), auroc=S.auroc(s[msk][yy], s[msk][~yy]),
                             lo=float(np.quantile(bt, .025)), hi=float(np.quantile(bt, .975)),
                             auroc_detstrat=S.stratified_auroc(s[msk], yy, st[msk]))
                    if s_in is not None:
                        bi = K.wauroc_fixed(s_in[msk], yy, W)
                        r.update(incontext_auroc=S.auroc(s_in[msk][yy], s_in[msk][~yy]),
                                 incontext_auroc_detstrat=S.stratified_auroc(s_in[msk], yy, st[msk]),
                                 transfer_minus_incontext=r["auroc"] - S.auroc(s_in[msk][yy], s_in[msk][~yy]),
                                 diff_lo=float(np.quantile(bt - bi, .025)), diff_hi=float(np.quantile(bt - bi, .975)))
                    rows.append(r)
            print(centre, L, C, [(r["test"], r["subset"], round(r["auroc"], 3)) for r in rows[-12:] if r["layer"] == L], flush=True)
    df = pd.DataFrame(rows)
    # layer-mean summary rows
    g = df.groupby(["centred", "train", "test", "subset"], as_index=False)[[c for c in ["auroc", "auroc_detstrat", "incontext_auroc",
                                                                             "incontext_auroc_detstrat", "transfer_minus_incontext"] if c in df]].mean()
    g["layer"] = "mean"
    pd.concat([df, g], ignore_index=True).to_csv(K.HERE / "transfer.csv", index=False)


if __name__ == "__main__":
    main()
