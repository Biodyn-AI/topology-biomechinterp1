"""E4c reference scorers on the IDENTICAL test folds (added by the resuming agent; no training).

Unsupervised H3 scorers evaluated with the same fold metric (e4c_lib.fold_metric -> stats3.StratifiedPairs)
on every outer test fold used by the probes, so that "does a supervised probe beat the primary scorer on
held-out pairs?" is a paired, same-fold comparison:
  ref / scgpt  / cosine    / layer 0..11 : signed cosine of the V-centred full-512 embeddings (H3 primary scorer)
  ref / coexpr / abs_coexpr / -1         : |Pearson co-expression| on the extraction cells (expr_extraction.npy)
usage: python e4c_ref.py <context>     -> partial/ref/<context>.csv
"""
import sys
import numpy as np, pandas as pd
import e4c_lib as L

ctx = sys.argv[1]
c = L.Ctx(ctx)
scores = {}
for l in range(12):
    E = c.layer(l).astype(np.float64)
    En = E / np.linalg.norm(E, axis=1, keepdims=True)       # 512-d, no exact-tie issue (k > 1)
    scores[("scgpt", "cosine", l)] = np.einsum("ij,ij->i", En[c.tf], En[c.tg])
scores[("coexpr", "abs_coexpr", -1)] = c.coexpr_features()[:, 0].astype(np.float64)
rows = []
for reg, reps in (("C3", range(5)), ("C2", range(5)), ("C1", range(5))):
    for rep in reps:
        for k, (tr, te) in enumerate(L.make_folds(c, reg, rep)):
            for (src, feat, l), s in scores.items():
                auc, mt, ntf, npos, nneg = L.fold_metric(c.tf[te], c.tg[te], c.y[te], s[te])
                rows.append(dict(context=ctx, model="ref", source=src, feature=feat, layer=l, regime=reg, repeat=rep,
                                 fold=k, label_col=0, auc=auc, metric_type=mt, n_tf_used=ntf, n_pos=npos, n_neg=nneg,
                                 n_train=tr.size, n_test=te.size))
    L.log(f"{ctx} ref {reg} done")
pd.DataFrame(rows).to_csv(L.HERE / "partial" / "ref" / f"{ctx}.csv", index=False)
