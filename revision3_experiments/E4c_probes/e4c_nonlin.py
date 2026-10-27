"""E4c M2 (HistGradientBoosting) / M3 (sklearn MLP, 2 hidden layers) edge probes on layers L0, L5, L11.

usage: python e4c_nonlin.py <context> <hgb|mlp> --phase true|null [--features concat,hadamard]
Budget-reduced design (see SUMMARY.md):
  true : regimes C3 and C2, repeat 0 (5 folds each), true labels; probes (--features), target-only and
         co-expression controls on the same folds. TF-only is not fitted: under the within-TF metric any
         TF-only score is constant within a TF and its AUROC is 0.5 by construction (pooled-fallback folds are
         flagged in the aggregation).
  null : C3 repeat 0, rewired draws 0-9 on identical folds (co-expression first, then probes). Hyperparameters are fixed a priori; the only in-fold tuning is early stopping on an
internal validation split of the training fold (n boosting iterations / n epochs), identical for true and
rewired labels. One checkpoint file per (source, feature, layer, regime, fold, label column).
"""
import sys, time, json
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
import e4c_lib as L

ctx_name, model = sys.argv[1], sys.argv[2]
PHASE = sys.argv[sys.argv.index("--phase") + 1]
FEATS = (sys.argv[sys.argv.index("--features") + 1].split(",") if "--features" in sys.argv else ["concat", "hadamard"])
LAYERS = [0, 5, 11]
N_NULL = 10
OUT = L.HERE / "partial" / model / ctx_name; OUT.mkdir(parents=True, exist_ok=True)
LOG = L.HERE / "logs" / f"{model}_{ctx_name}_{PHASE}.log"
c = L.Ctx(ctx_name)
z = np.load(L.HERE / "rewired" / f"{ctx_name}.npz")
Yr = np.unpackbits(z["Y"], axis=1, count=int(z["n"])).astype(bool)[:N_NULL]
YALL = np.column_stack([c.y, Yr.T])
MNAME = {"hgb": "M2", "mlp": "M3"}[model]


def make_model(seed):
    if model == "hgb":
        return HistGradientBoostingClassifier(max_iter=150, learning_rate=0.15, max_leaf_nodes=31, min_samples_leaf=40,
                                              l2_regularization=1.0, max_bins=63, early_stopping=True,
                                              validation_fraction=0.15, n_iter_no_change=8, random_state=seed)
    return MLPClassifier(hidden_layer_sizes=(128, 32), alpha=1e-3, batch_size=512, learning_rate_init=1e-3,
                         early_stopping=True, validation_fraction=0.1, n_iter_no_change=5, max_iter=60,
                         random_state=seed)


def run(source, feature, layer, regime, cols, Fall):
    folds = L.make_folds(c, regime, 0)
    for k, (tr, te) in enumerate(folds):
        for j in cols:
            fn = OUT / f"{source}_{feature}_L{layer}_{regime}_r0_f{k}_c{j}.json"
            if fn.exists(): continue
            t0 = time.time()
            F = Fall(tr) if callable(Fall) else Fall
            Xtr, Xte = F[tr], F[te]
            if model == "mlp":
                mu = Xtr.mean(0); sd = Xtr.std(0); sd[sd == 0] = 1
                Xtr = ((Xtr - mu) / sd).astype(np.float32); Xte = ((Xte - mu) / sd).astype(np.float32)
            y = YALL[:, j]
            m = make_model(seed=L.SEED % 100000 + 7 * k + 1000 * j).fit(Xtr, y[tr])
            s = m.predict_proba(Xte)[:, 1].astype(np.float64)
            auc, mt, ntf, npos, nneg = L.fold_metric(c.tf[te], c.tg[te], y[te], s)
            row = dict(context=ctx_name, model=MNAME, source=source, feature=feature, layer=layer, regime=regime,
                       repeat=0, fold=k, label_col=j, auc=auc, metric_type=mt, n_tf_used=ntf, n_pos=npos, n_neg=nneg,
                       n_train=int(tr.size), n_test=int(te.size), n_iter=int(m.n_iter_), sec=time.time() - t0)
            fn.write_text(json.dumps(row))
        L.log(f"{ctx_name} {MNAME} {source} {feature} L{layer} {regime} fold{k} cols{cols[0]}-{cols[-1]} done", LOG)


L.log(f"start {ctx_name} {MNAME} features={FEATS}", LOG)
Fco = c.coexpr_features()
TRUE, NULL = [0], list(range(1, N_NULL + 1))
Es = {l: c.layer(l) for l in LAYERS}
def feat_fn(l, f):
    if f == "bilinear":
        return lambda tr: L.pair_features("bilinear", Es[l], c.tf, c.tg, train_genes=np.unique(np.r_[c.tf[tr], c.tg[tr]]))
    return L.pair_features(f, Es[l], c.tf, c.tg)
if PHASE == "true":
    for reg in ("C3", "C2"):
        for l in LAYERS:
            for f in FEATS: run("scgpt", f, l, reg, TRUE, feat_fn(l, f))
            run("scgpt", "tg_only", l, reg, TRUE, L.pair_features("tg_only", Es[l], c.tf, c.tg))
        run("coexpr", "coexpr", -1, reg, TRUE, Fco)
else:
    run("coexpr", "coexpr", -1, "C3", NULL, Fco)
    for f in FEATS:
        for l in LAYERS: run("scgpt", f, l, "C3", NULL, feat_fn(l, f))
L.log(f"done {ctx_name} {MNAME}", LOG)
