"""E4c M2/M3 with in-fold hyperparameter tuning (small grids) — added by the resuming agent.

The pre-specified design asks for M2/M3 hyperparameters "tuned ONLY inside training folds (small grids)".
e4c_nonlin.py (previous agent) fixed them a priori (only early stopping is in-fold). This script adds the
tuned version on the SAME outer folds (C3 and C2, repeat 0; layers L0, L5, L11; features concat, hadamard,
the target-only control and the co-expression control), true labels only:
  for every outer training fold: regime-matched inner split (e4c_lib.inner_split, as for M1); each grid
  configuration fitted on the inner-training part, scored on the inner-validation part with the within-TF
  AUROC rule (e4c_lib.fast_within_tf_auroc; checked against stats3.StratifiedPairs); the best configuration
  is refitted on the whole outer training fold and scored on the outer test fold (e4c_lib.fold_metric).
Grids: HGB {lr 0.05, leaves 15, min_leaf 80} / {0.15, 31, 40} (= the fixed default) / {0.15, 63, 20};
       MLP alpha {1e-4, 1e-2} x hidden {(128, 32), (256, 64)}.
Rewired-label nulls are NOT re-run with tuning (budget); the tuned-vs-fixed comparison on true labels shows
whether the fixed-hyperparameter null comparison is affected by under-tuning.
usage: python e4c_nonlin_tuned.py <context> <hgb|mlp>   -> partial/<model>_tuned/<ctx>/*.json
"""
import sys, time, json
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
import e4c_lib as L

ctx_name, model = sys.argv[1], sys.argv[2]
LAYERS = [0, 5, 11]
OUT = L.HERE / "partial" / f"{model}_tuned" / ctx_name; OUT.mkdir(parents=True, exist_ok=True)
LOG = L.HERE / "logs" / f"{model}_tuned_{ctx_name}.log"
c = L.Ctx(ctx_name)
MNAME = {"hgb": "M2_tuned", "mlp": "M3_tuned"}[model]
GRID = ([dict(learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=80),
         dict(learning_rate=0.15, max_leaf_nodes=31, min_samples_leaf=40),
         dict(learning_rate=0.15, max_leaf_nodes=63, min_samples_leaf=20)] if model == "hgb" else
        [dict(alpha=a, hidden_layer_sizes=h) for a in (1e-4, 1e-2) for h in ((128, 32), (256, 64))])


def make_model(cfg, seed):
    if model == "hgb":
        return HistGradientBoostingClassifier(max_iter=150, l2_regularization=1.0, max_bins=63, early_stopping=True,
                                              validation_fraction=0.15, n_iter_no_change=8, random_state=seed, **cfg)
    return MLPClassifier(batch_size=512, learning_rate_init=1e-3, early_stopping=True, validation_fraction=0.1,
                         n_iter_no_change=5, max_iter=60, random_state=seed, **cfg)


def prep(Xa, Xb):
    if model != "mlp": return Xa, Xb
    mu = Xa.mean(0); sd = Xa.std(0); sd[sd == 0] = 1
    return ((Xa - mu) / sd).astype(np.float32), ((Xb - mu) / sd).astype(np.float32)


def run(source, feature, layer, regime, Fall):
    folds = L.make_folds(c, regime, 0)
    for k, (tr, te) in enumerate(folds):
        fn = OUT / f"{source}_{feature}_L{layer}_{regime}_r0_f{k}_c0.json"
        if fn.exists(): continue
        t0 = time.time()
        F = Fall(tr) if callable(Fall) else Fall
        y = c.y
        seed = L.SEED % 100000 + 7 * k
        itr, iva = L.inner_split(c, regime, tr, seed=k + 17 * {"C2": 2, "C3": 3}[regime])
        a_, b_ = prep(F[tr][itr], F[tr][iva]); sc = []
        for cfg in GRID:
            m = make_model(cfg, seed).fit(a_, y[tr][itr])
            sc.append(m.predict_proba(b_)[:, 1])
        Sv = np.column_stack(sc); Yv = np.tile(y[tr][iva][:, None], (1, len(GRID)))
        fast = L.fast_within_tf_auroc(c.tf[tr][iva], Yv, Sv); L.check_fast(c.tf[tr][iva], c.tg[tr][iva], Yv, Sv, fast)
        best = int(np.nanargmax(np.nan_to_num(fast, nan=-1)))
        Xtr, Xte = prep(F[tr], F[te])
        m = make_model(GRID[best], seed).fit(Xtr, y[tr])
        s = m.predict_proba(Xte)[:, 1].astype(np.float64)
        auc, mt, ntf, npos, nneg = L.fold_metric(c.tf[te], c.tg[te], y[te], s)
        row = dict(context=ctx_name, model=MNAME, source=source, feature=feature, layer=layer, regime=regime,
                   repeat=0, fold=k, label_col=0, auc=auc, metric_type=mt, n_tf_used=ntf, n_pos=npos, n_neg=nneg,
                   n_train=int(tr.size), n_test=int(te.size), n_iter=int(m.n_iter_), best_cfg=best,
                   inner_auc=[float(x) for x in fast], sec=time.time() - t0)
        fn.write_text(json.dumps(row))
    L.log(f"{ctx_name} {MNAME} {source} {feature} L{layer} {regime} done", LOG)


L.log(f"start {ctx_name} {MNAME}", LOG)
Fco = c.coexpr_features()
for reg in ("C3", "C2"):
    for l in LAYERS:
        E = c.layer(l)
        for f in ("concat", "hadamard", "tg_only"): run("scgpt", f, l, reg, L.pair_features(f, E, c.tf, c.tg))
    run("coexpr", "coexpr", -1, reg, Fco)
L.log(f"done {ctx_name} {MNAME}", LOG)
