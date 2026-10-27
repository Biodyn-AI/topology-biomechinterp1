"""E4c M1: L2 logistic edge probes, all 12 layers, regimes C1/C2/C3 (5 folds x repeats), true labels plus 100
rewired-label refits on identical folds, node-only and external controls, ordered direction test (C3 concat).

usage: python e4c_m1.py <context> --phase true|null|null_c2|c1extra
  true    : true labels only; C3 x 5 repeats, C2 x 5 repeats, C1 x repeat 0 (all sources, 12 layers)
  null    : true + 100 rewired label columns, C3 repeat 0 (identical folds), all sources, 12 layers
  null_c2 : the same for C2 repeat 0 (run only if the compute budget allows)
  c1extra : true labels, C1 repeats 1-4 (only if budget allows)
  combo   : sensitivity beyond the pre-specified grid: true labels, C3 x 5 repeats + C2 repeat 0, features
            [x_j, x_i*x_j] (tg_hadamard) and [x_j, vec(p_i p_j^T)] (tg_bilinear): does an interaction block add
            to the target-only block on identical folds?
  ext_node: true labels, node-only controls for the external sources (ESM2 / static TF-only, target-only) and
            the co-expression marginal-only feature set (mean, detection of each gene), regimes as 'true' 
Checkpoints: partial/m1/<ctx>/<source>_<feature>_L<layer>_<regime>_r<repeat>[_null].csv (resumable).
"""
import sys, time, json
import numpy as np, pandas as pd
import e4c_lib as L

ctx_name = sys.argv[1]
phase = sys.argv[sys.argv.index("--phase") + 1]
OUT = L.HERE / "partial" / "m1" / ctx_name; OUT.mkdir(parents=True, exist_ok=True)
LOG = L.HERE / "logs" / f"m1_{ctx_name}_{phase}.log"; LOG.parent.mkdir(exist_ok=True)
c = L.Ctx(ctx_name)
NODE_ONLY = ctx_name == "immune"                    # immune: sanity check of node-only baselines only
if NODE_ONLY:
    YALL = c.y[:, None]
else:
    z = np.load(L.HERE / "rewired" / f"{ctx_name}.npz")
    Yr = np.unpackbits(z["Y"], axis=1, count=int(z["n"])).astype(bool)
    YALL = np.column_stack([c.y, Yr.T])            # col 0 = true labels, 1..100 = rewired
UNITS_REG = {"true": [("C3", r) for r in range(5)] + [("C2", r) for r in range(5)] + [("C1", 0)],
             "null": [("C3", 0)], "null_c2": [("C2", 0)], "c1extra": [("C1", r) for r in range(1, 5)],
             "combo": [("C3", r) for r in range(5)] + [("C2", 0)],
             "ext_node": [("C3", r) for r in range(5)] + [("C2", r) for r in range(5)] + [("C1", 0)]}[phase]
WITH_NULL = phase in ("null", "null_c2")
if NODE_ONLY: assert phase in ("true", "c1extra")

# ---------------------------------------------------------------- direction-test pairs (C3 concat only)
tr_edges = np.load(L.CACHE / ctx_name / "edges_trrust.npy")
tset = {(int(a), int(b)) for a, b in tr_edges}
lam = c.is_tf_lambert
DIR = np.array([(a, b) for a, b in sorted(tset) if lam[a] and lam[b] and (b, a) not in tset], dtype=np.int64).reshape(-1, 2)


def run_unit(source, feature, layer, regime, rep, Egene, Fpair, extra, ncols):
    """One (source, feature, layer, regime, repeat): 5 outer folds."""
    fn = OUT / f"{source}_{feature}_L{layer}_{regime}_r{rep}{'_null' if ncols > 1 else ''}.csv"
    if fn.exists(): return
    t0 = time.time(); rows = []; drows = []
    folds = L.make_folds(c, regime, rep)
    Y = YALL[:, :ncols]
    Egene_arg = Egene
    for k, (tr, te) in enumerate(folds):
        trgenes = np.unique(np.r_[c.tf[tr], c.tg[tr]])
        Egene = Egene_arg(trgenes) if callable(Egene_arg) else Egene_arg
        # ---- design matrices (standardisation from training pairs only)
        if feature in ("concat", "tf_only", "tg_only"):
            blocks, idxs = [], []
            if feature in ("concat", "tf_only"): blocks.append(L.gene_block(Egene, c.tf[tr])); idxs.append(c.tf)
            if feature in ("concat", "tg_only"): blocks.append(L.gene_block(Egene, c.tg[tr])); idxs.append(c.tg)
            if extra is not None:
                if feature in ("concat", "tf_only"): blocks.append(L.gene_block(extra, c.tf[tr])); idxs.append(c.tf)
                if feature in ("concat", "tg_only"): blocks.append(L.gene_block(extra, c.tg[tr])); idxs.append(c.tg)
            Dall = L.Gather(blocks, idxs, c.G)
            Dtr, Dte = Dall.rows(tr), Dall.rows(te)
        else:
            if Fpair is not None: F = Fpair
            else: F = L.pair_features(feature, Egene, c.tf, c.tg, train_genes=trgenes, extra=extra)
            mu = F[tr].mean(0); sd = F[tr].std(0); sd[sd == 0] = 1
            Dtr = L.Dense(((F[tr] - mu) / sd).astype(np.float32)); Dte = L.Dense(((F[te] - mu) / sd).astype(np.float32))
        Ytr = Y[tr]
        # ---- inner, regime-matched split: choose alpha per label column
        itr, iva = L.inner_split(c, regime, tr, seed=1000 * rep + k + 17 * {"C1": 1, "C2": 2, "C3": 3}[regime])
        nA = L.ALPHAS.size
        Wi, bi, _ = L.fit_logistic_batched(Dtr.rows(itr), np.tile(Ytr[itr], (1, nA)), np.repeat(L.ALPHAS, ncols),
                                           max_iter=15, tol=1e-5)
        Sv = L.predict_linear(Dtr.rows(iva), Wi, bi)
        Yv = np.tile(Ytr[iva], (1, nA))
        fast = L.fast_within_tf_auroc(c.tf[tr][iva], Yv, Sv)
        L.check_fast(c.tf[tr][iva], c.tg[tr][iva], Yv, Sv, fast)
        fa = np.nan_to_num(fast.reshape(nA, ncols), nan=-1)
        best = L.ALPHAS[fa.argmax(0)]
        # ---- final fit on the whole training fold, score the test fold
        W, b, nit = L.fit_logistic_batched(Dtr, Ytr, best, max_iter=40, tol=1e-7)
        St = L.predict_linear(Dte, W, b)
        for j in range(ncols):
            auc, mt, ntf, npos, nneg = L.fold_metric(c.tf[te], c.tg[te], Y[te, j], St[:, j].astype(np.float64))
            rows.append(dict(context=ctx_name, model="M1", source=source, feature=feature, layer=layer, regime=regime,
                             repeat=rep, fold=k, label_col=j, auc=auc, metric_type=mt, n_tf_used=ntf, n_pos=npos,
                             n_neg=nneg, alpha=best[j], n_train=tr.size, n_test=te.size, n_iter=nit))
        # ---- ordered direction test: TRRUST i->j (both Lambert TFs, j->i undocumented), both genes unseen
        if source == "scgpt" and feature == "concat" and regime == "C3" and DIR.size:
            ok = ~np.isin(DIR[:, 0], trgenes) & ~np.isin(DIR[:, 1], trgenes)
            if ok.any():
                a, bb = DIR[ok, 0], DIR[ok, 1]
                fwd = Dall.__class__(Dall.B, [a, bb], c.G); rev = Dall.__class__(Dall.B, [bb, a], c.G)
                d = (L.predict_linear(fwd, W[:, :1], b[:1]) - L.predict_linear(rev, W[:, :1], b[:1]))[:, 0]
                for x, y_, dd in zip(a, bb, d):
                    drows.append(dict(context=ctx_name, layer=layer, repeat=rep, fold=k, i=int(x), j=int(y_),
                                      gene_i=c.genes[x], gene_j=c.genes[y_], diff=float(dd)))
    pd.DataFrame(rows).to_csv(fn.with_suffix(".tmp"), index=False); fn.with_suffix(".tmp").rename(fn)
    if drows and ncols == 1: pd.DataFrame(drows).to_csv(OUT / f"direction_L{layer}_r{rep}.csv", index=False)
    r0 = pd.DataFrame(rows).query("label_col == 0").auc.mean()
    L.log(f"{ctx_name} {source} {feature} L{layer} {regime} r{rep}: true-label mean AUROC {r0:.4f} "
          f"({time.time() - t0:.0f}s)", LOG)


L.log(f"start {ctx_name} phase={phase} n_pairs={c.n} n_pos={int(c.y.sum())} label_cols={YALL.shape[1]} "
      f"direction_pairs={len(DIR)}", LOG)
NC = YALL.shape[1] if WITH_NULL else 1
if phase == "combo":
    for layer in range(12):
        E = c.layer(layer)
        for reg, rep in UNITS_REG:
            for feat in ("tg_hadamard", "tg_bilinear"):
                run_unit("scgpt", feat, layer, reg, rep, E, None, None, 1)
    L.log(f"done {ctx_name} phase={phase}", LOG); sys.exit(0)
if phase == "ext_node":
    Fco = c.coexpr_features()
    for reg, rep in UNITS_REG:
        run_unit("coexpr", "coexpr_marg", -1, reg, rep, None, np.ascontiguousarray(Fco[:, 2:]), None, 1)
    for src in ("static", "esm2"):
        Mx, have = c.external_matrix(src); h3g = L.h3_genes(c)
        extra = None if have[h3g].all() else (~have).astype(np.float32)[:, None]
        if src == "static":
            E = np.where(have[:, None], np.nan_to_num(Mx), 0).astype(np.float64)
            E[have] -= E[have].mean(0); Eg = E.astype(np.float32)
        else:
            M0 = np.nan_to_num(Mx).astype(np.float64)
            _cache = {}
            def Eg(trgenes, M0=M0, have=have, _cache=_cache):
                key = trgenes.tobytes()
                if key not in _cache:                 # same training genes -> same PCA (cached across features)
                    if len(_cache) > 12: _cache.clear()
                    g = trgenes[have[trgenes]]
                    mu, C = L.pca_basis(M0, g, 128)
                    E = np.zeros((c.G, 128), np.float32); E[have] = ((M0[have] - mu) @ C.T).astype(np.float32)
                    _cache[key] = E
                return _cache[key]
        for reg, rep in UNITS_REG:
            for feat in ("tf_only", "tg_only"):
                run_unit(src, feat, -1, reg, rep, Eg, None, extra, 1)
    L.log(f"done {ctx_name} phase={phase}", LOG); sys.exit(0)
# ---------------------------------------------------------------- layer-free controls first (cheap)
if not NODE_ONLY:
    Fco = c.coexpr_features()
    for reg, rep in UNITS_REG:
        run_unit("coexpr", "coexpr", -1, reg, rep, None, Fco, None, NC)
    for src in ("static", "esm2"):
        Mx, have = c.external_matrix(src)
        h3g = L.h3_genes(c)
        L.log(f"{src}: coverage of h3 genes {have[h3g].mean():.4f} ({(~have[h3g]).sum()} missing -> zero vector + indicator)", LOG)
        extra = None if have[h3g].all() else (~have).astype(np.float32)[:, None]
        if src == "static":                       # full 512-d static token embedding, centred over V
            E = np.where(have[:, None], np.nan_to_num(Mx), 0).astype(np.float64)
            E[have] -= E[have].mean(0); Eg = E.astype(np.float32)
        else:                                     # ESM2: PCA-128 fitted on each fold's training genes
            M0 = np.nan_to_num(Mx).astype(np.float64)
            _cache = {}
            def Eg(trgenes, M0=M0, have=have, _cache=_cache):
                key = trgenes.tobytes()
                if key not in _cache:                 # same training genes -> same PCA (cached across features)
                    if len(_cache) > 12: _cache.clear()
                    g = trgenes[have[trgenes]]
                    mu, C = L.pca_basis(M0, g, 128)
                    E = np.zeros((c.G, 128), np.float32); E[have] = ((M0[have] - mu) @ C.T).astype(np.float32)
                    _cache[key] = E
                return _cache[key]
        for reg, rep in UNITS_REG:
            for feat in L.FEATURES:
                run_unit(src, feat, -1, reg, rep, Eg, None, extra, NC)
        del Mx
# ---------------------------------------------------------------- scGPT layers
for layer in range(12):
    E = c.layer(layer)
    dense = {} if NODE_ONLY else {f: L.pair_features(f, E, c.tf, c.tg) for f in ("hadamard", "absdiff")}
    for reg, rep in UNITS_REG:
        if not WITH_NULL:
            for feat in ("tf_only", "tg_only"):
                run_unit("scgpt", feat, layer, reg, rep, E, None, None, 1)
        if NODE_ONLY: continue
        for feat in L.FEATURES:
            run_unit("scgpt", feat, layer, reg, rep, E, dense.get(feat), None, NC)
    del dense
L.log(f"done {ctx_name} phase={phase}", LOG)
