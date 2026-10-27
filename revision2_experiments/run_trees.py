"""GENIE3 and GRNBoost2 implemented directly (arboreto/dask failed on this stack).

GENIE3 (Huynh-Thu et al. 2010): for each target, regress its expression on all
candidate regulators with an ensemble of extremely randomised trees; the edge
weight is the regulator's feature importance. GRNBoost2 (Moerman et al. 2019) is
the same recipe with gradient-boosted stumps. Both are the field-standard
comparators the manuscript says it "did not run here".
"""
import sys, json, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor
from joblib import Parallel, delayed


def main():
    rng = np.random.default_rng(20260819); t0 = time.time()
    def log(*a): print(*a, flush=True)
    genes = C.h5ad_genes(); emb = C.load_emb(); nz0 = C.nonzero_mask(emb, 0)
    sub = np.array([i for i in C.curated_genes(genes) if nz0[i]])
    names = [genes[i] for i in sub]; loc = {g: k for k, g in enumerate(names)}
    t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF","target","mode","pmid"])
    tfs_all = set(t.TF.astype(str))
    tr = sorted({(a,b) for a,b in zip(t.TF.astype(str), t.target.astype(str))
                 if a in loc and b in loc and a != b})
    import h5py
    with h5py.File(C.H5AD, "r") as f:
        gs = [g.decode() if isinstance(g,bytes) else str(g) for g in f["var"]["_index"][:]]
        gi = {g:i for i,g in enumerate(gs)}
        keep = [g for g in names if g in gi]; cols = [gi[g] for g in keep]
        X = f["X"]
        if isinstance(X, h5py.Group):
            data, indices, indptr = X["data"][:], X["indices"][:], X["indptr"][:]
            ncell = indptr.shape[0]-1
            sel = np.sort(rng.choice(ncell, min(3000, ncell), replace=False))
            cpos = {c:k for k,c in enumerate(cols)}; cset = set(cols)
            E = np.zeros((sel.size, len(cols)), dtype=np.float32)
            for r, ci in enumerate(sel):
                a,b = indptr[ci], indptr[ci+1]
                for j,v in zip(indices[a:b], data[a:b]):
                    if j in cset: E[r, cpos[j]] = v
        else:
            ncell = X.shape[0]
            sel = np.sort(rng.choice(ncell, min(3000, ncell), replace=False))
            E = np.asarray(X[sel,:][:,cols], dtype=np.float32)
    kloc = {g:k for k,g in enumerate(keep)}
    reg_names = [g for g in keep if g in tfs_all]
    reg_cols = [kloc[g] for g in reg_names]
    log(f"[setup] expression {E.shape}; {len(reg_names)} regulators, {len(keep)} targets")

    def one_target(ti, model):
        y = E[:, ti]
        use = [c for c in reg_cols if c != ti]
        if y.std() < 1e-8: return ti, use, np.zeros(len(use))
        Xr = E[:, use]
        m = (ExtraTreesRegressor(n_estimators=200, max_features="sqrt", random_state=42, n_jobs=1)
             if model == "GENIE3" else
             GradientBoostingRegressor(n_estimators=200, max_depth=3, learning_rate=0.05,
                                       subsample=0.9, random_state=42))
        m.fit(Xr, y)
        imp = m.feature_importances_
        s = imp.sum()
        return ti, use, (imp / s if s > 0 else imp)

    nets = {}
    for model in ("GENIE3", "GRNBoost2"):
        log(f"[{model}] fitting {len(keep)} targets ...")
        out = Parallel(n_jobs=8, backend="loky")(
            delayed(one_target)(ti, model) for ti in range(len(keep)))
        d = {}
        for ti, use, imp in out:
            tgt = keep[ti]
            for c, w in zip(use, imp):
                d[(keep[c], tgt)] = float(w)
        nets[model] = d
        log(f"[{model}] {len(d)} scored edges  ({time.time()-t0:.0f}s)")

    pos_list = sorted(tr); perm = rng.permutation(len(pos_list))
    n_test = max(1, int(0.20*len(pos_list)))
    held_pos = [pos_list[i] for i in perm[:n_test]]
    train_pos = {pos_list[i] for i in perm[n_test:]}; hp = set(held_pos)
    held_tfs = sorted({a for a,_ in held_pos})
    cand = [(tf,g) for tf in held_tfs for g in names
            if g != tf and (tf,g) not in train_pos and (tf,g) not in hp]
    held_neg = [cand[i] for i in rng.choice(len(cand), min(len(held_pos)*3, len(cand)), replace=False)]
    log(f"[split] {len(held_pos)} positives / {len(held_neg)} negatives")

    Zs = (E - E.mean(0)) / (E.std(0) + 1e-9)
    def coexp(a,b):
        if a not in kloc or b not in kloc: return 0.0
        return float(np.dot(Zs[:,kloc[a]], Zs[:,kloc[b]]) / Zs.shape[0])
    U0 = C.layer_us(emb, 0)[sub]; full0 = np.asarray(emb[0], dtype=np.float64)[sub]
    def cs(M):
        Mn = M/(np.linalg.norm(M,axis=1,keepdims=True)+1e-12)
        return lambda a,b: float(Mn[loc[a]] @ Mn[loc[b]])
    scorers = {
        "GRNBoost2": lambda a,b: nets["GRNBoost2"].get((a,b), 0.0),
        "GENIE3":    lambda a,b: nets["GENIE3"].get((a,b), 0.0),
        "scGPT full-512 (L0)": cs(full0),
        "scGPT spectral SV2-SV7 (L0)": cs(U0[:,1:7]),
        "scGPT spectral SV5-SV7 (L0)": cs(U0[:,4:7]),
        "|co-expression|": lambda a,b: abs(coexp(a,b)),
        "random": lambda a,b: float(rng.random()),
    }
    B = 2000; rows=[]; cache={}
    for nm, fn in scorers.items():
        p = np.array([fn(a,b) for a,b in held_pos]); q = np.array([fn(a,b) for a,b in held_neg])
        cache[nm]=(p,q); a = C.auroc(p,q)
        bs = np.array([C.auroc(p[rng.integers(0,p.size,p.size)], q[rng.integers(0,q.size,q.size)]) for _ in range(B)])
        lo,hi = np.percentile(bs,[2.5,97.5])
        rows.append(dict(method=nm, auroc=a, ci_lo=lo, ci_hi=hi))
        log(f"    {nm:30s} AUROC={a:.3f}  95% CI [{lo:.3f}, {hi:.3f}]")
    log("\n[paired bootstrap deltas vs |co-expression|]")
    deltas=[]
    p2,q2 = cache["|co-expression|"]
    for nm in scorers:
        if nm=="|co-expression|": continue
        p1,q1 = cache[nm]; d=np.empty(B)
        for i in range(B):
            ip=rng.integers(0,p1.size,p1.size); iq=rng.integers(0,q1.size,q1.size)
            d[i]=C.auroc(p1[ip],q1[iq])-C.auroc(p2[ip],q2[iq])
        lo,hi=np.percentile(d,[2.5,97.5]); sep = lo>0 or hi<0
        deltas.append(dict(method=nm, delta=float(d.mean()), ci_lo=float(lo), ci_hi=float(hi), separates=bool(sep)))
        log(f"    {nm:30s} {d.mean():+.3f}  [{lo:+.3f}, {hi:+.3f}]  {'separates' if sep else 'does NOT separate'} from 0")
    pd.DataFrame(rows).to_csv(C.OUT/"GRN_benchmark_final.csv", index=False)
    json.dump(dict(methods=rows, deltas=deltas, n_pos=len(held_pos), n_neg=len(held_neg)),
              open(C.OUT/"GRN_benchmark_final.json","w"), indent=1, default=float)
    log(f"\n[done] {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
