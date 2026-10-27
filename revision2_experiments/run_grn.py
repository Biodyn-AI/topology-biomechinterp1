"""Held-out GRN benchmark with the standard baselines round 1 never ran.

Dask is used with a threaded scheduler and everything is inside a __main__
guard, so worker processes cannot re-import and re-execute the module.
"""
import sys, json, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C


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
    log(f"[setup] n={len(names)}  TRRUST(dedup)={len(tr)}")

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
                a, b = indptr[ci], indptr[ci+1]
                for j, v in zip(indices[a:b], data[a:b]):
                    if j in cset: E[r, cpos[j]] = v
        else:
            ncell = X.shape[0]
            sel = np.sort(rng.choice(ncell, min(3000, ncell), replace=False))
            E = np.asarray(X[sel,:][:,cols], dtype=np.float32)
    kloc = {g:k for k,g in enumerate(keep)}
    Z = (E - E.mean(0)) / (E.std(0) + 1e-9)
    def coexp(a,b):
        if a not in kloc or b not in kloc: return 0.0
        return float(np.dot(Z[:,kloc[a]], Z[:,kloc[b]]) / Z.shape[0])
    log(f"[setup] expression {E.shape}")

    pos_list = sorted(tr); perm = rng.permutation(len(pos_list))
    n_test = max(1, int(0.20*len(pos_list)))
    held_pos = [pos_list[i] for i in perm[:n_test]]
    train_pos = {pos_list[i] for i in perm[n_test:]}
    hp = set(held_pos)
    held_tfs = sorted({a for a,_ in held_pos})
    cand = [(tf,g) for tf in held_tfs for g in names
            if g != tf and (tf,g) not in train_pos and (tf,g) not in hp]
    held_neg = [cand[i] for i in rng.choice(len(cand), min(len(held_pos)*3, len(cand)), replace=False)]
    log(f"[split] held-out {len(held_pos)} positives / {len(held_neg)} negatives   (round 1: 57 / 124)")

    scorers = {}
    U0 = C.layer_us(emb, 0)[sub]
    full0 = np.asarray(emb[0], dtype=np.float64)[sub]
    def cos_scorer(M):
        Mn = M/(np.linalg.norm(M,axis=1,keepdims=True)+1e-12)
        return lambda a,b: float(Mn[loc[a]] @ Mn[loc[b]])
    scorers["scGPT full-512 (L0)"] = cos_scorer(full0)
    scorers["scGPT spectral SV5-SV7 (L0)"] = cos_scorer(U0[:,4:7])
    scorers["scGPT spectral SV2-SV7 (L0)"] = cos_scorer(U0[:,1:7])
    scorers["|co-expression|"] = lambda a,b: abs(coexp(a,b))
    scorers["co-expression (signed)"] = lambda a,b: coexp(a,b)
    scorers["random"] = lambda a,b: float(rng.random())

    try:
        from arboreto.algo import grnboost2, genie3
        df = pd.DataFrame(E, columns=keep)
        tfn = [g for g in keep if g in tfs_all]
        for nm, fn in [("GRNBoost2", grnboost2), ("GENIE3", genie3)]:
            log(f"[baseline] {nm}: {df.shape[1]} genes x {df.shape[0]} cells, {len(tfn)} regulators ...")
            net = fn(expression_data=df, tf_names=tfn, seed=42, verbose=False)
            d = {(r.TF, r.target): float(r.importance) for r in net.itertuples()}
            scorers[nm] = (lambda dd: (lambda a,b: dd.get((a,b), 0.0)))(d)
            log(f"           {len(d)} scored edges  ({time.time()-t0:.0f}s)")
    except Exception as e:
        log(f"[baseline] arboreto FAILED: {type(e).__name__}: {e}")

    def ev(fn):
        return (np.array([fn(a,b) for a,b in held_pos]),
                np.array([fn(a,b) for a,b in held_neg]))
    B = 2000; rows = []; cache = {}
    for nm, fn in scorers.items():
        p, q = ev(fn); cache[nm] = (p,q)
        a = C.auroc(p,q)
        bs = np.array([C.auroc(p[rng.integers(0,p.size,p.size)], q[rng.integers(0,q.size,q.size)])
                       for _ in range(B)])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        rows.append(dict(method=nm, auroc=a, ci_lo=lo, ci_hi=hi))
        log(f"    {nm:30s} AUROC={a:.3f}  95% CI [{lo:.3f}, {hi:.3f}]")

    log("\n[paired bootstrap deltas, 95% CI]")
    deltas = []
    ref = "|co-expression|"
    for nm in scorers:
        if nm == ref: continue
        p1,q1 = cache[nm]; p2,q2 = cache[ref]
        d = np.empty(B)
        for i in range(B):
            ip = rng.integers(0,p1.size,p1.size); iq = rng.integers(0,q1.size,q1.size)
            d[i] = C.auroc(p1[ip],q1[iq]) - C.auroc(p2[ip],q2[iq])
        lo,hi = np.percentile(d,[2.5,97.5])
        sep = "separates" if (lo>0 or hi<0) else "does NOT separate"
        deltas.append(dict(method=nm, delta=float(d.mean()), ci_lo=float(lo), ci_hi=float(hi), separates=(lo>0 or hi<0)))
        log(f"    {nm:30s} vs {ref}: {d.mean():+.3f}  [{lo:+.3f}, {hi:+.3f}]  {sep} from 0")

    pd.DataFrame(rows).to_csv(C.OUT/"GRN_benchmark_round2.csv", index=False)
    json.dump(dict(methods=rows, deltas=deltas,
                   n_pos=len(held_pos), n_neg=len(held_neg)),
              open(C.OUT/"GRN_benchmark_round2.json","w"), indent=1, default=float)
    log(f"\n[done] {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
