"""Complete the SV5-SV7 confound residualisation across all 12 layers.

The earlier controls run was killed at layer 7, leaving 8 of 12 layers, and the
figure silently fell back to a stale CSV from an aborted run. This finishes it.
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
    tfs = set(t.TF.astype(str))
    tr = sorted({(a,b) for a,b in zip(t.TF.astype(str), t.target.astype(str))
                 if a in loc and b in loc and a != b})
    obj = json.loads(C.STRING_CACHE.read_text())
    st_all = [(p["g1"],p["g2"],p["score"]) for p in obj["pairs"]]
    gene2go = C.load_gene2go()

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
    Z = (E - E.mean(0)) / (E.std(0) + 1e-9)
    def coexp(a,b):
        if a not in kloc or b not in kloc: return 0.0
        return float(np.dot(Z[:,kloc[a]], Z[:,kloc[b]]) / Z.shape[0])
    deg = {g:0 for g in names}
    for g1,g2,s in st_all:
        if s>=0.7:
            if g1 in deg: deg[g1]+=1
            if g2 in deg: deg[g2]+=1
    ngo = {g: len(gene2go.get(g,set())) for g in names}
    mexp = {g: float(E[:,kloc[g]].mean()) if g in kloc else 0.0 for g in names}
    neg = []
    trset = set(tr)
    for a in {x for x,_ in tr}:
        for b in names:
            if b != a and (a,b) not in trset: neg.append((a,b))
    neg = [neg[i] for i in rng.choice(len(neg), min(4*len(tr), len(neg)), replace=False)]
    allp = tr + neg
    y = np.array([1]*len(tr) + [0]*len(neg))
    Xc = np.column_stack([
        [coexp(a,b) for a,b in allp],
        [np.log1p(mexp[a]) for a,b in allp], [np.log1p(mexp[b]) for a,b in allp],
        [deg[a]+deg[b] for a,b in allp], [ngo[a]+ngo[b] for a,b in allp]]).astype(np.float64)
    X1 = np.column_stack([np.ones(len(allp)), Xc])
    def rb(x, lab):
        r = pd.Series(x).rank().to_numpy()
        return float(2*(r[lab==1].mean() - r[lab==0].mean())/lab.size)
    log(f"[setup] {len(names)} genes, {len(tr)} positive pairs, {len(neg)} negatives")

    for tag, cols_ in [("SV5_SV7", slice(4,7)), ("SV2_SV7", slice(1,7))]:
        rows = []
        for L in range(12):
            U = C.layer_us(emb, L)[sub][:, cols_]
            Un = U/(np.linalg.norm(U,axis=1,keepdims=True)+1e-12)
            cos = np.array([float(Un[loc[a]] @ Un[loc[b]]) for a,b in allp])
            beta,*_ = np.linalg.lstsq(X1, cos, rcond=None)
            resid = cos - X1@beta
            obs = rb(resid, y)
            null = np.array([rb(resid, rng.permutation(y)) for _ in range(10000)])
            p,_,_ = C.emp_p(obs, null, tail="two-sided")
            rows.append(dict(layer=L, r_rb=obs, p=p))
            log(f"    [{tag}] L{L:<2} r_rb={obs:+.4f}  p={p:.2e}")
        df = pd.DataFrame(rows)
        df.to_csv(C.OUT/f"M1_confounds_{tag}.csv", index=False)
        Zs, pZ = C.stouffer(df.p.tolist())
        nsig = int((df.p < 0.05).sum())
        log(f"    [{tag}] {nsig}/12 layers p<0.05 | Stouffer Z={Zs:.2f} p={pZ:.3e}  -> {C.OUT.name}/M1_confounds_{tag}.csv")
    log(f"[done] {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
