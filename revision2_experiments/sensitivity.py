"""Basis-dependence of the spectral-collapse claim, and subset-SVD sensitivity."""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C
from scipy.stats import spearmanr

genes = C.h5ad_genes(); emb = C.load_emb(); nz0 = C.nonzero_mask(emb, 0)
cur209 = C.curated_genes(genes); cur187 = [i for i in cur209 if nz0[i]]

def eff_rank(X):
    s = np.linalg.svd(X - X.mean(0, keepdims=True), compute_uv=False)
    p = s / s.sum(); p = p[p > 0]
    return float(np.exp(-(p*np.log(p)).sum()))
def part_ratio(X):
    s = np.linalg.svd(X - X.mean(0, keepdims=True), compute_uv=False)
    v = s**2
    return float(v.sum()**2 / (v**2).sum())

print("Effective rank by basis   (manuscript claims 23.6 -> 1.6, a 14.4x drop)")
print(f"{'basis':28s} {'L0':>8s} {'L11':>8s} {'fold':>7s}")
for lbl, idx in [("full 4,803 genes", None), ("curated 209", cur209), ("curated 187 (L0 nonzero)", cur187)]:
    er = [eff_rank(np.asarray(emb[L], dtype=np.float64) if idx is None
                   else np.asarray(emb[L], dtype=np.float64)[idx]) for L in (0, 11)]
    print(f"{lbl:28s} {er[0]:8.2f} {er[1]:8.2f} {er[0]/er[1]:6.1f}x")
print()
print("Participation ratio       (manuscript claims 58 -> 9.5, a 6.1x drop, on '195-gene submatrix')")
for lbl, idx in [("full 4,803 genes", None), ("curated 209", cur209), ("curated 187 (L0 nonzero)", cur187)]:
    pr = [part_ratio(np.asarray(emb[L], dtype=np.float64) if idx is None
                     else np.asarray(emb[L], dtype=np.float64)[idx]) for L in (0, 11)]
    print(f"{lbl:28s} {pr[0]:8.2f} {pr[1]:8.2f} {pr[0]/pr[1]:6.1f}x")

print("\nSV1 variance fraction     (manuscript: 53.7% at L0 -> 93.4% at L11, full vocabulary)")
for lbl, idx in [("full 4,803 genes", None), ("curated 187", cur187)]:
    for L in (0, 11):
        X = np.asarray(emb[L], dtype=np.float64) if idx is None else np.asarray(emb[L], dtype=np.float64)[idx]
        s = np.linalg.svd(X - X.mean(0, keepdims=True), compute_uv=False)
        print(f"  {lbl:22s} L{L:<2} SV1 var fraction = {100*s[0]**2/(s**2).sum():5.1f}%")

# --- subset-SVD sensitivity for the headline co-pole and AUROC -------------
print("\nSubset-SVD sensitivity (loop convention) vs full-SVD (canonical), 187-gene set, N=20,000")
sub = np.array(cur187); names=[genes[i] for i in sub]; loc={g:k for k,g in enumerate(names)}
t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF","target","mode","pmid"])
tfs=set(t.TF.astype(str)); tgts=set(t.target.astype(str))-tfs
tr=sorted({(a,b) for a,b in zip(t.TF.astype(str),t.target.astype(str)) if a in loc and b in loc and a!=b})
obj=json.loads(C.STRING_CACHE.read_text())
st=sorted({(p["g1"],p["g2"]) for p in obj["pairs"] if p["score"]>=0.7 and p["g1"] in loc and p["g2"] in loc and p["g1"]!=p["g2"]})
tri=np.array([loc[a] for a,b in tr]); trj=np.array([loc[b] for a,b in tr])
sti=np.array([loc[a] for a,b in st]); stj=np.array([loc[b] for a,b in st])
tf_pos=np.array([k for k,g in enumerate(names) if g in tfs]); tgt_pos=np.array([k for k,g in enumerate(names) if g in tgts])
rng=np.random.default_rng(7)
def subset_us(L):
    X=np.asarray(emb[L],dtype=np.float64)[sub]; X=X-X.mean(0,keepdims=True)
    U,S,_=np.linalg.svd(X,full_matrices=False); return U*S
print(f"{'test':26s} {'full-SVD':>22s} {'subset-SVD':>22s}")
for tag, comp, pi, pj in [("SV2 TRRUST co-pole",1,tri,trj), ("SV2 STRING co-pole",1,sti,stj), ("SV3 STRING co-pole",2,sti,stj)]:
    out=[]
    for mode in ("full","subset"):
        best=(None,1.0,None)
        for L in range(12):
            U = C.layer_us(emb,L)[sub] if mode=="full" else subset_us(L)
            lab=C.pole_labels(U[:,comp]); obs=C.copole_rate(lab,pi,pj)
            p,_,_=C.emp_p(obs,C.copole_null(lab,pi,pj,20000,rng))
            if p<best[1]: best=(L,p,obs)
        out.append(f"L{best[0]} obs={best[2]:.3f} p={best[1]:.1e}")
    print(f"{tag:26s} {out[0]:>22s} {out[1]:>22s}")
aur=[]
for mode in ("full","subset"):
    vals=[C.tf_target_auroc(C.layer_us(emb,L)[sub] if mode=="full" else subset_us(L),
                            slice(1,7),tf_pos,tgt_pos,exclude_self=True) for L in range(12)]
    aur.append(f"mean={np.mean(vals):.3f} max={np.max(vals):.3f}")
print(f"{'TF/target AUROC (LOO)':26s} {aur[0]:>22s} {aur[1]:>22s}")
