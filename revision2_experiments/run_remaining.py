"""Round-2: the controls round 1 skipped or timed out on. See CONVENTION.md."""
import sys, json, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C
from scipy.stats import spearmanr

rng = np.random.default_rng(20260819); t0 = time.time()
OUT = C.OUT; res = {}
def log(*a): print(*a, flush=True)

# ---------- shared setup ----------------------------------------------------
genes = C.h5ad_genes(); emb = C.load_emb(); nz0 = C.nonzero_mask(emb, 0)
sub = np.array([i for i in C.curated_genes(genes) if nz0[i]])
names = [genes[i] for i in sub]; loc = {g: k for k, g in enumerate(names)}; n = len(names)
t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF","target","mode","pmid"])
tfs_all = set(t.TF.astype(str)); tgts_all = set(t.target.astype(str)) - tfs_all
tr = sorted({(a,b) for a,b in zip(t.TF.astype(str), t.target.astype(str)) if a in loc and b in loc and a!=b})
obj = json.loads(C.STRING_CACHE.read_text())
st_all = [(p["g1"],p["g2"],p["score"]) for p in obj["pairs"]]
gene2go = C.load_gene2go()
tf_pos = np.array([k for k,g in enumerate(names) if g in tfs_all])
tgt_pos = np.array([k for k,g in enumerate(names) if g in tgts_all])
log(f"[setup] n={n}  TRRUST={len(tr)}  TF={tf_pos.size}  {time.time()-t0:.0f}s")

def expression(h5, gene_list, max_cells=3000):
    """Return (cells x genes) dense log1p matrix for gene_list."""
    import h5py, scipy.sparse as sp
    with h5py.File(h5,"r") as f:
        gs = [g.decode() if isinstance(g,bytes) else str(g) for g in f["var"]["_index"][:]]
        gi = {g:i for i,g in enumerate(gs)}
        cols = [gi[g] for g in gene_list if g in gi]
        keep = [g for g in gene_list if g in gi]
        X = f["X"]
        if isinstance(X, h5py.Group):                       # CSR
            data,indices,indptr = X["data"], X["indices"], X["indptr"]
            ncell = indptr.shape[0]-1
            sel = np.sort(rng.choice(ncell, min(max_cells,ncell), replace=False))
            M = np.zeros((sel.size, len(cols)), dtype=np.float32)
            cpos = {c:k for k,c in enumerate(cols)}
            for r,ci in enumerate(sel):
                a,b = indptr[ci], indptr[ci+1]
                idx = indices[a:b]; val = data[a:b]
                m = np.isin(idx, cols)
                for j,v in zip(idx[m], val[m]): M[r, cpos[j]] = v
        else:
            ncell = X.shape[0]
            sel = np.sort(rng.choice(ncell, min(max_cells,ncell), replace=False))
            M = np.asarray(X[sel,:][:,cols], dtype=np.float32)
    return M, keep

# ============ M1: confound residualisation, all 12 layers ===================
log("\n[M1] multi-feature confound residualisation, all 12 layers")
E, kept = expression(C.H5AD, names)
kloc = {g:k for k,g in enumerate(kept)}
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
pos_pairs = tr
neg_pairs = []
tf_names = [names[k] for k in tf_pos]
for a in {x for x,_ in tr}:
    for b in names:
        if b!=a and (a,b) not in set(tr): neg_pairs.append((a,b))
neg_pairs = [neg_pairs[i] for i in rng.choice(len(neg_pairs), min(4*len(tr),len(neg_pairs)), replace=False)]
allp = pos_pairs + neg_pairs
y_lab = np.array([1]*len(pos_pairs) + [0]*len(neg_pairs))
Xc = np.column_stack([
    [coexp(a,b) for a,b in allp],
    [np.log1p(mexp[a]) for a,b in allp], [np.log1p(mexp[b]) for a,b in allp],
    [deg[a]+deg[b] for a,b in allp], [ngo[a]+ngo[b] for a,b in allp]]).astype(np.float64)
def rank_biserial(x, lab):
    r = pd.Series(x).rank().to_numpy(); n1,n0 = lab.sum(), (1-lab).sum()
    return float(2*(r[lab==1].mean() - r[lab==0].mean())/(n1+n0))
m1 = []
for L in range(12):
    U = C.layer_us(emb, L)[sub][:, 1:7]
    Un = U/(np.linalg.norm(U,axis=1,keepdims=True)+1e-12)
    cos = np.array([float(Un[loc[a]] @ Un[loc[b]]) for a,b in allp])
    X1 = np.column_stack([np.ones(len(allp)), Xc])
    beta,*_ = np.linalg.lstsq(X1, cos, rcond=None)
    resid = cos - X1@beta
    obs = rank_biserial(resid, y_lab)
    null = np.empty(10000)
    for i in range(10000):
        null[i] = rank_biserial(resid, rng.permutation(y_lab))
    p,_,_ = C.emp_p(obs, null, tail="two-sided")
    m1.append((L, obs, p)); log(f"    L{L:<2} r_rb={obs:+.4f}  p={p:.2e}")
Zs, pZ = C.stouffer([x[2] for x in m1])
log(f"    Stouffer across 12 layers: Z={Zs:.2f} p={pZ:.3e}   (round 1: L0 only, r_rb=0.061, p=0.080)")
res["M1"] = dict(per_layer=[dict(layer=l, r_rb=o, p=p) for l,o,p in m1], stouffer_Z=Zs, stouffer_p=pZ)
pd.DataFrame(m1, columns=["layer","r_rb","p"]).to_csv(OUT/"M1_confounds_all_layers.csv", index=False)

# ============ M2: matched ESM2 control ======================================
log("\n[M2] matched ESM2 control")
z = np.load(C.ESM2_NPZ, allow_pickle=True)
eg = [str(x) for x in z["gene_names"]]; eemb = z["embeddings"]
egi = {g:i for i,g in enumerate(eg)}
shared = [g for g in names if g in egi]
log(f"    shared genes (curated {n} ∩ ESM2): {len(shared)}")
sh_loc = {g:k for k,g in enumerate(shared)}
etf = np.array([k for k,g in enumerate(shared) if g in tfs_all])
etg = np.array([k for k,g in enumerate(shared) if g in tgts_all])
Xe = eemb[[egi[g] for g in shared]].astype(np.float64)
def us_of(X):
    Xc_ = X - X.mean(0, keepdims=True)
    U,S,_ = np.linalg.svd(Xc_, full_matrices=False); return U*S
esm_full = C.tf_target_auroc(us_of(Xe), slice(1,7), etf, etg)
P = rng.normal(size=(Xe.shape[1], 1024))/np.sqrt(1024)
esm_proj = C.tf_target_auroc(us_of(Xe@P), slice(1,7), etf, etg)
sc_sub = np.array([loc[g] for g in shared])
scg = np.asarray(emb[11], dtype=np.float64)[sub][sc_sub]
sc_full = C.tf_target_auroc(us_of(scg), slice(1,7), etf, etg)
P2 = rng.normal(size=(scg.shape[1], 102))/np.sqrt(102)     # same 5x reduction ratio
sc_proj = C.tf_target_auroc(us_of(scg@P2), slice(1,7), etf, etg)
log(f"    ESM2  un-projected 5120d : AUROC={esm_full:.3f}   (manuscript claims 0.889, on a DIFFERENT gene set after projection)")
log(f"    ESM2  projected -> 1024d : AUROC={esm_proj:.3f}")
log(f"    scGPT un-projected  512d : AUROC={sc_full:.3f}   [matched gene set]")
log(f"    scGPT projected  -> 102d : AUROC={sc_proj:.3f}   [MATCHED CONTROL, same 5x reduction]")
def sv_struct(X):
    U = us_of(X); lab = C.pole_labels(U[:,1])
    sp_ = [(a,b) for a,b,s in st_all if s>=0.7 and a in sh_loc and b in sh_loc and a!=b]
    sp_ = sorted(set(sp_))
    if not sp_: return float('nan'), float('nan')
    pi=np.array([sh_loc[a] for a,b in sp_]); pj=np.array([sh_loc[b] for a,b in sp_])
    obs=C.copole_rate(lab,pi,pj); nul=C.copole_null(lab,pi,pj,20000,rng)
    return obs, float((obs-nul.mean())/(nul.std()+1e-12))
for lbl, X in [("ESM2 5120d",Xe), ("ESM2 ->1024d",Xe@P), ("scGPT 512d",scg), ("scGPT ->102d",scg@P2)]:
    o,zz = sv_struct(X); log(f"    SV2 STRING co-pole  {lbl:14s} rate={o:.3f}  z={zz:+.2f}")
res["M2"] = dict(n_shared=len(shared), esm2_unprojected=esm_full, esm2_projected=esm_proj,
                 scgpt_matched=sc_full, scgpt_projected=sc_proj)

# ============ M3: kidney with matched TF set ================================
log("\n[M3] kidney replication with a matched TF set")
try:
    kg = C.h5ad_genes(C.H5AD_KIDNEY); kemb = np.load(C.EMB_KIDNEY, mmap_mode="r")
    log(f"    kidney matrix {kemb.shape}, gene axis {len(kg)}")
    kloc2 = {g:i for i,g in enumerate(kg)}
    knz = C.nonzero_mask(kemb, 0)
    imm_tf = {names[k] for k in tf_pos}
    kcur = [i for i,g in enumerate(kg) if g in set(names) and knz[i]]
    knames = [kg[i] for i in kcur]; kl = {g:k for k,g in enumerate(knames)}
    ktf_matched = np.array([k for k,g in enumerate(knames) if g in imm_tf])
    ktg = np.array([k for k,g in enumerate(knames) if g in tgts_all])
    ktf_all = np.array([k for k,g in enumerate(knames) if g in tfs_all])
    log(f"    kidney ∩ immune curated set: {len(knames)} genes | matched TF n={ktf_matched.size} | all-TRRUST TF n={ktf_all.size}")
    kk = np.array(kcur)
    kres = []
    for L in range(min(4, kemb.shape[0])):
        Uk = C.layer_us(kemb, L)[kk]
        a_m = C.tf_target_auroc(Uk, slice(1,7), ktf_matched, ktg) if ktf_matched.size>4 else float('nan')
        a_a = C.tf_target_auroc(Uk, slice(1,7), ktf_all, ktg) if ktf_all.size>4 else float('nan')
        kres.append((L,a_m,a_a)); log(f"    L{L}: matched-TF AUROC={a_m:.3f} | all-TRRUST-TF AUROC={a_a:.3f}")
    res["M3"] = dict(n_genes=len(knames), n_tf_matched=int(ktf_matched.size),
                     per_layer=[dict(layer=l, auroc_matched=m, auroc_all=a) for l,m,a in kres])
except Exception as e:
    log(f"    SKIPPED: {type(e).__name__}: {e}"); res["M3"] = dict(error=str(e))

# ============ P0.3 + M4: GRN benchmark with real baselines + bootstrap ======
log("\n[P0.3+M4] held-out GRN benchmark: embeddings vs GENIE3 / GRNBoost2 / co-expression")
pos_list = sorted(tr); perm = rng.permutation(len(pos_list))
n_test = max(1, int(0.20*len(pos_list)))
held_pos = [pos_list[i] for i in perm[:n_test]]
train_pos = {pos_list[i] for i in perm[n_test:]}
held_tfs = sorted({a for a,_ in held_pos})
held_neg = [(tf,g) for tf in held_tfs for g in names
            if g!=tf and (tf,g) not in train_pos and (tf,g) not in set(held_pos)]
held_neg = [held_neg[i] for i in rng.choice(len(held_neg), min(len(held_pos)*3, len(held_neg)), replace=False)]
log(f"    held-out: {len(held_pos)} positives, {len(held_neg)} negatives  (round 1: 57 / 124)")
scorers = {}
U0 = C.layer_us(emb,0)[sub]
full0 = np.asarray(emb[0], dtype=np.float64)[sub]
def cos_scorer(M):
    Mn = M/(np.linalg.norm(M,axis=1,keepdims=True)+1e-12)
    return lambda a,b: float(Mn[loc[a]] @ Mn[loc[b]])
scorers["scGPT_full512_L0"] = cos_scorer(full0)
scorers["scGPT_spectral_SV5_7_L0"] = cos_scorer(U0[:,4:7])
scorers["coexpression_abs"] = lambda a,b: abs(coexp(a,b))
scorers["coexpression_signed"] = lambda a,b: coexp(a,b)
try:
    from arboreto.algo import grnboost2, genie3
    from distributed import Client, LocalCluster
    df = pd.DataFrame(E, columns=kept)
    cluster = LocalCluster(n_workers=6, threads_per_worker=1, memory_limit="3GB",
                           processes=True, dashboard_address=None)
    client = Client(cluster)
    for nm, fn in [("GRNBoost2", grnboost2), ("GENIE3", genie3)]:
        log(f"    running {nm} on {df.shape[1]} genes x {df.shape[0]} cells ...")
        net = fn(expression_data=df, tf_names=[g for g in kept if g in tfs_all],
                 client_or_address=client, seed=42, verbose=False)
        d = {(r.TF, r.target): r.importance for r in net.itertuples()}
        scorers[nm] = (lambda dd: (lambda a,b: float(dd.get((a,b), 0.0))))(d)
        log(f"      {nm} produced {len(d)} scored edges")
    client.close(); cluster.close()
except Exception as e:
    log(f"    arboreto FAILED: {type(e).__name__}: {e}")
def evaluate(fn):
    p = np.array([fn(a,b) for a,b in held_pos]); q = np.array([fn(a,b) for a,b in held_neg])
    return C.auroc(p,q), p, q
B = 2000; grn_rows = []
for nm, fn in scorers.items():
    a, p, q = evaluate(fn)
    bs = np.empty(B)
    for i in range(B):
        bs[i] = C.auroc(p[rng.integers(0,p.size,p.size)], q[rng.integers(0,q.size,q.size)])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    grn_rows.append(dict(method=nm, auroc=a, ci_lo=lo, ci_hi=hi))
    log(f"    {nm:26s} AUROC={a:.3f}  95% CI [{lo:.3f}, {hi:.3f}]")
base = dict((r["method"], r) for r in grn_rows)
if "scGPT_full512_L0" in scorers and "coexpression_abs" in scorers:
    _,p1,q1 = evaluate(scorers["scGPT_full512_L0"]); _,p2,q2 = evaluate(scorers["coexpression_abs"])
    d = np.empty(B)
    for i in range(B):
        ip = rng.integers(0,p1.size,p1.size); iq = rng.integers(0,q1.size,q1.size)
        d[i] = C.auroc(p1[ip],q1[iq]) - C.auroc(p2[ip],q2[iq])
    lo,hi = np.percentile(d,[2.5,97.5])
    log(f"    paired delta scGPT_full512 - |coexpression| = {d.mean():+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]"
        f"  -> {'separates from 0' if lo>0 else 'DOES NOT separate from 0'}")
    res["grn_delta"] = dict(mean=float(d.mean()), ci_lo=float(lo), ci_hi=float(hi))
pd.DataFrame(grn_rows).to_csv(OUT/"GRN_benchmark_round2.csv", index=False)
res["GRN"] = grn_rows
json.dump(res, open(OUT/"remaining_results.json","w"), indent=1, default=float)
log(f"\n[done] {time.time()-t0:.0f}s")
