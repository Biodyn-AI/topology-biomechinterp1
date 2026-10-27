"""Test whether the manuscript's headline numbers come from the cycle4_immune
(4,941-gene) matrix rather than the cycle1_main (4,803-gene) matrix that the
deposited revision code and Methods 4.1 describe."""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C

EMB4 = C.REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main/layer_gene_embeddings.npy"
H5AD4 = C.REPO / "single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad"

genes = C.h5ad_genes(H5AD4)
emb = np.load(EMB4, mmap_mode="r")
print(f"matrix: {emb.shape}   gene axis: {len(genes)}")

e = pd.read_csv(C.REPO / "subproject_38_geometric_residual_stream_interpretability"
                / "implementation/outputs/cycle4_immune_main/cycle1_edge_dataset.tsv", sep="\t")
named = set(e["source"].astype(str)) | set(e["target"].astype(str))
cur = sorted([i for i, g in enumerate(genes) if g in named])
print(f"curated ∩ vocab      : {len(cur)}    (manuscript: 209)")
nz0 = C.nonzero_mask(emb, 0)
cur_nz = [i for i in cur if nz0[i]]
print(f"nonzero @L0          : {len(cur_nz)}    (manuscript: 195, i.e. 14 OOV)  -> dropped {len(cur)-len(cur_nz)}")

sub = np.array(cur_nz); names = [genes[i] for i in sub]; loc = {g: k for k, g in enumerate(names)}
t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF","target","mode","pmid"])
tfs = set(t.TF.astype(str)); tgts = set(t.target.astype(str)) - tfs
obj = json.loads(C.STRING_CACHE.read_text())
st = sorted({(p["g1"],p["g2"]) for p in obj["pairs"] if p["score"]>=0.7 and p["g1"] in loc and p["g2"] in loc and p["g1"]!=p["g2"]})
tr = sorted({(a,b) for a,b in zip(t.TF.astype(str), t.target.astype(str)) if a in loc and b in loc and a!=b})
print(f"TRRUST pairs in set  : {len(tr)}    (manuscript: 326)")
print(f"STRING>=0.7 in set   : {len(st)}   (manuscript: 1,022)")
tf_pos = np.array([k for k,g in enumerate(names) if g in tfs]); tgt_pos = np.array([k for k,g in enumerate(names) if g in tgts])
print(f"TF / target counts   : {tf_pos.size} / {tgt_pos.size}     (manuscript TF n=51)")

US11 = C.layer_us(emb, 11)[sub]
lab2 = C.pole_labels(US11[:,1]); rng = np.random.default_rng(0)
pi=np.array([loc[a] for a,b in tr]); pj=np.array([loc[b] for a,b in tr])
print(f"SV2 TRRUST co-pole L11: obs={C.copole_rate(lab2,pi,pj):.3f} null={C.copole_null(lab2,pi,pj,2000,rng).mean():.3f}   (manuscript 0.206 vs 0.122)")
gene2go = C.load_gene2go()
top,bot,K = C.poles(US11[:,0])
def fp(pole,term):
    a=sum(1 for k in range(len(names)) if k in pole and term in gene2go.get(names[k],set()))
    c=sum(1 for k in range(len(names)) if k not in pole and term in gene2go.get(names[k],set()))
    return C.fisher_2x2(a,len(pole)-a,c,(len(names)-len(pole))-c)
ot,pt=fp(top,C.GO_EXTRACELLULAR); ob,pb=fp(bot,C.GO_EXTRACELLULAR)
print(f"SV1 extracellular L11: OR={max(ot,ob):.2f} p={min(pt,pb):.2e}  K={K}   (manuscript OR 6.37, p 2.6e-4)")
r1=[];loo=[]
for L in range(12):
    U=C.layer_us(emb,L)[sub]
    r1.append(C.tf_target_auroc(U,slice(1,7),tf_pos,tgt_pos,exclude_self=False))
    loo.append(C.tf_target_auroc(U,slice(1,7),tf_pos,tgt_pos,exclude_self=True))
print(f"TF/target AUROC r1   : mean={np.mean(r1):.3f} max={np.max(r1):.3f}@L{int(np.argmax(r1))}   (manuscript mean 0.744, max 0.789@L3)")
print(f"TF/target AUROC LOO  : mean={np.mean(loo):.3f} max={np.max(loo):.3f}@L{int(np.argmax(loo))}    self-inclusion bias={np.mean(r1)-np.mean(loo):+.3f}")
