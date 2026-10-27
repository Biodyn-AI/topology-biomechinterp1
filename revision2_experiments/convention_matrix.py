"""Which analysis convention reproduces the manuscript's headline numbers?

The manuscript reports, for the immune curated set at L11:
   195 genes | 326 TRRUST pairs | 1,022 STRING>=0.7 pairs | TF n=51
   SV1 extracellular OR = 6.37 (p=2.6e-4) | SV2 TRRUST co-pole 0.206 vs null 0.122
   TF/target AUROC (SV2-SV7) mean 0.744, max 0.789 @L3

Three conventions are in play in the deposited material:
   A  loop (iter_0010): SVD on the 209-gene SUBmatrix, K=52 fixed,
      null = resample random gene pairs
   B  revision E3:      SVD on the FULL 4,803 matrix, K=0.27n within subset,
      null = gene-label shuffle
   C  Methods 4.4 text: K=52, "gene-label shuffle null, N=200-1,000"
"""
import sys, json, itertools, numpy as np, pandas as pd
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C

t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF","target","mode","pmid"])
tfs_all = set(t.TF.astype(str)); tgts_all = set(t.target.astype(str)) - tfs_all
tr_raw = list(zip(t.TF.astype(str), t.target.astype(str)))
obj = json.loads(C.STRING_CACHE.read_text())
st_raw = [(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= 0.7]
gene2go = C.load_gene2go()

MATS = {
 "cycle1_main(4803)": (C.EMB_IMMUNE, C.H5AD, C.REPO/"subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"),
 "cycle4_immune(4941)": (C.REPO/"subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main/layer_gene_embeddings.npy",
                         C.REPO/"single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad",
                         C.REPO/"subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main/cycle1_edge_dataset.tsv"),
}
print(f"{'matrix':22s} {'nzfilt':7s} {'svd':7s} {'dedup':6s} | {'nG':>4s} {'nTR':>5s} {'nST':>5s} {'nTF':>4s} | {'OR1':>6s} {'copole':>7s} {'AUROC':>6s}")
print("-"*104)
for mname,(embp,h5,tsv) in MATS.items():
    genes = C.h5ad_genes(h5); emb = np.load(embp, mmap_mode="r")
    e = pd.read_csv(tsv, sep="\t"); named = set(e["source"].astype(str))|set(e["target"].astype(str))
    base = [i for i,g in enumerate(genes) if g in named]
    nz0 = C.nonzero_mask(emb,0)
    for nzfilt in (False,True):
        sub = np.array([i for i in base if (nz0[i] or not nzfilt)])
        names=[genes[i] for i in sub]; loc={g:k for k,g in enumerate(names)}
        for dedup in (True,False):
            trp=[(a,b) for a,b in tr_raw if a in loc and b in loc and a!=b]
            stp=[(a,b) for a,b in st_raw if a in loc and b in loc and a!=b]
            if dedup: trp=sorted(set(trp)); stp=sorted(set(stp))
            for svd in ("subset","full"):
                US = C.layer_us(emb,11,gene_idx=sub) if svd=="full" else None
                if svd=="subset":
                    X=np.asarray(emb[11],dtype=np.float64)[sub]; X=X-X.mean(0,keepdims=True)
                    U,S,_=np.linalg.svd(X,full_matrices=False); US=U*S
                K=52
                order=np.argsort(US[:,1]); lab=np.zeros(len(sub),dtype=np.int8)
                lab[order[-K:]]=1; lab[order[:K]]=-1
                pi=np.array([loc[a] for a,b in trp]); pj=np.array([loc[b] for a,b in trp])
                cop=C.copole_rate(lab,pi,pj) if pi.size else float('nan')
                o1=np.argsort(US[:,0]); top=set(o1[-K:].tolist()); bot=set(o1[:K].tolist())
                def fp(pole):
                    a=sum(1 for k in range(len(names)) if k in pole and C.GO_EXTRACELLULAR in gene2go.get(names[k],set()))
                    c=sum(1 for k in range(len(names)) if k not in pole and C.GO_EXTRACELLULAR in gene2go.get(names[k],set()))
                    return C.fisher_2x2(a,len(pole)-a,c,(len(names)-len(pole))-c)
                ot,pt=fp(top); ob,pb=fp(bot)
                tf_pos=np.array([k for k,g in enumerate(names) if g in tfs_all]); tgt_pos=np.array([k for k,g in enumerate(names) if g in tgts_all])
                au=C.tf_target_auroc(US,slice(1,7),tf_pos,tgt_pos,exclude_self=False) if tf_pos.size>4 else float('nan')
                print(f"{mname:22s} {str(nzfilt):7s} {svd:7s} {str(dedup):6s} | {len(sub):4d} {len(trp):5d} {len(stp):5d} {tf_pos.size:4d} | {max(ot,ob):6.2f} {cop:7.3f} {au:6.3f}")
print("-"*104)
print(f"{'MANUSCRIPT TARGET':22s} {'':7s} {'':7s} {'':6s} | {195:4d} {326:5d} {1022:5d} {51:4d} | {6.37:6.2f} {0.206:7.3f} {0.744:6.3f}   (AUROC is the 12-layer mean)")
