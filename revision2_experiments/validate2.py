import sys, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C
import pandas as pd

genes = C.h5ad_genes(); n2i = {g: i for i, g in enumerate(genes)}
emb = C.load_emb()
cur_idx = C.curated_genes(genes)
print(f"curated ∩ vocab      : {len(cur_idx)}   (manuscript: 209)")
for L in (0, 11):
    nz = C.nonzero_mask(emb, L)
    print(f"  nonzero @L{L:<2}       : {sum(1 for i in cur_idx if nz[i])}  (manuscript: 195 after 14 OOV at input layer)")

t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF","target","mode","pmid"])
tfs_all = set(t.TF.astype(str)); targets_all = set(t.target.astype(str)) - tfs_all
import json
obj = json.loads(C.STRING_CACHE.read_text())
string_pairs = [(p["g1"], p["g2"]) for p in obj["pairs"] if p["score"] >= 0.7]
gene2go = C.load_gene2go()

US11 = C.layer_us(emb, 11)
sub = np.array(cur_idx); sub_names = [genes[i] for i in sub]
loc = {g: k for k, g in enumerate(sub_names)}

lab2 = C.pole_labels(US11[sub, 1])
tr = sorted({(a,b) for a,b in zip(t.TF.astype(str), t.target.astype(str)) if a in loc and b in loc and a != b})
st = sorted({(a,b) for a,b in string_pairs if a in loc and b in loc and a != b})
print(f"TRRUST pairs in set  : {len(tr)}   (manuscript: 326)")
print(f"STRING>=0.7 in set   : {len(st)}   (manuscript: 1,022)")
pi=np.array([loc[a] for a,b in tr]); pj=np.array([loc[b] for a,b in tr])
si=np.array([loc[a] for a,b in st]); sj=np.array([loc[b] for a,b in st])
rng=np.random.default_rng(0)
print(f"SV2 TRRUST co-pole L11: obs={C.copole_rate(lab2,pi,pj):.3f} null={C.copole_null(lab2,pi,pj,2000,rng).mean():.3f}  (manuscript 0.206 vs 0.122)")
print(f"SV2 STRING co-pole L11: obs={C.copole_rate(lab2,si,sj):.3f} null={C.copole_null(lab2,si,sj,2000,rng).mean():.3f}")

tf_pos=np.array([k for k,g in enumerate(sub_names) if g in tfs_all])
tgt_pos=np.array([k for k,g in enumerate(sub_names) if g in targets_all])
print(f"TF / target counts   : {tf_pos.size} / {tgt_pos.size}   (manuscript TF n=51)")
r1=[];loo=[]
for L in range(12):
    U=C.layer_us(emb,L)[sub]
    r1.append(C.tf_target_auroc(U,slice(1,7),tf_pos,tgt_pos,exclude_self=False))
    loo.append(C.tf_target_auroc(U,slice(1,7),tf_pos,tgt_pos,exclude_self=True))
print(f"TF/target AUROC r1   : mean={np.mean(r1):.3f} max={np.max(r1):.3f}@L{int(np.argmax(r1))}  (manuscript mean 0.744 max 0.789@L3)")
print(f"TF/target AUROC LOO  : mean={np.mean(loo):.3f} max={np.max(loo):.3f}@L{int(np.argmax(loo))}   bias={np.mean(r1)-np.mean(loo):+.3f}")

proj1=US11[sub,0]; top,bot,K=C.poles(proj1)
def fp(pole,term):
    a=sum(1 for k in range(len(sub_names)) if k in pole and term in gene2go.get(sub_names[k],set()))
    c=sum(1 for k in range(len(sub_names)) if k not in pole and term in gene2go.get(sub_names[k],set()))
    return C.fisher_2x2(a,len(pole)-a,c,(len(sub_names)-len(pole))-c)
ot,pt=fp(top,C.GO_EXTRACELLULAR); ob,pb=fp(bot,C.GO_EXTRACELLULAR)
print(f"SV1 extracellular L11: OR={max(ot,ob):.2f} p={min(pt,pb):.2e}  K={K}   (manuscript OR 6.37 p 2.6e-4)")
