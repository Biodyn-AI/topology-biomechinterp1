"""Reproduce the round-1 headline numbers with the round-2 loader.

If these match the manuscript, the re-implementation is faithful and any later
difference is attributable to the change under test (permutation count,
self-similarity handling, etc.) rather than to a different pipeline.
"""
import sys, json
import numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
import common as C

rng = np.random.default_rng(0)
genes = C.h5ad_genes()
n2i = {g: i for i, g in enumerate(genes)}
print(f"gene axis            : {len(genes)} genes  (manuscript: 4,803)")

curated = [l.strip() for l in open(C.REPO / "subproject_38_geometric_residual_stream_interpretability"
                                  / "implementation/outputs/cycle1_main/gene_list.txt") if l.strip()]
print(f"curated gene list    : {len(curated)}  (manuscript: 209)")

emb = C.load_emb()
print(f"embeddings           : {emb.shape}")

L = 11
nz = C.nonzero_mask(emb, L)
cur_idx = [n2i[g] for g in curated if g in n2i]
cur_in = [i for i in cur_idx if nz[i]]
print(f"curated in-vocab @L11: {len(cur_in)}  (manuscript: 195 after dropping 14 zero-norm)")

pairs, tfs, targets_only, signed = C.load_trrust()
cur_names = [genes[i] for i in cur_in]
cur_set = set(cur_names)
tr_pairs = [(a, b) for a, b in pairs if a in cur_set and b in cur_set and a != b]
tr_pairs = sorted(set(tr_pairs))
print(f"TRRUST pairs in set  : {len(tr_pairs)}  (manuscript: 326)")

string = C.load_string(min_score=0.7)
st_pairs = sorted({(a, b) for a, b, s in string if a in cur_set and b in cur_set and a != b})
print(f"STRING>=0.7 in set   : {len(st_pairs)}  (manuscript: 1,022)")

n_tf = sum(1 for g in cur_names if g in tfs)
n_tgt = sum(1 for g in cur_names if g in targets_only)
print(f"curated TF / target  : {n_tf} / {n_tgt}  (manuscript: TF n=51)")

# --- headline statistics, round-1 semantics --------------------------------
US = C.layer_us(emb, L, gene_idx=cur_in)
loc = {g: k for k, g in enumerate(cur_names)}

lab2 = C.pole_labels(US[:, 1])
pi = np.array([loc[a] for a, b in tr_pairs]); pj = np.array([loc[b] for a, b in tr_pairs])
obs_tr = C.copole_rate(lab2, pi, pj)
null_tr = C.copole_null(lab2, pi, pj, 2000, rng)
print(f"\nSV2 TRRUST co-pole L11: obs={obs_tr:.3f}  null_mean={null_tr.mean():.3f}   (manuscript: 0.206 vs 0.122)")

si = np.array([loc[a] for a, b in st_pairs]); sj = np.array([loc[b] for a, b in st_pairs])
obs_st = C.copole_rate(lab2, si, sj)
null_st = C.copole_null(lab2, si, sj, 2000, rng)
print(f"SV2 STRING co-pole L11: obs={obs_st:.3f}  null_mean={null_st.mean():.3f}   (manuscript mean over layers: 0.226 vs 0.124)")

tf_pos = np.array([k for k, g in enumerate(cur_names) if g in tfs])
tgt_pos = np.array([k for k, g in enumerate(cur_names) if g in targets_only])
aurocs_r1, aurocs_loo = [], []
for l in range(12):
    U = C.layer_us(emb, l, gene_idx=cur_in)
    aurocs_r1.append(C.tf_target_auroc(U, slice(1, 7), tf_pos, tgt_pos, exclude_self=False))
    aurocs_loo.append(C.tf_target_auroc(U, slice(1, 7), tf_pos, tgt_pos, exclude_self=True))
print(f"\nTF/target AUROC SV2-SV7 (round-1 estimator, self-similarity INCLUDED):")
print(f"   mean={np.mean(aurocs_r1):.3f}  max={np.max(aurocs_r1):.3f} @L{int(np.argmax(aurocs_r1))}   (manuscript: mean 0.744, max 0.789 @L3)")
print(f"TF/target AUROC SV2-SV7 (leave-one-out, self-similarity EXCLUDED):")
print(f"   mean={np.mean(aurocs_loo):.3f}  max={np.max(aurocs_loo):.3f} @L{int(np.argmax(aurocs_loo))}")
print(f"   bias from self-inclusion: {np.mean(aurocs_r1)-np.mean(aurocs_loo):+.3f} AUROC")

gene2go = C.load_gene2go()
proj1 = US[:, 0]
top, bot, K = C.poles(proj1)
def fisher_pole(pole, term):
    a = sum(1 for k in range(len(cur_names)) if k in pole and term in gene2go.get(cur_names[k], set()))
    b = len(pole) - a
    c = sum(1 for k in range(len(cur_names)) if k not in pole and term in gene2go.get(cur_names[k], set()))
    d = (len(cur_names) - len(pole)) - c
    return C.fisher_2x2(a, b, c, d)
or_t, p_t = fisher_pole(top, C.GO_EXTRACELLULAR); or_b, p_b = fisher_pole(bot, C.GO_EXTRACELLULAR)
print(f"\nSV1 extracellular L11 : OR={max(or_t,or_b):.2f}  p={min(p_t,p_b):.2e}   (manuscript: OR 6.37, p 2.6e-4)  K={K}")
