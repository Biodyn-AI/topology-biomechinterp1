"""Verbatim re-run of the SUBMITTED cross-model code (revision2_experiments/crossmodel.py).

Only change: the output CSV is written here (legacy_submitted_repro_raw.csv) instead of into
revision2_experiments/, which this experiment must not modify. The estimator, gene set, RNG seed and
model order are untouched, so the numbers must equal revision2_experiments/crossmodel_round2.csv.
"""
import sys, json, numpy as np, pandas as pd
from pathlib import Path
R2 = Path(__file__).resolve().parents[2] / "revision2_experiments"
sys.path.insert(0, str(R2))
import common as C
OUTCSV = Path(__file__).resolve().parent / "legacy_submitted_repro_raw.csv"
rng = np.random.default_rng(20260819)
genes = C.h5ad_genes(); emb = C.load_emb(); nz0 = C.nonzero_mask(emb, 0)
sub = np.array([i for i in C.curated_genes(genes) if nz0[i]])
names = [genes[i] for i in sub]
t = pd.read_csv(C.TRRUST_PATH, sep="\t", header=None, names=["TF","target","mode","pmid"])
tfs = set(t.TF.astype(str)); tgts = set(t.target.astype(str)) - tfs
obj = json.loads(C.STRING_CACHE.read_text())

def load_npz(p):
    z = np.load(p, allow_pickle=True)
    k = "gene_names" if "gene_names" in z else "genes"
    e = "embeddings" if "embeddings" in z else "emb"
    return [str(x) for x in z[k]], z[e]

models = {"scGPT L11": None}
for lbl, p in [("Geneformer V1-10M", C.GF_V1), ("Geneformer V2-104M", C.GF_V2), ("ESM2-3B", C.ESM2_NPZ)]:
    try: models[lbl] = load_npz(p)
    except Exception as ex: print(f"  skip {lbl}: {ex}")

common = set(names)
for lbl, mv in models.items():
    if mv is not None: common &= set(mv[0])
common = sorted(common)
print(f"curated 209-gene file -> in kidney_A axis -> nonzero at L0: {len(names)}")
print(f"matched gene set present in all {len(models)} models: {len(common)}")
loc = {g: k for k, g in enumerate(common)}
tf_pos = np.array([k for k, g in enumerate(common) if g in tfs])
tgt_pos = np.array([k for k, g in enumerate(common) if g in tgts])
st = sorted({(p["g1"],p["g2"]) for p in obj["pairs"] if p["score"]>=0.7
             and p["g1"] in loc and p["g2"] in loc and p["g1"]!=p["g2"]})
pi = np.array([loc[a] for a,b in st]); pj = np.array([loc[b] for a,b in st])
print(f"TF={tf_pos.size} target={tgt_pos.size} STRING pairs={len(st)}")

def us_of(X):
    Xc = np.asarray(X, dtype=np.float64); Xc = Xc - Xc.mean(0, keepdims=True)
    U,S,_ = np.linalg.svd(Xc, full_matrices=False); return U*S
rows = []
for lbl, mv in models.items():
    if lbl == "scGPT L11":
        X = np.asarray(emb[11], dtype=np.float64)[[list(genes).index(g) for g in common]]
    else:
        gn, ge = mv; gi = {g:i for i,g in enumerate(gn)}
        X = ge[[gi[g] for g in common]]
    U = us_of(X)
    au = C.tf_target_auroc(U, slice(1,7), tf_pos, tgt_pos)
    lab = C.pole_labels(U[:,1]); obs = C.copole_rate(lab, pi, pj)
    nul = C.copole_null(lab, pi, pj, 20000, rng)
    p,_,_ = C.emp_p(obs, nul)
    z = (obs - nul.mean())/(nul.std()+1e-12)
    rows.append(dict(model=lbl, dim=X.shape[1], tf_target_auroc=au,
                     ppi_copole=obs, ppi_z=z, ppi_p=p, n_genes=len(common), n_tf=tf_pos.size,
                     n_target=tgt_pos.size, n_string_functional_pairs=len(st)))
    print(f"  {lbl:20s} d={X.shape[1]:5d}  TF/target AUROC={au:.3f}  STRING co-pole={obs:.3f} z={z:+.2f} p={p:.1e}")
pd.DataFrame(rows).to_csv(OUTCSV, index=False)
(Path(__file__).resolve().parent / "legacy_genes.json").write_text(json.dumps(dict(
    curated_nonzero_kidneyA=names, legacy186=common)))
