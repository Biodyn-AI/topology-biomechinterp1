"""Independent verifier (A6 / E9): SV1 extracellular OR (iter_0006), legacy->observed crosswalk on kidney_A,
effective-rank folds, and the iter_0062 edge AUROC on immune. Own code; only raw inputs (+ gene2go pickle).
"""
import csv, json, pickle
from pathlib import Path
from collections import defaultdict
import numpy as np, h5py
from scipy.stats import fisher_exact, spearmanr

P41 = Path(__file__).resolve().parents[3]
REPO = P41.parent
S38 = REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs"
CYC1, CYC4 = S38 / "cycle1_main", S38 / "cycle4_immune_main"
out = {}

def svd_c(X):
    X = np.asarray(X, dtype=np.float64); X = X - X.mean(0, keepdims=True)
    return np.linalg.svd(X, full_matrices=False)

def er(S):
    p = S**2 / (S**2).sum(); p = p[p > 0]; return float(np.exp(-(p*np.log(p)).sum()))

def named(run):
    g2i = {}
    for r in csv.DictReader(open(run / "cycle1_edge_dataset.tsv"), delimiter="\t"):
        g2i.setdefault(r["source"], int(r["source_idx"])); g2i.setdefault(r["target"], int(r["target_idx"]))
    return g2i

# ================= kidney_A =================
E = np.load(CYC1 / "layer_gene_embeddings.npy", mmap_mode="r")
nz = np.linalg.norm(np.asarray(E[0], dtype=np.float64), axis=1) > 0
obs = np.flatnonzero(nz)
g2i = named(CYC1)

# --- iter_0006 H02: full 4,803-row SVD at L11, SV1 signed loading, top vs bottom quartile, GO:0005615 ---
gene2go = pickle.load(open(REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl", "rb"))
def go_or(loading_by_gene, term="GO:0005615"):
    genes = list(g2i.keys())
    gg = {g: gene2go.get(g, set()) for g in genes if g in gene2go}
    t2g = defaultdict(set)
    for g, ts in gg.items():
        for t in ts: t2g[t].add(g)
    order = sorted(genes, key=lambda g: -loading_by_gene[g])
    n25 = max(5, len(order) // 4); top, bot = set(order[:n25]), set(order[-n25:])
    tg = t2g[term]; a = len(tg & top); c = len(tg & bot)
    OR, p = fisher_exact([[a, n25 - a], [c, n25 - c]], alternative="greater")
    both = []
    for sgn in (1, -1):   # sign of an SVD vector is arbitrary: report the best pole too
        o = sorted(genes, key=lambda g: -sgn*loading_by_gene[g]); tp = set(o[:n25]); bt = set(o[-n25:])
        aa = len(tg & tp); cc = len(tg & bt)
        both.append(fisher_exact([[aa, n25-aa], [cc, n25-cc]], alternative="greater"))
    return dict(a=a, c=c, n25=n25, OR=float(OR), p=float(p), best_pole_p=float(min(b[1] for b in both)),
                term_size_in_named=len(tg))
U, S, Vt = svd_c(E[11])
out["kidA_L11_full_SV1_extracellular"] = go_or({g: U[i, 0] for g, i in g2i.items()})
out["kidA_L11_full_SV1_corr_with_observed_indicator"] = float(np.corrcoef(U[:, 0], nz.astype(float))[0, 1])
out["kidA_L11_full_SV1_varfrac"] = float(S[0]**2 / (S**2).sum())
# observed-gene SVD: genes with zero rows have no coordinate -> project them with observed mean centring
Xo = np.asarray(E[11], dtype=np.float64); mu = Xo[obs].mean(0)
Uo, So, Vto = np.linalg.svd(Xo[obs] - mu, full_matrices=False)
proj_all = (Xo - mu) @ Vto[0]
out["kidA_L11_observed_SV1_extracellular_projall"] = go_or({g: proj_all[i] for g, i in g2i.items()})
obs_named = {g: i for g, i in g2i.items() if nz[i]}
g2i_save = g2i; g2i = obs_named
out["kidA_L11_observed_SV1_extracellular_observed_named_only"] = go_or({g: proj_all[i] for g, i in obs_named.items()})
g2i = g2i_save

# --- crosswalk: legacy (full 4,803 incl zero rows) SVk vs observed SVD, layers 0..11 ---
cw = []; match_km1 = 0; tot = 0; cosv = []; sv1best = []
ang = {}
for l in range(12):
    _, _, VL = svd_c(E[l]); _, _, VO = svd_c(np.asarray(E[l])[obs])
    C = np.abs(VL[:12] @ VO[:13].T)
    sv1best.append(float(C[0].max()))
    for k in range(1, 12):           # legacy SV2..SV12 (0-based 1..11)
        tot += 1; b = int(C[k].argmax()); match_km1 += int(b == k - 1); cosv.append(float(C[k, k-1]))
    pa = lambda A, B: float(np.mean(np.linalg.svd(A @ B.T, compute_uv=False)**2))
    ang[l] = dict(SV2_7_vs_obs1_6=pa(VL[1:7], VO[0:6]), SV2_7_vs_obs2_7=pa(VL[1:7], VO[1:7]),
                  SV2_4_vs_obs1_3=pa(VL[1:4], VO[0:3]), SV5_7_vs_obs4_6=pa(VL[4:7], VO[3:6]))
out["crosswalk_frac_legacy_k_best_is_observed_k-1"] = match_km1 / tot
out["crosswalk_mean_abs_cos_legacy_k_vs_obs_k-1"] = float(np.mean(cosv))
out["crosswalk_legacy_SV1_mean_best_abs_cos"] = float(np.mean(sv1best))
out["crosswalk_band_mean_cos2_layer_avg"] = {k: float(np.mean([ang[l][k] for l in range(12)])) for k in ang[0]}

# --- effective rank ---
erf = [er(svd_c(E[l])[1]) for l in (0, 11)]; ero = [er(svd_c(np.asarray(E[l])[obs])[1]) for l in (0, 11)]
out["kidA_ER_full_L0_L11_fold"] = [erf[0], erf[1], erf[0]/erf[1]]
out["kidA_ER_observed_L0_L11_fold"] = [ero[0], ero[1], ero[0]/ero[1]]

# ================= immune =================
EI = np.load(CYC4 / "layer_gene_embeddings.npy", mmap_mode="r")
nzI = np.linalg.norm(np.asarray(EI[0], dtype=np.float64), axis=1) > 0; obsI = np.flatnonzero(nzI)
out["immune_G_nobs"] = [int(EI.shape[1]), int(nzI.sum())]
eri = [er(svd_c(np.asarray(EI[l])[obsI])[1]) for l in (0, 11)]
out["immune_ER_observed_L0_L11_fold"] = [eri[0], eri[1], eri[0]/eri[1]]
rs = np.random.default_rng(0); folds = []; sv1 = []
for s in range(20):
    sel = np.random.default_rng(1000 + s).choice(EI.shape[1], 2000, replace=False)
    a = svd_c(np.asarray(EI[0])[sel])[1]; b = svd_c(np.asarray(EI[11])[sel])[1]
    folds.append(er(a)/er(b)); sv1.append((a[0]**2/(a**2).sum(), b[0]**2/(b**2).sum()))
out["immune_2000subsample_ER_fold_mean_range"] = [float(np.mean(folds)), float(min(folds)), float(max(folds))]
out["immune_2000subsample_SV1frac_L0_L11_mean"] = [float(np.mean([x[0] for x in sv1])), float(np.mean([x[1] for x in sv1]))]

# --- iter_0062 edge AUROC on immune (SV5-7, observed SVD, cosine of U coordinates), both-observed pairs ---
rows = list(csv.DictReader(open(CYC4 / "cycle1_edge_dataset.tsv"), delimiter="\t"))
pos = [(int(r["source_idx"]), int(r["target_idx"])) for r in rows if r["label"] == "1"]
neg = [(int(r["source_idx"]), int(r["target_idx"])) for r in rows if r["label"] == "0"]
posO = [p for p in pos if nzI[p[0]] and nzI[p[1]]]; negO = [p for p in neg if nzI[p[0]] and nzI[p[1]]]
out["immune_edge_rows"] = dict(n_rows=len(rows), n_pos=len(pos), n_neg=len(neg), n_pos_both_obs=len(posO), n_neg_both_obs=len(negO))
rp = {r: i for i, r in enumerate(obsI)}
def auc(a, b):
    a = np.asarray(a); b = np.asarray(b)
    from scipy.stats import rankdata
    r = rankdata(np.r_[a, b]); return float((r[:len(a)].sum() - len(a)*(len(a)+1)/2) / (len(a)*len(b)))
def cosrow(M, P):
    A = M[[rp[i] for i, _ in P]]; B = M[[rp[j] for _, j in P]]
    return (A*B).sum(1) / (np.linalg.norm(A, axis=1)*np.linalg.norm(B, axis=1))
edge = []
for l in range(12):
    Uo, So, _ = svd_c(np.asarray(EI[l])[obsI])
    M = Uo[:, 4:7]
    edge.append(auc(cosrow(M, posO), cosrow(M, negO)))
out["immune_edge_auroc_sv57_U_by_layer"] = edge
out["immune_edge_auroc_depth_spearman"] = list(map(float, spearmanr(range(12), edge)))

json.dump(out, open(Path(__file__).with_suffix(".json"), "w"), indent=1)
print(json.dumps(out, indent=1))
