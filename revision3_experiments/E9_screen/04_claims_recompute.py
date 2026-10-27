"""E9 / A6 step 5 -- recompute preprint headline numbers under each claim's OWN convention on its OWN matrix.

Each block reproduces what the screen iteration actually executed (gene mapping, gene set, SVD convention,
estimator), as read from that iteration's script / method string, and then -- where informative -- the same
estimator under the CONVENTION3 observed-gene convention (labelled 'observed_variant'; these are NOT the
round-3 analyses, only a sanity contrast). Output: claims_recompute.json (all numbers used in
E9_claims_ledger.csv).
"""
from __future__ import annotations
import json, pickle, time
from collections import defaultdict
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, spearmanr, mannwhitneyu, linregress
import e9_lib as L

OUT = L.HERE
R = {}
rng_master = np.random.default_rng(20261026)


def js(x):
    if isinstance(x, dict): return {str(k): js(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [js(v) for v in x]
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, np.ndarray): return js(x.tolist())
    return x


# =====================================================================================================
# kidney_A (cycle1_main) blocks
# =====================================================================================================
embK, genesK = L.X.load_context("kidney_A")
nzK = L.nonzero_rows(embK)
true209 = L.named_genes_true()                    # gene -> true row
mis209 = L.named_genes_misaligned()               # gene -> row 0..208
names209 = sorted(true209)

# ---- C: counts behind '209 named / 14 OOV / 195 in-vocab' ----
zero_true = sorted(g for g in names209 if not nzK[true209[g]])
zero_mis = sorted(g for g in names209 if not nzK[mis209[g]])
R["gene_counts_kidney_A"] = dict(
    G=int(embK.shape[1]), n_observed=int(nzK.sum()), n_named=len(names209),
    n_named_zero_true_mapping=len(zero_true), named_zero_true_mapping=zero_true,
    n_named_observed_true_mapping=len(names209) - len(zero_true),
    n_named_zero_misaligned_mapping=len(zero_mis), named_zero_misaligned_mapping=zero_mis,
    n_named_nonzero_misaligned_mapping=len(names209) - len(zero_mis),
    misaligned_rows_range=[min(mis209.values()), max(mis209.values())],
    n_labels_on_correct_row_under_misaligned_mapping=int(sum(true209[g] == mis209[g] for g in names209)),
    gene_list_txt_lines=sum(1 for _ in open(L.CYC1 / "gene_list.txt")),
    gene_list_txt_nonblank=sum(1 for l in open(L.CYC1 / "gene_list.txt") if l.strip()))

# ---- STRING / TRRUST pair counts on the 209 set ----
sc = json.load(open(L.ITER / "iter_0012/string_ppi_named_genes.json"))
nameset_up = {g.upper() for g in names209}
def string_pairs(thr):
    seen = set(); out = []
    for e in sc:
        a, b, s = e["preferredName_A"].upper(), e["preferredName_B"].upper(), e["score"]
        if a in nameset_up and b in nameset_up and s >= thr:
            k = tuple(sorted((a, b)))
            if k not in seen: seen.add(k); out.append(k)
    return out
str07 = string_pairs(0.7)
c15 = json.load(open(L.ITER / "iter_0015/string_ppi_score04_cache.json"))["pairs"]
str04 = [(p["g1"], p["g2"], p["score"]) for p in c15]
loc = L.X.load_string_genomewide(400)
loc07 = [k for k, s in loc.items() if s >= 700 and k[0] in true209 and k[1] in true209]
loc04 = [k for k, s in loc.items() if s >= 400 and k[0] in true209 and k[1] in true209]
obs_true = {g for g in names209 if nzK[true209[g]]}
tr_rows = L.trrust_raw_rows()
tr_rows_209 = [(a.upper(), b.upper()) for a, b, m in tr_rows if a.upper() in nameset_up and b.upper() in nameset_up and a.upper() != b.upper()]
tr_pairs_209 = sorted(set(tr_rows_209))
R["pair_counts_209"] = dict(
    string07_api_cache_dedup=len(str07), string04_iter15_cache=len(str04),
    string07_local_v12_genomewide=len(loc07), string04_local_v12_genomewide=len(loc04),
    string07_api_both_observed_true=int(sum(a in {g.upper() for g in obs_true} and b in {g.upper() for g in obs_true} for a, b in str07)),
    trrust_rows_209_iter10_style=len(tr_rows_209), trrust_distinct_pairs_209=len(tr_pairs_209),
    trrust_distinct_pairs_both_observed_true=int(sum(a in {g.upper() for g in obs_true} and b in {g.upper() for g in obs_true} for a, b in tr_pairs_209)),
    edge_dataset_positives=int(sum(r["label"] for r in L.edge_dataset(L.CYC1))))

# ---- 209-subset SVD (true rows): SV1 variance fraction, PR, adjacent-layer SV1 Spearman, co-pole ----
idx_t = np.array([true209[g] for g in names209])             # sorted gene order (as iter_0010)
rows_m = np.array([mis209[g] for g in names209])
keep_m = nzK[rows_m]                                          # '195 in-vocab' = misaligned nonzero rows
blk = {"true209": [], "mis195": [], "true187_observed": []}
sv1_proj = []
cop = {k: [] for k in ("SV1", "SV2", "SV3")}
cop_tr = None
names_up = np.array([g.upper() for g in names209])
pos = {g: i for i, g in enumerate(names_up)}
str07_idx = np.array([(pos[a], pos[b]) for a, b in str07])
for l in range(12):
    E = np.asarray(embK[l], dtype=np.float64)
    U, Sv, Vt = L.svd_rows(E[idx_t])
    blk["true209"].append(dict(layer=l, sv1_frac=float(Sv[0] ** 2 / (Sv ** 2).sum()), pr=L.pr(Sv), er_sq=L.er_sq(Sv)))
    sv1_proj.append(U[:, 0] * Sv[0])
    for k, nm in enumerate(("SV1", "SV2", "SV3")):
        pr_ = U[:, k] * Sv[k]
        o = np.argsort(-pr_); top = set(o[:52]); bot = set(o[-52:])
        cop[nm].append(float(np.mean([(a in top and b in top) or (a in bot and b in bot) for a, b in str07_idx])))
    if l == 11:   # iter_0010 H02: raw TRRUST rows, U[:,1] argsort, top = last 52, bottom = first 52
        o = np.argsort(U[:, 1]); top = set(names_up[o[-52:]]); bot = set(names_up[o[:52]])
        cop_tr = float(np.mean([(a in top and b in top) or (a in bot and b in bot) for a, b in tr_rows_209]))
        o1 = np.argsort(U[:, 0]); top1 = set(names_up[o1[-52:]]); bot1 = set(names_up[o1[:52]])
        cop_tr_sv1 = float(np.mean([(a in top1 and b in top1) or (a in bot1 and b in bot1) for a, b in tr_rows_209]))
    _, Sm, _ = L.svd_rows(E[rows_m[keep_m]])
    blk["mis195"].append(dict(layer=l, sv1_frac=float(Sm[0] ** 2 / (Sm ** 2).sum()), pr=L.pr(Sm), er_sq=L.er_sq(Sm)))
    _, So, _ = L.svd_rows(E[idx_t[nzK[idx_t]]])
    blk["true187_observed"].append(dict(layer=l, sv1_frac=float(So[0] ** 2 / (So ** 2).sum()), pr=L.pr(So), er_sq=L.er_sq(So)))
sv1_proj = np.array(sv1_proj)
adj_raw = [float(spearmanr(sv1_proj[i], sv1_proj[i + 1])[0]) for i in range(11)]
R["subset_spectra_kidney_A"] = dict(per_layer=blk,
    sv1_adjacent_layer_spearman_raw_sign=adj_raw, sv1_adjacent_mean_raw=float(np.mean(adj_raw)),
    sv1_adjacent_mean_abs=float(np.mean(np.abs(adj_raw))))
R["copole_kidney_A_209"] = dict(
    note="209 true rows, SVD per layer, proj = U[:,k]*S[k], poles = top/bottom 52, STRING>=0.7 pairs from the iter_0012 API cache",
    n_string07_pairs=len(str07), per_layer=cop, mean_over_layers={k: float(np.mean(v)) for k, v in cop.items()},
    trrust_iter10_L11_sv2=cop_tr, trrust_iter10_L11_sv1=cop_tr_sv1, n_trrust_rows=len(tr_rows_209))

# ---- full-matrix (zero rows) spectra on kidney_A: what the preprint says the 14.4-fold drop was computed on ----
fullK = []
for l in (0, 11):
    _, Sv, _ = L.svd_rows(np.asarray(embK[l], dtype=np.float64))
    fullK.append(dict(layer=l, er_sq=L.er_sq(Sv), er_s=L.S.effective_rank(Sv), sv1_frac=float(Sv[0] ** 2 / (Sv ** 2).sum())))
obsK = []
for l in (0, 11):
    _, Sv, _ = L.S.svd_observed(embK, l, np.flatnonzero(nzK))
    obsK.append(dict(layer=l, er_sq=L.er_sq(Sv), er_s=L.S.effective_rank(Sv), sv1_frac=float(Sv[0] ** 2 / (Sv ** 2).sum())))
R["full_spectra_kidney_A"] = dict(legacy_full_zero=fullK, observed_variant=obsK,
                                  fold_er_sq_legacy=fullK[0]["er_sq"] / fullK[1]["er_sq"])

# ---- SV1 GO enrichment (iter_0006 H02 / iter_0007 H03): full 4,803 SVD at L11, 209 true rows, quartiles ----
g2go = pickle.load(open("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/single_cell_mechinterp/data/perturb/gene2go_all.pkl", "rb"))
term2genes = defaultdict(set)
for g in names209:
    for t in g2go.get(g, set()): term2genes[t].add(g)
valid = [t for t, gs in term2genes.items() if 5 <= len(gs) <= 200]
def go_poles(order_genes, n25):
    top_q = set(order_genes[:n25]); bot_q = set(order_genes[-n25:]); res = {}
    for t in valid:
        tg = term2genes[t]; a = len(tg & top_q); c = len(tg & bot_q)
        if a + c == 0: continue
        res[t] = dict(top=fisher_exact([[a, len(top_q) - a], [c, len(bot_q) - c]], alternative="greater"),
                      bot=fisher_exact([[c, len(bot_q) - c], [a, len(top_q) - a]], alternative="greater"), a=a, c=c)
    return res
U, Sv, Vt = L.svd_rows(np.asarray(embK[11], dtype=np.float64))
proj = {g: U[true209[g], 0] for g in names209}
order = sorted(names209, key=lambda g: -proj[g])
res = go_poles(order, max(5, len(names209) // 4))
best_top = min(res.items(), key=lambda kv: kv[1]["top"][1])
best_bot = min(res.items(), key=lambda kv: kv[1]["bot"][1])
go_obs = {}
# observed-variant: same quartile test on the observed-gene SVD SV1 (sign fixed by stats3), 187 observed named genes
USo, So, Vto = L.S.svd_observed(embK, 11, np.flatnonzero(nzK))
rowpos = {r: i for i, r in enumerate(np.flatnonzero(nzK))}
obs_named = [g for g in names209 if nzK[true209[g]]]
projo = {g: USo[rowpos[true209[g]], 0] for g in obs_named}
ordero = sorted(obs_named, key=lambda g: -projo[g])
# zero-row genes under the full-matrix SVD sit at one point; how many of the 52+52 pole genes are zero rows?
top52, bot52 = set(order[:52]), set(order[-52:])
R["sv1_go_kidney_A_L11"] = dict(
    n_valid_terms=len(valid), n_tested_terms=len(res),
    GO0005615_top=dict(OR=res["GO:0005615"]["top"][0], p=res["GO:0005615"]["top"][1], a=res["GO:0005615"]["a"], c=res["GO:0005615"]["c"]) if "GO:0005615" in res else None,
    GO0005829_bottom=dict(OR=res["GO:0005829"]["bot"][0], p=res["GO:0005829"]["bot"][1]) if "GO:0005829" in res else None,
    best_top_term=[best_top[0], best_top[1]["top"][0], best_top[1]["top"][1]],
    best_bot_term=[best_bot[0], best_bot[1]["bot"][0], best_bot[1]["bot"][1]],
    n_zero_row_genes_in_top52=int(sum(g in zero_true for g in top52)), n_zero_row_genes_in_bottom52=int(sum(g in zero_true for g in bot52)),
    sv1_full_matrix_corr_with_observed_indicator=float(np.corrcoef(U[:, 0], nzK.astype(float))[0, 1]))
n25o = max(5, len(obs_named) // 4)
top_o, bot_o = set(ordero[:n25o]), set(ordero[-n25o:])
def fisher_term(t, A, B):
    tg = term2genes[t]; a = len(tg & A); c = len(tg & B)
    return dict(OR=fisher_exact([[a, len(A) - a], [c, len(B) - c]], alternative="greater")[0],
                p=fisher_exact([[a, len(A) - a], [c, len(B) - c]], alternative="greater")[1], a=a, c=c)
R["sv1_go_kidney_A_L11"]["observed_variant"] = dict(
    n_genes=len(obs_named), n25=n25o,
    GO0005615_pos_pole=fisher_term("GO:0005615", top_o, bot_o), GO0005615_neg_pole=fisher_term("GO:0005615", bot_o, top_o))

# ---- cell-type marker AUROC (iter_0023 H01): true mapping (blank-preserving gene_list), unit-sphere L2 ----
MK = {"T_cell": ["CD3G", "CD8A", "CD8B", "CD6", "RUNX3", "FOXP3", "PRF1", "GZMB", "CD247"],
      "B_cell": ["CD19", "MS4A1", "CD79A", "SPIB", "PAX5", "BLNK"],
      "macrophage": ["ALOX5", "PTPRC", "CD68", "CSF1R", "MRC1", "ITGAM"],
      "NK_cell": ["NCAM1", "KLRB1", "KLRD1", "NKG7", "GNLY"],
      "fibroblast": ["DCN", "COL1A1", "VIM", "FAP", "PDGFRA"],
      "epithelial": ["EPCAM", "KRT18", "KRT19", "CLDN4", "MUC1"],
      "endothelial": ["KDR", "ICAM1", "PECAM1", "CDH5", "ENG"],
      "myeloid": ["CD14", "LYZ", "FCGR3A", "S100A8", "S100A9"]}
def celltype_auroc(mapping, drop_zero=False, dtype=np.float64):
    gl = sorted(g for g in mapping if (not drop_zero or nzK[mapping[g]]))
    loc = {g: i for i, g in enumerate(gl)}
    sets = {ct: [loc[g] for g in m if g in loc] for ct, m in MK.items()}
    sets = {ct: v for ct, v in sets.items() if len(v) >= 2}
    w = [(a, b) for v in sets.values() for i, a in enumerate(v) for b in v[i + 1:]]
    cts = list(sets); cr = [(a, b) for i in range(len(cts)) for j in range(i + 1, len(cts)) for a in sets[cts[i]] for b in sets[cts[j]]]
    idx = np.array([mapping[g] for g in gl]); out = []
    for l in range(12):
        e = np.asarray(embK[l], dtype=dtype)[idx]; n = np.linalg.norm(e, axis=1, keepdims=True); n[n == 0] = 1; ue = e / n
        wa, wb = np.array(w).T; ca, cb = np.array(cr).T
        dw = np.sqrt(np.sum((ue[wa] - ue[wb]) ** 2, axis=1)); dc = np.sqrt(np.sum((ue[ca] - ue[cb]) ** 2, axis=1))
        u, p = mannwhitneyu(dw, dc, alternative="less"); out.append(1 - u / (len(dw) * len(dc)))
    lg = [np.log(a / (1 - a + 1e-6) + 1e-6) for a in out]
    sl = linregress(range(12), lg)
    return dict(cell_types={ct: [gl[i] for i in v] for ct, v in sets.items()}, n_within=len(w), n_cross=len(cr),
                auroc_per_layer=out, auroc_mean=float(np.mean(out)), slope_logit=float(sl.slope), r=float(sl.rvalue), p=float(sl.pvalue),
                zero_row_markers=[gl[i] for v in sets.values() for i in v if not nzK[mapping[gl[i]]]])
R["celltype_kidney_A"] = dict(own_convention_true_mapping=celltype_auroc(true209, dtype=np.float32),
                              same_in_float64=celltype_auroc(true209),
                              observed_variant_drop_zero_rows=celltype_auroc(true209, drop_zero=True))

# ---- GC-plasma angle (iter_0038 H02) and B-cell precision@10 (iter_0037 H01): misaligned vs true mapping ----
from sklearn.neighbors import NearestNeighbors
def invocab(mapping):
    gl = sorted(g for g in mapping if nzK[mapping[g]])
    return gl, {g: i for i, g in enumerate(gl)}, np.array([mapping[g] for g in gl])
def angle_series(mapping):
    gl, loc, idx = invocab(mapping)
    bc = [g for g in ["MS4A1", "CD19", "CD79A", "BLK", "PRDM1"] if g in loc]
    gc = [g for g in ["BATF", "SPIB", "BACH2"] if g in loc]; pl = [g for g in ["IRF4", "PRDM1"] if g in loc]
    ang = []
    for l in range(12):
        E = np.asarray(embK[l], dtype=np.float64)[idx]
        c_b = E[[loc[g] for g in bc]].mean(0); c_g = E[[loc[g] for g in gc]].mean(0); c_p = E[[loc[g] for g in pl]].mean(0)
        u = c_g - c_b; v = c_p - c_b
        ang.append(float(np.degrees(np.arccos(np.clip(u @ v / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-10), -1, 1)))))
    return dict(n_invocab=len(gl), bcell=bc, gc=gc, plasma=pl, angle_deg=ang)
def precision_series(mapping, layers=(0, 2, 5, 8, 11), n_boot=500):
    gl, loc, idx = invocab(mapping)
    mk = [g for g in ["MS4A1", "CD19", "CD79A", "CD79B", "BLK", "BANK1", "PRDM1"] if g in loc]
    li = np.array([loc[g] for g in mk]); out = {}
    for l in layers:
        E = np.asarray(embK[l], dtype=np.float64)[idx]
        nn = NearestNeighbors(n_neighbors=11, metric="euclidean", algorithm="ball_tree").fit(E)
        def prec(ix):
            s = set(ix.tolist()); _, ind = nn.kneighbors(E[ix]); return float(np.mean([sum(n in s for n in r[1:]) / 10 for r in ind]))
        obs = prec(li); rb = np.random.default_rng(l * 100)
        null = np.array([prec(rb.choice(len(gl), size=len(li), replace=False)) for _ in range(n_boot)])
        out[f"L{l}"] = dict(obs=obs, null_mean=float(null.mean()), z=float((obs - null.mean()) / (null.std() + 1e-10)))
    return dict(n_invocab=len(gl), markers=mk, layers=out)
R["bcell_kidney_A"] = dict(angle_own_convention_misaligned=angle_series(mis209), angle_true_mapping=angle_series(true209),
                           precision_own_convention_misaligned=precision_series(mis209), precision_true_mapping=precision_series(true209))

# =====================================================================================================
# immune (cycle4_immune_main) blocks
# =====================================================================================================
embI, genesI = L.X.load_context("immune")
nzI = L.nonzero_rows(embI); obsI = np.flatnonzero(nzI)
edI = L.edge_dataset(L.CYC4)
R["gene_counts_immune"] = dict(G=int(embI.shape[1]), n_observed=int(nzI.sum()), n_zero=int((~nzI).sum()))

# ---- iter_0046 H02: ER / SV1 fraction on 2,000 random genes of 4,941 (zero rows included) ----
def spectra_subsample(seed, l, n=2000, shuffle=False):
    rs = np.random.default_rng(seed); sel = rs.choice(embI.shape[1], size=n, replace=False)
    E = np.asarray(embI[l], dtype=np.float64)[sel]
    if shuffle:
        E = np.column_stack([rs.permutation(E[:, j]) for j in range(E.shape[1])])
    _, Sv, _ = L.svd_rows(E)
    return L.er_sq(Sv), float(Sv[0] ** 2 / (Sv ** 2).sum()), float((Sv[:10] ** 2).sum() / (Sv ** 2).sum()), int(nzI[sel].sum())
sub = {l: np.array([spectra_subsample(s, l) for s in range(50)]) for l in (0, 11)}
seed42 = {l: spectra_subsample(42, l) for l in (0, 11)}
shuf = np.array([spectra_subsample(s, 11, shuffle=True)[0] for s in range(10)])
fullI = {}
for l in (0, 11):
    _, Sv, _ = L.svd_rows(np.asarray(embI[l], dtype=np.float64))
    fullI[l] = (L.er_sq(Sv), float(Sv[0] ** 2 / (Sv ** 2).sum()))
obsIv = {}
for l in (0, 11):
    _, Sv, _ = L.S.svd_observed(embI, l, obsI)
    obsIv[l] = dict(er_sq=L.er_sq(Sv), er_s=L.S.effective_rank(Sv), sv1_frac=float(Sv[0] ** 2 / (Sv ** 2).sum()))
R["spectra_immune_iter46"] = dict(
    subsample50={f"L{l}": dict(er_mean=float(v[:, 0].mean()), er_min=float(v[:, 0].min()), er_max=float(v[:, 0].max()),
                                sv1_mean=float(v[:, 1].mean()), sv1_min=float(v[:, 1].min()), sv1_max=float(v[:, 1].max()),
                                top10_mean=float(v[:, 2].mean()), n_nonzero_in_sample_mean=float(v[:, 3].mean())) for l, v in sub.items()},
    subsample_seed42={f"L{l}": dict(er_sq=v[0], sv1_frac=v[1], top10=v[2]) for l, v in seed42.items()},
    fold_er_subsample_mean=float(sub[0][:, 0].mean() / sub[11][:, 0].mean()),
    feature_shuffle_L11_er_mean=float(shuf.mean()),
    full4941_zero_rows={f"L{l}": dict(er_sq=v[0], sv1_frac=v[1]) for l, v in fullI.items()},
    fold_er_full4941=float(fullI[0][0] / fullI[11][0]),
    observed_variant={f"L{l}": v for l, v in obsIv.items()},
    fold_er_sq_observed=float(obsIv[0]["er_sq"] / obsIv[11]["er_sq"]))

# ---- iter_0045 H02: TwoNN on 2,000 random genes of 4,941 (estimator copied from iter_0042) ----
from scipy.spatial.distance import cdist
def twonn_id(Xm):
    d = cdist(Xm, Xm); np.fill_diagonal(d, np.inf); sd = np.sort(d, axis=1)
    mu = sd[:, 1] / (sd[:, 0] + 1e-12); mu = mu[mu > 1]
    if len(mu) < 3: return float("nan")
    ms = np.sort(mu); F = np.arange(1, len(ms) + 1) / len(ms)
    lm = np.log(ms); l1 = -np.log(1 - F + 1e-12)
    return float(np.sum(lm * l1) / np.sum(lm ** 2))
tw = {}
for seed in (42, 0, 1):
    rs = np.random.default_rng(seed); sel = rs.choice(embI.shape[1], size=2000, replace=False)
    tw[seed] = {f"L{l}": twonn_id(np.asarray(embI[l], dtype=np.float64)[sel]) for l in (0, 11)}
R["twonn_immune_iter45"] = dict(note="estimator = iter_0042 twonn_id (OLS through origin on -log(1-F) vs log mu, mu>1 only); "
                                    "sampling seed of iter_0045 unknown (no script survives)", by_seed=tw)

# ---- labels on the immune 2,039: 295 observed named genes, TF = edge-dataset source, target-only = rest ----
named_I = {}
for r in edI:
    named_I.setdefault(r["source"], r["source_idx"]); named_I.setdefault(r["target"], r["target_idx"])
sources = {r["source"] for r in edI}
obs_named_I = sorted(g for g in named_I if nzI[named_I[g]])
tf_I = [g for g in obs_named_I if g in sources]; tgt_I = [g for g in obs_named_I if g not in sources]
rowposI = {r: i for i, r in enumerate(obsI)}
pos_pairs = [(r["source_idx"], r["target_idx"]) for r in edI if r["label"] == 1 and nzI[r["source_idx"]] and nzI[r["target_idx"]]]
neg_pairs_obs = [(r["source_idx"], r["target_idx"]) for r in edI if r["label"] == 0 and nzI[r["source_idx"]] and nzI[r["target_idx"]]]
all_other_rows = [(r["source_idx"], r["target_idx"]) for r in edI if not (r["label"] == 1 and nzI[r["source_idx"]] and nzI[r["target_idx"]])]
tr_modes = defaultdict(set)
for a, b, m in L.trrust_raw_rows(): tr_modes[(a, b)].add(m)
pos_named = [(genesI[a], genesI[b]) for a, b in pos_pairs]
def mode_class(p):
    m = tr_modes.get(p, set())
    if m == {"Activation"}: return "activation_only"
    if m == {"Repression"}: return "repression_only"
    return "mixed_or_unknown"
mc = pd.Series([mode_class(p) for p in pos_named]).value_counts().to_dict()
def mode_class2(p):   # 'any Activation row' / 'any Repression row' / neither (alternative reading of 270/141/178)
    m = tr_modes.get(p, set())
    if "Activation" in m and "Repression" not in m: return "activation_no_repression"
    if "Repression" in m and "Activation" not in m: return "repression_no_activation"
    if "Activation" in m and "Repression" in m: return "both"
    return "unknown_only"
mc2 = pd.Series([mode_class2(p) for p in pos_named]).value_counts().to_dict()
mc3 = dict(any_activation=int(sum("Activation" in tr_modes.get(p, set()) for p in pos_named)),
           any_repression=int(sum("Repression" in tr_modes.get(p, set()) for p in pos_named)))
R["labels_immune"] = dict(n_named=len(named_I), n_named_observed=len(obs_named_I), n_tf=len(tf_I), n_target_only=len(tgt_I),
                          n_pos_both_observed=len(pos_pairs), n_neg_both_observed=len(neg_pairs_obs),
                          n_rows_other_than_valid_pos=len(all_other_rows), n_edge_rows=len(edI),
                          n_pos_rows=int(sum(r["label"] for r in edI)), trrust_mode_split_of_valid_pos=mc,
                          trrust_mode_split_alt=mc2, trrust_mode_any=mc3)

# ---- iter_0056 H01 / iter_0055 H02: TF vs target-only AUROC from SV2-4, SV5-7, joint SV2-7 ----
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
ytf = np.r_[np.ones(len(tf_I)), np.zeros(len(tgt_I))]
lab_rows = np.array([rowposI[named_I[g]] for g in tf_I + tgt_I])
def cv_auc(F, y, seed=42):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    m = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=1000))
    fold = []
    for tr, te in skf.split(F, y):
        m.fit(F[tr], y[tr]); fold.append(roc_auc_score(y[te], m.predict_proba(F[te])[:, 1]))
    oof = cross_val_predict(m, F, y, cv=skf, method="predict_proba")[:, 1]
    return float(np.mean(fold)), float(roc_auc_score(y, oof))
tfres = []
for l in range(12):
    US, Sv, Vt = L.S.svd_observed(embI, l, obsI)
    F24 = US[lab_rows, 1:4]; F57 = US[lab_rows, 4:7]; F27 = US[lab_rows, 1:7]
    r = dict(layer=l)
    for nm, F in (("sv24", F24), ("sv57", F57), ("joint", F27)):
        r[nm + "_foldmean"], r[nm + "_oof"] = cv_auc(F, ytf)
    tfres.append(r)
tfdf = pd.DataFrame(tfres)
# 100-permutation label-shuffle null for the joint model at L3 (the preprint's max layer) -- floor check only
US, Sv, Vt = L.S.svd_observed(embI, 3, obsI); F27 = US[lab_rows, 1:7]
rp = np.random.default_rng(7); nullj = np.array([cv_auc(F27, rp.permutation(ytf))[0] for _ in range(100)])
R["tf_target_immune_iter56"] = dict(per_layer=tfres, joint_foldmean_mean=float(tfdf.joint_foldmean.mean()),
                                    joint_foldmean_max=float(tfdf.joint_foldmean.max()), joint_foldmean_argmax=int(tfdf.joint_foldmean.idxmax()),
                                    joint_oof_mean=float(tfdf.joint_oof.mean()),
                                    joint_beats_both_layers=int(((tfdf.joint_foldmean > tfdf.sv24_foldmean) & (tfdf.joint_foldmean > tfdf.sv57_foldmean)).sum()),
                                    sv57_L0_3_mean=float(tfdf.sv57_foldmean[:4].mean()), sv24_L4_8_mean=float(tfdf.sv24_foldmean[4:9].mean()),
                                    L3_null100_mean=float(nullj.mean()), L3_null100_max=float(nullj.max()),
                                    L3_perm_p=float((1 + (nullj >= tfdf.joint_foldmean[3]).sum()) / 101),
                                    p_floor_100=1 / 101)
# cross-seed (iter_0057 H01): same model on immune_seed43 / immune_seed44 with their own observed sets / labels
cs = {}
for ctx, run in (("immune_seed43", "cycle4_immune_seed43"), ("immune_seed44", "cycle4_immune_seed44")):
    e2, g2 = L.X.load_context(ctx); nz2 = L.nonzero_rows(e2); ob2 = np.flatnonzero(nz2); rp2 = {r: i for i, r in enumerate(ob2)}
    ed2 = L.edge_dataset(L.X.S38 / run); nm2 = {}
    for r in ed2: nm2.setdefault(r["source"], r["source_idx"]); nm2.setdefault(r["target"], r["target_idx"])
    src2 = {r["source"] for r in ed2}; on2 = sorted(g for g in nm2 if nz2[nm2[g]])
    t2 = [g for g in on2 if g in src2]; n2 = [g for g in on2 if g not in src2]
    y2 = np.r_[np.ones(len(t2)), np.zeros(len(n2))]; lr2 = np.array([rp2[nm2[g]] for g in t2 + n2])
    a = [cv_auc(L.S.svd_observed(e2, l, ob2)[0][lr2, 1:7], y2)[0] for l in range(12)]
    cs[ctx] = dict(n_observed=int(nz2.sum()), n_tf=len(t2), n_target_only=len(n2), joint_mean=float(np.mean(a)),
                   joint_min=float(np.min(a)), joint_max=float(np.max(a)))
R["tf_target_cross_seed_iter57"] = cs

# ---- iter_0062 H01/H02: edge-level cosine AUROC in SV5-7 / SV2-4 ----
def cos_rows(A, B):
    return np.sum(A * B, 1) / (np.linalg.norm(A, axis=1) * np.linalg.norm(B, axis=1) + 1e-12)
P = np.array([(rowposI[a], rowposI[b]) for a, b in pos_pairs]); Nn = np.array([(rowposI[a], rowposI[b]) for a, b in neg_pairs_obs])
edge = []
for l in range(12):
    US, Sv, Vt = L.S.svd_observed(embI, l, obsI)
    Xc_all = np.asarray(embI[l], dtype=np.float64) - np.asarray(embI[l], dtype=np.float64)[obsI].mean(0)
    r = dict(layer=l)
    for nm, sl in (("sv57", slice(4, 7)), ("sv24", slice(1, 4))):
        for coord, M in (("US", US), ("U", US / Sv)):
            sp = cos_rows(M[P[:, 0], sl], M[P[:, 1], sl]); sn = cos_rows(M[Nn[:, 0], sl], M[Nn[:, 1], sl])
            r[f"{nm}_{coord}_neg1523"] = L.S.auroc(sp, sn)
        # variant: all 4,941 rows projected (centred by the observed mean), negatives = every other edge row (2,351)
        Pall = Xc_all @ Vt[sl].T
        pa = np.array(pos_pairs); oa = np.array(all_other_rows)
        r[f"{nm}_projall_neg2351"] = L.S.auroc(cos_rows(Pall[pa[:, 0]], Pall[pa[:, 1]]), cos_rows(Pall[oa[:, 0]], Pall[oa[:, 1]]))
    edge.append(r)
ed = pd.DataFrame(edge)
R["edge_immune_iter62"] = dict(per_layer=edge,
    rho_layer_sv57_US_neg1523=list(spearmanr(range(12), ed.sv57_US_neg1523)),
    rho_layer_sv57_U_neg1523=list(spearmanr(range(12), ed.sv57_U_neg1523)),
    sv24_U_neg1523_mean=float(ed.sv24_U_neg1523.mean()),
    rho_layer_sv57_projall_neg2351=list(spearmanr(range(12), ed.sv57_projall_neg2351)),
    sv24_US_neg1523_mean=float(ed.sv24_US_neg1523.mean()), sv24_projall_neg2351_mean=float(ed.sv24_projall_neg2351.mean()))

# ---- iter_0051 H01 / iter_0052 H03: distance residualised on full-embedding cosine (the screen's 'co-expression') ----
def rbc_resid(l, sl, coord="U"):
    US, Sv, Vt = L.S.svd_observed(embI, l, obsI); M = US / Sv if coord == "U" else US
    E = np.asarray(embI[l], dtype=np.float64)[obsI]
    allp = np.vstack([P, Nn]); y = np.r_[np.ones(len(P)), np.zeros(len(Nn))]
    d = np.linalg.norm(M[allp[:, 0], sl] - M[allp[:, 1], sl], axis=1)
    c = cos_rows(E[allp[:, 0]], E[allp[:, 1]])
    A = np.column_stack([np.ones_like(c), c]); beta = np.linalg.lstsq(A, d, rcond=None)[0]; res = d - A @ beta
    u, p = mannwhitneyu(res[y == 1], res[y == 0], alternative="less")
    return float(1 - 2 * u / ((y == 1).sum() * (y == 0).sum())), float(p)
rr = {}
for coord in ("U", "US"):
    rr[coord] = dict(sv24=[rbc_resid(l, slice(1, 4), coord) for l in range(12)], sv57=[rbc_resid(l, slice(4, 7), coord) for l in range(12)])
    rr[coord]["sv24_n_layers_p05"] = int(sum(p < 0.05 for _, p in rr[coord]["sv24"]))
    rr[coord]["sv57_n_layers_p05"] = int(sum(p < 0.05 for _, p in rr[coord]["sv57"]))
R["residualised_immune_iter51_52"] = dict(note="rbc = 1 - 2U/(n1 n0), one-sided pos<neg; covariate = cosine of the full 512-d embeddings "
                                               "(the screen's 'co-expression'), NOT expression correlation across cells", by_coord=rr)

# ---- iter_0042 H03: BATF / BACH2 rank to B-cell centroid over all 4,941 rows (zero rows included) ----
gi = {g: i for i, g in enumerate(genesI)}
bidx = [gi[g] for g in ["MS4A1", "CD19", "CD79A", "BLK", "PRDM1"] if g in gi]
gc = {}
for g in ("BATF", "BACH2"):
    rk, dp, rk_obs = [], [], []
    for l in range(12):
        E = np.asarray(embI[l], dtype=np.float64); c = E[bidx].mean(0); d = np.linalg.norm(E - c, axis=1)
        rk.append(int((d < d[gi[g]]).sum())); dp.append(float(np.linalg.norm(E[gi[g]] - E[gi["PAX5"]])))
        rk_obs.append(int((d[nzI] < d[gi[g]]).sum()))
    gc[g] = dict(rank_all_rows=rk, rank_observed_rows_only=rk_obs, dist_to_PAX5=dp, rho_rank_layer=list(spearmanr(range(12), rk)),
                 rho_rank_layer_observed=list(spearmanr(range(12), rk_obs)))
pax = []
for l in range(12):
    E = np.asarray(embI[l], dtype=np.float64); c = E[bidx].mean(0); d = np.linalg.norm(E - c, axis=1); pax.append(int((d < d[gi["PAX5"]]).sum()))
gc["PAX5_rank_all_rows"] = pax
gc["bcell_markers_present"] = [genesI[i] for i in bidx]
gc["bcell_markers_observed"] = [genesI[i] for i in bidx if nzI[i]]
R["gc_immune_iter42"] = gc

# ---- confidence gradient p-value (iter_0015 H02 'p=1.4e-24' vs preprint 'p=0.017') ----
R["spearman_n5_exact"] = dict(rho=1.0, p_two_sided_exact=2 / 120, p_one_sided_exact=1 / 120,
                              scipy_p=float(spearmanr([1, 2, 3, 4, 5], [1.00, 1.48, 2.09, 3.18, 5.02])[1]))

# ---- local STRING file used by round 3 (size check) ----
import gzip
with gzip.open(L.X.STRING_LINKS, "rt") as f: nlines = sum(1 for _ in f) - 1
R["string_local_file"] = dict(path=str(L.X.STRING_LINKS), n_links=nlines,
    example_api_vs_local=[(a, b, s_api, loc.get(tuple(sorted((a, b))))) for a, b, s_api in
                          [("DCN", "COL1A1", 0.923), ("DCN", "TNF", 0.9), ("DCN", "TLR2", 0.99)]])

# ---- screen size ----
recs = pd.read_csv(OUT / "audit_409_records.csv")
R["screen_size"] = dict(records_iter_1_63=int((recs.iter_num <= 63).sum()),
                        iterations_with_json_1_63=int(recs[recs.iter_num <= 63].iteration.nunique()),
                        iteration_dirs_1_63=int(sum((L.ITER / f"iter_{i:04d}").is_dir() for i in range(1, 64))),
                        json_families_1_63=sorted(recs[recs.iter_num <= 63].family.unique().tolist()))
(OUT / "claims_recompute.json").write_text(json.dumps(js(R), indent=1))
print("done")
