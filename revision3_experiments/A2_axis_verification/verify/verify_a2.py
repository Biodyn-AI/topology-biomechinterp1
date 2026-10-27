"""Independent verifier for A2 (does not import a2_lib or a2_axis_verification).

Uses anndata.read_h5ad + adata[keep_idx].copy() + per-row toarray(), i.e. the same calls the extractor
(run_layerwise_geometry_audit.py) and ScGPTDataset.__getitem__ make, instead of A2's raw-h5py CSR reader.
Writes verify/verify_results.json.
"""
import json, sys
from pathlib import Path
import numpy as np, scipy.sparse as sp, anndata as ad
from scipy.stats import hypergeom
HERE = Path(__file__).resolve().parent
R3 = HERE.parents[1]
sys.path.insert(0, str(R3))
import contexts as X
CACHE = R3 / "cache"
D38 = X.REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/data"
# my own provenance table (seed, max_cells, extraction file, subset column/value)
SUBF = {"immune_TSP14": (D38/"immune_donor_subsets/TSP14.h5ad", "donor_id", "TSP14"),
        "immune_TSP21": (D38/"immune_donor_subsets/TSP21.h5ad", "donor_id", "TSP21"),
        "immune_TSP25": (D38/"immune_donor_subsets/TSP25.h5ad", "donor_id", "TSP25"),
        "immune_Bcell": (D38/"immune_celltype_subsets/b_cell.h5ad", "cell_type", "B cell"),
        "immune_CD4T": (D38/"immune_celltype_subsets/cd4_positive_alpha_beta_t_cell.h5ad", "cell_type", "CD4-positive, alpha-beta T cell"),
        "immune_CD8T": (D38/"immune_celltype_subsets/cd8_positive_alpha_beta_t_cell.h5ad", "cell_type", "CD8-positive, alpha-beta T cell")}
SEED = {k: 42 for k in X.CONTEXTS}; SEED.update(immune_seed43=43, immune_seed44=44, kidney_B=67,
        kidney_ckpt_whole_human=87, kidney_ckpt_kidney=87, kidney_ckpt_brain=87)
NC = {k: (192 if (k == "kidney_B" or k.startswith("kidney_ckpt")) else 256) for k in X.CONTEXTS}
MG = 512

def tokens(values):
    mask = values > 0
    gi = np.arange(values.size, dtype=np.int64)[mask]; ev = values[mask]
    if gi.size > 1:
        o = np.argsort(-ev); gi = gi[o]
    return gi[:MG] if gi.size > MG else gi

def row_tokens(adata):
    out = []
    for i in range(adata.n_obs):
        row = adata.X[i]
        v = row.toarray().ravel() if sp.issparse(row) else np.asarray(row).ravel()
        out.append(tokens(v))
    return out

_cache = {}
def load(p):
    p = str(p)
    if p not in _cache: _cache[p] = ad.read_h5ad(p)
    return _cache[p]

res = {"contexts": {}}
drawn_names = {}
for ctx, c in X.CONTEXTS.items():
    ext_path = SUBF[ctx][0] if ctx in SUBF else c["h5ad"]
    A = load(ext_path)
    rng = np.random.default_rng(SEED[ctx])
    keep = np.sort(rng.choice(A.n_obs, size=NC[ctx], replace=False)) if A.n_obs > NC[ctx] else np.arange(A.n_obs)
    sub = A[keep].copy()
    toks = row_tokens(sub)
    cnt = np.zeros(sub.n_vars, dtype=np.int64)
    for t in toks: np.add.at(cnt, t, 1)
    ref = np.load(c["emb"].parent / "layer_gene_embedding_count.npy")
    emb = np.load(c["emb"], mmap_mode="r")
    Actx = load(c["h5ad"])
    var_same = list(sub.var_names.astype(str)) == list(Actx.var_names.astype(str))
    # cache checks
    d = CACHE / ctx
    obs_idx = np.load(d / "obs_idx.npy"); obs_genes = json.loads((d / "obs_genes.json").read_text())
    nz = np.flatnonzero(np.abs(np.asarray(emb[0], dtype=np.float64)).sum(1) > 0)
    names_ctx = np.asarray(Actx.obs_names.astype(str))
    pos = {n: i for i, n in enumerate(names_ctx)}
    my_cells = np.array([pos[n] for n in np.asarray(sub.obs_names.astype(str))], dtype=np.int64)
    a2_cells = np.load(d / "extraction_cells.npy")
    a2_nobs = np.load(d / "nobs_reproduced.npy")
    a2_expr = np.load(d / "expr_extraction.npy")
    Xc = Actx.X[my_cells]; Xc = Xc.toarray() if sp.issparse(Xc) else np.asarray(Xc)
    my_expr = Xc[:, obs_idx].astype(np.float32)
    # token mask
    col = np.full(sub.n_vars, -1); col[obs_idx] = np.arange(obs_idx.size)
    M = np.zeros((len(toks), obs_idx.size), dtype=bool)
    for k, t in enumerate(toks): M[k, col[t]] = True
    a2_mask = np.load(HERE.parent / "token_masks" / f"{ctx}.npy")
    info = dict(seed=SEED[ctx], n_cells=int(keep.size), n_obs_file=int(A.n_obs),
                count_exact_layer0=bool(np.array_equal(cnt, ref[0])),
                count_exact_all_layers=bool(all(np.array_equal(cnt, ref[l]) for l in range(ref.shape[0]))),
                n_genes_differ=int((cnt != ref[0]).sum()),
                ext_var_equals_ctx_var=var_same, emb_G=int(emb.shape[1]), G=int(sub.n_vars),
                obs_idx_equals_nonzero_emb_rows=bool(np.array_equal(nz, obs_idx)),
                obs_genes_equal_var_at_obs_idx=[str(g) for g in Actx.var_names[obs_idx]] == obs_genes,
                a2_extraction_cells_equal=bool(np.array_equal(a2_cells, my_cells)),
                a2_nobs_equal=bool(np.array_equal(a2_nobs, cnt[obs_idx])),
                a2_expr_equal=bool(a2_expr.dtype == np.float32 and np.array_equal(a2_expr, my_expr)),
                a2_token_mask_equal=bool(np.array_equal(a2_mask, M)),
                frac_detected_entries_tokenised=float(M.sum() / (my_expr > 0).sum()),
                n_detected_never_tokenised=int(((sub.X.toarray() if sp.issparse(sub.X) else sub.X) > 0).any(0).sum() - obs_idx.size))
    for colname in ("tissue_in_publication", "donor_id", "cell_type", "broad_cell_class", "tissue"):
        if colname in Actx.obs:
            info[colname] = {str(k): int(v) for k, v in Actx.obs[colname].iloc[my_cells].astype(str).value_counts().items()}
    if ctx in SUBF:
        _, cn, val = SUBF[ctx]
        par = names_ctx[np.asarray(Actx.obs[cn].astype(str) == val)]
        info["subset_equals_parent_rows_in_order"] = bool(np.array_equal(par, np.asarray(A.obs_names.astype(str))))
        info["subset_unique_values"] = sorted(map(str, A.obs[cn].astype(str).unique()))
    drawn_names[ctx] = set(names_ctx[my_cells])
    res["contexts"][ctx] = info
    print(ctx, json.dumps(info), flush=True)

# seed search for the checkpoint runs, using precomputed per-row tokens of the whole file
KB = load(X.H5_KIDNEY_B)
allt = row_tokens(KB)
ref = np.load(X.CONTEXTS["kidney_ckpt_whole_human"]["emb"].parent / "layer_gene_embedding_count.npy")[0]
refs = {k: np.load(X.CONTEXTS[k]["emb"].parent / "layer_gene_embedding_count.npy")[0] for k in ("kidney_ckpt_kidney", "kidney_ckpt_brain", "kidney_B")}
res["ckpt_count_vectors_identical"] = bool(all(np.array_equal(ref, refs[k]) for k in ("kidney_ckpt_kidney", "kidney_ckpt_brain")))
hits = {}
for mc in (192, 256):
    h = []
    for s in range(0, 2001):
        kk = np.sort(np.random.default_rng(s).choice(KB.n_obs, mc, replace=False))
        cc = np.bincount(np.concatenate([allt[i] for i in kk]), minlength=KB.n_vars)
        if np.array_equal(cc, ref): h.append(s)
    hits[mc] = h
res["ckpt_seed_search_0_2000"] = hits
print("ckpt seeds", hits, flush=True)

# overlaps
def ov(a, b, popb):
    k = len(drawn_names[a] & popb); n = len(drawn_names[b]); sh = len(drawn_names[a] & drawn_names[b])
    return dict(a=a, b=b, n_shared=sh, K=k, expected=k * n / len(popb), p=float(hypergeom.sf(sh - 1, len(popb), k, n)))
imm_pop = set(load(X.H5_IMMUNE).obs_names.astype(str))
res["overlaps"] = []
for b in ["immune_seed43", "immune_seed44", "immune_TSP14", "immune_TSP21", "immune_TSP25", "immune_Bcell", "immune_CD4T", "immune_CD8T"]:
    popb = set(load(SUBF[b][0]).obs_names.astype(str)) if b in SUBF else imm_pop
    res["overlaps"].append(ov("immune", b, popb))
res["overlaps"].append(ov("immune_seed43", "immune_seed44", imm_pop))
res["overlaps"].append(ov("kidney_A", "kidney_B", set(KB.obs_names.astype(str))))
res["overlaps"].append(ov("kidney_B", "kidney_ckpt_whole_human", set(KB.obs_names.astype(str))))
res["overlaps"].append(ov("lung", "lung_external", set(load(X.H5_EXT_LUNG).obs_names.astype(str))))
KA = load(X.H5_KIDNEY_A)
res["kidneyB_file_subset_of_kidneyA_file"] = bool(set(KB.obs_names.astype(str)) <= set(KA.obs_names.astype(str)))
res["kidneyA_file_donors"] = {str(k): int(v) for k, v in KA.obs["donor_id"].astype(str).value_counts().items()}
res["kidneyA_file_tissue_in_publication"] = {str(k): int(v) for k, v in KA.obs["tissue_in_publication"].astype(str).value_counts().items()}
for r in res["overlaps"]: print(r)
(HERE / "verify_results.json").write_text(json.dumps(res, indent=1))
print("ALL count exact:", all(v["count_exact_all_layers"] for v in res["contexts"].values()))
