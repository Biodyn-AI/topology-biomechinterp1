"""Independent re-computation of E11 part-1 headline numbers from raw inputs.

Does NOT import stats3, e11_common or the cache gene lists for the computation: observed mask from embedding norms,
gene axis from the h5ad var index, own centred SVD, own AUROC (scipy rankdata), scipy fisher_exact, own max-T
(Westfall-Young single-step, z-standardised by each test's own null, as CONVENTION3 / the spec requires),
own permutation seed. The cache obs_genes/obs_idx are only *compared* against.
Output: verify/verify_part1.json
"""
import json, pickle, sys
from pathlib import Path
import numpy as np, h5py
from scipy.stats import rankdata, fisher_exact, spearmanr

HERE = Path(__file__).resolve().parent; R3 = HERE.parents[1]; REPO = R3.parents[1]
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
H5 = {"lung": REPO / "single_cell_mechinterp/outputs/invariant_causal_edges/lung/processed.h5ad",
      "lung_external": REPO / "single_cell_mechinterp/outputs/invariant_causal_edges/external_lung/processed.h5ad",
      "kidney_B": REPO / "single_cell_mechinterp/outputs/invariant_causal_edges/kidney/processed.h5ad",
      "immune": REPO / "single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad",
      "kidney_A": REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad"}
G2GO = pickle.load(open(REPO / "single_cell_mechinterp/data/perturb/gene2go_all.pkl", "rb"))
EXT, CYT = "GO:0005615", "GO:0005829"
N_PERM = 2000


def genes_of(p):
    with h5py.File(p, "r") as f:
        return [g.decode() if isinstance(g, bytes) else str(g) for g in f["var"]["_index"][:]]


def coords(ctx):
    E = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    genes = genes_of(H5[ctx])
    obs = np.linalg.norm(np.asarray(E[0], np.float64), axis=1) > 0
    idx = np.flatnonzero(obs)
    C = np.empty((idx.size, 144))
    for L in range(12):
        Xl = np.asarray(E[L], np.float64)[idx]; Xl -= Xl.mean(0)
        U, s, _ = np.linalg.svd(Xl, full_matrices=False)
        C[:, L * 12:(L + 1) * 12] = U[:, :12] * s[:12]
    return C, [genes[i] for i in idx], idx


def auc_cols(R, y):
    n1 = y.sum(); n0 = y.size - n1
    return (R[y].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)


out = {}
coord_cache = {}
for ctx in ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]:
    C, g, idx = coords(ctx)
    cache_genes = json.loads((R3 / "cache" / ctx / "obs_genes.json").read_text())
    cache_idx = np.load(R3 / "cache" / ctx / "obs_idx.npy")
    y_ext = np.array([EXT in G2GO.get(x, ()) for x in g]); y_cyt = np.array([CYT in G2GO.get(x, ()) for x in g])
    R = rankdata(C, axis=0)
    a_ext = auc_cols(R, y_ext); a_cyt = auc_cols(R, y_cyt)
    sf_ext = 0.5 + np.abs(a_ext - 0.5); sf_cyt = 0.5 + np.abs(a_cyt - 0.5)
    best = int(np.argmax(sf_ext)); L, k = divmod(best, 12)
    # Fisher at K = floor(0.27 n), scipy sample OR, two-sided
    n = len(g); K = int(0.27 * n); o = np.argsort(C[:, best], kind="mergesort")
    top, bot = o[-K:], o[:K]
    fis = {}
    for nm, pole in (("top", top), ("bottom", bot)):
        a = int(y_ext[pole].sum()); b = K - a; c = int(y_ext.sum()) - a; d = n - K - c
        OR, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        fis[nm] = dict(a=a, OR_sample=float(OR), OR_haldane=float((a + .5) * (d + .5) / ((b + .5) * (c + .5))), p=float(p))
    # max-T over 144 cells x 2 primary terms, AUROC statistic, own seed
    rng = np.random.default_rng(777)
    tobs = np.r_[sf_ext, sf_cyt]
    null = np.empty((N_PERM, 288)); Y = np.c_[y_ext, y_cyt]
    for b in range(N_PERM):
        Yp = Y[rng.permutation(n)]
        null[b] = np.r_[0.5 + np.abs(auc_cols(R, Yp[:, 0]) - .5), 0.5 + np.abs(auc_cols(R, Yp[:, 1]) - .5)]
    mu, sd = null.mean(0), null.std(0, ddof=1)
    z = (tobs - mu) / sd; zmax = ((null - mu) / sd).max(1)
    padj = lambda j: float((np.sum(zmax >= z[j]) + 1) / (N_PERM + 1))
    l11 = 11 * 12 + 0
    out[ctx] = dict(
        n_V=n, cache_genes_match=bool(cache_genes == g), cache_idx_match=bool(np.array_equal(cache_idx, idx)),
        n_ext=int(y_ext.sum()), n_cyt=int(y_cyt.sum()),
        best_cell=f"L{L} SV{k + 1}", best_ext_auroc_signfree=float(sf_ext[best]), best_ext_auroc_signed=float(a_ext[best]),
        best_cyt_auroc_signed_same_cell=float(a_cyt[best]), best_z=float(z[best]), best_padj=padj(best),
        maxnull_z_q95=float(np.quantile(zmax, .95)), fisher_best_cell=fis,
        L11_SV1_ext_auroc_signfree=float(sf_ext[l11]), L11_SV1_ext_padj=padj(l11),
        L0_SV1_spearman_lognobs=None, n_cells_padj_lt_0p05=int(sum(padj(j) < 0.05 for j in range(288))))
    coord_cache[ctx] = dict(zip(g, C[:, best]))
    print(ctx, json.dumps(out[ctx], default=float), flush=True)

# lung vs other contexts: Pearson r of best-cell coordinates on matched genes
for ctx in ["lung_external", "kidney_B", "immune"]:
    a, b = coord_cache["lung"], coord_cache[ctx]; cm = sorted(set(a) & set(b))
    out[ctx]["r_with_lung_best_cell"] = float(np.corrcoef([a[x] for x in cm], [b[x] for x in cm])[0, 1])
    out[ctx]["n_matched"] = len(cm)
(HERE / "verify_part1.json").write_text(json.dumps(out, indent=1, default=float))
print(json.dumps({c: out[c].get("r_with_lung_best_cell") for c in out}))
