"""A2 supplement: hypergeometric tail for observed cell overlaps; kidney_A vs kidney_B X agreement detail."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import hypergeom
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE))
import contexts as X, a2_lib as L
ov = pd.read_csv(HERE / "overlaps.csv")
res = json.loads((HERE / "results.json").read_text())
CACHE = HERE.parent / "cache"
# population size of b
popsize = {c: v["n_obs_extraction_h5ad"] for c, v in res["contexts"].items()}
p = []
for r in ov.itertuples():
    N = popsize[r.b]; K = r.n_a_in_b_population; n = r.n_b
    p.append(float(hypergeom.sf(r.n_shared - 1, N, K, n)) if K > 0 else float("nan"))
ov["p_hypergeom_ge_observed"] = p
ov.to_csv(HERE / "overlaps.csv", index=False)
# kidney_A vs kidney_B X values on shared cells x shared genes (all 8000 shared cells would be slow; use kidney_B drawn cells)
gA = L.var_names(X.H5_KIDNEY_A); gB = L.var_names(X.H5_KIDNEY_B)
sh = sorted(set(gA) & set(gB)); pa = {g: i for i, g in enumerate(gA)}; pb = {g: i for i, g in enumerate(gB)}
iA = np.array([pa[g] for g in sh]); iB = np.array([pb[g] for g in sh])
posA = {n: i for i, n in enumerate(L.obs_names(X.H5_KIDNEY_A))}; namesB = L.obs_names(X.H5_KIDNEY_B)
cellsB = np.load(CACHE / "kidney_B" / "extraction_cells.npy")
RA, RB = L.H5Rows(X.H5_KIDNEY_A), L.H5Rows(X.H5_KIDNEY_B)
A = np.stack([RA.row_dense(posA[namesB[r]])[iA] for r in cellsB]); B = np.stack([RB.row_dense(int(r))[iB] for r in cellsB])
d = np.abs(A - B)
out = dict(n_cells=int(cellsB.size), n_shared_genes=len(sh), frac_entries_exactly_equal=float((A == B).mean()),
           max_abs_diff=float(d.max()), median_abs_diff_nonzero=float(np.median(d[(A > 0) | (B > 0)])),
           same_nonzero_pattern=bool(np.array_equal(A > 0, B > 0)), pearson=float(np.corrcoef(A.ravel(), B.ravel())[0, 1]),
           max_A=float(A.max()), max_B=float(B.max()))
res["kidney_A_vs_B_X_on_kidney_B_drawn_cells"] = out
res["overlaps"] = ov.to_dict(orient="records")
(HERE / "results.json").write_text(json.dumps(res, indent=1, default=str))
print(ov[["a", "b", "n_shared", "expected_shared_if_independent", "p_hypergeom_ge_observed"]].to_string()); print(out)

# whole-file provenance of each h5ad (tissue / donor across ALL cells of the file, not only the drawn ones)
files = {"immune": X.H5_IMMUNE, "kidney_A": X.H5_KIDNEY_A, "kidney_B": X.H5_KIDNEY_B, "lung": X.H5_LUNG, "lung_external": X.H5_EXT_LUNG}
wf = {}
for k, h in files.items():
    wf[k] = {"path": str(h)}
    for col in ("tissue_in_publication", "donor_id", "compartment", "assay"):
        v = L.obs_col(h, col)
        wf[k][col] = None if v is None else {str(a): int(b) for a, b in pd.Series(v).value_counts().items()}
res["whole_file_obs"] = wf
(HERE / "results.json").write_text(json.dumps(res, indent=1, default=str))
for k, v in wf.items(): print(k, {c: v[c] for c in ("tissue_in_publication", "donor_id", "assay")})
