"""Per-context tokenisation mask: token_masks/<context>.npy, bool (extraction cells x observed genes),
True where the gene entered that cell's scGPT input (X > 0 AND within the cell's top-512 by expression).
Rows follow cache/<context>/extraction_cells.npy, columns cache/<context>/obs_genes.json.
Column sums equal cache/<context>/nobs_reproduced.npy (asserted)."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE))
import contexts as X, a2_lib as L
from a2_axis_verification import PROV, MAX_GENES
CACHE = HERE.parent / "cache"; OUT = HERE / "token_masks"; OUT.mkdir(exist_ok=True)
summ = {}
for ctx, c in X.CONTEXTS.items():
    d = CACHE / ctx; P = PROV[ctx]
    obs_idx = np.load(d / "obs_idx.npy"); col = np.full(len(L.var_names(c["h5ad"])), -1); col[obs_idx] = np.arange(obs_idx.size)
    R = L.H5Rows(P["h5"]); cells_ext = L.sample_cells(R.shape[0], P["max_cells"], P["seed"])
    M = np.zeros((cells_ext.size, obs_idx.size), dtype=bool)
    det_any = np.zeros(col.size, dtype=bool)
    for k, r in enumerate(cells_ext):
        v = R.row_dense(int(r)); det_any |= v > 0
        gi = L.tokenised_genes(v, MAX_GENES); M[k, col[gi]] = True
    assert (col[np.concatenate([L.tokenised_genes(R.row_dense(int(r)), MAX_GENES) for r in cells_ext])] >= 0).all()
    assert np.array_equal(M.sum(0), np.load(d / "nobs_reproduced.npy"))
    E = np.load(d / "expr_extraction.npy")
    assert not (M & ~(E > 0)).any()
    np.save(OUT / f"{ctx}.npy", M)
    summ[ctx] = dict(tokens=int(M.sum()), detected_entries=int((E > 0).sum()),
                     n_genes_detected_in_any_extraction_cell=int(det_any.sum()), n_observed=int(obs_idx.size),
                     n_genes_detected_but_never_tokenised=int(det_any.sum() - obs_idx.size),
                     frac_detected_entries_tokenised=float(M.sum() / (E > 0).sum()))
    print(ctx, summ[ctx], flush=True)
(HERE / "token_mask_summary.json").write_text(json.dumps(summ, indent=1))
