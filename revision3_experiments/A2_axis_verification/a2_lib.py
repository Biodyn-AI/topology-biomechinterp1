"""Model-free re-implementation of the scGPT extraction sampler + tokeniser.

Mirrors, operation for operation, what
subproject_38/.../scripts/run_layerwise_geometry_audit.py and
single_cell_mechinterp/src/data/scgpt_dataset.py do before the model is called:

  rng = np.random.default_rng(seed)
  keep_idx = np.sort(rng.choice(n_obs, size=max_cells, replace=False))   # only if n_obs > max_cells
  for each kept cell (in keep_idx order):
      values = X[row].toarray().ravel()          # float32, h5ad X (log-normalised)
      mask = values > 0                          # include_zero=False
      gene_indices = arange(G)[mask]; expr = values[mask]
      order = np.argsort(-expr)                  # sort_by_expression=True, numpy default kind
      gene_indices = gene_indices[order][:max_genes]   # cls_token_id=None -> max_genes unchanged
      count[gene_indices] += 1                   # np.add.at(layer_count[l], valid_idx, 1)

ScGPTDataset raises if any var name is missing from the vocab, so no vocab filtering happens:
every var gene with X > 0 is a candidate token.
"""
from __future__ import annotations
import numpy as np, h5py
import scipy.sparse as sp


class H5Rows:
    """Row access to an h5ad CSR X without loading the whole matrix."""
    def __init__(self, path):
        self.path = str(path)
        self.f = h5py.File(self.path, "r")
        Xg = self.f["X"]
        assert isinstance(Xg, h5py.Group) and Xg.attrs["encoding-type"] == "csr_matrix", path
        self.shape = tuple(int(v) for v in Xg.attrs["shape"])
        self.indptr = Xg["indptr"][:]
        self.data = Xg["data"]; self.indices = Xg["indices"]
        self.dtype = self.data.dtype

    def row_dense(self, r):
        a, b = int(self.indptr[r]), int(self.indptr[r + 1])
        idx = self.indices[a:b]; val = self.data[a:b]
        # identical semantics to anndata row slicing + .toarray() (duplicates summed)
        m = sp.csr_matrix((val, idx, np.array([0, b - a])), shape=(1, self.shape[1]))
        return m.toarray().ravel()

    def rows_dense(self, rows):
        out = np.zeros((len(rows), self.shape[1]), dtype=self.dtype)
        for k, r in enumerate(rows):
            out[k] = self.row_dense(int(r))
        return out

    def close(self):
        self.f.close()


def var_names(f_or_path):
    own = not isinstance(f_or_path, h5py.File)
    f = h5py.File(f_or_path, "r") if own else f_or_path
    raw = f["var"]["_index"][:]
    if own: f.close()
    return [g.decode() if isinstance(g, bytes) else str(g) for g in raw]


def obs_names(path):
    with h5py.File(path, "r") as f:
        key = f["obs"].attrs.get("_index", "_index")
        raw = f["obs"][key][:]
    return np.array([g.decode() if isinstance(g, bytes) else str(g) for g in raw])


def obs_col(path, col):
    with h5py.File(path, "r") as f:
        if col not in f["obs"]:
            return None
        g = f["obs"][col]
        if isinstance(g, h5py.Group) and "categories" in g:
            cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]])
            codes = g["codes"][:]
            return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "<NA>")
        v = g[:]
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def sample_cells(n_obs, max_cells, seed):
    if max_cells and n_obs > max_cells:
        rng = np.random.default_rng(seed)
        return np.sort(rng.choice(n_obs, size=max_cells, replace=False))
    return np.arange(n_obs)


def tokenised_genes(values, max_genes=512):
    """Gene indices of one cell that enter the scGPT input (exact ScGPTDataset logic)."""
    gene_idx_all = np.arange(values.shape[0], dtype=np.int64)
    mask = values > 0
    gene_indices = gene_idx_all[mask]
    expr_values = values[mask]
    if gene_indices.size > 1:
        order = np.argsort(-expr_values)
        gene_indices = gene_indices[order]
        expr_values = expr_values[order]
    if gene_indices.size > max_genes > 0:
        gene_indices = gene_indices[:max_genes]
    return gene_indices


def tie_at_cutoff(values, max_genes=512):
    """True if the max_genes-th and (max_genes+1)-th largest positive values are equal (truncation
    then depends on argsort tie-breaking)."""
    v = values[values > 0]
    if v.size <= max_genes:
        return False
    s = np.sort(v)[::-1]
    return bool(s[max_genes - 1] == s[max_genes])


def count_vector(rows_reader, cells, max_genes=512, cache=None):
    G = rows_reader.shape[1]
    cnt = np.zeros(G, dtype=np.int64)
    for r in cells:
        r = int(r)
        if cache is not None and r in cache:
            gi = cache[r]
        else:
            gi = tokenised_genes(rows_reader.row_dense(r), max_genes)
            if cache is not None: cache[r] = gi
        np.add.at(cnt, gi, 1)
    return cnt
