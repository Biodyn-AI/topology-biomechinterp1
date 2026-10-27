"""Build the shared per-context cache every round-3 experiment reads.

cache/<context>/
  obs_genes.json      observed genes (non-zero embedding), in the order of obs_idx
  obs_idx.npy         their row indices in the context's embedding matrix
  ncells.npy          number of the 256 sampled cells whose input contained each observed gene
  expr.npy            (cells x observed genes) float32, log-normalised X, <= 3,000 cells from the
                      matching cell population (tissue / donor / cell type), fixed seed
  edges_trrust.npy    (m x 2) TRRUST TF->target pairs, both ends observed (indices into obs_genes)
  edges_string.npy    (m x 2) STRING v12 combined >= 700, undirected, both ends observed
  string_score.npy    combined score of each STRING edge
  edges_dorothea_chip_new.npy  DoRothEA ChIP-seq A-C pairs NOT in TRRUST (independent evidence type)
  any_regulatory.npy  (m x 2) every pair in TRRUST U DoRothEA(all levels, both files): used to
                      exclude undocumented-but-reported pairs from conservative negative sets
  meta.json           counts and provenance
Run once: python prep_cache.py
"""
import json, sys
from pathlib import Path
import numpy as np, h5py, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import contexts as X
import stats3 as S

OUT = Path(__file__).parent / "cache"
SEED = 20261026
SUBSET = {   # obs filter for expression of the matching cell population
    "immune_TSP14": ("donor_id", ["TSP14"]), "immune_TSP21": ("donor_id", ["TSP21"]),
    "immune_TSP25": ("donor_id", ["TSP25"]),
    "immune_Bcell": ("cell_type", ["B cell"]),
    "immune_CD4T": ("cell_type", ["CD4-positive, alpha-beta T cell"]),
    "immune_CD8T": ("cell_type", ["CD8-positive, alpha-beta T cell"]),
}

def obs_col(f, col):
    g = f["obs"][col]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]])
        codes = g["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "")
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in g[:]])

def load_expr(h5, genes_all, obs_idx, ctx, rng, max_cells=3000):
    with h5py.File(h5, "r") as f:
        n_cells = int(f["X"].attrs["shape"][0]) if isinstance(f["X"], h5py.Group) else f["X"].shape[0]
        cells = np.arange(n_cells)
        if ctx in SUBSET:
            col, vals = SUBSET[ctx]; cells = np.flatnonzero(np.isin(obs_col(f, col), vals))
        sel = np.sort(rng.choice(cells, min(max_cells, cells.size), replace=False))
        Xg = f["X"]
        colpos = np.full(len(genes_all), -1, dtype=np.int64); colpos[obs_idx] = np.arange(obs_idx.size)
        M = np.zeros((sel.size, obs_idx.size), dtype=np.float32)
        if isinstance(Xg, h5py.Group):
            indptr = Xg["indptr"][:]; data = Xg["data"]; indices = Xg["indices"]
            for r, ci in enumerate(sel):
                a, b = indptr[ci], indptr[ci + 1]
                idx = indices[a:b]; val = data[a:b]; cp = colpos[idx]; k = cp >= 0
                M[r, cp[k]] = val[k]
        else:
            M = np.asarray(Xg[sel, :][:, obs_idx], dtype=np.float32)
    return M, int(cells.size)

def main():
    t, trrust = X.load_trrust()
    stg = X.load_string_genomewide(700)
    chip = set(X.load_dorothea(chipseq=True, levels=("A", "B", "C")))
    allreg = set(trrust) | set(X.load_dorothea(False, ("A", "B", "C", "D"))) | set(X.load_dorothea(True, ("A", "B", "C", "D")))
    trs = set(trrust)
    expr_shared = {}
    for ctx, c in X.CONTEXTS.items():
        d = OUT / ctx; d.mkdir(parents=True, exist_ok=True)
        emb, genes = X.load_context(ctx)
        obs = S.observed_mask(emb); obs_idx = np.flatnonzero(obs); og = [genes[i] for i in obs_idx]
        pos = {g: k for k, g in enumerate(og)}
        cnt_p = c["emb"].parent / "layer_gene_embedding_count.npy"
        ncell = np.load(cnt_p)[0][obs_idx] if cnt_p.exists() else np.full(obs_idx.size, -1)
        e_tr = np.array([(pos[a], pos[b]) for a, b in trrust if a in pos and b in pos], dtype=np.int64).reshape(-1, 2)
        st = [((pos[a], pos[b]), s) for (a, b), s in stg.items() if a in pos and b in pos]
        e_st = np.array([p for p, _ in st], dtype=np.int64).reshape(-1, 2)
        s_st = np.array([s for _, s in st], dtype=np.int32)
        e_ch = np.array([(pos[a], pos[b]) for a, b in sorted(chip - trs) if a in pos and b in pos], dtype=np.int64).reshape(-1, 2)
        e_any = np.array([(pos[a], pos[b]) for a, b in sorted(allreg) if a in pos and b in pos], dtype=np.int64).reshape(-1, 2)
        key = (str(c["h5ad"]), SUBSET.get(ctx))
        rng = np.random.default_rng(SEED)
        M, pop = load_expr(c["h5ad"], genes, obs_idx, ctx, rng)
        np.save(d / "obs_idx.npy", obs_idx); np.save(d / "ncells.npy", ncell); np.save(d / "expr.npy", M)
        np.save(d / "edges_trrust.npy", e_tr); np.save(d / "edges_string.npy", e_st); np.save(d / "string_score.npy", s_st)
        np.save(d / "edges_dorothea_chip_new.npy", e_ch); np.save(d / "any_regulatory.npy", e_any)
        (d / "obs_genes.json").write_text(json.dumps(og))
        meta = dict(context=ctx, tissue=c["tissue"], group=c["group"], note=c["note"], emb=str(c["emb"]),
                    h5ad=str(c["h5ad"]), expr_population=SUBSET.get(ctx, "whole h5ad"), population_size=pop,
                    G=len(genes), n_observed=int(obs_idx.size), n_expr_cells=int(M.shape[0]),
                    n_trrust=int(len(e_tr)), n_string700=int(len(e_st)), n_dorothea_chip_new=int(len(e_ch)),
                    n_any_regulatory=int(len(e_any)))
        (d / "meta.json").write_text(json.dumps(meta, indent=1))
        print(json.dumps({k: meta[k] for k in ("context", "n_observed", "n_expr_cells", "population_size",
                                               "n_trrust", "n_string700", "n_dorothea_chip_new")}), flush=True)

if __name__ == "__main__":
    main()
