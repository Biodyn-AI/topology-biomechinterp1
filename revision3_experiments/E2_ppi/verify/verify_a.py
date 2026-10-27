"""Independent verifier for E2_ppi, part A (deterministic recomputation from raw inputs).

Does NOT import e2_common / stats3 / contexts. Re-derives from raw files:
  * gene axis from the h5ad var index, V = rows with non-zero embedding (local embedding copy),
  * STRING v12 physical (string_links_v12.txt.gz; checked to be the 1,477,610-line physical file) >= 700
    mapped with string_protein_info preferred_name, both ends in V -> compared with cache edges_string.npy,
  * SV coordinates via eigen-decomposition of the V-centred Gram matrix (not np.linalg.svd),
  * assortativity as np.corrcoef over both edge orientations,
  * max-T (two-sided) recomputed from the stored MS null draws with my own code,
  * expr_extraction.npy column alignment spot-check against the h5ad rows (extraction_cells) x obs_idx.
Writes verify/verify_a.json and verify/verify_a_cells.csv.
"""
import gzip, json, sys
from pathlib import Path
import numpy as np, pandas as pd, h5py

HERE = Path(__file__).resolve().parent
E2 = HERE.parent
R3 = E2.parent
REPO = R3.parents[1]
CACHE = R3 / "cache"
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
H5 = {"lung": REPO / "single_cell_mechinterp/outputs/invariant_causal_edges/lung/processed.h5ad",
      "lung_external": REPO / "single_cell_mechinterp/outputs/invariant_causal_edges/external_lung/processed.h5ad",
      "kidney_A": REPO / "single_cell_mechinterp/outputs/tabula_sapiens_processed.h5ad",
      "kidney_B": REPO / "single_cell_mechinterp/outputs/invariant_causal_edges/kidney/processed.h5ad",
      "immune": REPO / "single_cell_mechinterp/outputs/tabula_sapiens_immune_subset_hpn_processed.h5ad"}
SDIR = REPO / "subproject_61_two_axes_revision/data"


def var_names(p):
    with h5py.File(p, "r") as f:
        v = f["var"]; key = v.attrs.get("_index", "_index")
        return [x.decode() if isinstance(x, bytes) else str(x) for x in v[key][:]]


def string_phys700():
    e2s = {}
    with gzip.open(SDIR / "string_protein_info_v12.txt.gz", "rt") as f:
        next(f)
        for line in f:
            p = line.split("\t", 2); e2s[p[0]] = p[1]
    out = {}
    nlines = 0
    with gzip.open(SDIR / "string_links_v12.txt.gz", "rt") as f:
        next(f)
        for line in f:
            nlines += 1
            a, b, s = line.split()
            if int(s) < 700: continue
            ga, gb = e2s.get(a), e2s.get(b)
            if ga is None or gb is None or ga == gb: continue
            out[frozenset((ga, gb))] = 1
    return out, nlines


def coords_eig(E, layer, obs):
    X = np.asarray(E[layer], np.float64)[obs]
    X = X - X.mean(0)
    w, U = np.linalg.eigh(X.T @ X)
    o = np.argsort(w)[::-1][:10]
    return X @ U[:, o]                 # = U_k S_k up to sign


def assort(x, a, b):
    return np.corrcoef(np.r_[x[a], x[b]], np.r_[x[b], x[a]])[0, 1]


def maxt_two_sided(obs, null):
    mu = null.mean(0); sd = null.std(0, ddof=1)
    zo = np.abs((obs - mu) / sd); zn = np.abs((null - mu) / sd)
    mx = zn.max(1); N = null.shape[0]
    padj = np.array([(np.count_nonzero(mx >= z) + 1) / (N + 1) for z in zo])
    praw = np.array([(np.count_nonzero(zn[:, j] >= zo[j]) + 1) / (N + 1) for j in range(zo.size)])
    return (obs - mu) / sd, praw, padj, mx


def main():
    res = {}; rows = []
    phys, nlines = string_phys700()
    res["string_file_lines_excl_header"] = nlines
    res["string_phys700_symbol_pairs"] = len(phys)
    N_PRIM = {"lung": 10000, "lung_external": 10000, "kidney_B": 2000, "immune": 2000, "kidney_A": 2000}
    for ctx in N_PRIM:
        r = {}
        genes_all = var_names(H5[ctx])
        E = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
        assert E.shape[1] == len(genes_all)
        nz0 = np.linalg.norm(np.asarray(E[0], np.float64), axis=1) > 0
        nz11 = np.linalg.norm(np.asarray(E[11], np.float64), axis=1) > 0
        obs = np.flatnonzero(nz0 & nz11)
        r["V_equal_all_layers"] = bool(np.array_equal(nz0, nz11))
        cg = json.loads((CACHE / ctx / "obs_genes.json").read_text())
        ci = np.load(CACHE / ctx / "obs_idx.npy")
        r["obs_idx_match"] = bool(np.array_equal(ci, obs))
        r["obs_genes_match_h5ad"] = bool([genes_all[i] for i in obs] == cg)
        r["n_V"] = int(obs.size)
        V = [genes_all[i] for i in obs]; pos = {g: k for k, g in enumerate(V)}
        if len(pos) != len(V): r["duplicate_symbols_in_V"] = len(V) - len(pos)
        mine = set()
        for p in phys:
            a, b = tuple(p)
            if a in pos and b in pos: mine.add((min(pos[a], pos[b]), max(pos[a], pos[b])))
        ce = np.load(CACHE / ctx / "edges_string.npy")
        cset = {(min(a, b), max(a, b)) for a, b in ce.tolist()}
        r["edges_mine"] = len(mine); r["edges_cache"] = len(cset)
        r["edges_sym_diff"] = len(mine ^ cset)
        ed = np.array(sorted(mine)); a, b = ed[:, 0], ed[:, 1]
        # coordinates + assortativity
        Rm = np.empty(120)
        for L in range(12):
            C = coords_eig(E, L, obs)
            for k in range(10):
                Rm[L * 10 + k] = assort(C[:, k], a, b)
        df = pd.read_csv(E2 / f"h1_primary_{ctx}.csv")
        r["max_abs_diff_r_vs_csv"] = float(np.abs(df.r.values - Rm).max())
        # max-T from the stored MS draws
        z = np.load(E2 / "_ckpt" / f"{ctx}_ms_{N_PRIM[ctx]}.npz")
        assert z["null"].shape == (N_PRIM[ctx], 120)
        zz, praw, padj, mx = maxt_two_sided(z["obs"], z["null"])
        r["N_null_draws"] = int(z["null"].shape[0])
        r["max_abs_diff_z"] = float(np.abs(zz - df.z.values).max())
        r["max_abs_diff_padj"] = float(np.abs(padj - df.p_adj_maxT.values).max())
        r["max_abs_diff_praw"] = float(np.abs(praw - df.p_raw.values).max())
        r["family_min_padj_mine"] = float(padj.min())
        r["n_cells_padj_lt_05_mine"] = int((padj < 0.05).sum())
        r["maxT_null_q95"] = float(np.quantile(mx, 0.95)); r["maxT_null_max"] = float(mx.max())
        r["min_abs_z_obs"] = float(np.abs(zz).min())
        t = int(np.argmax(np.abs(zz)))
        r["top_cell"] = f"L{t // 10} SV{t % 10 + 1}"; r["top_r_mine"] = float(Rm[t]); r["top_z_mine"] = float(zz[t])
        r["median_r_mine"] = float(np.median(Rm))
        # null preserves degree? spot check: degrees of stored chain are not stored; check obs r equals mine
        r["stored_obs_vs_mine"] = float(np.abs(z["obs"] - Rm).max())
        for j in range(120):
            rows.append(dict(context=ctx, layer=j // 10, sv=j % 10 + 1, r_mine=Rm[j], r_csv=df.r.values[j],
                             z_mine=zz[j], z_csv=df.z.values[j], padj_mine=padj[j], padj_csv=df.p_adj_maxT.values[j]))
        # expression alignment spot check (extraction cells x obs_idx) on 20 random genes
        cells = np.load(CACHE / ctx / "extraction_cells.npy")
        X_ext = np.load(CACHE / ctx / "expr_extraction.npy", mmap_mode="r")
        rng = np.random.default_rng(7); gs = np.sort(rng.choice(obs.size, 20, replace=False))
        with h5py.File(H5[ctx], "r") as f:
            Xg = f["X"]
            if isinstance(Xg, h5py.Group):
                ip = Xg["indptr"][:]; ind = Xg["indices"]; dat = Xg["data"]
                M = np.zeros((cells.size, gs.size), np.float32); colmap = {int(obs[g]): k for k, g in enumerate(gs)}
                for rI, cI in enumerate(cells):
                    s, e = ip[cI], ip[cI + 1]; ii = ind[s:e]; vv = dat[s:e]
                    for c_, v_ in zip(ii, vv):
                        k = colmap.get(int(c_))
                        if k is not None: M[rI, k] = v_
            else:
                M = np.asarray(Xg[np.sort(cells), :][:, obs[gs]], np.float32)
        r["expr_extraction_alignment_maxdiff"] = float(np.abs(M - np.asarray(X_ext[:, gs])).max())
        res[ctx] = r; print(ctx, json.dumps(r), flush=True)
    pd.DataFrame(rows).to_csv(HERE / "verify_a_cells.csv", index=False)
    (HERE / "verify_a.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
