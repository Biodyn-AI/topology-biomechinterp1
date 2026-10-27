"""Verifier check 1 (independent code): input alignment and the direction-test pair set from raw files.

(a) expr_extraction.npy columns = obs_genes (by NAME via anndata var_names, not via obs_idx), rows = extraction cells;
(b) T2 expr.npy gene axis agrees with T1 (per-gene mean correlation; a misaligned axis would give r ~ 0);
(c) obs_genes == h5ad var_names[obs_idx]; embedding rows obs_idx are non-zero, other rows zero (V definition);
(d) direction pairs rebuilt from raw TRRUST / DoRothEA TSVs with pandas (not via cache any_regulatory.npy);
(e) regulator / target sets rebuilt from h3_pairs + Lambert file directly.
Writes verify/v1_inputs.json.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

V = Path(__file__).resolve().parent; E = V.parent; R3 = E.parent; CACHE = R3 / "cache"
sys.path.insert(0, str(R3))
import contexts as X  # only for the h5ad / embedding paths
import anndata as ad

NET = X.NET
LAMBERT = X.REPO / "subproject_61_two_axes_revision/data/humantfs_v1.01.txt"
out = {}
tr = pd.read_csv(NET / "trrust_human.tsv", sep="\t", header=None, names=["TF", "target", "mode", "pmid"]).astype(str)
trp = {(a, b) for a, b in zip(tr.TF, tr.target) if a != b}
dm = pd.read_csv(NET / "dorothea_human.tsv", sep="\t"); dc = pd.read_csv(NET / "dorothea_chipseq_human.tsv", sep="\t")
reg_any = trp | set(zip(dm.source.astype(str), dm.target.astype(str))) | set(zip(dc.source.astype(str), dc.target.astype(str)))
lam = pd.read_csv(LAMBERT, sep="\t")
lam_tf = set(lam.loc[lam["Is TF?"].astype(str).str.strip() == "Yes", "HGNC symbol"].astype(str))

for ctx in ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]:
    d = CACHE / ctx; o = {}
    genes = json.loads((d / "obs_genes.json").read_text()); obs_idx = np.load(d / "obs_idx.npy")
    h5 = X.CONTEXTS[ctx]["h5ad"]
    a = ad.read_h5ad(h5, backed="r"); vn = np.array(a.var_names.astype(str))
    o["obs_genes_eq_varnames_obs_idx"] = bool(list(vn[obs_idx]) == genes)
    emb = np.load(Path.home() / f"biomechinterp-local/subproject41_r3/emb/{ctx}.npy", mmap_mode="r")
    nz = np.linalg.norm(np.asarray(emb[0], np.float64), axis=1) > 0
    o["V_equals_nonzero_rows_L0"] = bool(np.array_equal(np.flatnonzero(nz), obs_idx))
    # (a) rebuild expression of the extraction cells by gene NAME
    cells = np.load(d / "extraction_cells.npy")
    col = {g: k for k, g in enumerate(vn)}; cols = np.array([col[g] for g in genes])
    Xs = a.X[np.sort(cells)] if not np.all(np.diff(cells) > 0) else a.X[cells]
    Xs = Xs.toarray() if hasattr(Xs, "toarray") else np.asarray(Xs)
    if not np.all(np.diff(cells) > 0):
        Xs = Xs[np.argsort(np.argsort(cells))]
    mine = Xs[:, cols].astype(np.float32)
    t1 = np.load(d / "expr_extraction.npy")
    o["expr_extraction_max_abs_diff_vs_h5ad_by_name"] = float(np.abs(mine - t1).max())
    # (b) T2 axis agreement
    t2 = np.load(d / "expr.npy", mmap_mode="r")
    m1 = t1.mean(0); m2 = np.asarray(t2).mean(0)
    o["T1_T2_gene_mean_pearson"] = float(np.corrcoef(m1, m2)[0, 1])
    o["T1_T2_gene_mean_pearson_shuffled_axis"] = float(np.corrcoef(m1, np.random.default_rng(0).permutation(m2))[0, 1])
    o["T1_cells"] = int(t1.shape[0]); o["T2_cells"] = int(t2.shape[0])
    # (e) regulators / targets
    h = np.load(d / "h3_pairs.npz")
    regs_mine = sorted(set(h["tf"].tolist()) | {k for k, g in enumerate(genes) if g in lam_tf})
    o["n_regulators_mine"] = len(regs_mine)
    rset = set(regs_mine); pos = {g: k for k, g in enumerate(genes)}
    # (d) direction pairs from raw files
    dirp = sorted((pos[a_], pos[b_]) for a_, b_ in trp if a_ in pos and b_ in pos and pos[a_] in rset and pos[b_] in rset
                  and (b_, a_) not in reg_any)
    o["n_direction_pairs_mine"] = len(dirp)
    tg_all = np.unique(np.concatenate([h["tg"], np.array(dirp, dtype=np.int64).ravel()]))
    o["n_targets_mine"] = int(tg_all.size); o["n_h3_targets_mine"] = int(np.unique(h["tg"]).size)
    # compare with what E4b fitted / tested
    z = np.load(E / "partial" / f"{ctx}_T1_genie3_s0.npz")
    o["regs_equal_E4b"] = bool(np.array_equal(z["regs"], np.array(regs_mine)))
    o["targets_equal_E4b"] = bool(np.array_equal(z["targets"], tg_all))
    sys.path.insert(0, str(E)); import e4b_common as C
    o["dir_pairs_equal_E4b"] = bool(np.array_equal(C.load_sets(ctx)["dir_pairs"], np.array(dirp).reshape(-1, 2)))
    # h3 pair sanity: no TF == target, positives are TRRUST, candidates are in no regulatory set either way
    tfs = [genes[i] for i in h["tf"]]; tgs = [genes[i] for i in h["tg"]]; y = h["y"].astype(bool)
    o["h3_self_pairs"] = int(sum(a_ == b_ for a_, b_ in zip(tfs, tgs)))
    o["h3_pos_all_in_trrust"] = bool(all((a_, b_) in trp for a_, b_, yy in zip(tfs, tgs, y) if yy))
    o["h3_neg_in_any_reg_either_dir"] = int(sum(((a_, b_) in reg_any) or ((b_, a_) in reg_any) for a_, b_, yy in zip(tfs, tgs, y) if not yy))
    o["h3_n_pairs"] = int(y.size); o["h3_n_pos"] = int(y.sum()); o["h3_n_tf"] = int(np.unique(h["tf"]).size)
    out[ctx] = o; print(ctx, json.dumps(o), flush=True)
(V / "v1_inputs.json").write_text(json.dumps(out, indent=1))
