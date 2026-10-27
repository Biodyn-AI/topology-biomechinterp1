"""A2 - gene-axis verification for every context in contexts.CONTEXTS.

For each context: reproduce the extraction's cell draw and tokenisation without the model,
compare the per-gene count vector with layer_gene_embedding_count.npy[0] (exact equality is the pass
criterion), verify the gene axis independently via the run's own cycle1_edge_dataset.tsv
(name -> index written by the extractor from adata.var_names), save the reproduced cells and their
expression to cache/<context>/, and tabulate tissue / donor / cell-type truth plus cell overlaps.

Writes (A2_axis_verification/): results.json, per_context.csv, obs_tables.csv, overlaps.csv,
a2_axis_verification.log (stdout). Writes to cache/<context>/ only:
extraction_cells.npy, nobs_reproduced.npy, expr_extraction.npy.
"""
from __future__ import annotations
import json, sys, platform
from pathlib import Path
import numpy as np, pandas as pd, h5py

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE))
import contexts as X
import a2_lib as L

CACHE = HERE.parent / "cache"
S38_ROOT = X.REPO / "subproject_38_geometric_residual_stream_interpretability"
S38_DATA = S38_ROOT / "implementation/data"
IT = X.S41

# ---------------------------------------------------------------------------------------------
# Extraction provenance. "documented" = where the seed / max_cells / input h5ad come from.
# max_genes = 512 everywhere (script default; every documented command passes --max-genes 512 or omits it).
# ---------------------------------------------------------------------------------------------
README38 = "subproject_38/implementation/README.md (Reproduction commands)"
PROV = {
    "immune": dict(h5=X.H5_IMMUNE, seed=42, max_cells=256, sub=None,
                   doc=README38 + ": 'Cross-domain (immune) main run' --seed 42 --max-cells 256"),
    "immune_seed43": dict(h5=X.H5_IMMUNE, seed=43, max_cells=256, sub=None,
                          doc="research_log.md Cycle 4 'Replicates: cycle4_immune_seed43' (seed from dir name; no command recorded)"),
    "immune_seed44": dict(h5=X.H5_IMMUNE, seed=44, max_cells=256, sub=None,
                          doc="research_log.md Cycle 4 'Replicates: cycle4_immune_seed44' (seed from dir name; no command recorded)"),
    "immune_TSP14": dict(h5=S38_DATA / "immune_donor_subsets/TSP14.h5ad", seed=42, max_cells=256, sub=("donor_id", "TSP14"),
                         doc="research_log.md Cycle 10 donor-specific h5ad subsets; seed 42 from dir name"),
    "immune_TSP21": dict(h5=S38_DATA / "immune_donor_subsets/TSP21.h5ad", seed=42, max_cells=256, sub=("donor_id", "TSP21"),
                         doc="research_log.md Cycle 10 donor-specific h5ad subsets; seed 42 from dir name"),
    "immune_TSP25": dict(h5=S38_DATA / "immune_donor_subsets/TSP25.h5ad", seed=42, max_cells=256, sub=("donor_id", "TSP25"),
                         doc="research_log.md Cycle 10 donor-specific h5ad subsets; seed 42 from dir name"),
    "immune_Bcell": dict(h5=S38_DATA / "immune_celltype_subsets/b_cell.h5ad", seed=42, max_cells=256, sub=("cell_type", "B cell"),
                         doc="research_log.md Cycle 11 make_obs_subsets.py; seed = script default 42 (no command recorded)"),
    "immune_CD4T": dict(h5=S38_DATA / "immune_celltype_subsets/cd4_positive_alpha_beta_t_cell.h5ad", seed=42, max_cells=256,
                        sub=("cell_type", "CD4-positive, alpha-beta T cell"),
                        doc="research_log.md Cycle 11 make_obs_subsets.py; seed = script default 42 (no command recorded)"),
    "immune_CD8T": dict(h5=S38_DATA / "immune_celltype_subsets/cd8_positive_alpha_beta_t_cell.h5ad", seed=42, max_cells=256,
                        sub=("cell_type", "CD8-positive, alpha-beta T cell"),
                        doc="research_log.md Cycle 11 make_obs_subsets.py; seed = script default 42 (no command recorded)"),
    "kidney_A": dict(h5=X.H5_KIDNEY_A, seed=42, max_cells=256, sub=None,
                     doc=README38 + ": 'Cycle 1 (seed 42, max_genes=512)', default --processed-h5ad tabula_sapiens_processed.h5ad"),
    "kidney_B": dict(h5=X.H5_KIDNEY_B, seed=67, max_cells=192, sub=None,
                     doc="iter_0067/executor_iteration_report.md Command Trace: --max-cells 192 --seed 67; log 'cells=192'"),
    "kidney_ckpt_whole_human": dict(h5=X.H5_KIDNEY_B, seed=87, max_cells=192, sub=None,
                                    doc="iter_0087 report command omits --seed/--max-cells; log says cells=192; seed 87 found by search (see results.json)"),
    "kidney_ckpt_kidney": dict(h5=X.H5_KIDNEY_B, seed=87, max_cells=192, sub=None,
                               doc="as kidney_ckpt_whole_human"),
    "kidney_ckpt_brain": dict(h5=X.H5_KIDNEY_B, seed=87, max_cells=192, sub=None,
                              doc="as kidney_ckpt_whole_human"),
    "lung": dict(h5=X.H5_LUNG, seed=42, max_cells=256, sub=None,
                 doc=README38 + ": 'Cross-domain (lung) main run' --seed 42 --max-cells 256"),
    "lung_external": dict(h5=X.H5_EXT_LUNG, seed=42, max_cells=256, sub=None,
                          doc="research_log.md Cycle 7 cycle7_external_lung_main (seed42 per Cycle 11 'seed42/43/44' naming; no command recorded)"),
}
MAX_GENES = 512
OBS_COLS = ["tissue", "tissue_in_publication", "compartment", "broad_cell_class", "donor_id", "cell_type"]
CACHE_NEW = ["extraction_cells.npy", "nobs_reproduced.npy", "expr_extraction.npy"]


def log(*a):
    print(*a, flush=True)


def edge_axis_check(run_dir, genes):
    """The extractor wrote source_idx/target_idx from adata.var_names; check against the CONTEXTS h5ad var."""
    p = run_dir / "cycle1_edge_dataset.tsv"
    if not p.exists():
        return dict(edge_file=False)
    e = pd.read_csv(p, sep="\t")
    pos = {g: i for i, g in enumerate(genes)}
    names = pd.concat([e.source, e.target]).astype(str).to_numpy()
    idx = pd.concat([e.source_idx, e.target_idx]).to_numpy()
    ok = np.array([pos.get(n, -1) == i for n, i in zip(names, idx)])
    return dict(edge_file=True, n_edge_rows=int(len(e)), n_gene_index_pairs=int(ok.size),
                n_distinct_genes=int(len(set(names))), n_index_mismatch=int((~ok).sum()))


def search_seed(R, ref, max_cells, seeds):
    cache = {}
    for s in seeds:
        c = L.count_vector(R, L.sample_cells(R.shape[0], max_cells, s), MAX_GENES, cache)
        if (c == ref).all():
            return int(s)
    return None


def main():
    log("python", sys.executable, "numpy", np.__version__, "h5py", h5py.__version__, platform.platform())
    res = dict(numpy=np.__version__, max_genes=MAX_GENES, contexts={}, overlaps={}, seed_search={})
    rows, obs_rows = [], []
    immune_obs = L.obs_names(X.H5_IMMUNE)
    imm_pos = {n: i for i, n in enumerate(immune_obs)}
    immune_genes = L.var_names(X.H5_IMMUNE)
    cells_ctx, obsnames_ctx = {}, {}

    for ctx, c in X.CONTEXTS.items():
        P = PROV[ctx]
        run_dir = c["emb"].parent
        ref_all = np.load(run_dir / "layer_gene_embedding_count.npy")
        ref = ref_all[0].astype(np.int64)
        layers_equal = bool(all((ref_all[l] == ref_all[0]).all() for l in range(ref_all.shape[0])))
        genes_ctx = L.var_names(c["h5ad"])                       # the assumed axis
        genes_ext = L.var_names(P["h5"])                         # the axis the extractor actually saw
        emb = np.load(c["emb"], mmap_mode="r")
        axis_same = genes_ext == genes_ctx
        R = L.H5Rows(P["h5"])
        n_obs_ext = R.shape[0]

        # ---- reproduce draw + tokenisation -------------------------------------------------
        cells_ext = L.sample_cells(n_obs_ext, P["max_cells"], P["seed"])
        cnt = L.count_vector(R, cells_ext, MAX_GENES)
        exact = bool((cnt == ref).all())
        # ties at the max_genes cutoff (where argsort tie-breaking could matter)
        nnz = np.zeros(cells_ext.size, dtype=np.int64); ties = 0
        for k, r in enumerate(cells_ext):
            v = R.row_dense(int(r)); nnz[k] = int((v > 0).sum()); ties += L.tie_at_cutoff(v, MAX_GENES)

        # ---- seed uniqueness / search for the checkpoint runs (command trace incomplete) ----
        if ctx.startswith("kidney_ckpt"):
            if "kidney_ckpt" not in res["seed_search"]:
                seeds = list(range(0, 2001))
                hits = []
                cache = {}
                for s in seeds:
                    cc = L.count_vector(R, L.sample_cells(n_obs_ext, P["max_cells"], s), MAX_GENES, cache)
                    if (cc == ref).all(): hits.append(s)
                # also: does 256 cells with any small seed match? (checks the max_cells assumption)
                hits256 = [s for s in (42, 67, 87) if (L.count_vector(R, L.sample_cells(n_obs_ext, 256, s), MAX_GENES) == ref).all()]
                res["seed_search"]["kidney_ckpt"] = dict(max_cells=P["max_cells"], seeds_tested="0..2000",
                                                         seeds_matching=hits, max_cells_256_seeds_42_67_87_matching=hits256)
                log("seed search kidney_ckpt:", res["seed_search"]["kidney_ckpt"])

        # ---- map extraction rows to rows of the CONTEXTS h5ad ------------------------------
        obs_ext = L.obs_names(P["h5"])
        if P["sub"] is None:
            assert Path(P["h5"]) == Path(c["h5ad"])
            cells = cells_ext.copy(); sub_check = None
        else:
            col, val = P["sub"]
            parent_col = L.obs_col(c["h5ad"], col)
            mask_rows = np.flatnonzero(parent_col == val)
            sub_names_match = bool(np.array_equal(immune_obs[mask_rows], obs_ext))
            cells = np.array([imm_pos[n] for n in obs_ext[cells_ext]], dtype=np.int64)
            # X rows identical between subset file and parent file for the drawn cells?
            RP = L.H5Rows(c["h5ad"])
            x_equal = bool(all(np.array_equal(R.row_dense(int(a)), RP.row_dense(int(b))) for a, b in zip(cells_ext, cells)))
            RP.close()
            sub_vals = np.unique(L.obs_col(P["h5"], col))
            sub_check = dict(subset_file=str(P["h5"]), obs_column=col, value=val, subset_n_obs=int(n_obs_ext),
                             parent_n_with_value=int(mask_rows.size), subset_equals_parent_mask_in_order=sub_names_match,
                             subset_unique_values=sub_vals.tolist(), drawn_rows_X_identical_to_parent=x_equal)
        cells_ctx[ctx] = cells
        obsnames_ctx[ctx] = (str(c["h5ad"]), L.obs_names(c["h5ad"])[cells])

        # ---- gene-axis checks --------------------------------------------------------------
        nonzero_emb = np.abs(np.asarray(emb[0])).sum(1) > 0
        observed_eq_count = bool(np.array_equal(nonzero_emb, ref > 0))
        edge_chk = edge_axis_check(run_dir, genes_ctx)

        # ---- cache outputs -----------------------------------------------------------------
        d = CACHE / ctx
        obs_idx = np.load(d / "obs_idx.npy"); obs_genes = json.loads((d / "obs_genes.json").read_text())
        assert [genes_ctx[i] for i in obs_idx] == obs_genes
        ncells_cached = np.load(d / "ncells.npy")
        nobs_rep = cnt[obs_idx].astype(np.int64)
        RC = L.H5Rows(c["h5ad"])
        expr = RC.rows_dense(cells)[:, obs_idx].astype(np.float32)
        RC.close()
        for fn in CACHE_NEW:
            p = d / fn
            if p.exists():
                log(f"[note] {p} exists, overwriting (file created by this task)")
        np.save(d / "extraction_cells.npy", cells.astype(np.int64))
        np.save(d / "nobs_reproduced.npy", nobs_rep)
        np.save(d / "expr_extraction.npy", expr)
        # consistency: genes present (expr>0) in cells >= tokenised count; equal unless truncated
        det = (expr > 0).sum(0)

        # ---- tissue / donor tables ---------------------------------------------------------
        obs_tab = {}
        for col in OBS_COLS:
            v = L.obs_col(c["h5ad"], col)
            if v is None:
                obs_tab[col] = None; continue
            vc = pd.Series(v[cells]).value_counts()
            obs_tab[col] = {str(k): int(n) for k, n in vc.items()}
            for k, n in vc.items():
                obs_rows.append(dict(context=ctx, column=col, value=str(k), n_cells=int(n), frac=float(n / cells.size)))

        info = dict(
            emb=str(c["emb"]), h5ad_context=str(c["h5ad"]), h5ad_extraction=str(P["h5"]), provenance=P["doc"],
            seed=P["seed"], max_cells=P["max_cells"], max_genes=MAX_GENES, n_obs_extraction_h5ad=int(n_obs_ext),
            n_cells=int(cells.size), subsetting=sub_check, G=len(genes_ctx), emb_G=int(emb.shape[1]),
            extraction_axis_equals_context_axis=bool(axis_same), count_layers_identical=layers_equal,
            count_exact_match=exact, n_genes_differ=int((cnt != ref).sum()), count_sum_reproduced=int(cnt.sum()),
            count_sum_saved=int(ref.sum()), n_observed=int((ref > 0).sum()),
            observed_mask_equals_count_gt0=observed_eq_count,
            nobs_reproduced_equals_cache_ncells=bool(np.array_equal(nobs_rep, ncells_cached)),
            cells_with_nnz_gt_max_genes=int((nnz > MAX_GENES).sum()), cells_with_tie_at_cutoff=int(ties),
            nnz_median=float(np.median(nnz)), nnz_max=int(nnz.max()),
            detected_ge_tokenised_all=bool((det >= nobs_rep).all()), n_genes_detected_gt_tokenised=int((det > nobs_rep).sum()),
            edge_axis_check=edge_chk, obs_tables=obs_tab,
        )
        if not exact:
            from scipy.stats import spearmanr
            info["spearman"] = float(spearmanr(cnt, ref).correlation); info["frac_equal"] = float((cnt == ref).mean())
        res["contexts"][ctx] = info
        rows.append(dict(context=ctx, seed=P["seed"], max_cells=P["max_cells"], n_cells=int(cells.size),
                         extraction_h5ad=Path(P["h5"]).name, subset=None if P["sub"] is None else f"{P['sub'][0]}=={P['sub'][1]}",
                         G=len(genes_ctx), axis_same=bool(axis_same), count_exact=exact, n_genes_differ=info["n_genes_differ"],
                         n_observed=info["n_observed"], edge_index_mismatch=edge_chk.get("n_index_mismatch"),
                         edge_gene_index_pairs=edge_chk.get("n_gene_index_pairs"),
                         ties_at_cutoff=int(ties), cells_trunc=info["cells_with_nnz_gt_max_genes"],
                         nobs_eq_cache=info["nobs_reproduced_equals_cache_ncells"]))
        log(json.dumps(rows[-1]))
        R.close()

    # ---- overlaps -----------------------------------------------------------------------------
    ov = []
    def overlap(a, b, how):
        ha, na = obsnames_ctx[a]; hb, nb = obsnames_ctx[b]
        sa, sb = set(na), set(nb)
        inter = sa & sb
        r = dict(a=a, b=b, how=how, same_h5ad=ha == hb, n_a=len(sa), n_b=len(sb), n_shared=len(inter),
                 jaccard=len(inter) / len(sa | sb))
        ov.append(r); return r
    for b in ["immune_seed43", "immune_seed44", "immune_TSP14", "immune_TSP21", "immune_TSP25",
              "immune_Bcell", "immune_CD4T", "immune_CD8T"]:
        overlap("immune", b, "obs_names (immune h5ad)")
    overlap("immune_seed43", "immune_seed44", "obs_names (immune h5ad)")
    overlap("kidney_A", "kidney_B", "obs_names (different h5ads)")
    overlap("kidney_A", "kidney_ckpt_whole_human", "obs_names (different h5ads)")
    overlap("kidney_B", "kidney_ckpt_whole_human", "obs_names (kidney_B h5ad)")
    overlap("lung", "lung_external", "obs_names (different h5ads)")
    # whole-file obs overlap kidney_A h5ad vs kidney_B h5ad (is kidney_B's file a subset of TS kidney?)
    nA, nB = set(L.obs_names(X.H5_KIDNEY_A)), set(L.obs_names(X.H5_KIDNEY_B))
    ojA = L.obs_col(X.H5_KIDNEY_A, "observation_joinid"); ojB = L.obs_col(X.H5_KIDNEY_B, "observation_joinid")
    res["kidney_h5ad_file_overlap"] = dict(n_kidney_A_file=len(nA), n_kidney_B_file=len(nB), n_shared_obs_names=len(nA & nB),
                                           n_shared_observation_joinid=len(set(ojA) & set(ojB)),
                                           kidney_A_B_var_shared=len(set(L.var_names(X.H5_KIDNEY_A)) & set(L.var_names(X.H5_KIDNEY_B))))
    # cross-check the obs_names overlap with observation_joinid for the drawn cells
    ojA_d = set(ojA[cells_ctx["kidney_A"]]); ojB_d = set(ojB[cells_ctx["kidney_B"]])
    res["kidney_A_vs_B_drawn_observation_joinid_shared"] = len(ojA_d & ojB_d)
    # expected overlap under two independent uniform draws: |A_draw in B's population| * n_B / |B population|
    pop = {}
    for ctx, P in PROV.items():
        hname = str(X.CONTEXTS[ctx]["h5ad"])
        if P["sub"] is None:
            pop[ctx] = (hname, set(L.obs_names(P["h5"])))
        else:
            pop[ctx] = (hname, set(L.obs_names(P["h5"])))
    for r in ov:
        a, b = r["a"], r["b"]
        na = set(obsnames_ctx[a][1]); popb = pop[b][1]
        r["n_a_in_b_population"] = len(na & popb)
        r["expected_shared_if_independent"] = len(na & popb) * r["n_b"] / len(popb)
    # same cell, same gene: are kidney_A and kidney_B X values identical (same normalisation)?
    gA = L.var_names(X.H5_KIDNEY_A); gB = L.var_names(X.H5_KIDNEY_B)
    shared_g = sorted(set(gA) & set(gB)); iA = np.array([gA.index(g) for g in shared_g]); iB = np.array([gB.index(g) for g in shared_g])
    posA = {n: i for i, n in enumerate(L.obs_names(X.H5_KIDNEY_A))}
    namesB = L.obs_names(X.H5_KIDNEY_B)
    probe = cells_ctx["kidney_B"][:32]
    RA, RB = L.H5Rows(X.H5_KIDNEY_A), L.H5Rows(X.H5_KIDNEY_B)
    maxdiff = max(float(np.abs(RA.row_dense(posA[namesB[r]])[iA] - RB.row_dense(int(r))[iB]).max()) for r in probe)
    RA.close(); RB.close()
    res["kidney_A_vs_B_X_shared_genes_maxabsdiff_32cells"] = maxdiff
    # gene_list.txt in the kidney_A run dir (written 2 days after the extraction): what is it?
    gl = (X.CONTEXTS["kidney_A"]["emb"].parent / "gene_list.txt").read_text().split("\n")
    nonempty = [(i, g) for i, g in enumerate(gl) if g.strip()]
    gAset = set(gA)
    res["kidney_A_gene_list_txt"] = dict(n_lines=len(gl), n_nonempty=len(nonempty),
                                         n_nonempty_in_kidney_A_var=sum(g in gAset for _, g in nonempty),
                                         n_nonempty_line_index_equals_var_position=sum(i < len(gA) and gA[i] == g for i, g in nonempty))
    res["overlaps"] = ov
    pd.DataFrame(ov).to_csv(HERE / "overlaps.csv", index=False)
    pd.DataFrame(rows).to_csv(HERE / "per_context.csv", index=False)
    pd.DataFrame(obs_rows).to_csv(HERE / "obs_tables.csv", index=False)
    res["all_exact"] = bool(all(v["count_exact_match"] for v in res["contexts"].values()))
    (HERE / "results.json").write_text(json.dumps(res, indent=1, default=str))
    log("all_exact", res["all_exact"])
    for r in ov: log(json.dumps(r))
    log(json.dumps(res["kidney_h5ad_file_overlap"]), res["kidney_A_vs_B_drawn_observation_joinid_shared"])


if __name__ == "__main__":
    main()
