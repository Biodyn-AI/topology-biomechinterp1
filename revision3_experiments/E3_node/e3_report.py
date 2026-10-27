"""Add headline numbers (read from the CSV/JSON outputs; no new statistics) to results.json."""
import json
import numpy as np, pandas as pd
import e3_lib as L
H = L.HERE
res = json.loads((H / "results.json").read_text())
assert res["conv_hash"] == res["frozen_hash"], "CONVENTION3 hash mismatch"
hl = {}
for c in ("lung", "lung_external", "kidney_B", "immune", "kidney_A"):
    t = pd.read_csv(H / f"h2_primary_{c}.csv"); lay = t[t.row_type == "H2_layer"]; b = t[t.row_type == "baseline"].set_index("baseline")
    lg = pd.read_csv(H / "logistic_adjustment.csv"); lg = lg[(lg.context == c)]
    rv = pd.read_csv(H / "restricted_V.csv"); rv = rv[(rv.context == c) & (rv.row_type == "H2_layer")]
    bp = pd.read_csv(H / "baselines_paired.csv"); bp = bp[bp.context == c]
    tt = pd.read_csv(H / "trrust_task.csv"); tt = tt[tt.context == c]; ttl = tt[tt.row_type == "layer"]; ttb = tt[tt.row_type == "baseline"].set_index("baseline")
    hl[c] = dict(
        n_V=int(lay.n_genes.iloc[0]), n_tf=int(lay.n_tf.iloc[0]),
        strat_auroc_range=[float(lay.strat_auroc.min()), float(lay.strat_auroc.max())],
        best_layer=int(lay.layer.iloc[int(np.argmax(lay.strat_auroc.values))]),
        best=dict(zip(["strat_auroc", "lo", "hi", "auroc"], map(float, lay.iloc[int(np.argmax(lay.strat_auroc.values))][["strat_auroc", "strat_lo", "strat_hi", "auroc"]]))),
        min_lo=float(lay.strat_lo.min()), z_range=[float(lay.z.min()), float(lay.z.max())], max_p_maxT=float(lay.p_maxT.max()),
        baselines={k: dict(auroc=float(b.loc[k, "auroc"]), strat_auroc=float(b.loc[k, "strat_auroc"]),
                           strat_ci=[float(b.loc[k, "strat_lo"]), float(b.loc[k, "strat_hi"])]) for k in b.index},
        static_minus_layer_delta={s: dict(range=[float(g.delta.min()), float(g.delta.max())], min_lo=float(g.delta_lo.min()))
                                  for s, g in bp.groupby("static_label")},
        logistic_extraction=dict(beta_range=[float(lg[lg.expr_covariate == "extraction"].beta_score_per_sd.min()),
                                             float(lg[lg.expr_covariate == "extraction"].beta_score_per_sd.max())],
                                 beta_lo_min=float(lg[lg.expr_covariate == "extraction"].beta_lo.min()),
                                 LR_p_max=float(lg[lg.expr_covariate == "extraction"].LR_p.max())),
        logistic_population=dict(beta_range=[float(lg[lg.expr_covariate == "population"].beta_score_per_sd.min()),
                                             float(lg[lg.expr_covariate == "population"].beta_score_per_sd.max())],
                                 LR_p_max=float(lg[lg.expr_covariate == "population"].LR_p.max())),
        restricted={int(k): dict(n=int(g.n_genes.iloc[0]), n_tf=int(g.n_tf.iloc[0]), range=[float(g.strat_auroc.min()), float(g.strat_auroc.max())],
                                 max_p_maxT=float(g.p_maxT.max())) for k, g in rv.groupby("min_nobs")},
        trrust_task=dict(n=int(ttl.n_genes.iloc[0]), n_tf=int(ttl.n_tf.iloc[0]), range=[float(ttl.strat_auroc.min()), float(ttl.strat_auroc.max())],
                         max_p_maxT=float(ttl.p_maxT.max()), log_nobs_strat=float(ttb.loc["log_nobs", "strat_auroc"]),
                         log_nobs_auroc=float(ttb.loc["log_nobs", "auroc"]), log_expr_strat=float(ttb.loc["log_mean_expr_extraction", "strat_auroc"])))
res["headline"] = hl
cf = pd.read_csv(H / "crossfit.csv"); res["crossfit"] = cf.to_dict(orient="records")
bm = pd.read_csv(H / "bandscan_maxT.csv"); res["bandscan_maxT"] = bm[bm.context.isin(["lung", "lung_external", "immune", "kidney_A", "kidney_B"])].to_dict(orient="records")
res["ledger_0744_file"] = "ledger_0744.json"
res["labels"] = dict(static_esm2="ESM2-15B (UCE per-protein; 5,120-d; file subproject_60.../uce_esm2_human_gene_embeddings.npz) -- not ESM2-3B",
                     static_scgpt="scGPT static gene-token embedding (subproject_53 scgpt_gene_embeddings.npz)")
res["errata_item3"] = ("not applicable: E3 reports node statistics only; every interval is a gene bootstrap (B=1,000) and no "
                       "vertex-subsampling interval is computed, so the delete-d rescaling of CONVENTION3_ERRATA item 3 does not arise")
res["permutation_resolution"] = dict(primary_N=10000, primary_floor=1 / 10001, bandscan_N=2000, bandscan_floor=1 / 2001,
                                     restricted_and_static_N=2000, ledger_screen_procedure_N=1000, ledger_loo_N=10000,
                                     assert_resolution_called=True)
res["notes"] = ["bootstrap CIs resample genes with the LOO scores held fixed (TF centroid not re-estimated per resample)",
                "band-scan 'lda' variant is fitted in-sample (refit in every null draw, so p-values are valid but AUROCs are optimistic); "
                "the cross-fit (cos variant) is the selection-unbiased estimate",
                "kidney_ckpt_* contexts share cells and genes (n_V=4923); kidney_A/kidney_B share a donor; immune/lung share donors"]
(H / "results.json").write_text(json.dumps(res, indent=1, default=float))
print(json.dumps(hl["lung"], indent=0)[:3000]); print(json.dumps(hl["lung_external"], indent=0)[:3000])
