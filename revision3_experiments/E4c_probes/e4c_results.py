"""Write results.json from the aggregated CSVs (run e4c_aggregate.py first). Bookkeeping only."""
import json
import numpy as np, pandas as pd
import e4c_lib as L
H = L.HERE
P = pd.read_csv(H / "probe_results.csv"); R = pd.read_csv(H / "relational_signal.csv")
N = pd.read_csv(H / "null_results.csv"); D = pd.read_csv(H / "direction_test.csv")
RW = pd.read_csv(H / "rewire_stats.csv")
CTX = ["lung", "lung_external", "kidney_B"]


def f(x, nd=4):
    return None if pd.isna(x) else round(float(x), nd)


def layer_summary(sub):
    """min / max over layers and counts of layers whose interval excludes 0 (both directions)."""
    o = dict(n_layers=int(sub.layer.nunique()),
             probe_auc_min=f(sub.probe_auc.min()), probe_auc_max=f(sub.probe_auc.max()),
             tg_only_auc_min=f(sub.tg_only_auc.min()), tg_only_auc_max=f(sub.tg_only_auc.max()),
             delta_node_min=f(sub.delta_node.min()), delta_node_max=f(sub.delta_node.max()),
             n_layers_delta_node_lo_gt0=int((sub.delta_node_lo > 0).sum()),
             n_layers_delta_node_hi_lt0=int((sub.delta_node_hi < 0).sum()))
    if "delta_cosine" in sub:
        o.update(delta_cosine_min=f(sub.delta_cosine.min()), delta_cosine_max=f(sub.delta_cosine.max()),
                 n_layers_delta_cosine_lo_gt0=int((sub.delta_cosine_lo > 0).sum()),
                 delta_abs_coexpr_min=f(sub.delta_abs_coexpr.min()), delta_abs_coexpr_max=f(sub.delta_abs_coexpr.max()),
                 n_layers_delta_abs_coexpr_lo_gt0=int((sub.delta_abs_coexpr_lo > 0).sum()))
    if "relational_signal" in sub and sub.relational_signal.notna().any():
        s = sub[sub.relational_signal.notna()]
        o.update(rewired_auc_min=f(s.rewired_auc.min()), rewired_auc_max=f(s.rewired_auc.max()),
                 relational_signal_min=f(s.relational_signal.min()), relational_signal_max=f(s.relational_signal.max()),
                 n_layers_relational_lo_gt0=int((s.relational_signal_lo > 0).sum()),
                 n_layers_relational_hi_lt0=int((s.relational_signal_hi < 0).sum()),
                 n_layers_p_emp_rewired_lt05=int((s.p_emp_rewired < 0.05).sum()),
                 n_layers_p_maxT_lt05=int((s.p_maxT_layers < 0.05).sum()) if s.p_maxT_layers.notna().any() else None,
                 null_N=int(s.null_N.iloc[0]), rewired_claim_allowed=bool(s.rewired_claim_allowed.iloc[0]))
    return o


out = dict(experiment="E4c_probes", role="sensitivity analysis (CONVENTION3 §5, E4 extras: supervised probes); no pass/fail",
           convention3_sha256=L.conv_hash(), frozen_txt=(L.R3 / "FROZEN.txt").read_text().strip(),
           contexts=CTX, immune="node-only (TF-only / target-only) sanity check only",
           esm2_label="ESM2-15B (UCE per-protein), 5,120-d, PCA-128 fitted on training-fold genes",
           static_label="scGPT static token embedding (subproject_53 npz), 512-d, centred over V",
           metric="within-TF stratified AUROC (stats3.StratifiedPairs) over TFs with >=3 positives and >=3 candidate non-edges in the test fold, else pooled AUROC (counts per row in probe_results.csv)",
           interval="Nadeau-Bengio corrected resampled t over outer folds (var = (1/J + n_test/n_train) s^2, df = J-1), 95% two-sided",
           errata_item3="not applicable: no network-statistic interval with vertex subsampling is computed in E4c",
           rewiring=RW.to_dict("records"))
# design actually run
P["n_folds"] = P.n_folds.astype(int)
design = P.groupby(["model", "regime"]).agg(n_rows=("auc_mean", "size"), max_repeats=("n_repeats", "max")).reset_index()
out["design_rows"] = design.to_dict("records")
out["pooled_fallback_folds_total"] = int(P.n_pooled_folds.sum()); out["within_tf_folds_total"] = int(P.n_within_tf_folds.sum())
# M1 summaries per context x regime x source x feature
m1 = {}
for (ctx, reg, src, feat), g in R[R.model == "M1"].groupby(["context", "regime", "source", "feature"]):
    m1.setdefault(ctx, {}).setdefault(reg, {})[f"{src}/{feat}"] = layer_summary(g)
out["M1_relational_by_context_regime"] = m1
# nonlinear (fixed and tuned)
nl = R[R.model.isin(["M2", "M3", "M2_tuned", "M3_tuned"])]
cols = ["context", "model", "source", "feature", "layer", "regime", "probe_auc", "tg_only_auc", "delta_node", "delta_node_lo",
        "delta_node_hi", "rewired_auc", "relational_signal", "relational_signal_lo", "relational_signal_hi", "p_emp_rewired",
        "delta_cosine", "delta_abs_coexpr", "notes"]
out["nonlinear_rows"] = [{k: (f(v) if isinstance(v, float) else v) for k, v in r.items()} for r in nl[[c for c in cols if c in nl]].to_dict("records")]
# tuned vs fixed on identical folds
PF = pd.read_csv(H / "probe_folds.csv")
tv = []
for base in ("M2", "M3"):
    a = PF[PF.model == base].set_index(["context", "source", "feature", "layer", "regime", "repeat", "fold"]).auc
    b = PF[PF.model == base + "_tuned"].set_index(["context", "source", "feature", "layer", "regime", "repeat", "fold"]).auc
    j = pd.concat([a.rename("fixed"), b.rename("tuned")], axis=1, join="inner").reset_index()
    for k, g in j.groupby(["context", "source", "feature", "layer", "regime"]):
        tv.append(dict(zip(["context", "source", "feature", "layer", "regime"], k), model=base, n_folds=len(g),
                       fixed=f(g.fixed.mean()), tuned=f(g.tuned.mean()), tuned_minus_fixed=f((g.tuned - g.fixed).mean())))
TV = pd.DataFrame(tv); TV.to_csv(H / "tuned_vs_fixed.csv", index=False)
out["tuned_vs_fixed_summary"] = (TV.groupby(["model", "regime"]).tuned_minus_fixed.describe()[["count", "mean", "min", "max"]]
                                 .round(4).reset_index().to_dict("records") if len(TV) else [])
# direction test
out["direction_test"] = {ctx: dict(n_pairs=int(g.n_pairs.iloc[0]), frac_forward_min=f(g.frac_forward_higher.min()),
                                   frac_forward_max=f(g.frac_forward_higher.max()),
                                   n_layers_ci_lo_gt05=int((g.ci_lo > 0.5).sum()), n_layers_ci_hi_lt05=int((g.ci_hi < 0.5).sum()),
                                   n_strict=int(g.n_strict.iloc[0]) if "n_strict" in g else None,
                                   frac_forward_strict_min=f(g.frac_forward_strict.min()) if "n_strict" in g else None,
                                   frac_forward_strict_max=f(g.frac_forward_strict.max()) if "n_strict" in g else None,
                                   outdeg_baseline_frac_i_gt_j=f(g.frac_outdeg_i_gt_j.iloc[0]),
                                   indeg_baseline_frac_j_gt_i=f(g.frac_indeg_j_gt_i.iloc[0]),
                                   probe_frac_when_outdeg_favours_j_min=f(g.probe_frac_when_outdeg_favours_j.min()),
                                   probe_frac_when_outdeg_favours_j_max=f(g.probe_frac_when_outdeg_favours_j.max()),
                                   n_outdeg_favours_j=int(g.n_outdeg_favours_j.iloc[0]))
                         for ctx, g in D.groupby("context")}
# immune node-only
imm = P[(P.context == "immune")]
out["immune_node_only"] = imm.groupby(["feature", "regime"]).auc_mean.agg(["min", "max"]).round(4).reset_index().to_dict("records")
# reference scorers on the same folds
ref = P[P.model == "ref"]
out["reference_scorers_same_folds"] = ref.groupby(["context", "feature", "regime"]).auc_mean.agg(["min", "max"]).round(4).reset_index().to_dict("records")
out["null_consistency_max_abs_diff_true_vs_null_run"] = f(N.max_abs_diff_vs_true_run.max(), 6)
out["reductions_vs_specification"] = [
    "C1 random-pair regime: 1 repeat (5 folds) instead of 5 (allowed by the specification's compute clause); C2 and C3: 5 x 5 for M1",
    "M1 rewired-label null (100 degree-preserving rewirings): run on C3 repeat 0 only (5 folds); not on C2, C1 or C3 repeats 1-4",
    "M2 (HistGradientBoosting) / M3 (MLP): layers L0, L5, L11; regimes C3 and C2 repeat 0 only; features concat and hadamard only (no absdiff / bilinear); scGPT and co-expression sources only (no ESM2 / static)",
    "M2/M3 TF-only control not fitted: under the within-TF metric a TF-only score is constant within TF, AUROC = 0.5 by construction",
    "M2/M3 rewired null: 10 rewirings, C3 repeat 0; M2 null for concat only; N = 10 has resolution floor 1/11 = 0.091 > 0.05, so no significance claim (stats3.assert_resolution fails); effect sizes and z only",
    "M2/M3 main runs: hyperparameters fixed a priori, only early stopping in-fold; tuned variant (small grids, regime-matched inner split) added on true labels for the same folds (M2_tuned / M3_tuned); nulls not re-run with tuning",
    "co-expression node baseline for M2/M3 rows = M1 marginal-only (mean, detection) model",
    "run-to-run SD = SD over the 5 CV repetitions (fold-assignment seeds) for M1; model seeds fixed; M2/M3 have one repetition (no run-to-run SD)",
    "relational-signal intervals: Nadeau-Bengio corrected resampled t (no paired cluster CI)"]
out["additions_by_resuming_agent"] = [
    "e4c_ref.py: unsupervised H3 scorers (signed cosine full-512 per layer; |co-expression|) on the identical outer test folds -> delta_cosine / delta_abs_coexpr columns in relational_signal.csv",
    "e4c_nonlin_tuned.py: M2/M3 with in-fold small-grid tuning (tuned_vs_fixed.csv)",
    "direction test strict subset: j->i absent from TRRUST and every DoRothEA level (columns *_strict in direction_test.csv)",
    "auc_sd_repeats column in probe_results.csv; results.json writer (e4c_results.py)"]
json.dump(out, open(H / "results.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("results.json written")
