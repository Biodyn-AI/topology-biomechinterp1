"""E4c verifier part 4: pull every number quoted in the implementer's SUMMARY.txt from the CSVs (bookkeeping)."""
from pathlib import Path
import pandas as pd, numpy as np
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 400)
H = Path(__file__).resolve().parent.parent
P = pd.read_csv(H / "probe_results.csv"); R = pd.read_csv(H / "relational_signal.csv")
N = pd.read_csv(H / "null_results.csv"); D = pd.read_csv(H / "direction_test.csv"); TV = pd.read_csv(H / "tuned_vs_fixed.csv")
CT = ["lung", "lung_external", "kidney_B"]
print("## 1. best C3 probe means per context (M1, scgpt; non-node features)")
c3 = P[(P.regime == "C3") & (P.model == "M1") & ~P.feature.isin(["tf_only", "tg_only", "coexpr_marg"])]
for ctx in CT:
    g = c3[c3.context == ctx].sort_values("auc_mean", ascending=False)
    print(ctx, g[["source", "feature", "layer", "auc_mean", "auc_ci_lo", "auc_ci_hi", "auc_sd_repeats"]].head(4).round(4).to_string(index=False))
    for src in ("static", "esm2", "coexpr"):
        gg = g[g.source == src].head(1); print("  best", src, gg[["feature", "auc_mean", "auc_ci_hi"]].round(4).values.tolist())
print("max auc_sd_repeats M1 C3:", P[(P.model == "M1") & (P.regime == "C3")].auc_sd_repeats.max().round(4))
print("upper NB of best scgpt probes per ctx:", c3[c3.source == "scgpt"].groupby("context").auc_ci_hi.max().round(4).to_dict())
print("\n## 2. target-only / tf-only by regime (M1 scgpt)")
node = P[(P.model == "M1") & (P.source == "scgpt") & P.feature.isin(["tf_only", "tg_only"]) & P.context.isin(CT)]
print(node.groupby(["feature", "regime"]).auc_mean.agg(["min", "max"]).round(4))
print(node.groupby(["context", "feature", "regime"]).auc_mean.agg(["min", "max"]).round(4))
m1 = R[(R.model == "M1") & (R.source == "scgpt")]
def lay(sub, col="delta_node"):
    return dict(min=round(sub[col].min(), 4), max=round(sub[col].max(), 4), n_lo_gt0=int((sub[col + "_lo"] > 0).sum()),
                n_hi_lt0=int((sub[col + "_hi"] < 0).sum()), n=len(sub))
print("\n## concat/hadamard/absdiff/bilinear/tg_* minus node, per ctx/regime")
for feat in ("concat", "hadamard", "absdiff", "bilinear", "tg_bilinear", "tg_hadamard"):
    for reg in ("C3", "C2", "C1"):
        for ctx in CT:
            s = m1[(m1.feature == feat) & (m1.regime == reg) & (m1.context == ctx)]
            if len(s): print(feat, reg, ctx, lay(s))
print("\n## 3. lung bilinear C3 per layer")
for ctx in CT:
    s = m1[(m1.feature == "bilinear") & (m1.regime == "C3") & (m1.context == ctx)].sort_values("layer")
    print(s[["context", "layer", "probe_auc", "tg_only_auc", "delta_node", "delta_node_lo", "delta_node_hi", "rewired_auc", "p_emp_rewired",
             "p_maxT_layers", "relational_signal", "relational_signal_lo", "relational_signal_hi", "delta_cosine", "delta_cosine_lo",
             "delta_abs_coexpr_lo", "cosine_auc"]].round(4).to_string(index=False))
    print(" n p_emp<.05", int((s.p_emp_rewired < .05).sum()), "n maxT<.05", int((s.p_maxT_layers < .05).sum()),
          "n rel lo>0", int((s.relational_signal_lo > 0).sum()), "n dcos lo>0", int((s.delta_cosine_lo > 0).sum()),
          "dcos range", round(s.delta_cosine.min(), 4), round(s.delta_cosine.max(), 4), "n dcoexpr lo>0", int((s.delta_abs_coexpr_lo > 0).sum()),
          "cos range", round(s.cosine_auc.min(), 4), round(s.cosine_auc.max(), 4))
print("\n## tg_bilinear C3 per layer")
for ctx in CT:
    s = m1[(m1.feature == "tg_bilinear") & (m1.regime == "C3") & (m1.context == ctx)].sort_values("layer")
    print(ctx, s[["layer", "probe_auc", "delta_node", "delta_node_lo"]].round(4).values.tolist())
print("\n## 4. static / esm2 bilinear C3")
print(R[(R.model == "M1") & R.source.isin(["static", "esm2"]) & (R.regime == "C3")][["context", "source", "feature", "probe_auc", "tg_only_auc",
      "delta_node", "delta_node_lo", "delta_node_hi", "p_emp_rewired", "rewired_auc"]].round(4).to_string(index=False))
print("\n## 6. nonlinear")
nl = R[R.model.isin(["M2", "M3", "M2_tuned", "M3_tuned"]) & (R.source == "scgpt")]
print("C3 probe auc range", nl[nl.regime == "C3"].probe_auc.min().round(4), nl[nl.regime == "C3"].probe_auc.max().round(4))
print("any C3 delta_node_lo>0:", int((nl[nl.regime == "C3"].delta_node_lo > 0).sum()), " C2:", int((nl[nl.regime == "C2"].delta_node_lo > 0).sum()))
print(nl[(nl.delta_node_hi < 0)][["context", "model", "feature", "layer", "regime", "delta_node", "delta_node_lo", "delta_node_hi"]].round(4).to_string(index=False))
print(nl[nl.regime == "C2"][["context", "model", "feature", "layer", "probe_auc", "delta_node", "delta_node_lo"]].round(4).to_string(index=False))
print("tuned-fixed", TV.tuned_minus_fixed.min(), TV.tuned_minus_fixed.max(), TV.groupby("model").tuned_minus_fixed.mean().round(4).to_dict())
nn = N[N.model.isin(["M2", "M3"])]
print(nn[["context", "model", "source", "feature", "layer", "obs_auc_mean", "null_mean", "z", "p_emp"]].round(3).to_string(index=False))
print("\n## 7. direction")
print(D[["context", "layer", "n_pairs", "frac_forward_higher", "ci_lo", "ci_hi", "frac_forward_strict", "strict_ci_lo", "n_strict",
         "frac_outdeg_i_gt_j", "frac_indeg_j_gt_i", "frac_nobs_i_gt_j"]].round(3).to_string(index=False))
print("\n## 8. immune")
print(P[P.context == "immune"][["feature", "layer", "regime", "auc_mean", "n_pooled_folds", "mean_n_tf_used", "mean_n_pos"]].groupby(["feature", "regime"]).agg(
    lo=("auc_mean", "min"), hi=("auc_mean", "max"), pooled=("n_pooled_folds", "sum"), ntf=("mean_n_tf_used", "mean"), npos=("mean_n_pos", "mean")).round(3))
print("\n## fold sizes C3 M1 scgpt bilinear")
print(P[(P.model == "M1") & (P.regime == "C3") & (P.feature == "bilinear") & (P.source == "scgpt")].groupby("context")[["mean_n_tf_used", "mean_n_pos", "mean_n_test", "mean_n_train", "n_pooled_folds"]].mean().round(1))
print("pooled folds by context:", P.groupby("context").n_pooled_folds.sum().to_dict())
print("null max_abs_diff>0.01:", int((N.max_abs_diff_vs_true_run > 0.01).sum()), "of", int(N.max_abs_diff_vs_true_run.notna().sum()), "M1 rows", int((N.model=='M1').sum()), "median", N.max_abs_diff_vs_true_run.median())
