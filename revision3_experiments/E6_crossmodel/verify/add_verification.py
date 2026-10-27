"""Verifier: table of stored-vs-independently-recomputed numbers (verify/reproduced_numbers.csv) and a 'verification' block
in ../results.json. Run after e6_assemble.py (which rewrites results.json)."""
import json
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; E6 = H.parent
v = json.loads((H / "verify_results.json").read_text())
A = pd.read_csv(E6 / "tierA_results.csv"); A = A[A.estimator == "full"]
D = pd.read_csv(E6 / "node_detstrat_results.csv"); S = pd.read_csv(E6 / "selection_bias.csv"); Lr = pd.read_csv(E6 / "legacy_reproduction.csv")
def a(e, t): return float(A[(A.embedding == e) & (A.task == t)].value.iloc[0])
def d(s, e): return float(D[(D.set == s) & (D.embedding == e)].value.iloc[0])
def sb(m, c): x = S[(S.row_type == "subset_null") & (S.variant == "legacy186") & (S.model == m)]; return float(x[c].iloc[0])
lr = Lr[Lr.section == "submitted Table 1 models"].set_index("embedding")
R = [
 ("Tier A size", 17693, v["vocab"]["tierA"], "own vocabulary intersection of raw npz (all-zero rows dropped)"),
 ("3 submitted models vocab intersection", 18907, v["vocab"]["three_submitted"], "same"),
 ("legacy set size", 186, v["legacy"]["n_legacy"], "209 named genes -> observed in kidney_A -> 3 vocabularies; identical gene list"),
 ("Tier A Lambert AUROC ESM2-15B", a("esm2_uce", "node_lambert"), v["tierA_node"]["esm2_uce"]["lambert_auc"], "own LOO cosine + sklearn AUROC"),
 ("Tier A Lambert AUROC scGPT static", a("scgpt_static", "node_lambert"), v["tierA_node"]["scgpt_static"]["lambert_auc"], "same"),
 ("Tier A Lambert detstrat ESM2-15B", d("tierA", "esm2_uce"), v["tierA_node"]["esm2_uce"]["lambert_detstrat_cache_deciles"], "own stratified AUROC, cache deciles"),
 ("Tier A Lambert detstrat scGPT static", d("tierA", "scgpt_static"), v["tierA_node"]["scgpt_static"]["lambert_detstrat_cache_deciles"], "same"),
 ("Tier A Lambert detstrat ESM2-15B (own rank deciles)", d("tierA", "esm2_uce"), v["tierA_node"]["esm2_uce"]["lambert_detstrat_own_deciles"], "own decile convention"),
 ("Tier B lung contextual L1 detstrat", d("tierB_lung", "lung@L1"), v["tierB_lung_node"]["lung_L1_detstrat"], "own, from local lung.npy"),
 ("Tier B lung contextual L11 detstrat", d("tierB_lung", "lung@L11"), v["tierB_lung_node"]["lung_L11_detstrat"], "same"),
 ("Tier B lung ESM2-15B detstrat", d("tierB_lung", "esm2_uce"), v["tierB_lung_node"]["esm2_detstrat"], "same"),
 ("legacy est(i) TRRUST scGPT L11", float(lr.loc["kidney_A@L11", "est_i_tf_target_auroc"]), v["legacy"]["kidney_A@L11"]["est_i"], "own"),
 ("legacy est(i) TRRUST GF V2-104M", float(lr.loc["gf_v2_104m", "est_i_tf_target_auroc"]), v["legacy"]["gf_v2_104m"]["est_i"], "own"),
 ("legacy est(i) TRRUST ESM2-15B", float(lr.loc["esm2_uce", "est_i_tf_target_auroc"]), v["legacy"]["esm2_uce"]["est_i"], "own"),
 ("legacy subset-SVD SV2-SV7 scGPT L11 (submitted 0.650)", float(lr.loc["kidney_A@L11", "submitted_tf_target_auroc"]), v["legacy"]["kidney_A@L11"]["subsetSVD_SV2_SV7"], "own SVD"),
 ("legacy subset-SVD SV2-SV7 GF V1-10M (submitted 0.713)", float(lr.loc["gf_v1_10m", "submitted_tf_target_auroc"]), v["legacy"]["gf_v1_10m"]["subsetSVD_SV2_SV7"], "own SVD"),
 ("legacy subset-SVD SV2-SV7 GF V2-104M (submitted 0.837)", float(lr.loc["gf_v2_104m", "submitted_tf_target_auroc"]), v["legacy"]["gf_v2_104m"]["subsetSVD_SV2_SV7"], "own SVD"),
 ("legacy subset-SVD SV2-SV7 ESM2-15B (submitted 0.924)", float(lr.loc["esm2_uce", "submitted_tf_target_auroc"]), v["legacy"]["esm2_uce"]["subsetSVD_SV2_SV7"], "own SVD"),
 ("Tier A STRING-physical pairs", 73430, v["tierA_ppi_edge"]["n_pairs"], "own parse of STRING v12 files; identical pair set"),
 ("Tier A PPI AUROC ESM2-15B (E6 non-pairs)", a("esm2_uce", "ppi"), v["tierA_ppi_edge"]["esm2_uce"]["ppi_auc_E6_nonpairs"], "own cosines + sklearn"),
 ("Tier A PPI AUROC ESM2-15B (own independent rewiring)", a("esm2_uce", "ppi"), v["tierA_ppi_edge"]["esm2_uce"]["ppi_auc_own_rewiring"], "own double-edge swap, 10m accepted swaps"),
 ("Tier A PPI AUROC scGPT static (E6 non-pairs)", a("scgpt_static", "ppi"), v["tierA_ppi_edge"]["scgpt_static"]["ppi_auc_E6_nonpairs"], "own"),
 ("Tier A PPI AUROC scGPT static (own independent rewiring)", a("scgpt_static", "ppi"), v["tierA_ppi_edge"]["scgpt_static"]["ppi_auc_own_rewiring"], "own rewiring"),
 ("Tier A PPI delete-d half-width scGPT static (upper side)", float(A[(A.embedding == "scgpt_static") & (A.task == "ppi")].ci_hi_deleted.iloc[0]) - a("scgpt_static", "ppi"),
  v["tierA_ppi_edge"]["scgpt_static"]["deleted_d_halfwidth_approx"], "own B=200 subsampling, 1.96*SD*sqrt(m/(n-m))"),
 ("Tier A edge within-TF AUROC ESM2-15B", a("esm2_uce", "edge"), v["tierA_ppi_edge"]["esm2_uce"]["edge_within_tf_auc"], "own per-TF sklearn AUROC on stored pairs"),
 ("Tier A edge within-TF AUROC scGPT static", a("scgpt_static", "edge"), v["tierA_ppi_edge"]["scgpt_static"]["edge_within_tf_auc"], "same"),
 ("selection legacy186 scGPT static legacy value", sb("scgpt_static", "legacy_value"), v["selection_bias"]["scgpt_static"]["legacy"], "own"),
 ("selection legacy186 scGPT static null mean", sb("scgpt_static", "null_mean"), v["selection_bias"]["scgpt_static"]["null_mean"], "own matched draws, N=2,000, different seed"),
 ("selection legacy186 scGPT static percentile", sb("scgpt_static", "legacy_percentile"), v["selection_bias"]["scgpt_static"]["percentile"], "same"),
 ("selection legacy186 ESM2-15B null mean", sb("esm2_uce", "null_mean"), v["selection_bias"]["esm2_uce"]["null_mean"], "same"),
 ("selection legacy186 ESM2-15B percentile", sb("esm2_uce", "legacy_percentile"), v["selection_bias"]["esm2_uce"]["percentile"], "same"),
]
df = pd.DataFrame(R, columns=["quantity", "stored", "recomputed", "how"]); df["abs_diff"] = (df.stored - df.recomputed).abs()
df.to_csv(H / "reproduced_numbers.csv", index=False); print(df.round(4).to_string())
cs = pd.read_csv(H / "centring_sensitivity.csv"); cd = pd.read_csv(H / "centring_delta.csv")
res = json.loads((E6 / "results.json").read_text())
res["verification"] = dict(
    verifier="independent adversarial verifier, 2026-10-26 (post-reboot)",
    reproduced_numbers=df.to_dict("records"),
    edge_pair_construction_checks=v["tierA_ppi_edge"]["edge_checks"],
    trrust_centring_sensitivity_tierA=cs.to_dict("records"),
    trrust_centring_paired_deltas=cd.to_dict("records"),
    corrections=[
        "Estimator (iii) wording: the best band within SV1-SV12 does NOT generally reach the full-dimension value (e.g. Geneformer "
        "V2-104M Lambert 0.630 vs 0.836, V2-316M 0.633 vs 0.817, V1-10M 0.571 vs 0.742); it does for the protein-sequence models.",
        "Tier A TRRUST TF/target-only AUROC depends on the centring set for expression-trained models: centring over all 17,693 Tier A "
        "genes (E6 convention) vs over the 2,753 TRRUST-labelled genes changes scGPT static 0.684 -> 0.789, AIDO 0.579 -> 0.690, tGPT "
        "0.520 -> 0.639; protein models change <= 0.007. ESM2-15B minus scGPT static: +0.186 [0.160, 0.212] vs +0.088 [0.066, 0.111]. "
        "In Tier B lung the convention changes little (ESM2-15B minus scGPT static +0.078 under both).",
        "pairwise_deltas_vs_scgpt.csv: new column ci_excludes_0_deleted for PPI rows (ci_excludes_0 used the sqrt(m/n) interval).",
        "Minor wording/number fixes in SUMMARY.txt (see VERIFICATION section)."])
(E6 / "results.json").write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
print("results.json updated")
