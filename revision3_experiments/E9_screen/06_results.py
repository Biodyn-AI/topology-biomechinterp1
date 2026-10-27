"""E9: collect headline numbers from this folder's outputs into results.json (no new computation)."""
import json
from pathlib import Path
import pandas as pd

H = Path(__file__).resolve().parent
A = json.loads((H / "audit_409.json").read_text())
R = json.loads((H / "claims_recompute.json").read_text())
T = json.loads((H / "screen_tables_summary.json").read_text())
M = json.loads((H / "legacy_band_mapping.json").read_text())
L = pd.read_csv(H / "E9_claims_ledger.csv")
S2 = pd.read_csv(H / "S2_hypotheses.csv")
cw = pd.read_csv(H / "crosswalk.csv")
v = cw[(cw.kind == "vector") & (cw.context == "kidney_A") & (cw.legacy_convention == "legacy_full_zero")]
shift = v[v.legacy_sv != "SV1"].apply(lambda r: int(r.best_observed_sv[2:]) == int(r.legacy_sv[2:]) - 1, axis=1)
out = dict(
    records=dict(all=A["n_records_all"], json_files=int(S2.iteration.nunique()),
                 aborted_duplicate_run_records=int(S2.aborted_duplicate_run.sum()),
                 in_183=T["n_183"], json_families_all=len(T["json_families_all"]), json_families_183=len(T["json_families_183"])),
    matrices_183=dict(kidney_A=T["n_183_kidney_A"], immune=T["n_183_immune"], misaligned_gene_rows=T["n_183_misaligned"],
                      spectral=T["n_183_spectral"], by_svd_convention=T["n_183_by_svd_convention"]),
    audit_183_rederived=A["rederive_183"],
    audit_all=A["all_records"], audit_all_excl_aborted=A["all_records_excluding_aborted_duplicate_run"],
    headline_survivors_all_records=dict(bonferroni=A["headline_survivors_all_records_bonferroni"],
                                        bonferroni_exact_gradient_p=A["headline_survivors_all_records_bonferroni_exact_gradient_p"],
                                        bonferroni_denominator_409=A["headline_survivors_denominator_409_bonferroni"],
                                        bonferroni_denominator_409_exact_gradient_p=A["headline_survivors_denominator_409_bonferroni_exact_gradient_p"],
                                        bh=A["headline_survivors_all_records_bh"]),
    ledger=dict(n_claims=len(L), reproduces_counts=L.reproduces.value_counts().to_dict()),
    crosswalk_kidney_A_legacy_full_zero=dict(
        frac_layer_sv_pairs_with_best_match_k_minus_1=float(shift.mean()),
        mean_abs_cos_best_SV2_12=float(v[v.legacy_sv != "SV1"].abs_cos_best.mean()),
        mean_abs_cos_best_SV1=float(v[v.legacy_sv == "SV1"].abs_cos_best.mean())),
    legacy_band_mapping=[{k: b[k] for k in ("legacy_band", "context", "legacy_convention", "mapped_observed_band",
                                            "mean_cos2_layer_avg", "identity_band_mean_cos2")} for b in M["bands"]],
    key_recomputations=dict(
        sv1_extracellular_OR=R["sv1_go_kidney_A_L11"]["GO0005615_top"], sv1_zero_rows_in_top52=R["sv1_go_kidney_A_L11"]["n_zero_row_genes_in_top52"],
        sv1_extracellular_observed_variant=R["sv1_go_kidney_A_L11"]["observed_variant"]["GO0005615_pos_pole"],
        named_genes=dict(true_observed=R["gene_counts_kidney_A"]["n_named_observed_true_mapping"], true_zero=R["gene_counts_kidney_A"]["n_named_zero_true_mapping"],
                         misaligned_zero=R["gene_counts_kidney_A"]["n_named_zero_misaligned_mapping"]),
        immune_labels=R["labels_immune"], tf_target_joint_mean=R["tf_target_immune_iter56"]["joint_foldmean_mean"],
        edge_auroc_L0=R["edge_immune_iter62"]["per_layer"][0]["sv57_U_neg1523"],
        celltype_slope=R["celltype_kidney_A"]["own_convention_true_mapping"]["slope_logit"],
        string_local_links=R["string_local_file"]["n_links"]),
)
(H / "results.json").write_text(json.dumps(out, indent=1, default=str))
print(json.dumps(out, indent=1, default=str)[:4000])
