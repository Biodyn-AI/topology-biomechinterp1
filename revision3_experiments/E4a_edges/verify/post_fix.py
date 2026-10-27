"""After fix_ties.py / fix_nobs_strata.py: compile the audit tables and the verifier block of results.json.
Run AFTER `python e4a_aggregate.py` has regenerated the CSVs from the corrected work/ files.
Writes verify/cx_tie_audit.csv, verify/old_vs_new_summary.csv, verify/old_vs_new_rank1.csv and adds results.json['verifier']."""
import json
from pathlib import Path
import pandas as pd, numpy as np

V = Path(__file__).resolve().parent; E4 = V.parent; OLD = V / "pre_verify_outputs"
rows = [json.loads(p.read_text()) for p in sorted((V / "fixed").glob("*_audit.json")) if not p.name.endswith("_chip_audit.json") and not p.name.startswith("._")]
aud = pd.DataFrame(rows)
for p in sorted((V / "fixed").glob("*_chip_audit.json")):
    if p.name.startswith("._"): continue
    a = json.loads(p.read_text()); i = aud.index[aud.context == a["context"]][0]
    for k in ("chip_cx_ext_old", "chip_cx_ext_new", "chip_cx_3k_old", "chip_cx_3k_new"): aud.loc[i, k + "_bootstrap_rerun"] = a[k]
aud.to_csv(V / "cx_tie_audit.csv", index=False)
o = pd.read_csv(OLD / "h3_summary.csv"); n = pd.read_csv(E4 / "h3_summary.csv")
cols = ["auroc_cx_ext", "auroc_cx_3k", "delta_cxext_min", "delta_cxext_max", "max_delta_cxext_ub95", "max_delta_cx3k_ub95", "verdict", "verdict_3k_comparator",
        "auroc_min", "auroc_max", "max_ub95", "min_p_maxT"]
m = o[["context"] + cols].merge(n[["context"] + cols], on="context", suffixes=("_old", "_new"))
m.to_csv(V / "old_vs_new_summary.csv", index=False)
so = pd.read_csv(OLD / "sens_scorers.csv"); sn = pd.read_csv(E4 / "sens_scorers.csv")
r1 = so[so.scorer == "rank-1 cosine"][["context", "layer", "auroc", "ci_lo", "ci_hi"]].merge(
    sn[sn.scorer == "rank-1 cosine"][["context", "layer", "auroc", "ci_lo", "ci_hi"]], on=["context", "layer"], suffixes=("_old", "_new"))
r1.to_csv(V / "old_vs_new_rank1.csv", index=False)
res = json.loads((E4 / "results.json").read_text())
p2 = json.loads((V / "verify_part2.json").read_text()); p3 = json.loads((V / "verify_part3.json").read_text())
res["verifier"] = dict(
    corrections=[
        "|co-expression| comparator: float-noise tie breaking replaced by exact ties (|r| rounded to 12 decimals; 10/11/13 give identical AUROCs). "
        "Recomputed with the original seeds/weight streams: rows_coexp, the PRIMARY rows' delta columns (all 16 contexts), |cos| and rank-1 rows, "
        "comparator rows of sens_macro, both TF-detection strata in every context (incl. the contextual analysis), ChIP rows in lung and lung_external. "
        "Rows flagged delta_cx_basis='original float-tie-broken ...' (rank-k>=2, bands, mapped bands, PMID/mode strata, sens_matching, sens_legacy_pooled, "
        "ChIP outside the lung pair) keep the original comparator; their comparator point shift is in verify/cx_tie_audit.csv.",
        "rank-1 cosine: exact sign(x_tf)*sign(x_tg) (stats3.l2norm_rows epsilon pitfall).",
        "Added strata with the CONVENTION3 section 2 detection definition (n_obs = cells whose scGPT input contained the TF): "
        "'TF n_obs (scGPT input) >= / < 10% of extraction cells [CONVENTION3 detection]'.",
        "power.json not re-run: it uses the comparator only through a_cx, rho and tau_cx (changes <= 0.0011 in AUROC)."],
    verdicts_changed=[c for c, a, b in zip(m.context, m.verdict_old, m.verdict_new) if a != b] +
                     [c + " (3k)" for c, a, b in zip(m.context, m.verdict_3k_comparator_old, m.verdict_3k_comparator_new) if a != b],
    independent_recomputation=dict(part1="verify/verify_results.json", part2="verify/verify_part2.json", part3="verify/verify_part3.json"),
    kidney_ckpt_kidney_L2_delta_ub95=dict(replay=p3["replay_seed"]["delta_ub95"], indep1=p3["indep_seed_1"]["delta_ub95"], indep2=p3["indep_seed_2"]["delta_ub95"]),
    power_size_own_simulation={k: p2[k] for k in ("G_lung", "G_kidney_B")})
(E4 / "results.json").write_text(json.dumps(res, indent=1, default=float))
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(m[["context", "auroc_cx_ext_old", "auroc_cx_ext_new", "auroc_cx_3k_old", "auroc_cx_3k_new", "max_delta_cxext_ub95_old", "max_delta_cxext_ub95_new",
         "max_delta_cx3k_ub95_old", "max_delta_cx3k_ub95_new"]].round(4).to_string())
print(m[["context", "verdict_old", "verdict_new", "verdict_3k_comparator_new"]].to_string())
print(r1[r1.context.isin(["lung", "lung_external", "immune"])].round(4).to_string())
print("verdicts changed:", res["verifier"]["verdicts_changed"])
