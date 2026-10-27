"""E4c verifier part 9: collect the verifier's numbers into verify/verification.json and add a pointer block to ../results.json."""
import json
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; E = H.parent; R3 = E.parent
h = (R3 / "FROZEN.txt").read_text().split()[2]
NB = pd.read_csv(H / "v6_nb_counts.csv"); DC = pd.read_csv(H / "v6_direction_cluster.csv")
DID = pd.read_csv(H / "v1_rewired_did.csv"); IND = pd.read_csv(H / "v7_independent_vs_implementer.csv")
S8 = pd.read_csv(H / "v8_did_sources.csv"); V1 = json.load(open(H / "v1_summary.json"))
sel = NB[(NB.source == "scgpt") & NB.feature.isin(["bilinear", "tg_bilinear", "concat"])]
out = dict(
    convention3_sha256=h, verifier="independent adversarial verifier (v1-v3 by a first verifier cut off by the reboot; v4-v9 by the second)",
    code_bugs_changing_results="none found",
    alignment_folds_rewiring=dict(c3_max_gene_overlap=V1["c3_max_gene_overlap"], c2_max_tf_overlap=V1["c2_max_tf_overlap"],
                                  ref_cosine_max_absdiff=V1["ref_cosine_max_absdiff"], null_recompute_max_absdiff_p=V1["null_recompute_max_absdiff_p"]),
    independent_refits=IND.to_dict("records"),
    nb_ratio_sensitivity_counts=sel.to_dict("records"),
    direction_cluster_bootstrap=DC.to_dict("records"),
    did_scgpt=DID[DID.test == "DiD_bilinear_minus_concat"].to_dict("records"),
    did_layer_free_sources=S8.to_dict("records"))
json.dump(out, open(H / "verification.json", "w"), indent=1, default=float)
R = json.load(open(E / "results.json"))
R["verification"] = dict(file="verify/verification.json", convention3_sha256=h, code_bugs_changing_results="none found",
    corrections=["direction-test CIs: gene-cluster bootstrap replaces the independence binomial CI (verify/v6_direction_cluster.csv)",
                 "Nadeau-Bengio ratio: gene-level 0.25 reported beside the pair-level ~0.0625 for C3 (verify/v6_nb_counts.csv)",
                 "rewired null is not relation-clean for additive probes; DiD bilinear-minus-concat added (verify/v1_rewired_did.csv, v8_did_sources.csv)",
                 "3 (not 2) M1 null rows differ from the stand-alone true-label run by > 0.01",
                 "tuned-minus-fixed means +0.007/-0.002 are C3 only", "full list C1-C11 in SUMMARY.txt section V"])
json.dump(R, open(E / "results.json", "w"), indent=1, default=float)
print("written")
