"""Collect the verifier's checks into verify_summary.json (plus derived numbers quoted in the corrected SUMMARY)."""
import json, glob
from pathlib import Path
import numpy as np, pandas as pd
V = Path(__file__).resolve().parent; E5 = V.parent; R3 = E5.parent
out = dict(part1=json.loads((V / "verify_part1.json").read_text()),
           continuity_C=json.loads((V / "verify_continuity_C.json").read_text()),
           resid_boot=json.loads((V / "verify_resid_boot_kidney_A_0_kidney_B_0.json").read_text()), clogit={})
for f in sorted(V.glob("verify_clogit_*.json")):
    txt = f.read_text().strip()
    r = json.loads(txt.splitlines()[0]) if f.name.endswith("_point.json") else json.loads(txt)
    out["clogit"][r["tag"]] = {k: r[k] for k in r if k not in ("beta_adj_verify", "beta_adj_e5")}
ra = pd.read_csv(E5 / "residual_auroc.csv")
ra["w_ratio"] = (ra.auroc_resid_tfonly_hi - ra.auroc_resid_tfonly_lo) / (ra.auroc_resid_hi - ra.auroc_resid_lo)
out["tfonly_to_twoway_width_resid_auroc"] = {c: [float(g.w_ratio.min()), float(g.w_ratio.max())] for c, g in ra.groupby("context")}
e4 = {}
for c in ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]:
    e = pd.read_csv(R3 / f"E4a_edges/h3_primary_{c}.csv"); e = e[e.scorer.str.contains("PRIMARY")]
    r = (e.tfonly_ci_hi - e.tfonly_ci_lo) / (e.ci_hi - e.ci_lo)
    e4[c] = dict(tfonly_to_twoway_width=[float(r.min()), float(r.max())], null_sd_mean=float(e.null_sd.mean()), boot_sd_range=[float(e.boot_sd.min()), float(e.boot_sd.max())])
out["e4a_raw_auroc_interval_scale"] = e4
m = pd.read_csv(R3 / "E4a_edges/sens_macro.csv"); m = m[(m.context == "lung") & (m.scorer == "signed cosine full-512")]
p = pd.read_csv(R3 / "E4a_edges/h3_primary_lung.csv"); p = p[p.scorer.str.contains("PRIMARY")]
out["e4a_lung_macro_vs_pooled"] = dict(macro_range=[float(m.macro_auroc.min()), float(m.macro_auroc.max())],
                                       pooled_range=[float(p.auroc.min()), float(p.auroc.max())],
                                       L0_macro=float(m[m.layer == 0].macro_auroc.iloc[0]), L0_pooled=float(p[p.layer == 0].auroc.iloc[0]))
inf = pd.read_csv(E5 / "inference_scale.csv")
out["ratio_boot_sd_to_model_se_overall"] = [float(inf.ratio_boot_sd_to_model_se.min()), float(inf.ratio_boot_sd_to_model_se.max())]
p_lrt = {}
for c in ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]:
    d = pd.read_csv(E5 / f"clogit_{c}.csv"); d = d[d.layer != "joint_mean_over_layers"]; p_lrt[c] = [float(d.p_lrt_adj.min()), float(d.p_lrt_adj.max())]
out["p_lrt_adj_range"] = p_lrt
r = json.loads((E5 / "results.json").read_text())
out["reduced_model_beta"] = {c: v["reduced_model"]["beta"] for c, v in r["contexts"].items() if "reduced_model" in v}
(V / "verify_summary.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))
