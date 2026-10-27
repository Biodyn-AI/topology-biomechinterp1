"""E1 post-processing (no new random draws; reads the CSVs written by e1_compression.py).

* per_context_layer.csv: adds basic (reverse-percentile) bootstrap intervals for ER, PR and SV1 variance and
  a flag for whether the percentile interval covers the point estimate. The gene bootstrap (resampling with
  replacement) duplicates ~37% of rows, which biases ER downward by 1-7% while its spread is < 1%: neither the
  percentile nor the basic interval of ER covers the estimate, so NO bootstrap interval of ER is reported as a CI
  (diagnostic columns only). Sampling variability of ER is taken from the without-replacement draws in
  matched.csv instead. SV1-variance percentile intervals cover the estimate in every context-layer and are used.
* core_table.csv: compact main-text table for the five core contexts.
* detection_summary.csv: per context and layer, the SV1 Spearman with log n_obs and the SV1..10 direction
  that tracks detection most strongly.
* gap_summary.csv: per context, number of near-degenerate adjacent pairs among k = 1..11 x 12 layers and
  the pairs among SV1..SV4 (the directions a main-text claim could name).
* null_direction.json: the last singular direction (near-zero singular value) vs the uniform vector.
* results.json: updated with these summaries.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import stats3 as S  # noqa: E402

CORE = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"

per = pd.read_csv(HERE / "per_context_layer.csv")
for k in ("er", "pr", "sv1_var"):
    per[f"{k}_basic_lo"] = 2 * per[k] - per[f"{k}_boot_hi"]
    per[f"{k}_basic_hi"] = 2 * per[k] - per[f"{k}_boot_lo"]
    per[f"{k}_pct_ci_covers_estimate"] = (per[f"{k}_boot_lo"] <= per[k]) & (per[k] <= per[f"{k}_boot_hi"])
per.to_csv(HERE / "per_context_layer.csv", index=False)

tr = pd.read_csv(HERE / "depth_trend.csv").set_index("context")
sh = pd.read_csv(HERE / "shuffle.csv")
mt = pd.read_csv(HERE / "matched.csv")
rows = []
for c in CORE:
    p0 = per[(per.context == c) & (per.layer == 0)].iloc[0]
    p11 = per[(per.context == c) & (per.layer == 11)].iloc[0]
    m0 = mt[(mt.context == c) & (mt.layer == 0) & (mt.scheme == "pooled")].iloc[0]
    m11 = mt[(mt.context == c) & (mt.layer == 11) & (mt.scheme == "pooled")].iloc[0]
    s11 = sh[(sh.context == c) & (sh.layer == 11)].iloc[0]
    t = tr.loc[c]
    rows.append(dict(
        context=c, n_V=int(p0.n_V),
        ER_L0=p0.er, ER_L11=p11.er,
        ER_fold=p0.er / p11.er,
        PR_L0=p0.pr, PR_L11=p11.pr, PR_fold=p0.pr / p11.pr,
        SV1var_L0=p0.sv1_var, SV1var_L0_lo=p0.sv1_var_boot_lo, SV1var_L0_hi=p0.sv1_var_boot_hi,
        SV1var_L11=p11.sv1_var, SV1var_L11_lo=p11.sv1_var_boot_lo, SV1var_L11_hi=p11.sv1_var_boot_hi,
        cum10var_L0=p0.cum10_var, cum10var_L11=p11.cum10_var,
        ER_shuffle_L11_mean=s11.er_shuf_mean, ER_shuffle_L11_sd=s11.er_shuf_sd,
        ER_matched1300_L0=m0.er_mean, ER_matched1300_L0_sd=m0.er_sd,
        ER_matched1300_L11=m11.er_mean, ER_matched1300_L11_sd=m11.er_sd,
        ER_matched_fold=m0.er_mean / m11.er_mean,
        rho_ER_layer=t.spearman_er_layer, rho_PR_layer=t.spearman_pr_layer,
        ER_strictly_decreasing=bool(t.er_strictly_decreasing), ER_n_increases=int(t.er_n_increases),
        PR_strictly_decreasing=bool(t.pr_strictly_decreasing), PR_n_increases=int(t.pr_n_increases)))
core = pd.DataFrame(rows)
core.to_csv(HERE / "core_table.csv", index=False)

det = pd.read_csv(HERE / "detection_corr.csv")
ds = []
for (c, L), g in det.groupby(["context", "layer"], sort=False):
    top = g.loc[g.abs_spearman.idxmax()]
    sv1 = g[g.k == 1].iloc[0]
    ds.append(dict(context=c, layer=L, sv1_spearman=sv1.spearman, sv1_abs=sv1.abs_spearman,
                   sv1_var_frac=sv1.var_frac, top_k=int(top.k), top_abs_spearman=top.abs_spearman,
                   n_dirs_abs_ge_0p3=int((g.abs_spearman >= 0.3).sum()),
                   var_frac_on_dirs_abs_ge_0p3=float(g.loc[g.abs_spearman >= 0.3, "var_frac"].sum())))
ds = pd.DataFrame(ds)
ds.to_csv(HERE / "detection_summary.csv", index=False)

gaps = pd.read_csv(HERE / "gaps.csv")
gs = []
for c, g in gaps.groupby("context", sort=False):
    dg = g[g.near_degenerate]
    gs.append(dict(context=c, n_pairs=len(g), n_near_degenerate=len(dg), frac=len(dg) / len(g),
                   n_near_degenerate_k1to3=int((dg.k <= 3).sum()),
                   pairs_k1to3=";".join(f"L{int(r.layer)}:SV{int(r.k)}-SV{int(r.k) + 1}" for r in dg[dg.k <= 3]
                                        .itertuples())))
gs = pd.DataFrame(gs)
gs.to_csv(HERE / "gap_summary.csv", index=False)

# last singular direction ~ uniform vector (affine 511-dim constraint)
nd = {}
for c in CORE:
    emb = np.load(EMB / f"{c}.npy", mmap_mode="r")
    oi = np.load(R3 / "cache" / c / "obs_idx.npy")
    nd[c] = {}
    for L in (0, 11):
        _, Sv, Vt = S.svd_observed(emb, L, oi)
        nd[c][L] = dict(rel_last_sv=float(Sv[-1] / Sv[0]),
                        abs_cos_last_dir_uniform=float(abs(Vt[-1].sum()) / np.sqrt(Vt.shape[1])),
                        rel_second_last_sv=float(Sv[-2] / Sv[0]))
(HERE / "null_direction.json").write_text(json.dumps(nd, indent=1))

res = json.loads((HERE / "results.json").read_text())
res["core_table_csv"] = core.round(4).to_dict("records")
res["gap_summary"] = gs.to_dict("records")
res["detection_summary_core_L0_L11"] = ds[ds.context.isin(CORE) & ds.layer.isin([0, 11])].round(4).to_dict("records")
res["null_direction"] = nd
res["bootstrap_note"] = ("percentile CI of ER covers the estimate in "
                         f"{int(per.er_pct_ci_covers_estimate.sum())}/{len(per)} context-layers and the basic CI in 0; ER gene-bootstrap is biased (duplicated rows), so no bootstrap CI of ER is reported; ER sampling sd comes from matched.csv")
(HERE / "results.json").write_text(json.dumps(res, indent=1, default=float))
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(core.round(3).T.to_string())
print(gs.to_string())
print(ds[ds.context.isin(CORE) & ds.layer.isin([0, 11])].round(3).to_string())
print(json.dumps(nd))
print(res["bootstrap_note"])
