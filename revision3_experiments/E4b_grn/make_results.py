"""Assemble results.json (headline numbers, every one read from grn_summary.csv / grn_deltas.csv /
direction_test.csv) and print a compact table for SUMMARY.md."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import e4b_common as C  # noqa: E402

H = C.HERE
s = pd.read_csv(H / "grn_summary.csv", dtype={"seed": str})
d = pd.read_csv(H / "grn_deltas.csv")
r = pd.read_csv(H / "direction_test.csv", dtype={"seed": str})


def f(x, k=4):
    return None if pd.isna(x) else round(float(x), k)


head = {}
for ctx in C.CONTEXTS:
    sc = s[s.context == ctx]
    if sc.empty: continue
    out = {}
    for method in C.METHODS + ["abs_coexpr"]:
        for tier in ("T1", "T2"):
            row = sc[(sc.method == method) & (sc.tier == tier) & sc.seed.isin(["ens", "-"])]
            if row.empty: continue
            row = row.iloc[0]
            e = dict(auroc=f(row.auroc), ci=[f(row.ci_lo), f(row.ci_hi)], upper95=f(row.upper95),
                     tf_only_ci=[f(row.tf_only_ci_lo), f(row.tf_only_ci_hi)], perm_p=f(row.perm_p_one_sided, 5),
                     macro_auroc=f(row.macro_auroc), frac_zero=f(row.frac_zero), frac_tied_within_tf=f(row.frac_tied_within_tf),
                     frac_tied_pos_neg=f(row.frac_tied_pos_neg), n_cells=int(row.n_cells), n_pairs=int(row.n_pairs),
                     n_pos=int(row.n_pos), n_tfs=int(row.n_tfs))
            if method in C.METHODS:
                e["n_seeds"] = int(row.n_seeds); e["seeds"] = str(row.seeds)
                v = sc[(sc.method == method) & (sc.tier == tier) & (sc.seed == "seed_mean_sd")]
                if len(v):
                    v = v.iloc[0]
                    e["seed_mean_auroc"] = f(v.auroc); e["seed_sd_auroc"] = f(v.seed_sd, 5)
                    e["seed_range"] = [f(v.seed_min), f(v.seed_max)]
                dc = d[(d.context == ctx) & (d.method == method) & (d.tier == tier) & (d.comparison == "grn_minus_abs_coexpr_same_cells")]
                if len(dc):
                    dc = dc.iloc[0]
                    e["delta_vs_abs_coexpr"] = dict(delta=f(dc.delta), ci=[f(dc.ci_lo), f(dc.ci_hi)], upper95=f(dc.upper95))
                ds = d[(d.context == ctx) & (d.method == method) & (d.tier == tier) & (d.comparison == "grn_minus_scgpt_cos")]
                if len(ds):
                    e["delta_vs_scgpt_cos"] = dict(
                        per_layer=[dict(layer=int(x.comparator_layer), delta=f(x.delta), ci=[f(x.ci_lo), f(x.ci_hi)])
                                   for x in ds.itertuples()],
                        n_layers_ci_above_0=int((ds.ci_lo > 0).sum()), n_layers_ci_below_0=int((ds.ci_hi < 0).sum()),
                        delta_range=[f(ds.delta.min()), f(ds.delta.max())])
            out[f"{method}_{tier}"] = e
    sg = sc[sc.method == "scgpt_cos"]
    out["scgpt_cos"] = [dict(layer=int(x.layer), auroc=f(x.auroc), ci=[f(x.ci_lo), f(x.ci_hi)], upper95=f(x.upper95))
                        for x in sg.itertuples()]
    dt = d[(d.context == ctx) & (d.comparison == "T2_minus_T1_same_method")]
    out["T2_minus_T1"] = {x.method: dict(delta=f(x.delta), ci=[f(x.ci_lo), f(x.ci_hi)]) for x in dt.itertuples()}
    meta = json.loads((H / "eval" / f"{ctx}_meta.json").read_text())
    out["meta"] = meta
    head[ctx] = out

direction = {}
for x in r.itertuples():
    if x.seed not in ("ens", "-"): continue
    direction.setdefault(x.context, {})[f"{x.method}_{x.tier}"] = dict(
        n_pairs=int(x.n_pairs), n_untied=int(x.n_untied), frac_tied=f(x.frac_tied), frac_forward_gt=f(x.frac_forward_gt),
        sign_test_p=f(getattr(x, "sign_test_p_two_sided", np.nan), 5),
        cp95=[f(getattr(x, "cp_lo", np.nan)), f(getattr(x, "cp_hi", np.nan))],
        src_cluster95=[f(getattr(x, "src_cluster_lo", np.nan)), f(getattr(x, "src_cluster_hi", np.nan))])

# ---- run inventory (which GRN fits exist; wall time under load) -> grn_runs.csv
runs = []
for ctx in C.CONTEXTS:
    for tier in ("T1", "T2"):
        for method in C.METHODS:
            for seed in C.SEEDS[tier]:
                fz = H / "partial" / f"{ctx}_{tier}_{method}_s{seed}.npz"
                row = dict(context=ctx, tier=tier, method=method, seed=seed, done=fz.exists())
                if fz.exists():
                    z = np.load(fz)
                    row.update(n_cells=int(z["n_cells"]), n_targets=int(z["targets"].size), n_regulators=int(z["regs"].size),
                               n_targets_unfitted=int((z["fitted"] == 0).sum()), wall_seconds_this_process=float(z["seconds"]),
                               median_n_iter=(float(np.median(z["n_iter"])) if method == "grnboost2" else None))
                else:
                    row["reason_missing"] = "not completed (compute reduction; see results.json compute_plan_resume)"
                runs.append(row)
pd.DataFrame(runs).to_csv(H / "grn_runs.csv", index=False)

# ---- cross-check: our scGPT signed-cosine AUROCs vs E4a's primary H3 table (identical pairs, same estimator)
xchk = {}
for ctx in C.CONTEXTS:
    fa = C.R3 / "E4a_edges" / f"h3_primary_{ctx}.csv"
    sg = s[(s.context == ctx) & (s.method == "scgpt_cos")]
    if fa.exists() and len(sg):
        a = pd.read_csv(fa); a = a[a.scorer.str.contains("PRIMARY")].set_index("layer").auroc
        m = sg.set_index("layer").auroc
        xchk[ctx] = dict(max_abs_diff_point_auroc=float((m - a.reindex(m.index)).abs().max()), n_layers=int(len(m)))

res = dict(
    experiment="E4b_grn", convention3_sha256=C.convention_hash(),
    convention3_sha256_verified=__import__("hashlib").sha256((C.R3 / "CONVENTION3.md").read_bytes()).hexdigest() == C.convention_hash(),
    status="sensitivity analysis (CONVENTION3 §5, E4 extras: GENIE3 / GRNBoost2, direction test); no pass/fail",
    spec=C.__doc__,
    compute_plan_resume=("Resumed twice (session cut-off; machine reboot at ~22:47). A temporary plan to drop kidney_B T2 "
                         "(load average ~470) was cancelled at 19:18 when load fell; T2 seed 0 then finished for lung, "
                         "lung_external and kidney_B, and the optional T2 seed 1 jobs (queues E and F, the latter restarted "
                         "from batch checkpoints after the reboot) restored the full design. grn_runs.csv lists every "
                         "(context, tier, method, seed) fit and whether it is complete; every completed fit is used."),
    convention3_errata_item3=("Not applicable: E4b reports no vertex-subsampling interval. All edge-level intervals are the "
                              "two-way (TF x target) pigeonhole cluster bootstrap of CONVENTION3 §6; the direction-test intervals "
                              "are exact Clopper-Pearson and a source-regulator cluster bootstrap."),
    scgpt_crosscheck_vs_E4a=xchk,
    esm2_note="No ESM2 embedding is used in E4b.",
    symmetric_scorers_note=("Signed cosine and |co-expression| are symmetric (score(i,j) = score(j,i)); in the direction test "
                            "every pair is tied for them by construction (n_untied = 0), so they cannot carry direction."),
    compute_plan_original=("Initial plan: T1 5 seeds (0-4) in every context; T2 2 seeds (0,1) in immune and kidney_A and 1 seed (0) in lung, "
                  "lung_external, kidney_B (fallback pre-authorised by the task when T2 x 2 seeds x 5 contexts exceeds ~8 h; "
                  "measured: lung T2 = 35 min GENIE3 + 102 min GRNBoost2 per seed on an uncontended estimate, ~3x slower "
                  "under the machine's actual load). Seed 1 for the three large contexts was added later as optional work."),
    seed_design_final={f"{x['context']}_{x['tier']}_{x['method']}": [y["seed"] for y in runs if y["done"] and (y["context"], y["tier"], y["method"]) == (x["context"], x["tier"], x["method"])] for x in runs},
    bootstrap=dict(B=2000, type="two-way TF x target pigeonhole (stats3.StratifiedPairs.bootstrap); TF-only as check",
                   paired="identical weight draws for every score within a context"),
    permutation=dict(N=2000, type="within-TF label permutation (stats3.permute_within_tf), one-sided",
                     resolution="stats3.assert_resolution(2000, 0.05, 6) passed (floor 5.0e-4 < 0.05/6)"),
    direction_test_inference="exact two-sided binomial sign test on untied pairs, Clopper-Pearson 95% CI, plus a "
                             "source-regulator cluster bootstrap (B = 2,000); exact test, no permutation floor",
    headline=head, direction=direction,
    files=dict(runs="grn_runs.csv", summary="grn_summary.csv", deltas="grn_deltas.csv", direction="direction_test.csv",
               scores="grn_scores_<ctx>_<tier>_<seed>.npz", per_context="eval/<ctx>_*.csv"))
(H / "results.json").write_text(json.dumps(res, indent=1, default=str))
print("results.json written")

# compact table for SUMMARY.md
pd.set_option("display.width", 250)
rows = []
for ctx, o in head.items():
    for k, e in o.items():
        if not isinstance(e, dict) or "auroc" not in e: continue
        rows.append(dict(ctx=ctx, score=k, auroc=e["auroc"], ci=f"[{e['ci'][0]:.3f}, {e['ci'][1]:.3f}]", up95=e["upper95"],
                         perm_p=e["perm_p"], seedSD=e.get("seed_sd_auroc"), nseed=e.get("n_seeds"),
                         d_coex=(e.get("delta_vs_abs_coexpr") or {}).get("delta"),
                         d_coex_ci=str((e.get("delta_vs_abs_coexpr") or {}).get("ci")),
                         d_coex_up=(e.get("delta_vs_abs_coexpr") or {}).get("upper95"),
                         d_scgpt_range=str((e.get("delta_vs_scgpt_cos") or {}).get("delta_range")),
                         scgpt_above=(e.get("delta_vs_scgpt_cos") or {}).get("n_layers_ci_above_0"),
                         scgpt_below=(e.get("delta_vs_scgpt_cos") or {}).get("n_layers_ci_below_0"),
                         zero=e["frac_zero"], tiedPN=e["frac_tied_pos_neg"]))
print(pd.DataFrame(rows).to_string(index=False))
print(r[r.seed.isin(["ens"]) | r.method.isin(["abs_coexpr", "scgpt_cos_L0"])][
    ["context", "method", "tier", "seed", "n_pairs", "n_untied", "frac_forward_gt", "sign_test_p_two_sided", "cp_lo", "cp_hi",
     "src_cluster_lo", "src_cluster_hi"]].to_string(index=False))
