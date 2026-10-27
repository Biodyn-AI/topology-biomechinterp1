"""Assemble E3 tables, figures and results.json from partial/ outputs (no new statistics)."""
import json
import numpy as np, pandas as pd
import e3_lib as L

H = L.HERE; P = H / "partial"
ctxs = [c for c in L.CONTEXTS if (P / f"primary_{c}.json").exists()]
prim = {c: json.loads((P / f"primary_{c}.json").read_text()) for c in ctxs}
BASE_ORDER = ["log_nobs", "log_mean_expr_extraction", "log_mean_expr_population", "static_scgpt",
              "static_scgpt_uncentred", "static_esm2", "static_esm2_uncentred"]

# ------------------------------------------------------------ baselines.csv (+ paired)
b = pd.DataFrame([dict(context=c, **r) for c in ctxs for r in prim[c]["baselines"]])
b["order"] = b.baseline.map({k: i for i, k in enumerate(BASE_ORDER)})
b = b.sort_values(["order", "context"], key=lambda s: s.map({c: i for i, c in enumerate(L.CONTEXTS)}) if s.name == "context" else s)
LABEL = {"log_nobs": "log n_obs (detection)", "log_mean_expr_extraction": "log mean expression (extraction cells)",
         "log_mean_expr_population": "log mean expression (population)",
         "static_scgpt": "scGPT static gene-token embedding (V-centred)", "static_scgpt_uncentred": "scGPT static gene-token embedding (uncentred)",
         "static_esm2": "ESM2-15B (UCE per-protein), V-centred", "static_esm2_uncentred": "ESM2-15B (UCE per-protein), uncentred"}
b.insert(2, "model_label", b.baseline.map(LABEL))
b.drop(columns="order").to_csv(H / "baselines.csv", index=False)
bp = pd.concat([pd.read_csv(P / f"baselines_paired_{c}.csv") for c in ctxs if (P / f"baselines_paired_{c}.csv").exists()])
bp.insert(1, "static_label", bp.static.map({"scgpt": "scGPT static gene-token embedding", "esm2": "ESM2-15B (UCE per-protein)"}))
bp.to_csv(H / "baselines_paired.csv", index=False)

# ------------------------------------------------------------ h2_summary.csv (baseline columns first)
rows = []
for c in ctxs:
    d = prim[c]; lay = pd.DataFrame(d["layers"]); bb = {r["baseline"]: r for r in d["baselines"]}
    r = dict(context=c, tier={"lung": "primary", "lung_external": "replication"}.get(c, "sensitivity"),
             discovery_context=c in L.DISCOVERY)
    for nm in BASE_ORDER:
        if nm in bb:
            r[f"base_{nm}_auroc"] = bb[nm]["auroc"]; r[f"base_{nm}_strat"] = bb[nm]["strat_auroc"]
    lg = pd.DataFrame(d["logistic"]); lg = lg[lg.expr_covariate == "extraction"]
    best = int(lay.strat_auroc.idxmax())
    r.update(n_V=d["n_V"], n_tf=d["n_tf"], strat_auroc_min=lay.strat_auroc.min(), strat_auroc_max=lay.strat_auroc.max(),
             best_layer=best, best_strat_lo=lay.strat_lo[best], best_strat_hi=lay.strat_hi[best],
             best_unstrat_auroc=lay.auroc[best], best_z=lay.z[best], min_p_maxT=lay.p_maxT.min(),
             n_layers_maxT_lt_0p05=int((lay.p_maxT < 0.05).sum()), n_layers_ci_excl_0p5=int((lay.strat_lo > 0.5).sum()),
             null_mean_range=f"{min(x['null']['mean'] for x in d['layers']):.4f}-{max(x['null']['mean'] for x in d['layers']):.4f}",
             logit_beta_per_sd_min=lg.beta_score_per_sd.min(), logit_beta_per_sd_max=lg.beta_score_per_sd.max(),
             logit_LR_p_max=lg.LR_p.max(), N_perm=10000, conv_hash=d["conv_hash"])
    rows.append(r)
h2s = pd.DataFrame(rows); h2s.to_csv(H / "h2_summary.csv", index=False)
pl, pe = prim["lung"], prim["lung_external"]
lungP = min(x["p_maxT"] for x in pl["layers"]); extP = min(x["p_maxT"] for x in pe["layers"])
rep_layers = [l for l in range(12) if pl["layers"][l]["p_maxT"] < 0.05 and pe["layers"][l]["strat_lo"] > 0.5
              and pe["layers"][l]["strat_auroc"] > 0.5]
decision = dict(rule="H2 supported iff max-T p < 0.05 in lung AND lung_external (CONVENTION3 §4)",
                lung_min_p_maxT=lungP, lung_external_min_p_maxT=extP, supported=bool(lungP < 0.05 and extP < 0.05),
                replication_criterion_layers=rep_layers,
                replication_note="§1: same sign in lung_external and its interval (gene bootstrap 95%) excludes 0.5")

# ------------------------------------------------------------ logistic / restricted / trrust
pd.DataFrame([r for c in ctxs for r in prim[c]["logistic"]]).to_csv(H / "logistic_adjustment.csv", index=False)
rv = pd.DataFrame([r for c in ctxs for r in prim[c]["restricted"]]); rv.insert(2, "row_type", "H2_layer")
rb = pd.read_csv(P / "restricted_baselines.csv")
pd.concat([rb, rv], ignore_index=True).sort_values(["context", "min_nobs"], key=lambda s: s.map({c: i for i, c in enumerate(L.CONTEXTS)}) if s.name == "context" else s,
    kind="stable").to_csv(H / "restricted_V.csv", index=False)
tt = pd.DataFrame([r for c in ctxs for r in prim[c]["trrust_task"]])
tt["rt_order"] = (tt.row_type != "baseline").astype(int)
tt = tt.sort_values(["rt_order"], kind="stable").drop(columns="rt_order")
tt.to_csv(H / "trrust_task.csv", index=False)

# ------------------------------------------------------------ band scans: baseline rows prepended; maxT summary
bm = []
for task, suf in (("lambert", ""), ("trrust", "_trrust")):
    for c in L.CONTEXTS:
        f = H / f"bandscan_{c}{suf}.csv"; j = P / f"bandscan_{c}_{task}.json"
        if not f.exists() or not j.exists(): continue
        df = pd.read_csv(f)
        if "row_type" not in df.columns:
            if task == "lambert":
                base = [dict(row_type="baseline", context=c, task=task, band=LABEL.get(r["baseline"], r["baseline"]), strat_auroc_cos=r["strat_auroc"],
                             auroc_cos_unstrat=r["auroc"]) for r in prim[c]["baselines"]] if c in prim else []
                base += [dict(row_type="baseline", context=c, task=task, band="full512_LOO_cosine", layer=x["layer"],
                              strat_auroc_cos=x["strat_auroc"], auroc_cos_unstrat=x["auroc"], p_maxT_cos=x["p_maxT"])
                         for x in prim[c]["layers"]] if c in prim else []
            else:
                base = [dict(row_type="baseline", context=c, task=task, band=r.get("baseline") or "full512_LOO_cosine",
                             layer=r.get("layer"), strat_auroc_cos=r["strat_auroc"], auroc_cos_unstrat=r.get("auroc"),
                             p_maxT_cos=r.get("p_maxT")) for r in prim[c]["trrust_task"]] if c in prim else []
            df.insert(0, "row_type", "band_cell")
            df = pd.concat([pd.DataFrame(base), df], ignore_index=True)
            df.to_csv(f, index=False)
        js = json.loads(j.read_text())
        for v in ("cos", "lda"):
            s = js["summary"][v]
            bm.append(dict(context=c, task=task, variant=v, n=js["n"], n_tf=js["n_tf"], N_perm=js["N_perm"],
                           n_cells=s["n_cells"], best_z_layer=s["best_cell"]["layer"], best_z_band=s["best_cell"]["band"],
                           best_z_strat_auroc=s["best_cell"]["strat_auroc"], best_z=s["best_cell"]["z"],
                           best_z_p_maxT=s["best_cell"]["p_maxT"], best_auroc_layer=s["best_auroc_cell"]["layer"],
                           best_auroc_band=s["best_auroc_cell"]["band"], best_auroc=s["best_auroc_cell"]["strat_auroc"],
                           n_cells_maxT_lt_0p05=s["n_cells_maxT_lt_0p05"], null_max_auroc_q97_5=s["max_auroc_null"]["q"]["97.5"],
                           discovery_context=c in L.DISCOVERY))
bmt = pd.DataFrame(bm)
if len(bmt):
    # baseline rows first: the full-dimensional primary statistic of the same context
    base = [dict(context=c, task="lambert", variant="baseline_full512_primary", best_z_layer=int(np.argmax([x["z"] for x in prim[c]["layers"]])),
                 best_z_strat_auroc=max(x["strat_auroc"] for x in prim[c]["layers"]), best_z_p_maxT=min(x["p_maxT"] for x in prim[c]["layers"]),
                 best_auroc=max(x["strat_auroc"] for x in prim[c]["layers"]))
            for c in ctxs]
    pd.concat([pd.DataFrame(base), bmt], ignore_index=True).to_csv(H / "bandscan_maxT.csv", index=False)

# ------------------------------------------------------------ crossfit
cf = []
for task in ("lambert", "trrust"):
    for c in L.CONTEXTS:
        f = P / f"crossfit_{c}_{task}.csv"
        if not f.exists(): continue
        d = pd.read_csv(f); vc = d.chosen_band.astype(str) + "@L" + d.chosen_layer.astype(str)
        top = vc.value_counts().head(3)
        q = lambda x, a: float(np.quantile(x, a))
        cf.append(dict(context=c, task=task, n_splits=len(d), heldout_median=d.heldout_B.median(), heldout_lo=q(d.heldout_B, .025),
                       heldout_hi=q(d.heldout_B, .975), insample_median=d.insample_A.median(),
                       optimism_median=float((d.insample_A - d.heldout_B).median()),
                       full512_heldout_median=d.full_heldout_B.median(), full512_heldout_lo=q(d.full_heldout_B, .025),
                       full512_heldout_hi=q(d.full_heldout_B, .975), full512_insample_median=d.full_insample_A.median(),
                       band_minus_full512_median=float((d.heldout_B - d.full_heldout_B).median()),
                       band_minus_full512_lo=q(d.heldout_B - d.full_heldout_B, .025), band_minus_full512_hi=q(d.heldout_B - d.full_heldout_B, .975),
                       top_choices="; ".join(f"{k} ({v})" for k, v in top.items())))
pd.DataFrame(cf).to_csv(H / "crossfit.csv", index=False)

# ------------------------------------------------------------ legacy bands
lg = [pd.read_csv(f) for f in sorted(P.glob("legacy_*_*.csv"))]
if lg:
    lgd = pd.concat(lg, ignore_index=True)
    lgd.insert(lgd.columns.get_loc("screen_derived"), "label",
               np.where(lgd.discovery_contaminated, "screen-derived; discovery-contaminated", "screen-derived"))
    base = []
    for c in ctxs:
        for r in prim[c]["baselines"][:3]:
            base.append(dict(row_type="baseline", context=c, task="lambert", legacy_band=r["baseline"], strat_auroc=r["strat_auroc"],
                             p_raw=r["strat_p_raw"], N_perm=r["N_perm"]))
        for x in prim[c]["layers"]:
            base.append(dict(row_type="baseline", context=c, task="lambert", legacy_band="full512_LOO_cosine (primary H2)", layer=x["layer"],
                             strat_auroc=x["strat_auroc"], p_raw=x["p_raw"], N_perm=10000))
        for r in prim[c]["trrust_task"]:
            base.append(dict(row_type="baseline", context=c, task="trrust", legacy_band=r.get("baseline") or "full512_LOO_cosine",
                             layer=r.get("layer"), strat_auroc=r["strat_auroc"], p_raw=r["p_raw"], N_perm=r["N_perm"]))
    lgd.insert(0, "row_type", "legacy_band")
    pd.concat([pd.DataFrame(base), lgd], ignore_index=True).to_csv(H / "legacy_bands.csv", index=False)

# ------------------------------------------------------------ singular structure: E1 gaps + |Spearman| with log n_obs
sg = [pd.read_csv(f) for f in sorted(P.glob("singular_*.csv"))]
if sg:
    sgd = pd.concat(sg, ignore_index=True)
    g = pd.read_csv(L.R3 / "E1_compression/gaps.csv")[["context", "layer", "k", "ratio", "pct_lo", "pct_hi", "near_degenerate"]]
    g = g.rename(columns={"ratio": "gap_ratio_Sk_Sk1", "pct_lo": "gap_lo", "pct_hi": "gap_hi", "near_degenerate": "near_degenerate_with_next"})
    sgd.merge(g, on=["context", "layer", "k"], how="left").to_csv(H / "singular_structure.csv", index=False)

# ------------------------------------------------------------ heatmaps
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    show = [c for c in ("lung", "lung_external", "immune", "kidney_A") if (P / f"heatmap_cos_{c}.csv").exists()]
    if show:
        fig, axes = plt.subplots(1, len(show), figsize=(3.3 * len(show), 9), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, c in zip(axes, show):
            m = pd.read_csv(P / f"heatmap_cos_{c}.csv", index_col=0)
            im = ax.imshow(m.values, aspect="auto", cmap="Blues", vmin=0.5, vmax=max(0.75, float(np.nanmax(m.values))))
            ax.set_title(c + (" (discovery)" if c in L.DISCOVERY else ""), fontsize=9)
            ax.set_xticks(range(12)); ax.set_xticklabels(range(12), fontsize=6); ax.set_xlabel("layer", fontsize=8)
            ax.set_yticks(range(len(m.index))); ax.set_yticklabels(m.index, fontsize=4.5)
            bs = pd.read_csv(H / f"bandscan_{c}.csv"); bs = bs[bs.row_type == "band_cell"]
            sig = bs[bs.p_maxT_cos < 0.05]
            order = {b: i for i, b in enumerate(m.index)}
            ax.scatter(sig.layer, sig.band.map(order), s=2, c="k", marker="o", linewidths=0)
        fig.colorbar(im, ax=axes, shrink=0.4, label="detection-stratified AUROC (sign-invariant for width 1)")
        fig.suptitle("E3 band scan, Lambert TF vs other (dot: max-T p < 0.05 over 936 cells)", fontsize=9)
        fig.savefig(H / "bandscan_heatmaps.png", dpi=200, bbox_inches="tight")
except Exception as e:                      # figure is optional
    print("figure failed:", e)

# ------------------------------------------------------------ results.json
led = json.loads((H / "ledger_0744.json").read_text()) if (H / "ledger_0744.json").exists() else None
res = dict(experiment="E3_node", conv_hash=L.conv_hash(), frozen_hash=(L.R3 / "FROZEN.txt").read_text().split()[2],
           contexts_done=ctxs, H2_decision=decision,
           H2_primary={c: dict(strat_auroc_by_layer=[x["strat_auroc"] for x in prim[c]["layers"]],
                               p_maxT_by_layer=[x["p_maxT"] for x in prim[c]["layers"]],
                               ci_by_layer=[[x["strat_lo"], x["strat_hi"]] for x in prim[c]["layers"]]) for c in ctxs},
           baselines_file="baselines.csv", ledger_0744=None if led is None else dict(
               n_tf=led["n_tf"], n_target_only=led["n_target_only"],
               screen_procedure_mean=led["screen_procedure"]["mean_over_layers"],
               screen_procedure_perm_p_mean=led["screen_procedure"]["perm_null"]["p_mean"],
               loo_SV2_7_unstrat_mean=led["loo_cos_SV2_7"]["unstratified"]["mean_over_layers"],
               loo_SV2_7_strat_mean=led["loo_cos_SV2_7"]["stratified"]["mean_over_layers"],
               loo_512_unstrat_mean=led["loo_cos_full512"]["unstratified"]["mean_over_layers"],
               loo_512_strat_mean=led["loo_cos_full512"]["stratified"]["mean_over_layers"]),
           files=sorted(p.name for p in H.iterdir() if p.is_file()))
(H / "results.json").write_text(json.dumps(res, indent=1, default=float))
print(json.dumps(decision, indent=1))
