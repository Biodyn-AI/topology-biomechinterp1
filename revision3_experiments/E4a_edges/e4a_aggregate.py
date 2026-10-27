"""E4a aggregation: per-context primary tables, pre-declared wording rule, power, sensitivity CSVs, results.json.

Decision rule (CONVENTION3 §4 H3), applied per context with the primary comparator (|co-expression| on the
extraction cells) and, as a secondary column, the 3,000-cell tier:
  detected                    max-T p (12 layers, N = 10,000, tail 'greater') < 0.05 at >= 1 layer
  AUROC >= 0.60 excluded      every layer's one-sided 95% upper bound < 0.60 (intersection-union)
  not above co-expression     every layer's one-sided 95% upper bound of delta (layer - |co-expr|) < 0.05
  inconclusive                none of the above
Labels are not mutually exclusive (e.g. 'detected' and 'AUROC >= 0.60 excluded' = a small but non-zero signal).
Replication criterion (CONVENTION3 §1): lung_external has the same sign (AUROC - 0.5) as lung and its 95% CI
excludes 0.5 -- evaluated per layer.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; R3 = HERE.parent
sys.path.insert(0, str(R3))
import stats3 as S          # noqa: E402
import contexts as X        # noqa: E402

W = HERE / "work"
CTXS = list(X.CONTEXTS)
FROZEN = (R3 / "FROZEN.txt").read_text().split()[2]
PIS = [0.05, 0.1, 0.2]


def jl(p): return json.loads(Path(p).read_text())


def power_long(design, statistic, pw, pre, base):
    rows = []
    for test in ("permutation", "two_way_bootstrap", "tf_only_bootstrap", "iid_delong"):
        if f"{pre}auroc_mde_{test}" not in pw and f"{pre}delta_mde_{test}" not in pw: continue
        rows.append(dict(design=design, statistic=statistic, test=test,
                         mde_auroc_sim=pw.get(f"{pre}auroc_mde_{test}", np.nan), size_auroc_sim=pw.get(f"{pre}auroc_size_{test}", np.nan),
                         mde_delta_sim=pw.get(f"{pre}delta_mde_{test}", np.nan), size_delta_sim=pw.get(f"{pre}delta_size_{test}", np.nan),
                         sampling_sd_auroc_at_null=pw.get(f"{pre}auroc_sampling_sd_at_null"), mean_se_two_way_at_null=pw.get(f"{pre}auroc_mean_se_two_way_at_null"),
                         mean_se_tf_only_at_null=pw.get(f"{pre}auroc_mean_se_tf_only_at_null"), sd_perm_at_null=pw.get(f"{pre}auroc_sd_perm_at_null"),
                         mde_auroc_oracle_sd=pw.get(f"{pre}auroc_mde_oracle_sd"), mde_delta_oracle_sd=pw.get(f"{pre}delta_mde_oracle_sd"), **base))
    return rows


def main():
    done = [c for c in CTXS if (W / c / "DONE").exists()]
    missing = [c for c in CTXS if c not in done]
    summ, power_rows, sc_rows, macro_rows, leg_rows, short_rows, strata_rows, chip_rows, noise_rows, ctx_rows, extra_rows, mapped_rows, example_rows = ([] for _ in range(13))
    prim_all = {}
    for c in done:
        meta = jl(W / c / "meta.json"); perm = jl(W / c / "perm.json")
        rows = [jl(W / c / f"rows_L{l:02d}.json") for l in range(12)]
        rcx = jl(W / c / "rows_coexp.json")
        pw = jl(W / c / "power.json") if (W / c / "power.json").exists() else {}
        prim = pd.DataFrame([r[0] for r in rows])
        prim.insert(0, "context", c)
        prim["n_pos"] = meta["n_pos"]; prim["n_neg"] = meta["n_neg"]; prim["n_tf"] = meta["n_tf"]; prim["n_targets"] = meta["n_targets"]
        prim["p_raw"] = perm["p_raw"]; prim["p_maxT"] = perm["p_maxT"]; prim["p_maxT_two_sided"] = perm["p_maxT_two_sided"]
        prim["z_null_standardised"] = perm["z"]
        prim["null_mean"] = [s["mean"] for s in perm["null_summary"]]; prim["null_sd"] = [s["sd"] for s in perm["null_summary"]]
        prim["null_q95"] = [s["q95"] for s in perm["null_summary"]]; prim["perm_N"] = perm["N"]; prim["p_floor"] = perm["floor"]
        assert np.allclose(prim["auroc"], perm["obs"])
        prim["mde_auroc_bootstrap_wald_layer"] = 0.5 + 2.8015857 * prim["boot_sd"]
        prim["mde_auroc_sim_permutation"] = pw.get("auroc_mde_permutation", np.nan)
        prim["mde_auroc_sim_two_way_bootstrap"] = pw.get("auroc_mde_two_way_bootstrap", np.nan)
        prim["mde_auroc_hanley_mcneil"] = S.mde_auroc(meta["n_pos"], meta["n_neg"])
        prim["mde_delta_sim_two_way_bootstrap"] = pw.get("delta_mde_two_way_bootstrap", np.nan)
        prim["auroc_cx_ext"] = rcx[0]["auroc"]; prim["auroc_cx_3k"] = rcx[1]["auroc"]
        prim.to_csv(HERE / f"h3_primary_{c}.csv", index=False); prim_all[c] = prim
        det = bool(prim.p_maxT.min() < 0.05); excl = S.equivalence_excluded(prim.ub95, 0.60)
        nab = S.equivalence_excluded(prim.delta_cxext_ub95, 0.05); nab3 = S.equivalence_excluded(prim.delta_cx3k_ub95, 0.05)
        labs = [x for x, f in (("detected", det), ("AUROC>=0.60 excluded", excl), ("not above co-expression", nab)) if f] or ["inconclusive"]
        labs3 = [x for x, f in (("detected", det), ("AUROC>=0.60 excluded", excl), ("not above co-expression (3k tier)", nab3)) if f] or ["inconclusive"]
        best = prim.loc[prim.auroc.idxmax()]
        summ.append(dict(context=c, group=X.CONTEXTS[c]["group"], n_pos=meta["n_pos"], n_neg=meta["n_neg"], n_tf=meta["n_tf"],
                         n_targets=meta["n_targets"], frac_neg_fully_matched=meta["frac_neg_fully_matched"],
                         auroc_min=prim.auroc.min(), auroc_max=prim.auroc.max(), best_layer=int(best.layer),
                         best_auroc=best.auroc, best_ci_lo=best.ci_lo, best_ci_hi=best.ci_hi, best_ub95=best.ub95,
                         L0_auroc=prim.auroc[0], L0_ci_lo=prim.ci_lo[0], L0_ci_hi=prim.ci_hi[0],
                         L11_auroc=prim.auroc[11], L11_ci_lo=prim.ci_lo[11], L11_ci_hi=prim.ci_hi[11],
                         max_ub95=prim.ub95.max(), min_p_maxT=prim.p_maxT.min(), n_layers_detected=int((prim.p_maxT < 0.05).sum()),
                         layers_detected=";".join(str(int(l)) for l in prim.layer[prim.p_maxT < 0.05]),
                         auroc_cx_ext=rcx[0]["auroc"], cx_ext_ci_lo=rcx[0]["ci_lo"], cx_ext_ci_hi=rcx[0]["ci_hi"],
                         auroc_cx_3k=rcx[1]["auroc"], cx_3k_ci_lo=rcx[1]["ci_lo"], cx_3k_ci_hi=rcx[1]["ci_hi"],
                         delta_cxext_min=prim.delta_cxext.min(), delta_cxext_max=prim.delta_cxext.max(), max_delta_cxext_ub95=prim.delta_cxext_ub95.max(),
                         min_delta_cxext_ci_lo=prim.delta_cxext_ci_lo.min(), n_layers_delta_cxext_ci_above0=int((prim.delta_cxext_ci_lo > 0).sum()),
                         max_delta_cx3k_ub95=prim.delta_cx3k_ub95.max(), n_layers_delta_cx3k_ci_above0=int((prim.delta_cx3k_ci_lo > 0).sum()),
                         detected=det, auroc_ge_060_excluded=excl, not_above_coexpression=nab, not_above_coexpression_3k=nab3,
                         verdict=" + ".join(labs), verdict_3k_comparator=" + ".join(labs3),
                         mde_auroc_sim_permutation=pw.get("auroc_mde_permutation", np.nan),
                         mde_auroc_sim_two_way_bootstrap=pw.get("auroc_mde_two_way_bootstrap", np.nan),
                         mde_auroc_sim_tf_only_bootstrap=pw.get("auroc_mde_tf_only_bootstrap", np.nan),
                         mde_auroc_hanley_mcneil=S.mde_auroc(meta["n_pos"], meta["n_neg"]),
                         mde_auroc_bootstrap_wald_median=0.5 + 2.8015857 * prim.boot_sd.median(),
                         mde_delta_sim_two_way_bootstrap=pw.get("delta_mde_two_way_bootstrap", np.nan),
                         mde_delta_sim_tf_only_bootstrap=pw.get("delta_mde_tf_only_bootstrap", np.nan),
                         mde_delta_bootstrap_wald_median=2.8015857 * prim.delta_cxext_sd.median(),
                         macro_auroc_min=prim.macro_auroc.min(), macro_auroc_max=prim.macro_auroc.max()))
        if pw:
            power_rows += power_long(f"{c} h3_pairs", "within-TF stratified AUROC", pw, "", dict(
                n_pos=pw["n_pos"], n_neg=pw["n_neg"], n_tf=pw["n_tf"], n_targets=pw["n_targets"], rho_target_icc=pw["rho_target_icc"],
                rho_scgpt_vs_coexp=pw["rho_scgpt_vs_coexp"], nsim_per_grid=pw["nsim_per_grid"], B_sim=pw["B_sim"],
                mde_auroc_hanley_mcneil=pw["mde_auroc_hanley_mcneil"], mde_delta_hanley_mcneil=pw["mde_delta_hanley_mcneil"],
                observed_auroc_min=prim.auroc.min(), observed_auroc_max=prim.auroc.max(),
                observed_delta_cxext_min=prim.delta_cxext.min(), observed_delta_cxext_max=prim.delta_cxext.max()))
        # secondary scorers
        for l in range(12):
            for r in rows[l][1:]:
                sc_rows.append(dict(context=c, **r, label="hypothesis-generating (nominal)" if r["scorer"].startswith("band") else "sensitivity"))
        for r in rcx[2:]:
            short_rows.append(dict(context=c, pairs="h3_pairs (target-matched)", metric="within-TF stratified AUROC", **r))
        for r in jl(W / c / "sens_macro.json"): macro_rows.append(dict(context=c, **r))
        lg = jl(W / c / "sens_legacy_pooled.json")
        for r in lg["rows"]:
            row = dict(context=c, n_pos=lg["n_pos"], n_neg=lg["n_neg"], n_tf=lg["n_tf"], **r)
            (short_rows if r["scorer"].startswith("shortcut") else leg_rows).append(row if not r["scorer"].startswith("shortcut") else dict(row, pairs="h3_legacy_pairs (all-V non-edges)"))
        for r in jl(W / c / "sens_strata.json"): strata_rows.append(dict(context=c, **r))
        if (W / c / "sens_chip.json").exists():
            ch = jl(W / c / "sens_chip.json")
            for r in ch["rows"]: chip_rows.append(dict(context=c, n_pos=ch["n_pos"], n_neg=ch["n_neg"], n_tf=ch["n_tf"], **r))
        for _, r in prim.iterrows():
            for pi in PIS:
                noise_rows.append(dict(context=c, layer=int(r.layer), pi=pi, auroc_obs=r.auroc, ub95_obs=r.ub95,
                                       auroc_corrected=S.label_noise_corrected(r.auroc, pi), ub95_corrected=S.label_noise_corrected(r.ub95, pi),
                                       ub95_corrected_below_060=bool(S.label_noise_corrected(r.ub95, pi) < 0.60)))
        ex = jl(W / c / "sens_example.json")
        if ex.get("pair_in_benchmark"):
            for r in ex["per_layer"]:
                example_rows.append(dict(context=c, pair="PAX5->CD79A", label=ex["label"], n_candidates=ex["n_candidates"], n_pos=ex["n_pos"],
                                         n_neg=ex["n_neg"], PAX5_det_frac=ex["PAX5_det_frac"], CD79A_det_frac=ex["CD79A_det_frac"], **r))
        if (W / c / "sens_extra.json").exists():
            e = jl(W / c / "sens_extra.json")
            for r in e["rows"]: extra_rows.append(dict(context=c, **r))
            if e.get("mapped_bands", {}).get("status") == "ok":
                for r in e["mapped_bands"]["rows"]: mapped_rows.append(dict(context=c, label="hypothesis-generating (nominal)", **r))
    # contextual: immune cell-type contexts restricted to TFs detected in >= 10% of extraction cells
    for r in strata_rows:
        if r["context"] in ("immune_Bcell", "immune_CD4T", "immune_CD8T") and r.get("stratum", "").startswith("TF detected in >= 10%"):
            ctx_rows.append(r)
    # legacy power
    if (W / "power_legacy.json").exists():
        pl = jl(W / "power_legacy.json")
        base = dict(n_pos=pl["n_pos"], n_neg=pl["n_neg"], n_tf=pl["n_tf"], n_targets=pl["n_targets"], rho_target_icc=pl["rho_target_icc_within_tf"],
                    rho_tf_icc=pl["rho_tf_icc"], rho_scgpt_vs_coexp=pl["rho_full_vs_coexp"], nsim_per_grid=pl["nsim_per_grid"], B_sim=pl["B_sim"],
                    mde_auroc_hanley_mcneil=pl["mde_auroc_hanley_mcneil"], mde_delta_hanley_mcneil=pl["mde_delta_hanley_mcneil"],
                    observed_auroc_min=pl["observed_auroc_full512"], observed_auroc_max=pl["observed_auroc_full512"],
                    observed_delta_cxext_min=pl["observed_auroc_full512"] - pl["observed_auroc_coexp"],
                    observed_delta_cxext_max=pl["observed_auroc_full512"] - pl["observed_auroc_coexp"])
        power_rows += power_long("kidney_A legacy (run_trees.py split, 49/147)", "pooled AUROC (as submitted)", pl, "pooled_", base)
        power_rows += power_long("kidney_A legacy (run_trees.py split, 49/147)", "within-TF stratified AUROC", pl, "stratified_", base)
    out = {"h3_summary.csv": summ, "power.csv": power_rows, "sens_scorers.csv": sc_rows, "sens_macro.csv": macro_rows,
           "sens_legacy_pooled.csv": leg_rows, "sens_shortcut.csv": short_rows, "sens_strata.csv": strata_rows,
           "sens_chip.csv": chip_rows, "sens_label_noise.csv": noise_rows, "sens_contextual.csv": ctx_rows, "sens_example_pax5_cd79a.csv": example_rows,
           "sens_matching.csv": extra_rows, "sens_mapped_bands.csv": mapped_rows}
    for fn, rws in out.items():
        if rws: pd.DataFrame(rws).to_csv(HERE / fn, index=False)
    # replication (per layer)
    rep = []
    if "lung" in prim_all and "lung_external" in prim_all:
        a, b = prim_all["lung"], prim_all["lung_external"]
        for l in range(12):
            rep.append(dict(layer=l, lung_auroc=a.auroc[l], lung_ext_auroc=b.auroc[l], lung_ext_ci_lo=b.ci_lo[l], lung_ext_ci_hi=b.ci_hi[l],
                            same_sign=bool(np.sign(a.auroc[l] - 0.5) == np.sign(b.auroc[l] - 0.5)),
                            lung_ext_ci_excludes_05=bool(b.ci_lo[l] > 0.5 or b.ci_hi[l] < 0.5)))
        pd.DataFrame(rep).to_csv(HERE / "h3_replication_lung_external.csv", index=False)
    verif = {c: jl(W / c / "meta.json")["fast_vs_stats3"] for c in done}
    perm_verif = {c: jl(W / c / "perm.json")["fast_vs_stats3_max_absdiff"] for c in done}
    res = dict(experiment="E4a_edges", hypothesis="H3 edge relation (pre-registered primary) + edge sensitivity analyses",
               convention3_sha256=FROZEN, frozen_note=(R3 / "FROZEN.txt").read_text().strip(),
               contexts_done=done, contexts_missing=missing,
               primary_spec=dict(pairs="cache/<ctx>/h3_pairs.npz", scorer="signed cosine, V-centred full 512-d, per layer (12)",
                                 statistic="stats3.StratifiedPairs.auroc (within-TF stratified)",
                                 detection="within-TF label permutation applied jointly to 12 layers, N=10,000, stats3.max_t_adjust tail='greater'",
                                 intervals="two-way TF x target cluster bootstrap B=2,000 (stats3 weight stream), percentile 95% CI + one-sided 95% UB",
                                 comparator="|Pearson r| on expr_extraction.npy (primary) and expr.npy (3,000 cells), paired delta on same resamples",
                                 decision="detected: min p_maxT<0.05; AUROC>=0.60 excluded: all UB<0.60; not above co-expression: all delta UB<0.05"),
               implementation_notes=["Bootstrap and permutation use vectorised evaluators (e4a_lib.py) verified equal to stats3 "
                                     "(point, two-way and TF-only bootstrap replicates from the same seed, permuted-label AUROC); max abs diff recorded.",
                                     "Permutation random stream differs from StratifiedPairs.permute_within_tf (vectorised generator of the same within-TF permutation distribution).",
                                     "rank-k and band scorers, strata on large contexts and ChIP outside lung/lung_external use the first 1,000 of the 2,000 resamples / B=1,000.",
                                     "Strata: positives restricted, all candidate non-edges of each TF kept; TF detection = fraction of extraction cells with expression > 0."],
               verification_fast_vs_stats3=verif, verification_perm_fast_vs_stats3=perm_verif,
               verdicts={r["context"]: r["verdict"] for r in summ},
               verdicts_3k_comparator={r["context"]: r["verdict_3k_comparator"] for r in summ},
               replication_lung_external=rep,
               headline={r["context"]: {k: r[k] for k in ("auroc_min", "auroc_max", "best_layer", "best_ci_lo", "best_ci_hi", "max_ub95",
                                                          "min_p_maxT", "auroc_cx_ext", "max_delta_cxext_ub95", "mde_auroc_sim_permutation", "mde_auroc_sim_two_way_bootstrap",
                                                          "mde_delta_sim_two_way_bootstrap", "mde_auroc_hanley_mcneil")} for r in summ if r["context"] in ("lung", "lung_external", "kidney_A")},
               legacy_power=jl(W / "power_legacy.json") if (W / "power_legacy.json").exists() else None,
               legacy_split_meta=jl(W / "legacy_meta.json") if (W / "legacy_meta.json").exists() else None,
               files=sorted(p.name for p in HERE.glob("*.csv")))
    (HERE / "results.json").write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    print("done:", done, "missing:", missing)


if __name__ == "__main__":
    main()
