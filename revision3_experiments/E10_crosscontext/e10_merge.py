"""E10 final step: collect headline numbers from the CSVs into results.json (every number here is read from a
CSV / npz written by an E10 script; nothing is recomputed)."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e10_common as K
H = K.HERE


def f(x):
    try:
        x = float(x); return None if not np.isfinite(x) else round(x, 5)
    except Exception:
        return x


def main():
    res = dict(experiment="E10_crosscontext", kind="sensitivity (CONVENTION3 §5: effect sizes with intervals; no pass/fail)",
               convention3_sha256=K.convention_hash(), answers=["R2.7", "R1.7", "R2.8"],
               python=sys.executable, numpy=np.__version__, contexts=K.CTXS)
    # per-context layer-mean statistics
    pc = {}
    for c in K.CTXS:
        z = np.load(H / "partial" / f"stats_{c}.npz"); p = z["point"]; b = z["boot"]
        pc[c] = dict(meta=json.loads(str(z["meta"])),
                     **{str(s): dict(layer_mean=f(p[i].mean()), boot_se=f(np.nanstd(b[i].mean(0), ddof=1)),
                                     boot_mean=f(np.nanmean(b[i].mean(0))), per_layer=[f(v) for v in p[i]])
                        for i, s in enumerate(z["stats"])})
    res["per_context_statistics"] = pc
    # POST HOC (exploratory, defined after seeing varcomp): extraction regime = HVG-axis size of the context
    small = ["immune", "immune_seed43", "immune_seed44", "immune_TSP14", "immune_TSP21", "immune_TSP25",
             "immune_Bcell", "immune_CD4T", "immune_CD8T", "kidney_A"]           # 4,941 / 4,803-gene axes
    large = ["kidney_B", "kidney_ckpt_whole_human", "lung", "lung_external"]     # 8,143-8,229-gene axes, whole-human ckpt
    res["posthoc_extraction_regime"] = dict(
        note="exploratory grouping chosen after the variance decomposition; not a test",
        small_axis_contexts=small, large_axis_contexts=large,
        mean_small={st: f(np.mean([pc[c][st]["layer_mean"] for c in small])) for st in ["S1", "S2", "S2b_full", "S2b_band", "S3"]},
        mean_large={st: f(np.mean([pc[c][st]["layer_mean"] for c in large])) for st in ["S1", "S2", "S2b_full", "S2b_band", "S3"]},
        range_small={st: [f(min(pc[c][st]["layer_mean"] for c in small)), f(max(pc[c][st]["layer_mean"] for c in small))] for st in ["S1", "S2", "S2b_full", "S2b_band", "S3"]},
        range_large={st: [f(min(pc[c][st]["layer_mean"] for c in large)), f(max(pc[c][st]["layer_mean"] for c in large))] for st in ["S1", "S2", "S2b_full", "S2b_band", "S3"]},
        kidneyA_minus_kidneyB_fullV={st: f(pc["kidney_A"][st]["layer_mean"] - pc["kidney_B"][st]["layer_mean"]) for st in ["S1", "S2", "S2b_full", "S2b_band", "S3"]})
    v = pd.read_csv(H / "varcomp.csv")
    m = v[v.layer.astype(str) == "mean"]
    res["varcomp_layer_mean"] = {
        st: {r.group: {k: f(getattr(r, k)) for k in ("mean", "range", "tau", "tau2", "sampling_var", "I2", "p_Q",
                                                     "tau2_cond_lo", "tau2_cond_hi", "tau2_qprof_lo", "tau2_qprof_hi")}
             for r in m[(m.stat == st) & (m.group != "D_tissue_minus_within_pooled")].itertuples()}
        for st in m.stat.unique()}
    d = m[m.group == "D_tissue_minus_within_pooled"]
    res["tissue_dependence_decision"] = {
        r.stat: dict(tau2_tissue=f(r.tau2_tissue), tau2_within_pooled=f(r.tau2_within_pooled), tau2_within_max=f(r.tau2_within_max),
                     D=f(r.D), D_param_95ci=[f(r.D_param_lo), f(r.D_param_hi)], D_cond_95ci=[f(r.D_cond_lo), f(r.D_cond_hi)],
                     p_param=f(r.p_param_H0_tissue_eq_floor), claimed=bool(r.tissue_dependence_claimed))
        for r in d.itertuples()}
    dl = v[(v.group == "D_tissue_minus_within_pooled") & (v.layer.astype(str) != "mean")]
    res["tissue_dependence_per_layer_D_lower_gt0"] = {st: int((g.D_param_lo > 0).sum()) for st, g in dl.groupby("stat")}
    # kidney A vs B
    k = pd.read_csv(H / "kidneyA_vs_B.csv")
    km = k[k.layer.astype(str) == "mean"]
    res["kidneyA_vs_B_layer_mean"] = {r.statistic: dict(A=f(r.A), B=f(r.B), delta=f(r.delta), ci=[f(r.delta_lo), f(r.delta_hi)])
                                      for r in km.itertuples()}
    kp = k[k.layer.astype(str) == "profile"]
    res["kidneyA_vs_B_profile_mean_absdelta"] = {r.statistic: f(r.delta) for r in kp.itertuples()}
    br = k[k.statistic.str.startswith("legacy_bridge")]
    res["legacy_bridge_layer_mean"] = {s: f(g.A.mean()) for s, g in br.groupby("statistic")}
    res["legacy_bridge_kidneyB_L0_L3_submitted_band"] = [f(x) for x in
        br[br.statistic == "legacy_bridge_B_curated167_submitted_zero_row_SVD_band"].assign(L=lambda x: x.layer.astype(int)).sort_values("L").A.values[:4]]
    res["legacy_bridge_hanley_mcneil_se"] = dict(auc_0p55_n59_108=f(K.S.hanley_mcneil_se(0.55, 59, 108)),
                                                 auc_0p61_n65_122=f(K.S.hanley_mcneil_se(0.61, 65, 122)))
    res["kidney_mechanism"] = json.loads((H / "partial" / "kidney_mech.json").read_text())
    bd = pd.read_csv(H / "partial" / "kidney_basis_diag.csv")
    res["kidney_basis_diag_mean"] = {c: f(bd[c].mean()) for c in bd.columns if bd[c].dtype != object and c != "layer"}
    # detection
    de = pd.read_csv(H / "detection.csv", low_memory=False)
    a = de[de.part == "a"]
    res["detection_lognobs_auroc"] = {r.context: {} for r in a.itertuples()}
    for r in a.itertuples(): res["detection_lognobs_auroc"][r.context][r.test] = dict(value=f(r.value), ci=[f(r.lo), f(r.hi)])
    b = de[de.part == "b"]
    res["detection_interaction_layers_ci_excl0"] = {
        c: dict(positive=int(((g.test == "logit_score_x_lognobs") & (g.lo > 0)).sum()),
                negative=int(((g.test == "logit_score_x_lognobs") & (g.hi < 0)).sum()),
                mean_coef=f(g[g.test == "logit_score_x_lognobs"].value.mean()))
        for c, g in b.groupby("context")}
    cc = de[(de.part == "c") & (de.layer.astype(str) == "mean")]
    res["S2_restricted_layer_mean"] = {c: {r.test: dict(value=f(r.value), ci=[f(r.lo), f(r.hi)], n_tf=int(r.n_pos)) for r in g.itertuples()}
                                       for c, g in cc.groupby("context")}
    # lineage
    li = pd.read_csv(H / "lineage_tfs.csv", low_memory=False)
    sel = li[li.part == "selection"]
    res["lineage_sets"] = {t: sel[sel[f"in_set_{t}"] == True].sort_values(f"rank_{t}").symbol.tolist() for t in ("immune", "kidney", "lung")}
    ev = li[(li.part == "eval") & (li.layer.astype(str) == "mean")]
    res["lineage_eval_layer_mean"] = [dict(context=r.context, set=r.tissue_set, home=bool(r.home), n_observed=f(r.n_observed),
                                           median_nobs_set=f(r.median_nobs_set), auroc=f(r.auroc), ci_clipped=[f(r.lo_clip), f(r.hi_clip)],
                                           ci_reliable=bool(r.ci_reliable), random_subset_mean=f(r.random_subset_mean),
                                           random_subset_q975=f(r.random_subset_q975)) for r in ev.itertuples()]
    res["lineage_summary"] = {r.test: dict(rho=f(r.auroc), p=f(r.p), n=f(r.n_observed)) for r in li[li.part == "summary"].itertuples()}
    # transfer
    tr = pd.read_csv(H / "transfer.csv")
    tm = tr[tr.layer.astype(str) == "mean"]
    res["transfer_layer_mean"] = [{k: (f(getattr(r, k)) if k not in ("test", "subset", "centred") else getattr(r, k))
                                   for k in ("centred", "test", "subset", "auroc", "auroc_detstrat", "incontext_auroc",
                                             "incontext_auroc_detstrat", "transfer_minus_incontext") if hasattr(r, k)}
                                  for r in tm.itertuples()]
    tl = tr[(tr.layer.astype(str) != "mean") & (tr.centred == True) & (tr.subset == "novel")]
    res["transfer_novel_layers_transfer_below_incontext_ci"] = {t: dict(n_layers_ci_below0=int((g.diff_hi < 0).sum()), n_layers=int(len(g)),
                                                                       auroc_range=[f(g.auroc.min()), f(g.auroc.max())],
                                                                       lo_min=f(g.lo.min()))
                                                                for t, g in tl.groupby("test")}
    # stability
    sb = pd.read_csv(H / "stability.csv", low_memory=False)
    sm = sb[(sb.part == "knn_rsa") & (sb.layer.astype(str) == "mean")]
    res["stability_layer_mean"] = [dict(group=r.group, a=r.ctx_a, b=r.ctx_b, n_M=f(r.n_M), jaccard_TF=f(r.jaccard_TF),
                                        jaccard_nonTF=f(r.jaccard_nonTF), jaccard_TF_minus_nonTF_detadj=f(r.jaccard_TF_minus_nonTF_detadj),
                                        rsa_TF=f(r.rsa_TF), rsa_nonTF_sizematched=f(r.rsa_nonTF_sizematched), rsa_random1500=f(r.rsa_random1500))
                                   for r in sm.itertuples()]
    pm = sb[sb.part == "powermatch"]
    res["powermatch"] = [dict(ctx=r.ctx_a, stat=r.statistic, immune=f(r.immune_value), full=f(r.full_value), matched_mean=f(r.matched_mean),
                              matched_95=[f(r.matched_q025), f(r.matched_q975)], frac_matched_gt_immune=f(r.frac_matched_gt_immune))
                         for r in pm.itertuples()]
    res["outputs"] = ["varcomp.csv", "kidneyA_vs_B.csv", "detection.csv", "lineage_tfs.csv", "transfer.csv", "stability.csv",
                      "partial/stats_<ctx>.npz", "partial/kidney_basis_diag.csv", "partial/kidney_mech.json",
                      "partial/lineage_selection_full.csv"]
    (H / "results.json").write_text(json.dumps(res, indent=1, default=str))
    print("results.json written", K.convention_hash())


if __name__ == "__main__":
    main()
