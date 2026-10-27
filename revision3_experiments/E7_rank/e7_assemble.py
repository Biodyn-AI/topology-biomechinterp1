"""E7 assembly: delta.csv, any_k_test.csv, band/removal/static summaries, results.json, SUMMARY.txt."""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import e7lib as L

S = L.S
HERE = L.HERE
TASK_X = dict(PPI="STRING-physical partner recovery (pair-cosine AUROC vs degree-preserving rewired non-pairs)",
              NODE="Lambert TF identity (H2 detection-stratified LOO AUROC)",
              EDGE="TF->target edge recovery (H3 within-TF stratified AUROC)")
TASK_SHORT = dict(PPI="PPI-partner recovery", NODE="TF-identity recovery", EDGE="regulatory-edge recovery")


def fmt(x, nd=3):
    return f"{x:+.{nd}f}" if np.isfinite(x) else "nan"


def deleted_ci(est, lo, hi, n, m):
    """CONVENTION3_ERRATA.md item 3 (added by the E7 verifier): the pre-registered vertex-subsampling interval
    rescales deviations by sqrt(m/n); the manuscript uses the delete-d scaling sqrt(m/(n-m)). Both intervals are
    basic intervals est - quantile(scale * dev), so the delete-d interval is the pre-registered one with each
    half-width multiplied by sqrt(n/(n-m)) (exact; = sqrt(2) for even n)."""
    c = np.sqrt(n / (n - m))
    return float(est - (est - lo) * c), float(est + (hi - est) * c)


def main():
    kj = json.loads((HERE / "kstar.json").read_text()); kstar = {t: kj["tasks"][t]["k_star"] for t in L.TASKS}
    assert S.assert_resolution(L.N_PERM_ANYK, 0.05, len(L.TASKS) * len(L.EVAL_CTX))
    drows, arows, ver = [], [], {}
    for ctx in L.EVAL_CTX:
        meta = json.loads((L.WORK / f"meta_{ctx}.json").read_text()); n = meta["n_V"]; m = int(0.5 * n)
        for t in L.TASKS:
            z = np.load(L.WORK / f"{ctx}_delta_{t}.npz"); pts, dr = z["pts"], z["draws"]
            ver[f"{ctx}|{t}|delta"] = json.loads(str(z["verify"]))
            dd = dr[0] - dr[1]
            for lay in list(range(L.N_LAYERS)) + ["mean"]:
                if lay == "mean":
                    a, b = pts[0].mean(), pts[1].mean(); d = dd.mean(0)
                else:
                    a, b = pts[0, lay], pts[1, lay]; d = dd[lay]
                est = a - b
                lo_dd = hi_dd = np.nan
                if kstar[t] == 512:
                    lo = hi = 0.0
                elif t == "PPI":
                    lo, hi = L.vertex_ci_from_draws(est, d, n, m)
                    lo_dd, hi_dd = deleted_ci(est, lo, hi, n, m)
                else:
                    lo, hi = float(np.nanquantile(d, 0.025)), float(np.nanquantile(d, 0.975))
                lo_ms, hi_ms = (lo_dd, hi_dd) if t == "PPI" and kstar[t] != 512 else (lo, hi)
                drows.append(dict(ctx=ctx, task=t, layer=lay, k_star=kstar[t], metric_kstar=a, metric_full=b, delta=est,
                                  ci_lo=lo, ci_hi=hi, B=int(np.isfinite(d).sum()), ci_method=str(z["ci_method"]),
                                  ci_excludes_0_above=bool(lo > 0), ci_excludes_0_below=bool(hi < 0),
                                  ci_lo_deleted=lo_dd, ci_hi_deleted=hi_dd, ci_lo_ms=lo_ms, ci_hi_ms=hi_ms,
                                  ci_ms_excludes_0_above=bool(lo_ms > 0), ci_ms_excludes_0_below=bool(hi_ms < 0)))
            # any-k (descriptive)
            za = np.load(L.WORK / f"{ctx}_anyk_{t}.npz"); obs, null = za["obs"], za["null"]
            ver[f"{ctx}|{t}|anyk"] = float(za["verify"])
            iF = L.K_T.index(512); red = [i for i in range(len(L.K_T)) if i != iF]
            gap_obs = obs.mean(1)[red] - obs.mean(1)[iF]
            gap_null = null.mean(2)[:, red] - null.mean(2)[:, [iF]]
            T_obs = float(gap_obs.max()); T_null = gap_null.max(1)
            p, bcount, N = S.emp_p(T_obs, T_null)
            g2_obs = obs[red] - obs[[iF]]; g2_null = null[:, red, :] - null[:, [iF], :]
            T2_obs = float(g2_obs.max()); T2_null = g2_null.reshape(N, -1).max(1); p2, b2, _ = S.emp_p(T2_obs, T2_null)
            ia, la = np.unravel_index(np.argmax(g2_obs), g2_obs.shape)
            ns = L.summarise_null(T_null); ns2 = L.summarise_null(T2_null)
            arows.append(dict(ctx=ctx, task=t, statistic="max_k<512 [layer-mean metric(k) - layer-mean metric(512)]",
                              observed_max_gap=T_obs, argmax_k=L.K_T[red[int(np.argmax(gap_obs))]], null_p=p, null_b=bcount,
                              N=N, null_mean=ns["mean"], null_sd=ns["sd"], null_q95=ns["q95"], null_max=ns["max"],
                              observed_max_gap_any_layer=T2_obs, argmax_k_any_layer=L.K_T[red[ia]], argmax_layer=int(la),
                              null_p_any_layer=p2, null_q95_any_layer=ns2["q95"],
                              null=("within-detection-stratum Lambert-label permutation" if t == "NODE" else
                                    "within-TF label permutation" if t == "EDGE" else
                                    "STRING-pair vs rewired-non-pair label permutation"),
                              note="DESCRIPTIVE: the permutation removes the signal from every representation, so p "
                                   "measures the gap against label noise, not against H0: metric(k) <= metric(512)"))
    dl = pd.DataFrame(drows); dl.to_csv(HERE / "delta.csv", index=False)
    ak = pd.DataFrame(arows); ak.to_csv(HERE / "any_k_test.csv", index=False)

    # ---------------------------------------------------------------- answer sentences (frozen form)
    sentences = {}
    for t in L.TASKS:
        mrow = dl[(dl.task == t) & (dl.layer == "mean")]
        n_up = int(mrow.ci_ms_excludes_0_above.sum()); n_dn = int(mrow.ci_ms_excludes_0_below.sum())
        parts = "; ".join(f"{r.ctx} {fmt(r.delta)} [{r.ci_lo_ms:+.3f}, {r.ci_hi_ms:+.3f}]" for r in mrow.itertuples())
        verb = "improves" if n_up > 0 else "does not improve"
        if kstar[t] == 512:
            s = (f"Rank reduction does not improve {TASK_SHORT[t]} by Delta = 0 [0, 0] in 0 of 4 contexts: the selection "
                 f"context chose k* = 512 (the full representation), so Delta is identically zero; every truncation "
                 f"scored below full in immune_seed43.")
        else:
            s = (f"Rank reduction (k* = {kstar[t]} of 512) {verb} {TASK_SHORT[t]} by layer-mean Delta "
                 f"[95% CI] = {parts} in {n_up} of 4 contexts" + (f" (CI below 0 in {n_dn} of 4)" if n_dn else "") + ".")
            if t == "PPI":
                n_up0 = int(mrow.ci_excludes_0_above.sum()); n_dn0 = int(mrow.ci_excludes_0_below.sum())
                s += (f" [PPI CIs: delete-d vertex-subsampling scaling per CONVENTION3_ERRATA item 3; with the "
                      f"pre-registered sqrt(m/n) scaling: " + "; ".join(
                          f"{r.ctx} [{r.ci_lo:+.3f}, {r.ci_hi:+.3f}]" for r in mrow.itertuples()) +
                      f", CI above 0 in {n_up0} of 4, below 0 in {n_dn0} of 4.]")
        sentences[t] = s

    # ---------------------------------------------------------------- family R / N summaries (descriptive)
    cur = pd.concat([pd.read_csv(HERE / f"curves_{c}.csv") for c in [L.SELECTION_CTX] + L.EVAL_CTX + L.EXTRA_CTX])
    cm = cur[cur.layer == "mean"]
    fam_R = cm[cm.family == "R"][["ctx", "task", "rep", "value", "full_value", "diff_vs_full"]]
    fam_N = cm[cm.family.str.startswith("N_")][["ctx", "task", "family", "rep", "param", "value", "partner_rep",
                                                 "partner_value", "diff_vs_partner", "full_value", "diff_vs_full"]]
    fam_N = fam_N.copy()
    fam_N["screen_mapping_for_ctx"] = [
        ("identity (screen used observed basis)" if (c == "immune" and f == "N_nominal") else
         "legacy_band_mapping.json kidney_A legacy_full_zero" if (c == "kidney_A" and f == "N_kidneyA_rule") else
         "no screen mapping (index rule only)") for c, f in zip(fam_N.ctx, fam_N.family)]
    fam_N.to_csv(HERE / "bands_vs_topk.csv", index=False)
    fam_R.to_csv(HERE / "remove_top_j.csv", index=False)
    # curve shape: best top-k per eval ctx (descriptive, post hoc)
    best = []
    for (c, t), g in cm[cm.family == "T"].groupby(["ctx", "task"]):
        g = g.assign(k=g.param.astype(int)).sort_values("k")
        i = g.value.idxmax(); r = g.loc[i]
        best.append(dict(ctx=c, task=t, best_k_posthoc=int(r.k), best_value=float(r.value),
                         full_value=float(r.full_value), best_minus_full=float(r.value - r.full_value),
                         value_k10=float(g[g.k == 10].value.iloc[0]), value_k64=float(g[g.k == 64].value.iloc[0]),
                         value_k128=float(g[g.k == 128].value.iloc[0])))
    best = pd.DataFrame(best); best.to_csv(HERE / "posthoc_best_k.csv", index=False)

    st = pd.read_csv(HERE / "static_curves.csv") if (HERE / "static_curves.csv").exists() else None
    st_sum = []
    if st is not None:
        for (mdl, t), g in st[st.layer.astype(str).isin(["static", "mean"])].groupby(["model", "task"], sort=False):
            i = g.value.idxmax(); r = g.loc[i]
            st_sum.append(dict(model=mdl, task=t, dim_full=int(r.dim_full), full_value=float(r.full_value),
                               best_k=int(r.k), best_value=float(r.value), best_minus_full=float(r.value - r.full_value),
                               value_k10=float(g[g.k == 10].value.iloc[0]) if (g.k == 10).any() else np.nan,
                               value_k64=float(g[g.k == 64].value.iloc[0]) if (g.k == 64).any() else np.nan))
        st_sum = pd.DataFrame(st_sum); st_sum.to_csv(HERE / "static_summary.csv", index=False)

    metas = {c: json.loads((L.WORK / f"meta_{c}.json").read_text()) for c in [L.SELECTION_CTX] + L.EVAL_CTX + L.EXTRA_CTX}
    _, fs_sha = L.load_faststrat()
    res = dict(experiment="E7_rank", hypothesis="H4 (primary, pre-registered, CONVENTION3 §4)",
               convention3=L.convention_hash(), assembled=time.strftime("%Y-%m-%d %H:%M:%S"),
               kstar=kj, delta_layer_mean=dl[dl.layer == "mean"].to_dict("records"),
               any_k_test=ak.to_dict("records"), answers=sentences, verification=ver, context_meta=metas,
               settings=dict(K_T=L.K_T, J_R=L.J_R, bands_nominal=L.BANDS_NOMINAL, bands_kidneyA_rule=L.BANDS_KIDNEYA_RULE,
                             B_boot=L.B_BOOT, N_perm_anyk=L.N_PERM_ANYK, seed=L.SEED,
                             seeds=dict(ms_rewiring=L.SEED, node_boot=L.SEED + 101, edge_boot=L.SEED + 102,
                                        ppi_vertex=L.SEED + 103, node_perm=L.SEED + 201, edge_perm=L.SEED + 202,
                                        ppi_perm=L.SEED + 203),
                             e4a_lib_sha256=fs_sha, legacy_band_mapping=L.check_mapping(),
                             improves_rule="a context counts as 'improves' when the layer-mean Delta 95% CI lies above 0"),
               files=sorted(p.name for p in HERE.iterdir() if p.suffix in (".csv", ".json", ".txt") and not p.name.startswith("._")))
    if (HERE / "verify" / "audit.json").exists():
        res["adversarial_verification"] = json.loads((HERE / "verify" / "audit.json").read_text())
    (HERE / "results.json").write_text(json.dumps(res, indent=1, default=float))

    # ---------------------------------------------------------------- SUMMARY.txt
    Ls = []
    Ls.append("E7 - H4: does truncating to the top-k spectral directions improve recovery over the full representation?")
    Ls.append(f"Pre-registered primary (CONVENTION3 §4), sha256 {L.convention_hash()['actual']} (matches FROZEN.txt: "
              f"{L.convention_hash()['match']}).")
    Ls.append("")
    Ls.append("SELECTION (immune_seed43, layer-mean, family T, ties -> smallest k; kstar.json written before any "
              "evaluation context was computed):")
    for t in L.TASKS:
        v = kj["tasks"][t]
        Ls.append(f"  {t:4s} k* = {v['k_star']:3d}   layer-mean metric(k*) = {v['layer_mean_at_kstar']:.4f}, "
                  f"metric(512) = {v['layer_mean_full']:.4f}, selection gap {v['selection_gap']:+.4f}")
    Ls.append("")
    Ls.append("MAIN TABLE: layer-mean Delta = metric(k*) - metric(512) [95% CI], delta.csv (layer == 'mean')")
    Ls.append(f"  {'task':5s} {'k*':>4s}  " + "  ".join(f"{c:>26s}" for c in L.EVAL_CTX))
    for t in L.TASKS:
        mrow = dl[(dl.task == t) & (dl.layer == "mean")].set_index("ctx")
        Ls.append(f"  {t:5s} {kstar[t]:4d}  " + "  ".join(
            f"{mrow.loc[c,'delta']:+.4f} [{mrow.loc[c,'ci_lo_ms']:+.4f},{mrow.loc[c,'ci_hi_ms']:+.4f}]" for c in L.EVAL_CTX))
        if t == "PPI":
            Ls.append(f"  {'(PPI, pre-registered sqrt(m/n) scaling)':>11s} " + "  ".join(
                f"{'':7s} [{mrow.loc[c,'ci_lo']:+.4f},{mrow.loc[c,'ci_hi']:+.4f}]" for c in L.EVAL_CTX))
    Ls.append("  PPI CIs in the table use the delete-d vertex-subsampling scaling sqrt(m/(n-m)) (CONVENTION3_ERRATA item 3,"
              " manuscript convention; delta.csv ci_lo_deleted/ci_hi_deleted = ci_lo_ms/ci_hi_ms); the pre-registered"
              " sqrt(m/n) intervals are kept in delta.csv ci_lo/ci_hi. NODE/EDGE: ci_lo_ms = ci_lo.")
    Ls.append("  full-dimension layer-mean metric(512): " + "; ".join(
        f"{t}: " + ", ".join(f"{c} {dl[(dl.task==t)&(dl.layer=='mean')&(dl.ctx==c)].metric_full.iloc[0]:.3f}"
                             for c in L.EVAL_CTX) for t in L.TASKS))
    Ls.append("  per-layer CIs: delta.csv; CI methods: NODE gene bootstrap, EDGE two-way TF x target cluster bootstrap, "
              "PPI vertex subsampling (m = n/2) of the paired difference; B = 2,000.")
    Ls.append("  per-layer: number of layers (of 12) with CI above 0 / below 0 (PPI: delete-d; pre-registered in brackets):")
    for t in L.TASKS:
        g = dl[(dl.task == t) & (dl.layer != "mean")]
        Ls.append("    " + t + ": " + "; ".join(
            f"{c} {int(g[g.ctx==c].ci_ms_excludes_0_above.sum())}/{int(g[g.ctx==c].ci_ms_excludes_0_below.sum())}"
            + (f" [{int(g[g.ctx==c].ci_excludes_0_above.sum())}/{int(g[g.ctx==c].ci_excludes_0_below.sum())}]"
               if t == "PPI" else "")
            for c in L.EVAL_CTX))
    Ls.append("")
    Ls.append("ANSWERS (frozen form):")
    for t in L.TASKS: Ls.append(f"  {t}: {sentences[t]}")
    Ls.append("")
    Ls.append("DESCRIPTIVE joint 'any k beats full' (any_k_test.csv; N = 500 joint permutations across k and layers; "
              "statistic max over k<512 of the layer-mean gap):")
    for r in ak.itertuples():
        Ls.append(f"  {r.ctx:14s} {r.task:4s} observed max gap {r.observed_max_gap:+.4f} (k={r.argmax_k}), "
                  f"null p = {r.null_p:.3f} (null q95 {r.null_q95:+.4f});  any layer: {r.observed_max_gap_any_layer:+.4f} "
                  f"(k={r.argmax_k_any_layer}, L{r.argmax_layer}), p = {r.null_p_any_layer:.3f}")
    Ls.append("  Caveat: the label permutation erases the signal in every representation, so this p measures the gap "
              "against label noise only; it is not a test of H0 metric(k) <= metric(512).")
    Ls.append("")
    Ls.append("POST-HOC curve shape (posthoc_best_k.csv; best top-k per context, layer-mean; descriptive):")
    for r in best.itertuples():
        Ls.append(f"  {r.ctx:14s} {r.task:4s} best k={r.best_k_posthoc:3d} {r.best_value:.4f} vs full {r.full_value:.4f} "
                  f"({r.best_minus_full:+.4f}); k=10 {r.value_k10:.4f}, k=64 {r.value_k64:.4f}, k=128 {r.value_k128:.4f}")
    Ls.append("")
    Ls.append("FAMILY R remove-top-j (remove_top_j.csv; layer-mean diff vs full):")
    for (c, t), g in fam_R.groupby(["ctx", "task"], sort=False):
        Ls.append(f"  {c:14s} {t:4s} " + "  ".join(f"{r.rep}:{r.diff_vs_full:+.4f}" for r in g.itertuples()))
    Ls.append("")
    Ls.append("FAMILY N screen bands vs same-width top-k (bands_vs_topk.csv; layer-mean value, diff vs partner):")
    for (c, t), g in fam_N.groupby(["ctx", "task"], sort=False):
        Ls.append(f"  {c:14s} {t:4s} " + "  ".join(
            f"{'nom' if f=='N_nominal' else 'kA-rule'} {p}:{v:.3f}({d:+.3f} vs {pr})"
            for f, p, v, d, pr in zip(g.family, g.param, g.value, g.diff_vs_partner, g.partner_rep)))
    Ls.append("  kidney_A rule: SV1-3 and SV1-6 ARE top-3 / top-6 (diff 0 by construction); immune: mapping = identity.")
    sp_ = HERE / "sens_R_N_ci.csv"
    if sp_.exists():
        sv = pd.read_csv(sp_)
        # errata item 3: delete-d scaling for the PPI vertex-subsampling intervals (pre-registered columns kept)
        nV = {c: json.loads((L.WORK / f"meta_{c}.json").read_text())["n_V"] for c in sv.ctx.unique()}
        lo_dd, hi_dd = [], []
        for r in sv.itertuples():
            if r.task == "PPI":
                a_, b_ = deleted_ci(r.diff, r.ci_lo, r.ci_hi, nV[r.ctx], int(0.5 * nV[r.ctx]))
            else:
                a_, b_ = np.nan, np.nan
            lo_dd.append(a_); hi_dd.append(b_)
        sv["ci_lo_deleted"] = lo_dd; sv["ci_hi_deleted"] = hi_dd
        sv["ci_lo_ms"] = np.where(sv.task == "PPI", sv.ci_lo_deleted, sv.ci_lo)
        sv["ci_hi_ms"] = np.where(sv.task == "PPI", sv.ci_hi_deleted, sv.ci_hi)
        sv.to_csv(sp_, index=False)
        sm = sv[sv.layer.astype(str) == "mean"]
        res["sensitivity_R_N_layer_mean"] = sm.to_dict("records")
        (HERE / "results.json").write_text(json.dumps(res, indent=1, default=float))
        Ls.append("")
        Ls.append("SENSITIVITY with paired CIs (sens_R_N_ci.csv; same resamples as the primary; layer-mean diff [95% CI]; "
                  "PPI CIs delete-d scaled, columns ci_lo_ms/ci_hi_ms):")
        for (c, t), g in sm.groupby(["ctx", "task"], sort=False):
            Ls.append(f"  {c:14s} {t:4s} " + "  ".join(
                f"{r.rep}-{r.reference}:{r.diff:+.3f}[{r.ci_lo_ms:+.3f},{r.ci_hi_ms:+.3f}]" for r in g.itertuples()))
    if st is not None:
        Ls.append("")
        Ls.append("STATIC MODELS on the lung Tier B gene set (static_curves.csv / static_summary.csv; descriptive):")
        for r in st_sum.itertuples():
            Ls.append(f"  {r.model:22s} {r.task:4s} d={r.dim_full:5d} full {r.full_value:.4f}  best k={r.best_k:4d} "
                      f"{r.best_value:.4f} ({r.best_minus_full:+.4f}); k=10 {r.value_k10:.4f}")
    ip = HERE / "summary_interpretation.txt"          # hand-written interpretation block (numbers from the CSVs)
    if ip.exists():
        Ls.append(""); Ls += ip.read_text().rstrip("\n").split("\n")
    (HERE / "SUMMARY.txt").write_text("\n".join(Ls) + "\n")
    print("\n".join(Ls))


if __name__ == "__main__":
    main()
