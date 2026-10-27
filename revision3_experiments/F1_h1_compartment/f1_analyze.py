"""F1 analysis: reads _ckpt/<ctx>_*.npz from f1_run.py and E2's h1_primary_<ctx>.csv; computes max-T summaries
(E2 convention: two-sided, stats3.max_t_adjust via e2_common.null_summary) and vertex-subsampling intervals
(stats3.vertex_subsample_ci, B = 500; pre-registered sqrt(m/n) and delete-d sqrt(m/(n-m)) scalings).
Writes classes.csv, strat_null.csv, residual_assort.csv, within_between.csv, nobs_threshold.csv,
nobs_threshold_cells.csv, intervals.csv, results.json.
Usage: python f1_analyze.py --ctx lung,lung_external,kidney_B
"""
import argparse, json, time
import numpy as np
import pandas as pd
import f1_common as F

E, S = F.E, F.S
B_CI = 500
H = F.HERE


def cells():
    return [E.cell_name(c) for c in range(E.N_LAYERS * E.K_SV)]


def fam(obs, null):
    sm = E.null_summary(obs, null, tail="two-sided")
    return sm


def top_by_absz(z):
    return int(np.argsort(-np.abs(z), kind="mergesort")[0])


def ci_row(ctx, analysis, stat_name, col, x_fn, n, pi, pj, seed_parts, null_mean=np.nan, **extra):
    r = F.vs_ci_both(x_fn, n, pi, pj, B_CI, E.seed_for("F1ci", ctx, analysis, stat_name, col, *seed_parts))
    L, k = E.cell_name(col)
    row = dict(context=ctx, analysis=analysis, statistic=stat_name, layer=L, sv=k, null_mean=null_mean, **r, **extra)
    if np.isfinite(null_mean):
        for s in ("prereg", "deleted"):
            row[f"excess_lo_{s}"] = row[f"ci_lo_{s}"] - null_mean; row[f"excess_hi_{s}"] = row[f"ci_hi_{s}"] - null_mean
    return row


def analyse(ctx, res, rows):
    d = dict(np.load(F.CKPT / f"{ctx}_setup.npz", allow_pickle=True))
    Xc, R30, Rcl, cls = d["Xc"], d["R30"], d["Rcl"], d["cls"]
    n = int(d["n"]); ed = d["edges"].astype(np.int64); pi, pj = ed[:, 0], ed[:, 1]
    cl_ = cells(); Ls = [c[0] for c in cl_]; Ks = [c[1] for c in cl_]
    e2 = pd.read_csv(F.R3 / "E2_ppi" / f"h1_primary_{ctx}.csv")
    assert np.allclose(e2.r.values, E.assort_cols(Xc, pi, pj), atol=1e-12), "raw r differs from E2"
    top3 = [F.col_of(L, k) for L, k in F.E2_TOP[ctx]]
    out = dict(n_V=n, m=int(pi.size), e2_top3=[f"L{L} SV{k}" for L, k in F.E2_TOP[ctx]])

    # ------------------------------------------------------------- (a) classes
    ms = dict(np.load(F.CKPT / f"{ctx}_ms.npz"))
    obs_same = F.same_class_frac(cls, pi, pj); cfg = F.config_same_class_expectation(cls, pi, pj)
    sc = ms["same_class"]; z_same = (obs_same - sc.mean()) / sc.std(ddof=1)
    st = cls[np.concatenate([pi, pj])]; stub = np.bincount(st, minlength=8) / st.size
    within_counts = np.bincount(cls[pi][cls[pi] == cls[pj]], minlength=8)
    for ci, nm in enumerate(F.CLASS_ORDER):
        rows["classes"].append(dict(context=ctx, cls=nm, n_genes=int((cls == ci).sum()), frac_V=float((cls == ci).mean()),
                                    stub_share=float(stub[ci]), n_within_class_edges=int(within_counts[ci]),
                                    frac_edges_within_this_class=float(within_counts[ci] / pi.size),
                                    config_expected_frac_edges_within_this_class=float(stub[ci] ** 2)))
    rows["classes"].append(dict(context=ctx, cls="ALL_same_class", n_genes=n, frac_V=1.0, stub_share=1.0,
                                n_within_class_edges=int(within_counts.sum()), frac_edges_within_this_class=obs_same,
                                config_expected_frac_edges_within_this_class=cfg, ms_null_mean=float(sc.mean()),
                                ms_null_sd=float(sc.std(ddof=1)), ms_z=float(z_same), ms_N=int(sc.size),
                                ms_p_raw_greater=float(S.emp_p(obs_same, sc, "greater")[0])))
    out["classes"] = dict(sizes={nm: int((cls == ci).sum()) for ci, nm in enumerate(F.CLASS_ORDER)},
                          same_class_frac_obs=obs_same, same_class_frac_ms_mean=float(sc.mean()),
                          same_class_frac_ms_sd=float(sc.std(ddof=1)), same_class_frac_config=cfg, same_class_z=float(z_same),
                          n_any_go=int(d["n_any_annot"]), n_cc_annotated=int(d["n_cc_annot"]))

    # ------------------------------------------------------------- (b) stratified permutation + class-pair MS
    pm = dict(np.load(F.CKPT / f"{ctx}_perm.npz")); cp = dict(np.load(F.CKPT / f"{ctx}_cpms.npz"))
    r = ms["obs_raw"]; assert np.allclose(pm["obs"], r) and np.allclose(cp["obs"], r)
    s_pm, s_cp, s_ms = fam(r, pm["null"]), fam(r, cp["null"]), fam(r, ms["null_raw"])
    for c in range(r.size):
        rows["strat"].append(dict(
            context=ctx, layer=Ls[c], sv=Ks[c], r=r[c],
            strat_null_mean=s_pm["null_mean"][c], strat_null_sd=s_pm["null_sd"][c], strat_z=s_pm["z"][c],
            strat_p_raw=s_pm["p_raw"][c], strat_p_adj_maxT=s_pm["p_adj"][c], strat_excess_r=r[c] - s_pm["null_mean"][c],
            cpms_null_mean=s_cp["null_mean"][c], cpms_null_sd=s_cp["null_sd"][c], cpms_z=s_cp["z"][c],
            cpms_p_raw=s_cp["p_raw"][c], cpms_p_adj_maxT=s_cp["p_adj"][c], cpms_excess_r=r[c] - s_cp["null_mean"][c],
            f1ms_null_mean=s_ms["null_mean"][c], f1ms_null_sd=s_ms["null_sd"][c], f1ms_z=s_ms["z"][c],
            f1ms_p_adj_maxT=s_ms["p_adj"][c], f1ms_excess_r=r[c] - s_ms["null_mean"][c],
            e2ms_null_mean=e2.null_mean.iloc[c], e2ms_null_sd=e2.null_sd.iloc[c], e2ms_z=e2.z.iloc[c],
            e2ms_p_adj_maxT=e2.p_adj_maxT.iloc[c], e2_N=int(e2.N.iloc[c]),
            strat_excess_over_ms_excess=(r[c] - s_pm["null_mean"][c]) / (r[c] - s_ms["null_mean"][c]),
            cpms_excess_over_ms_excess=(r[c] - s_cp["null_mean"][c]) / (r[c] - s_ms["null_mean"][c]),
            N=int(pm["null"].shape[0])))
    df = pd.DataFrame([x for x in rows["strat"] if x["context"] == ctx])
    def fsum(zcol, pcol, mcol):
        t = top_by_absz(df[zcol].values)
        return dict(family_min_p_adj=float(df[pcol].min()), n_cells_p_adj_lt_05=int((df[pcol] < 0.05).sum()),
                    n_cells_p_adj_lt_05_positive=int(((df[pcol] < 0.05) & (df[zcol] > 0)).sum()),
                    top_cell=f"L{Ls[t]} SV{Ks[t]}", top_r=float(r[t]), top_null_mean=float(df[mcol].iloc[t]),
                    top_z=float(df[zcol].iloc[t]), median_z=float(df[zcol].median()),
                    median_excess_ratio=float(np.median((r - df[mcol].values) / (r - s_ms["null_mean"]))),
                    e2_top3={f"L{Ls[c]} SV{Ks[c]}": dict(r=float(r[c]), null_mean=float(df[mcol].iloc[c]),
                                                          z=float(df[zcol].iloc[c]), p_adj=float(df[pcol].iloc[c]),
                                                          excess_over_ms_excess=float((r[c] - df[mcol].iloc[c]) / (r[c] - s_ms["null_mean"][c])))
                             for c in top3}, _top=t)
    out["strat_perm"] = fsum("strat_z", "strat_p_adj_maxT", "strat_null_mean")
    out["classpair_ms"] = fsum("cpms_z", "cpms_p_adj_maxT", "cpms_null_mean")
    out["classpair_ms"].update(acc_mean=float(cp["acc"].mean()), completed_frac=float(cp["completed"].mean()),
                               same_class_frac_constant=bool(np.allclose(cp["same_class"], obs_same)))
    out["f1_ms_replica"] = fsum("f1ms_z", "f1ms_p_adj_maxT", "f1ms_null_mean")
    out["f1_ms_replica"]["acc_mean"] = float(ms["acc"].mean())
    out["f1_ms_replica"]["max_abs_null_mean_diff_vs_E2"] = float(np.abs(s_ms["null_mean"] - e2.null_mean.values).max())
    S.assert_resolution(int(pm["null"].shape[0]), 0.05, 1)

    # ------------------------------------------------------------- (c) residualised coordinates
    for nm, key, R in (("resid30cc", "r30", R30), ("residclass", "rcl", Rcl)):
        o = ms[f"obs_{key}"]; sm = fam(o, ms[f"null_{key}"]); r2 = d["r2_30"] if key == "r30" else d["r2_cl"]
        for c in range(o.size):
            rows["resid"].append(dict(context=ctx, design=nm, layer=Ls[c], sv=Ks[c], r_raw=r[c], r_resid=o[c],
                                      ratio_resid_over_raw=o[c] / r[c],
                                      excess_ratio=(o[c] - sm["null_mean"][c]) / (r[c] - s_ms["null_mean"][c]),
                                      R2_coord_on_design=r2[c], null_mean=sm["null_mean"][c], null_sd=sm["null_sd"][c],
                                      z=sm["z"][c], p_raw=sm["p_raw"][c], p_adj_maxT=sm["p_adj"][c], raw_z_f1ms=s_ms["z"][c],
                                      N=int(ms[f"null_{key}"].shape[0])))
        t = top_by_absz(sm["z"])
        out[nm] = dict(family_min_p_adj=float(sm["p_adj"].min()), n_cells_p_adj_lt_05=int((sm["p_adj"] < 0.05).sum()),
                       n_cells_sig_positive=int(((sm["p_adj"] < 0.05) & (sm["z"] > 0)).sum()),
                       top_cell=f"L{Ls[t]} SV{Ks[t]}", top_r=float(o[t]), top_z=float(sm["z"][t]),
                       top_raw_r_same_cell=float(r[t]),
                       median_ratio_resid_over_raw=float(np.median(o / r)),
                       iqr_ratio=[float(np.quantile(o / r, 0.25)), float(np.quantile(o / r, 0.75))],
                       median_R2_coord=float(np.median(r2)), median_resid_r=float(np.median(o)),
                       top_raw_cell_ratio=float(o[s_ms["z"].argmax()] / r[s_ms["z"].argmax()]),
                       e2_top3={f"L{Ls[c]} SV{Ks[c]}": dict(r_raw=float(r[c]), r_resid=float(o[c]), ratio=float(o[c] / r[c]),
                                                             R2=float(r2[c]), z=float(sm["z"][c]), p_adj=float(sm["p_adj"][c]))
                                for c in top3}, _top=t)
        if key == "r30":
            out[nm]["cc_terms"] = [dict(id=str(i), name=str(nn), n=int(k)) for i, nn, k in
                                   zip(d["cc_ids"], d["cc_names"], d["cc_counts"])]
    S.assert_resolution(int(ms["null_raw"].shape[0]), 0.05, 1)

    # ------------------------------------------------------------- (d) within vs between
    wb = dict(np.load(F.CKPT / f"{ctx}_wb.npz")); same = cls[pi] == cls[pj]
    out["within_between"] = {}
    for part in ("within", "between"):
        for ver in ("raw", "cc"):
            o = wb[f"{part}_obs_{ver}"]; sm = fam(o, wb[f"{part}_null_{ver}"])
            for c in range(o.size):
                rows["wb"].append(dict(context=ctx, edges=part, coords=("raw" if ver == "raw" else "class_centred"),
                                       m=int(wb[f"{part}_m"]), layer=Ls[c], sv=Ks[c], r=o[c], r_all_edges=r[c],
                                       null_mean=sm["null_mean"][c], null_sd=sm["null_sd"][c], z=sm["z"][c],
                                       p_raw=sm["p_raw"][c], p_adj_maxT=sm["p_adj"][c], N=int(wb[f"{part}_null_{ver}"].shape[0])))
            t = top_by_absz(sm["z"])
            out["within_between"][f"{part}_{ver}"] = dict(
                m=int(wb[f"{part}_m"]), acc_mean=float(wb[f"{part}_acc"].mean()),
                family_min_p_adj=float(sm["p_adj"].min()), n_cells_p_adj_lt_05=int((sm["p_adj"] < 0.05).sum()),
                n_cells_sig_positive=int(((sm["p_adj"] < 0.05) & (sm["z"] > 0)).sum()),
                top_cell=f"L{Ls[t]} SV{Ks[t]}", top_r=float(o[t]), top_z=float(sm["z"][t]), median_r=float(np.median(o)),
                median_ratio_to_all_edges=float(np.median(o / r)),
                e2_top3={f"L{Ls[c]} SV{Ks[c]}": dict(r=float(o[c]), z=float(sm["z"][c]), p_adj=float(sm["p_adj"][c]))
                         for c in top3}, _top=t)
    S.assert_resolution(int(wb["within_null_raw"].shape[0]), 0.05, 1)

    # ------------------------------------------------------------- (e) n_obs thresholds
    out["nobs"] = {}
    for t_ in (3, 5, 10):
        z = dict(np.load(F.CKPT / f"{ctx}_nobs{t_}.npz")); sm = fam(z["obs"], z["null"]); o = z["obs"]
        keep = z["keep"]; Xs = z["Xs"]
        # match each E2 top cell to the restricted SVk at the same layer with max |corr| of coordinates on kept genes
        match = {}
        for c in top3:
            L, k = E.cell_name(c)
            cs = [abs(np.corrcoef(Xc[keep, c], Xs[:, L * E.K_SV + j])[0, 1]) for j in range(E.K_SV)]
            j = int(np.argmax(cs)); cc_ = L * E.K_SV + j
            match[f"L{L} SV{k}"] = dict(matched=f"L{L} SV{j+1}", abs_corr=float(cs[j]), r=float(o[cc_]), z=float(sm["z"][cc_]),
                                        p_adj=float(sm["p_adj"][cc_]), _col=cc_)
        for c in range(o.size):
            rows["nobs_cells"].append(dict(context=ctx, nobs_min=t_, n_V=int(z["n"]), m=int(z["m"]), layer=Ls[c], sv=Ks[c],
                                           r=o[c], null_mean=sm["null_mean"][c], null_sd=sm["null_sd"][c], z=sm["z"][c],
                                           p_raw=sm["p_raw"][c], p_adj_maxT=sm["p_adj"][c], N=int(z["null"].shape[0])))
        tt = top_by_absz(sm["z"])
        out["nobs"][str(t_)] = dict(n_V=int(z["n"]), m=int(z["m"]), acc_mean=float(z["acc"].mean()),
                                    family_min_p_adj=float(sm["p_adj"].min()),
                                    n_cells_p_adj_lt_05=int((sm["p_adj"] < 0.05).sum()),
                                    n_cells_sig_positive=int(((sm["p_adj"] < 0.05) & (sm["z"] > 0)).sum()),
                                    top_cell=f"L{Ls[tt]} SV{Ks[tt]}", top_r=float(o[tt]), top_z=float(sm["z"][tt]),
                                    top_null_mean=float(sm["null_mean"][tt]), median_r=float(np.median(o)),
                                    max_r=float(o.max()), e2_top3_matched=match, _top=tt)
        S.assert_resolution(int(z["null"].shape[0]), 0.05, 1)

    # ------------------------------------------------------------- intervals
    def assort_fn(x):
        return lambda nodes, a, b: S.assortativity(x, a, b)
    ci = rows["ci"]
    raw_cols = sorted(set(top3 + [out["strat_perm"]["_top"]]))
    for c in raw_cols:
        ci.append(ci_row(ctx, "raw_all_edges", "r", c, assort_fn(Xc[:, c]), n, pi, pj, (),
                         null_mean=float(s_pm["null_mean"][c]), null_label="strat_perm_null_mean",
                         ms_null_mean=float(s_ms["null_mean"][c]), cpms_null_mean=float(s_cp["null_mean"][c])))
    for nm, R, key in (("resid30cc", R30, "r30"), ("residclass", Rcl, "rcl")):
        nmean = fam(ms[f"obs_{key}"], ms[f"null_{key}"])["null_mean"]
        for c in sorted(set(top3 + [out[nm]["_top"]])):
            ci.append(ci_row(ctx, nm, "r_resid", c, assort_fn(R[:, c]), n, pi, pj, (), null_mean=float(nmean[c])))
        for c in top3:
            xr, xx = R[:, c], Xc[:, c]
            ci.append(ci_row(ctx, nm, "ratio_resid_over_raw", c,
                             lambda nodes, a, b, xr=xr, xx=xx: S.assortativity(xr, a, b) / S.assortativity(xx, a, b),
                             n, pi, pj, ()))
    Xcc = F.class_centre(Xc, cls)
    for part, msk in (("within", same), ("between", ~same)):
        a0, b0 = pi[msk], pj[msk]
        for ver, M in (("raw", Xc), ("cc", Xcc)):
            key = f"{part}_{ver}"; nmean = fam(wb[f"{part}_obs_{ver}"], wb[f"{part}_null_{ver}"])["null_mean"]
            for c in sorted(set([top3[0], out["within_between"][key]["_top"]])):
                ci.append(ci_row(ctx, f"{part}_class_edges_{'raw' if ver == 'raw' else 'class_centred'}", "r", c,
                                 assort_fn(M[:, c]), n, a0, b0, (), null_mean=float(nmean[c])))
    for t_ in (3, 5, 10):
        z = dict(np.load(F.CKPT / f"{ctx}_nobs{t_}.npz")); Xs = z["Xs"]; nn = int(z["n"])
        nmean = fam(z["obs"], z["null"])["null_mean"]
        cs_ = sorted(set([out["nobs"][str(t_)]["_top"]] + [v["_col"] for v in out["nobs"][str(t_)]["e2_top3_matched"].values()]))
        for c in cs_:
            ci.append(ci_row(ctx, f"nobs_ge_{t_}", "r", c, assort_fn(Xs[:, c]), nn, z["a"], z["b"], (t_,),
                             null_mean=float(nmean[c])))
    res[ctx] = out


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if not k.startswith("_")}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def main(ctxs):
    t0 = time.time()
    h = E.frozen_hash_check()
    rows = {k: [] for k in ("classes", "strat", "resid", "wb", "nobs_cells", "ci")}
    res = {}
    for ctx in ctxs:
        analyse(ctx, res, rows); print(f"[{ctx}] analysed {time.time()-t0:.0f}s", flush=True)
    pd.DataFrame(rows["classes"]).to_csv(H / "classes.csv", index=False)
    pd.DataFrame(rows["strat"]).to_csv(H / "strat_null.csv", index=False)
    pd.DataFrame(rows["resid"]).to_csv(H / "residual_assort.csv", index=False)
    pd.DataFrame(rows["wb"]).to_csv(H / "within_between.csv", index=False)
    pd.DataFrame(rows["nobs_cells"]).to_csv(H / "nobs_threshold_cells.csv", index=False)
    pd.DataFrame(rows["ci"]).to_csv(H / "intervals.csv", index=False)
    summ = []
    for ctx in ctxs:
        for t_, v in res[ctx]["nobs"].items():
            cis = [x for x in rows["ci"] if x["context"] == ctx and x["analysis"] == f"nobs_ge_{t_}"
                   and f"L{x['layer']} SV{x['sv']}" == v["top_cell"]][0]
            summ.append(dict(context=ctx, nobs_min=int(t_), n_V=v["n_V"], m=v["m"], N=2000, family_min_p_adj=v["family_min_p_adj"],
                             n_cells_p_adj_lt_05=v["n_cells_p_adj_lt_05"], n_cells_sig_positive=v["n_cells_sig_positive"],
                             top_cell=v["top_cell"], top_r=v["top_r"], top_z=v["top_z"], top_null_mean=v["top_null_mean"],
                             top_r_ci_lo_prereg=cis["ci_lo_prereg"], top_r_ci_hi_prereg=cis["ci_hi_prereg"],
                             top_r_ci_lo_deleted=cis["ci_lo_deleted"], top_r_ci_hi_deleted=cis["ci_hi_deleted"],
                             median_r=v["median_r"], ms_acc_mean=v["acc_mean"],
                             **{f"e2top1_{k}": vv for k, vv in list(v["e2_top3_matched"].values())[0].items() if not k.startswith("_")}))
    pd.DataFrame(summ).to_csv(H / "nobs_threshold.csv", index=False)
    meta = dict(convention3_sha256=h, frozen_txt=(F.R3 / "FROZEN.txt").read_text().strip(), experiment="F1",
                role="sensitivity analysis (CONVENTION3 s5 / s2 V-restriction); no pass/fail",
                contexts=ctxs, N_main=2000, N_within_between=1000, B_ci=B_CI, seed_base=E.SEED,
                ms_chain="stats3.maslov_sneppen_chain burn 10m thin m",
                classpair_chain="f1_common.ms_classpair_chain: MS swaps accepted only if cls[v]==cls[y] or cls[u]==cls[x]; burn 10m thin m accepted swaps",
                tail="two-sided max-T (stats3.max_t_adjust), as E2",
                ci_note="vertex subsampling m=n/2; prereg sqrt(m/n) and delete-d sqrt(m/(n-m)) scalings from the same draws (ERRATA 3)",
                go_source=str(F.GENE2GO), go_obo=str(F.GO_OBO), classes_priority=F.CLASS_ORDER,
                results=strip(res), runtime_s=time.time() - t0)
    F.jdump(H / "results.json", meta)
    print("done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--ctx", default="lung,lung_external,kidney_B")
    main(ap.parse_args().ctx.split(","))
