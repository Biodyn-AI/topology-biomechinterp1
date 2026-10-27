"""E2 assembly: checkpoints -> CSVs + results.json. Every number in SUMMARY.md comes from a file written here
(or by e2_bridge.py / e2_prep.py)."""
import json
import numpy as np, pandas as pd
import e2_common as E

S = E.S
H = E.HERE


def fam_rows(ctx, z, tail="two-sided", extra=None):
    obs, null = z["obs"], z["null"]
    sm = E.null_summary(obs, null, tail=tail)
    rows = []
    for cidx in range(obs.size):
        L, k = E.cell_name(cidx)
        d = dict(context=ctx, layer=L, sv=k, r=obs[cidx], null_mean=sm["null_mean"][cidx], null_sd=sm["null_sd"][cidx],
                 z=sm["z"][cidx], p_raw=sm["p_raw"][cidx], p_adj_maxT=sm["p_adj"][cidx], null_q025=sm["q025"][cidx],
                 null_q975=sm["q975"][cidx], excess_r=obs[cidx] - sm["null_mean"][cidx],
                 n_V=int(z["n"]), m=int(z["m"]), N=int(null.shape[0]))
        if extra: d.update(extra)
        rows.append(d)
    df = pd.DataFrame(rows)
    df["rank_absz"] = df.z.abs().rank(ascending=False, method="first").astype(int)
    return df


def fam_summary(df, label):
    t = df.sort_values("rank_absz").iloc[0]
    return dict(context=df.context.iloc[0], null=label, n_V=int(t.n_V), m=int(t.m), N=int(t.N),
                family_min_p_adj=df.p_adj_maxT.min(), n_cells_p_adj_lt_05=int((df.p_adj_maxT < 0.05).sum()),
                n_cells_positive_r=int((df.r > 0).sum()), top_layer=int(t.layer), top_sv=int(t.sv), top_r=t.r,
                top_null_mean=t.null_mean, top_null_sd=t.null_sd, top_z=t.z, top_p_adj=t.p_adj_maxT,
                median_z=df.z.median(), max_r=df.r.max(), median_r=df.r.median(), p_floor=1 / (t.N + 1))


def main():
    hash_ = E.frozen_hash_check()
    res = dict(experiment="E2_ppi", hypothesis="H1", convention3_sha256=hash_, primary={}, sensitivity={})
    # ---------------------------------------------------------------- primary
    summ, diag = [], []
    prim = {}
    for ctx in E.PRIMARY_CTX:
        N = E.N_PRIMARY[ctx]
        z = np.load(E.CKPT / f"{ctx}_ms_{N}.npz")
        S.assert_resolution(N, 0.05, 1)            # family-wise max-T claim at alpha = 0.05
        df = fam_rows(ctx, z); prim[ctx] = df
        df.to_csv(H / f"h1_primary_{ctx}.csv", index=False)
        s = fam_summary(df, "maslov_sneppen")
        p2 = json.loads((E.CKPT / f"{ctx}_pass2.json").read_text())
        top = p2["top_cells"]
        assert E.cell_name(top[0]) == (s["top_layer"], s["top_sv"])
        vr = {v["cell"]: v for v in p2["vertex_ci"] if v["stat"] == "r"}
        s.update(top_r_ci_lo=vr[top[0]]["lo"], top_r_ci_hi=vr[top[0]]["hi"], ci_B=vr[top[0]]["B"],
                 top_excess_r_lo=vr[top[0]]["lo"] - s["top_null_mean"], top_excess_r_hi=vr[top[0]]["hi"] - s["top_null_mean"],
                 second=f"L{E.cell_name(top[1])[0]} SV{E.cell_name(top[1])[1]}",
                 second_r=float(df.r.iloc[top[1]]), second_r_ci=f"[{vr[top[1]]['lo']:.3f}, {vr[top[1]]['hi']:.3f}]",
                 third=f"L{E.cell_name(top[2])[0]} SV{E.cell_name(top[2])[1]}",
                 third_r=float(df.r.iloc[top[2]]), third_r_ci=f"[{vr[top[2]]['lo']:.3f}, {vr[top[2]]['hi']:.3f}]",
                 ms_acc_mean=float(z["acc"].mean()), role={"lung": "primary", "lung_external": "replication"}.get(ctx, "secondary"))
        summ.append(s)
        # chain diagnostics
        nul = z["null"]; lag = np.array([E.lag1(nul[:, j]) for j in range(nul.shape[1])])
        cp = pd.DataFrame(p2["copole"]); cpm = cp[(cp.null == "ms") & (cp.pole_frac == 0.27)]
        diag.append(dict(context=ctx, N=N, m=int(z["m"]), burn_swaps=10 * int(z["m"]), thin_swaps=int(z["m"]),
                         acc_mean=float(z["acc"].mean()), acc_min=float(z["acc"].min()), acc_max=float(z["acc"].max()),
                         lag1_top_cell=float(lag[top[0]]), lag1_L0_SV1=float(lag[0]), lag1_L11_SV1=float(lag[110]),
                         lag1_median_120=float(np.median(lag)), lag1_max_abs_120=float(np.abs(lag).max()),
                         copole027_ms_null_mean_top=float(cpm[cpm.cell == top[0]].null_mean.iloc[0]),
                         copole027_expected_configuration_top=float(cpm[cpm.cell == top[0]].analytic_null_mean.iloc[0]),
                         copole027_max_abs_diff_null_vs_expected=float((cpm.null_mean - cpm.analytic_null_mean).abs().max()),
                         copole027_max_rel_diff=float(((cpm.null_mean - cpm.analytic_null_mean).abs() / cpm.analytic_null_mean).max()),
                         copole027_cells_checked=int(cpm.shape[0]), replay_identical=True,
                         replay_acc_mean=p2["ms_acc_replay_mean"]))
    sdf = pd.DataFrame(summ)
    sdf.to_csv(H / "h1_summary.csv", index=False)
    pd.DataFrame(diag).to_csv(H / "chain_diagnostics.csv", index=False)
    pl = sdf.set_index("context").family_min_p_adj
    verdict = bool(pl["lung"] < 0.05 and pl["lung_external"] < 0.05)
    res["primary"] = dict(rule="H1 supported if family max-T p < 0.05 in lung AND in lung_external (CONVENTION3 s4)",
                          family_min_p_adj=pl.to_dict(), supported=verdict,
                          summary=json.loads(sdf.to_json(orient="records")))
    # replication criterion: same sign + interval excluding the null value (top cell of lung, matched in ext)
    # ---------------------------------------------------------------- direction match lung -> lung_external
    cl, cx = E.load_ctx("lung"), E.load_ctx("lung_external")
    Xl, Vl, _ = E.coords("lung", cl["obs_idx"], want_vt=True)
    Xx, Vx, _ = E.coords("lung_external", cx["obs_idx"], want_vt=True)
    posx = {g: i for i, g in enumerate(cx["genes"])}
    shared = [(i, posx[g]) for i, g in enumerate(cl["genes"]) if g in posx]
    si = np.array([a for a, _ in shared]); sx = np.array([b for _, b in shared])
    dm = []
    p2l = json.loads((E.CKPT / "lung_pass2.json").read_text())
    dfx = prim["lung_external"]
    for rank, cell in enumerate(p2l["top_cells"]):
        L, k = E.cell_name(cell); v = Vl[L, k - 1]
        cs = np.abs(Vx[L] @ v); j = int(np.argmax(cs))
        allc = np.abs(np.einsum("lkd,d->lk", Vx, v)); La, ja = np.unravel_index(np.argmax(allc), allc.shape)
        rowx = dfx[(dfx.layer == L) & (dfx.sv == j + 1)].iloc[0]
        gene_corr = float(np.corrcoef(Xl[si, cell], Xx[sx, L * 10 + j])[0, 1])
        dm.append(dict(lung_rank=rank + 1, lung_layer=L, lung_sv=k, lung_r=float(prim["lung"].r.iloc[cell]),
                       ext_best_same_layer_sv=j + 1, abs_cos_same_layer=float(cs[j]),
                       abs_cos_second_best_same_layer=float(np.sort(cs)[-2]),
                       ext_matched_r=float(rowx.r), ext_matched_z=float(rowx.z), ext_matched_p_adj=float(rowx.p_adj_maxT),
                       gene_space_corr_matched=gene_corr, n_shared_genes=int(si.size),
                       ext_best_any_layer=int(La), ext_best_any_sv=int(ja) + 1, abs_cos_any_layer=float(allc[La, ja]),
                       random_unit_vector_abs_cos_expected=float(np.sqrt(2 / (np.pi * 512)))))
    # replication criterion (CONVENTION3 s1): same sign in lung_external and interval excluding the null value,
    # for the lung_external direction matched to each lung top cell; vertex-subsampling CI, B = 500
    pix, pjx = cx["edges"][:, 0], cx["edges"][:, 1]
    for d in dm:
        col = d["lung_layer"] * 10 + d["ext_best_same_layer_sv"] - 1; x = Xx[:, col]
        est, lo, hi, _ = S.vertex_subsample_ci(lambda nodes, a, b: S.assortativity(x, a, b), cx["n"], pix, pjx, 500,
                                               np.random.default_rng(E.seed_for("lung_external", "repl", col)))
        nm = float(dfx.null_mean.iloc[col])
        d.update(ext_matched_r_ci_lo=lo, ext_matched_r_ci_hi=hi, ext_matched_null_mean=nm,
                 replicates=bool(np.sign(est) == np.sign(d["lung_r"]) and (lo > nm or hi < nm)))
    pd.DataFrame(dm).to_csv(H / "h1_direction_match.csv", index=False)
    res["primary"]["direction_match_lung_to_ext"] = dm
    res["primary"]["replication_criterion"] = dict(
        rule="lung_external statistic has the same sign as lung and its interval excludes the null value (CONVENTION3 s1); "
             "evaluated on the lung_external direction best matched by |cos| of right-singular vectors to each lung top-3 cell",
        all_top3_replicate=bool(all(d["replicates"] for d in dm)))
    # ---------------------------------------------------------------- descriptive: what the top cells are
    from scipy.stats import spearmanr
    desc, poles = [], []
    for ctx in E.PRIMARY_CTX:
        c_ = E.load_ctx(ctx); Xc_ = Xl if ctx == "lung" else (Xx if ctx == "lung_external" else E.coords(ctx, c_["obs_idx"]))
        p2 = json.loads((E.CKPT / f"{ctx}_pass2.json").read_text())
        deg = np.bincount(c_["edges"].ravel(), minlength=c_["n"]); lno = c_["labels"]["log_nobs"]
        for rank, cell in enumerate(p2["top_cells"]):
            x = Xc_[:, cell]; L, k = E.cell_name(cell)
            desc.append(dict(context=ctx, top_rank=rank + 1, layer=L, sv=k, r=float(prim[ctx].r.iloc[cell]),
                             spearman_coord_vs_log_nobs=float(spearmanr(x, lno)[0]),
                             spearman_abscoord_vs_log_nobs=float(spearmanr(np.abs(x), lno)[0]),
                             spearman_coord_vs_degree=float(spearmanr(x, deg)[0]),
                             frac_genes_with_edges=float((deg > 0).mean())))
            if rank == 0:
                o = np.argsort(x)
                for end, idx in (("low", o[:20]), ("high", o[::-1][:20])):
                    for i, g in enumerate(idx):
                        poles.append(dict(context=ctx, layer=L, sv=k, end=end, rank=i + 1, gene=c_["genes"][g],
                                          coord=float(x[g]), degree_withinV=int(deg[g]),
                                          det_decile=int(c_["labels"]["det_decile"][g]),
                                          is_tf_lambert=bool(c_["labels"]["is_tf_lambert"][g])))
    pd.DataFrame(desc).to_csv(H / "h1_top_cell_descriptors.csv", index=False)
    pd.DataFrame(poles).to_csv(H / "h1_top_cell_extremes.csv", index=False)
    # ---------------------------------------------------------------- (a)-(c) alternative nulls
    nsum = []
    for kind, fname in (("label", "sens_label"), ("deg_v", "sens_degree_withinV"), ("deg_g", "sens_degree_genomewide"),
                        ("det", "sens_detection")):
        allrows = []
        for ctx in E.PRIMARY_CTX:
            z = np.load(E.CKPT / f"{ctx}_{kind}_{E.N_PRIMARY[ctx]}.npz")
            df = fam_rows(ctx, z, extra=dict(null=kind)); allrows.append(df)
            s = fam_summary(df, kind)
            t = prim[ctx].sort_values("rank_absz").iloc[0]
            tr = df[(df.layer == t.layer) & (df.sv == t.sv)].iloc[0]
            s.update(ms_top_layer=int(t.layer), ms_top_sv=int(t.sv), ms_top_z_under_this_null=float(tr.z),
                     ms_top_null_mean_this=float(tr.null_mean), ms_top_excess_r_this=float(tr.excess_r),
                     ms_top_p_adj_this=float(tr.p_adj_maxT), ms_top_z_under_ms=float(t.z),
                     z_ratio_this_over_ms_median=float(np.median(df.z.values / prim[ctx].z.values)))
            nsum.append(s)
        pd.concat(allrows).to_csv(H / f"{fname}.csv", index=False)
    ms_s = [dict(fam_summary(prim[c], "maslov_sneppen"), ms_top_z_under_this_null=fam_summary(prim[c], "ms")["top_z"])
            for c in E.PRIMARY_CTX]
    nsdf = pd.DataFrame(ms_s + nsum)
    nsdf.to_csv(H / "sens_nulls_summary.csv", index=False)
    res["sensitivity"]["alt_nulls"] = json.loads(nsdf.to_json(orient="records"))
    # ---------------------------------------------------------------- (d) co-pole K/n, (e) hubs + jackknife
    cpr, hub, jk, jks = [], [], [], []
    for ctx in E.PRIMARY_CTX:
        p2 = json.loads((E.CKPT / f"{ctx}_pass2.json").read_text())
        cp = pd.DataFrame(p2["copole"]); cp.insert(0, "context", ctx)
        vc = pd.DataFrame([v for v in p2["vertex_ci"] if v["stat"] != "r"])
        vc["null"] = vc.stat.str.replace("copole_ratio_", "").str.replace("_analytic", "")
        cp = cp.merge(vc[["cell", "pole_frac", "null", "est", "lo", "hi", "B"]].rename(
            columns=dict(est="ratio_analytic_est", lo="ratio_ci_lo", hi="ratio_ci_hi", B="ci_B")),
            on=["cell", "pole_frac", "null"], how="left")
        cp["is_top_cell"] = cp.cell.isin(p2["top_cells"])
        cp["top_rank"] = cp.cell.map({c: i + 1 for i, c in enumerate(p2["top_cells"])})
        cpr.append(cp)
        h = pd.DataFrame(p2["hub"]); h.insert(0, "context", ctx); hub.append(h)
        j = pd.DataFrame(p2["jackknife"]); j.insert(0, "context", ctx)
        L, k = E.cell_name(p2["top_cells"][0]); j.insert(1, "layer", L); j.insert(2, "sv", k); jk.append(j)
        jks.append(dict(context=ctx, layer=L, sv=k, **p2["jackknife_summary"]))
    pd.concat(cpr).to_csv(H / "sens_copole_K.csv", index=False)
    pd.concat(hub).to_csv(H / "sens_hub_exclusion.csv", index=False)
    pd.concat(jk).to_csv(H / "sens_jackknife_top10.csv", index=False)
    pd.DataFrame(jks).to_csv(H / "sens_jackknife_summary.csv", index=False)
    # ---------------------------------------------------------------- (f) co-expression matched
    ce = []
    for ctx in E.PRIMARY_CTX:
        d = json.loads((E.CKPT / f"{ctx}_coexpr.json").read_text())
        for src, sd in d["sources"].items():
            for name, v in sd["scorers"].items():
                ce.append(dict(context=ctx, coexpr_source=src, n_cells=sd["n_cells"], scorer=name,
                               top_layer=d["top_cell"]["layer"], top_sv=d["top_cell"]["sv"], m=d["m"],
                               n_matched_neg=sd["n_matched_neg"], shortfall=sd["shortfall"],
                               mean_abs_coexpr_pos=sd["mean_abs_coexpr_pos"], mean_abs_coexpr_pool=sd["mean_abs_coexpr_pool"],
                               mean_abs_coexpr_matched=sd["mean_abs_coexpr_matched"], **v))
    pd.DataFrame(ce).to_csv(H / "sens_coexpr_matched.csv", index=False)
    # ---------------------------------------------------------------- (g) complexes
    cx_rows, cx_s = [], []
    for ctx in ("lung", "lung_external"):
        for tag, lab in (("cp", "ComplexPortal"), ("cpx", "ComplexPortal minus STRING700"),
                         ("hm", "hu.MAP 2.0"), ("hmx", "hu.MAP 2.0 minus STRING700")):
            z = np.load(E.CKPT / f"{ctx}_{tag}_ms_2000.npz")
            df = fam_rows(ctx, z, extra=dict(edge_set=lab)); cx_rows.append(df)
            s = fam_summary(df, "maslov_sneppen"); s["edge_set"] = lab; s["ms_acc_mean"] = float(z["acc"].mean())
            t = prim[ctx].sort_values("rank_absz").iloc[0]; tr = df[(df.layer == t.layer) & (df.sv == t.sv)].iloc[0]
            s.update(string_top_cell=f"L{int(t.layer)} SV{int(t.sv)}", r_at_string_top=float(tr.r), z_at_string_top=float(tr.z),
                     p_adj_at_string_top=float(tr.p_adj_maxT))
            cx_s.append(s)
    pd.concat(cx_rows).to_csv(H / "sens_complexes.csv", index=False)
    pd.DataFrame(cx_s).to_csv(H / "sens_complexes_summary.csv", index=False)
    res["sensitivity"]["complexes"] = cx_s
    # ---------------------------------------------------------------- (h) robustness contexts
    rb_rows, rb_s = [], []
    for ctx in E.ROBUST_CTX:
        for kind in ("ms", "label"):
            z = np.load(E.CKPT / f"{ctx}_{kind}_1000.npz")
            df = fam_rows(ctx, z, extra=dict(null=kind)); rb_rows.append(df)
            s = fam_summary(df, kind)
            if kind == "ms": s["ms_acc_mean"] = float(z["acc"].mean())
            rb_s.append(s)
    pd.concat(rb_rows).to_csv(H / "sens_robustness.csv", index=False)
    pd.DataFrame(rb_s).to_csv(H / "sens_robustness_summary.csv", index=False)
    res["sensitivity"]["robustness"] = rb_s
    # ---------------------------------------------------------------- (i) bridge
    br = pd.read_csv(H / "bridge_kidneyA.csv")
    res["sensitivity"]["bridge_kidneyA"] = json.loads(br.to_json(orient="records"))
    res["prep"] = json.loads((H / "prep_meta.json").read_text())
    res["files"] = sorted(p.name for p in H.glob("*.csv") if not p.name.startswith("._"))
    (H / "results.json").write_text(json.dumps(res, indent=1, default=float))
    print(sdf.to_string()); print(pd.DataFrame(diag).to_string()); print(nsdf.to_string())


if __name__ == "__main__":
    main()
