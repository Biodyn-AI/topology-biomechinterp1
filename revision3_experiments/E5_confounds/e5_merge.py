"""Merge E5 checkpoints into clogit_<ctx>.csv, residual_auroc.csv, continuity_hashseed_spread.csv, results.json."""
import sys, json, time
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import e5_lib as L
S = L.S
HERE = Path(__file__).parent; CK = HERE / "ckpt"
NAMES = ["score"] + L.COV_NAMES


def q(a, lv):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return float(np.quantile(a, lv)) if a.size else float("nan")


def main():
    frozen, actual, ok = L.convention_hash()
    assert ok, "CONVENTION3.md hash mismatch"
    res = dict(experiment="E5_confounds", convention3_sha256_frozen=frozen, convention3_sha256_actual=actual,
               convention3_hash_ok=ok, created=time.strftime("%Y-%m-%d %H:%M:%S"),
               spec=dict(model="conditional logit, stratum = TF (exact conditional likelihood; engine validated "
                               "against statsmodels ConditionalLogit, validation_statsmodels.csv)",
                         score="z-scored signed cosine of V-centred full-dimensional (512) embeddings",
                         covariates_estimated=L.COV_NAMES,
                         covariates_absorbed_by_stratum=L.TF_COV_NAMES,
                         not_adjusted="TRRUST target in-degree (label leakage; matched by construction)",
                         coexp="Pearson r over the extraction cells (expr_extraction.npy), entered as log(|r| + 1e-3)",
                         bootstrap=f"two-way (TF x target) pigeonhole, B = {L.B_BOOT}, refit each time, percentile CI; "
                                   "the same draws for every layer and for the residual AUROC",
                         permutation=f"within-TF label permutation, N = {L.N_PERM}, one permutation for all 12 layers; "
                                     "statistic = mean over layers of the adjusted score coefficient",
                         residual_auroc="score residualised on the covariates by OLS after within-TF centring "
                                        "(TF-level covariates vanish), within-TF stratified AUROC; r_rb = 2*AUROC - 1",
                         seeds=dict(base=L.SEED)),
               contexts={})
    ra_rows = []
    for ctx in L.CTXS:
        files = [CK / f"{ctx}_L{l:02d}.npz" for l in range(L.N_LAYERS)]
        if not all(f.exists() for f in files):
            print(ctx, "incomplete:", sum(f.exists() for f in files), "/12"); continue
        red = json.loads((CK / f"{ctx}_reduced.json").read_text())
        Z = [dict(np.load(f)) for f in files]
        rows = []
        perm_b = np.column_stack([z["perm_adj"][:, 0] for z in Z])          # N x 12
        perm_ok = np.all(np.column_stack([z["perm_conv"] for z in Z]), axis=1)
        obs_b = np.array([float(z["beta_adj"][0]) for z in Z])
        zobs, p_adj, p_raw = S.max_t_adjust(obs_b, perm_b[perm_ok], tail="greater")
        for l, z in enumerate(Z):
            ba = z["boot_adj"]; bu = z["boot_unadj"]; bc = z["boot_conv"]
            r = dict(context=ctx, layer=l, n_pairs=red["n_rows"], n_pos=red["n_pos"], n_tf=red["n_strata"],
                     n_targets=red["n_targets"],
                     beta_adj=float(z["beta_adj"][0]), se_adj=float(z["se_adj"][0]),
                     or_adj=float(np.exp(z["beta_adj"][0])),
                     or_adj_wald_lo=float(np.exp(z["beta_adj"][0] - 1.96 * z["se_adj"][0])),
                     or_adj_wald_hi=float(np.exp(z["beta_adj"][0] + 1.96 * z["se_adj"][0])),
                     or_adj_boot_lo=float(np.exp(q(ba[bc[:, 0], 0], 0.025))),
                     or_adj_boot_hi=float(np.exp(q(ba[bc[:, 0], 0], 0.975))),
                     boot_sd_beta_adj=float(np.nanstd(ba[bc[:, 0], 0], ddof=1)), n_boot_conv_adj=int(bc[:, 0].sum()),
                     lrt_adj=float(z["lrt_adj"]), p_lrt_adj=float(z["p_lrt_adj"]),
                     beta_unadj=float(z["beta_unadj"][0]), se_unadj=float(z["se_unadj"][0]),
                     or_unadj=float(np.exp(z["beta_unadj"][0])),
                     or_unadj_boot_lo=float(np.exp(q(bu[bc[:, 1]], 0.025))),
                     or_unadj_boot_hi=float(np.exp(q(bu[bc[:, 1]], 0.975))),
                     n_boot_conv_unadj=int(bc[:, 1].sum()),
                     lrt_unadj=float(z["lrt_unadj"]), p_lrt_unadj=float(z["p_lrt_unadj"]),
                     attenuation_frac=float(1 - z["beta_adj"][0] / z["beta_unadj"][0]),
                     boot_attenuation_lo=q(1 - ba[bc.all(1), 0] / bu[bc.all(1)], 0.025),
                     boot_attenuation_hi=q(1 - ba[bc.all(1), 0] / bu[bc.all(1)], 0.975),
                     perm_null_mean_beta=float(perm_b[perm_ok, l].mean()), perm_null_sd_beta=float(perm_b[perm_ok, l].std(ddof=1)),
                     p_perm_layer_onesided=float(p_raw[l]), p_perm_layer_maxT=float(p_adj[l]),
                     beta_sens_log1p=float(z["beta_sens"][0]), or_sens_log1p=float(np.exp(z["beta_sens"][0])),
                     p_lrt_sens_log1p=float(z["p_lrt_sens"]),
                     conv_adj=bool(z["conv_adj"]), conv_unadj=bool(z["conv_unadj"]),
                     refit_engine=str(z["refit_engine"]) if "refit_engine" in z else "clogit_fit")
            for j, nm in enumerate(L.COV_NAMES):
                r[f"beta_{nm}"] = float(z["beta_adj"][j + 1]); r[f"se_{nm}"] = float(z["se_adj"][j + 1])
                r[f"corr_score_{nm}"] = float(z["corr_score_covs"][j])
            rows.append(r)
            # residual AUROC
            br, bs, bt = z["boot_auc_raw"], z["boot_auc_res"], z["boot_auc_res_tfonly"]
            ra_rows.append(dict(context=ctx, layer=l, auroc_raw=float(z["auc_raw"]),
                                auroc_raw_lo=q(br, 0.025), auroc_raw_hi=q(br, 0.975),
                                auroc_resid=float(z["auc_res"]), auroc_resid_lo=q(bs, 0.025), auroc_resid_hi=q(bs, 0.975),
                                auroc_resid_upper95_onesided=q(bs, 0.95),
                                auroc_resid_tfonly_lo=q(bt, 0.025), auroc_resid_tfonly_hi=q(bt, 0.975),
                                r_rb_resid=float(2 * z["auc_res"] - 1), r_rb_resid_lo=2 * q(bs, 0.025) - 1,
                                r_rb_resid_hi=2 * q(bs, 0.975) - 1,
                                delta_resid_minus_raw=float(z["auc_res"] - z["auc_raw"]),
                                delta_lo=q(bs - br, 0.025), delta_hi=q(bs - br, 0.975),
                                n_pos=red["n_pos"], n_pairs=red["n_rows"], n_tf=red["n_strata"], n_boot=int(br.size)))
        df = pd.DataFrame(rows)
        # joint permutation over layers
        T_obs = float(obs_b.mean()); T_null = perm_b[perm_ok].mean(1); N = int(T_null.size)
        floor = S.assert_resolution(N, 0.05, 1)
        S.assert_resolution(N, 0.05, len(L.CTXS))
        b1 = int(np.sum(T_null >= T_obs)); b2 = int(np.sum(np.abs(T_null - np.median(T_null)) >= abs(T_obs - np.median(T_null))))
        joint = dict(statistic="mean over 12 layers of the adjusted z-score coefficient", T_obs=T_obs,
                     or_at_T=float(np.exp(T_obs)), N=N, n_perm_all_layers_converged=N, b_greater=b1,
                     p_onesided=(b1 + 1) / (N + 1), b_two_sided=b2, p_two_sided=(b2 + 1) / (N + 1), p_floor=floor,
                     null_mean=float(T_null.mean()), null_sd=float(T_null.std(ddof=1)),
                     null_q=dict(zip(["0.025", "0.5", "0.975", "0.99"], np.quantile(T_null, [0.025, 0.5, 0.975, 0.99]).tolist())),
                     z_vs_null=float((T_obs - T_null.mean()) / T_null.std(ddof=1)))
        # bootstrap CI of the mean-over-layers coefficient (same draws across layers)
        bmat = np.column_stack([z["boot_adj"][:, 0] for z in Z]); bok = np.all(np.column_stack([z["boot_conv"][:, 0] for z in Z]), 1)
        umat = np.column_stack([z["boot_unadj"] for z in Z]); uok = np.all(np.column_stack([z["boot_conv"][:, 1] for z in Z]), 1)
        joint.update(T_boot_lo=q(bmat[bok].mean(1), 0.025), T_boot_hi=q(bmat[bok].mean(1), 0.975),
                     T_unadj=float(np.mean([z["beta_unadj"][0] for z in Z])),
                     T_unadj_boot_lo=q(umat[uok].mean(1), 0.025), T_unadj_boot_hi=q(umat[uok].mean(1), 0.975))
        jr = dict(context=ctx, layer="joint_mean_over_layers", beta_adj=T_obs, or_adj=float(np.exp(T_obs)),
                  or_adj_boot_lo=float(np.exp(joint["T_boot_lo"])), or_adj_boot_hi=float(np.exp(joint["T_boot_hi"])),
                  beta_unadj=joint["T_unadj"], or_unadj=float(np.exp(joint["T_unadj"])),
                  or_unadj_boot_lo=float(np.exp(joint["T_unadj_boot_lo"])), or_unadj_boot_hi=float(np.exp(joint["T_unadj_boot_hi"])),
                  perm_null_mean_beta=joint["null_mean"], perm_null_sd_beta=joint["null_sd"],
                  p_perm_layer_onesided=joint["p_onesided"])
        df = pd.concat([df, pd.DataFrame([jr])], ignore_index=True)
        df.to_csv(HERE / f"clogit_{ctx}.csv", index=False)
        lay = df[df.layer != "joint_mean_over_layers"]
        res["contexts"][ctx] = dict(
            n_pairs=red["n_rows"], n_pos=red["n_pos"], n_tf=red["n_strata"], n_targets=red["n_targets"],
            joint_permutation=joint,
            or_adj_range=[float(lay.or_adj.min()), float(lay.or_adj.max())],
            or_unadj_range=[float(lay.or_unadj.min()), float(lay.or_unadj.max())],
            attenuation_range=[float(lay.attenuation_frac.min()), float(lay.attenuation_frac.max())],
            min_or_adj_boot_lo=float(lay.or_adj_boot_lo.min()),
            max_p_lrt_adj=float(lay.p_lrt_adj.max()), n_layers_maxT_p_lt_005=int((lay.p_perm_layer_maxT < 0.05).sum()),
            best_layer_or_adj=int(lay.or_adj.idxmax()),
            reduced_model=dict(beta=dict(zip(L.COV_NAMES, red["beta_reduced"])), se=dict(zip(L.COV_NAMES, red["se_reduced"]))),
            n_boot_nonconverged_max=int(L.B_BOOT - lay.n_boot_conv_adj.min()))
        print(ctx, json.dumps({k: res["contexts"][ctx][k] for k in ("or_adj_range", "or_unadj_range", "attenuation_range")}),
              "joint p", joint["p_onesided"], flush=True)
    if ra_rows:
        ra = pd.DataFrame(ra_rows); ra.to_csv(HERE / "residual_auroc.csv", index=False)
        for ctx, g in ra.groupby("context"):
            res["contexts"].setdefault(ctx, {})["residual_auroc"] = dict(
                auroc_raw_range=[float(g.auroc_raw.min()), float(g.auroc_raw.max())],
                auroc_resid_range=[float(g.auroc_resid.min()), float(g.auroc_resid.max())],
                r_rb_resid_range=[float(g.r_rb_resid.min()), float(g.r_rb_resid.max())],
                min_auroc_resid_lo=float(g.auroc_resid_lo.min()), max_auroc_resid_hi=float(g.auroc_resid_hi.max()),
                delta_resid_minus_raw_range=[float(g.delta_resid_minus_raw.min()), float(g.delta_resid_minus_raw.max())])
    # validation
    v = pd.read_csv(HERE / "validation_statsmodels.csv")
    res["engine_validation"] = dict(file="validation_statsmodels.csv", n_fits=int(len(v)),
                                    max_abs_ll_diff_same_beta=float(v.abs_ll_diff_same_beta.max()),
                                    max_abs_statsmodels_score_at_engine_beta=float(v.max_abs_score_statsmodels_at_engine_beta.max()),
                                    max_abs_beta_diff_vs_statsmodels_bfgs=float(v.max_abs_beta_diff.max()),
                                    engine_ll_ge_statsmodels_ll=bool((v.ll_engine >= v.ll_statsmodels_fit - 1e-9).all()))
    # continuity
    cp = HERE / "continuity_kidneyA.csv"
    if cp.exists():
        c = pd.read_csv(cp); res["continuity_kidneyA"] = dict(file="continuity_kidneyA.csv",
                                                             reproduction=json.loads((HERE / "continuity_reproduction_check.json").read_text()))
        summ = {}
        for (vname, band), g in c.groupby(["variant", "band"]):
            gl = g[g.layer.astype(str) != "joint_mean_over_layers"]
            d = dict(r_rb_range=[float(gl.r_rb.min()), float(gl.r_rb.max())])
            if "p_two_sided" in gl and gl.p_two_sided.notna().any():
                d["n_layers_p_lt_005"] = int((gl.p_two_sided < 0.05).sum())
            gj = g[g.layer.astype(str) == "joint_mean_over_layers"]
            if len(gj): d["joint_mean_r_rb"] = float(gj.r_rb.iloc[0]); d["joint_p_two_sided"] = float(gj.p_two_sided.iloc[0])
            if "r_rb_ci_lo" in gl and gl.r_rb_ci_lo.notna().any():
                d["r_rb_ci_lo_min"] = float(gl.r_rb_ci_lo.min()); d["r_rb_ci_hi_max"] = float(gl.r_rb_ci_hi.max())
                d["n_layers_ci_excludes_0"] = int(((gl.r_rb_ci_lo > 0) | (gl.r_rb_ci_hi < 0)).sum())
            summ[f"{vname}|{band}"] = d
        res["continuity_kidneyA"]["summary"] = summ
    sp_dir = HERE / "continuity_hashseed_spread"
    if sp_dir.exists():
        parts = [pd.read_csv(f) for f in sorted(sp_dir.glob("hashseed_*.csv"))]
        if parts:
            spd = pd.concat(parts, ignore_index=True)
            agg = spd.groupby(["variant", "band", "layer"]).r_rb.agg(["count", "mean", "std", "min", "max"]).reset_index()
            sub = pd.read_csv(HERE.parents[1] / "revision2_experiments/M1_confounds_SV2_SV7.csv").assign(band="SV2_SV7")
            sub = pd.concat([sub, pd.read_csv(HERE.parents[1] / "revision2_experiments/M1_confounds_SV5_SV7.csv").assign(band="SV5_SV7")])
            sub = sub.rename(columns={"r_rb": "r_rb_submitted"})[["band", "layer", "r_rb_submitted"]]
            agg = agg.merge(sub.assign(variant="A_submitted_zero_row"), on=["variant", "band", "layer"], how="left")
            agg.to_csv(HERE / "continuity_hashseed_spread.csv", index=False)
            spm = spd.groupby(["variant", "band", "hashseed"]).r_rb.mean().reset_index()
            res.setdefault("continuity_kidneyA", {})["hashseed_spread"] = dict(
                file="continuity_hashseed_spread.csv", n_hashseeds=int(spd.hashseed.nunique()),
                mean_over_layers_by_variant={f"{a}|{b}": [float(g.r_rb.min()), float(g.r_rb.max())]
                                             for (a, b), g in spm.groupby(["variant", "band"])})
    (HERE / "results.json").write_text(json.dumps(res, indent=1, default=float))
    print("results.json written")


if __name__ == "__main__":
    main()
