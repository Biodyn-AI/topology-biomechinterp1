"""E10 step 2: random-effects variance decomposition of S1, S2, S2b_full, S2b_band, S3 across contexts.

For a group of k contexts with estimates theta (k,) and joint gene-bootstrap covariance Sigma (k x k, from
partial/stats_<ctx>.npz; contexts that share genes have correlated errors, which the joint Poisson bootstrap
retains), the moment estimator of the between-context variance is
    tau2 = S2_between - samp,   samp = (tr Sigma - 1'Sigma 1 / k) / (k - 1)
(E[S2_between] = tau2 + samp exactly under theta_i = mu + u_i + e_i, u ~ (0, tau2), e ~ (0, Sigma)).
I2 = tau2+ / (tau2+ + samp). Intervals:
  * cond  : gene-bootstrap interval conditional on these contexts (centred replicate deviations);
  * qprof : Q-profile interval treating contexts as random draws, generalised Q with full Sigma;
  * Q test: generalised Q(0) ~ chi2_{k-1} (heterogeneity beyond gene sampling).
Pre-stated decision rule (fixed before varcomp results were viewed): tissue dependence is claimed for a
statistic only if D = tau2_tissue - tau2_within_pooled has a 95% lower bound > 0 under the parametric
random-effects bootstrap (R = 10,000), where within_pooled pools resample (immune x3; kidney_B vs
kidney_ckpt_whole_human), donor (TSP14/21/25) and HVG-axis (kidney_A vs kidney_B) groups with weights k_g - 1.
The layer-mean statistic is the decision unit; per-layer rows are descriptive.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import chi2
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e10_common as K
S = K.S

GROUPS = {
    "resample_immune": ["immune", "immune_seed43", "immune_seed44"],
    "resample_kidney": ["kidney_B", "kidney_ckpt_whole_human"],
    "donor": ["immune_TSP14", "immune_TSP21", "immune_TSP25"],
    "celltype": ["immune_Bcell", "immune_CD4T", "immune_CD8T"],
    "hvg_axis": ["kidney_A", "kidney_B"],
    "checkpoint": ["kidney_ckpt_whole_human", "kidney_ckpt_kidney", "kidney_ckpt_brain"],
    "tissue": ["immune", "kidney_B", "lung", "lung_external"],
    # sensitivity splits of the tissue group
    "tissue3_no_lung_external": ["immune", "kidney_B", "lung"],
    "study_lung_vs_lung_external": ["lung", "lung_external"],
}
WITHIN = ["resample_immune", "resample_kidney", "donor", "hvg_axis"]
R_PARAM = 10000


def load_all():
    P, Bt, meta = {}, {}, {}
    for c in K.CTXS:
        z = np.load(K.HERE / "partial" / f"stats_{c}.npz")
        P[c] = z["point"]; Bt[c] = z["boot"]; meta[c] = json.loads(str(z["meta"]))
        stats = [str(s) for s in z["stats"]]
    return P, Bt, meta, stats


def samp_var(Sig):
    k = Sig.shape[0]
    return (np.trace(Sig) - Sig.sum() / k) / (k - 1)


def tau2_mom(theta, Sig):
    return float(np.var(theta, ddof=1) - samp_var(Sig))


def gen_Q(theta, Sig, t2):
    k = theta.size; M = Sig + (t2 + 1e-12) * np.eye(k)
    Wi = np.linalg.pinv(M); one = np.ones(k)
    mu = (one @ Wi @ theta) / (one @ Wi @ one); r = theta - mu
    return float(r @ Wi @ r)


def qprofile(theta, Sig, level=0.95):
    k = theta.size; df = k - 1
    lo_c, hi_c = chi2.ppf(1 - (1 - level) / 2, df), chi2.ppf((1 - level) / 2, df)
    f = lambda t2, c: gen_Q(theta, Sig, t2) - c
    big = 10.0 * (np.var(theta, ddof=1) + np.trace(Sig)) + 1e-6
    lo = 0.0 if f(0.0, lo_c) <= 0 else brentq(f, 0.0, big, args=(lo_c,))
    hi = 0.0 if f(0.0, hi_c) <= 0 else brentq(f, 0.0, big * 1e5, args=(hi_c,))
    q0 = gen_Q(theta, Sig, 0.0)
    return lo, hi, q0, float(chi2.sf(q0, df))


def group_rows(stat, layer, g, ctxs, th, bt):
    k = len(ctxs); Sig = np.cov(bt); t2 = tau2_mom(th, Sig); sv = samp_var(Sig)
    t2b = np.var(bt, axis=0, ddof=1) - sv
    dev = t2b - t2b.mean()
    lo_c, hi_c = t2 - np.quantile(dev, 0.975), t2 - np.quantile(dev, 0.025)
    qlo, qhi, q0, pq = qprofile(th, Sig)
    t2p = max(t2, 0.0)
    return dict(stat=stat, layer=layer, group=g, k=k, contexts=";".join(ctxs),
                estimates=";".join(f"{x:.4f}" for x in th), mean=float(th.mean()),
                range=float(th.max() - th.min()), S2_between=float(np.var(th, ddof=1)),
                sampling_var=float(sv), tau2=t2, tau=float(np.sqrt(t2p)),
                tau2_cond_lo=float(lo_c), tau2_cond_hi=float(hi_c), tau2_qprof_lo=qlo, tau2_qprof_hi=qhi,
                I2=float(t2p / (t2p + sv)) if (t2p + sv) > 0 else np.nan, Q=q0, p_Q=pq,
                se_boot_mean=float(np.sqrt(np.diag(Sig)).mean())), Sig, t2b


def main():
    P, Bt, meta, stats = load_all()
    rng = np.random.default_rng(K.SEED + 20)
    S.assert_resolution(R_PARAM, 0.05, len(stats))       # decision family: 5 statistics (layer-mean)
    rows, drows = [], []
    for si, stat in enumerate(stats):
        for layer in list(range(K.N_LAYERS)) + ["mean"]:
            res = {}
            for g, ctxs in GROUPS.items():
                if layer == "mean":
                    th = np.array([P[c][si].mean() for c in ctxs]); bt = np.array([Bt[c][si].mean(0) for c in ctxs])
                else:
                    th = np.array([P[c][si, layer] for c in ctxs]); bt = np.array([Bt[c][si, layer] for c in ctxs])
                ok = np.all(np.isfinite(bt), axis=0); bt = bt[:, ok]
                r, Sig, t2b = group_rows(stat, layer, g, ctxs, th, bt)
                rows.append(r); res[g] = (th, Sig, r, t2b)
            # ---- D = tau2_tissue - pooled within-tissue floor
            dfw = {g: len(GROUPS[g]) - 1 for g in WITHIN}
            pooled = sum(dfw[g] * res[g][2]["tau2"] for g in WITHIN) / sum(dfw.values())
            maxw = max(res[g][2]["tau2"] for g in WITHIN)
            D = res["tissue"][2]["tau2"] - pooled
            # conditional (gene-bootstrap) interval: same replicates for every group (joint bootstrap)
            nb = min(res[g][3].size for g in WITHIN + ["tissue"])
            Db = res["tissue"][3][:nb] - sum(dfw[g] * res[g][3][:nb] for g in WITHIN) / sum(dfw.values())
            dev = Db - Db.mean()
            # parametric random-effects bootstrap (contexts as random draws) and test under H0: tau2_t = floor
            def sim(t2_tissue):
                out = np.empty(R_PARAM)
                draws = {}
                for g in WITHIN + ["tissue"]:
                    th, Sig, r, _ = res[g]; k = th.size
                    t2g = max(r["tau2"], 0.0) if g != "tissue" else max(t2_tissue, 0.0)
                    e = rng.multivariate_normal(np.zeros(k), Sig + 1e-14 * np.eye(k), R_PARAM, method="eigh")
                    u = rng.normal(0, np.sqrt(t2g), (R_PARAM, k))
                    thr = r["mean"] + u + e
                    draws[g] = np.var(thr, axis=1, ddof=1) - samp_var(Sig)
                out = draws["tissue"] - sum(dfw[g] * draws[g] for g in WITHIN) / sum(dfw.values())
                return out
            Dsim = sim(res["tissue"][2]["tau2"])
            ddev = Dsim - Dsim.mean()
            D0 = sim(max(pooled, 0.0))
            p_param = (np.sum(D0 >= D) + 1) / (R_PARAM + 1)
            drows.append(dict(stat=stat, layer=layer, group="D_tissue_minus_within_pooled", k=4,
                              contexts=";".join(GROUPS["tissue"]), tau2_tissue=res["tissue"][2]["tau2"],
                              tau2_within_pooled=pooled, tau2_within_max=maxw,
                              tau2_within_components=";".join(f"{g}={res[g][2]['tau2']:.3g}" for g in WITHIN),
                              D=D, D_cond_lo=float(D - np.quantile(dev, 0.975)), D_cond_hi=float(D - np.quantile(dev, 0.025)),
                              D_param_lo=float(D - np.quantile(ddev, 0.975)), D_param_hi=float(D - np.quantile(ddev, 0.025)),
                              p_param_H0_tissue_eq_floor=float(p_param), N_param=R_PARAM,
                              D_vs_max=res["tissue"][2]["tau2"] - maxw,
                              tissue_dependence_claimed=bool((D - np.quantile(ddev, 0.975)) > 0) if layer == "mean" else None))
    out = pd.DataFrame(rows + drows)
    out.to_csv(K.HERE / "varcomp.csv", index=False)
    m = out[(out.layer.astype(str) == "mean")]
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(m[m.group.isin(list(GROUPS))][["stat", "group", "mean", "range", "tau", "tau2", "sampling_var", "I2", "p_Q",
                                          "tau2_qprof_lo", "tau2_qprof_hi"]].to_string())
    print(m[m.group == "D_tissue_minus_within_pooled"][["stat", "tau2_tissue", "tau2_within_pooled", "tau2_within_max", "D",
                                                         "D_cond_lo", "D_param_lo", "D_param_hi", "p_param_H0_tissue_eq_floor",
                                                         "tissue_dependence_claimed"]].to_string())


if __name__ == "__main__":
    main()
