"""E4a — a priori power: minimum detectable effect (MDE) by cluster-respecting simulation.

For each benchmark design (every context's h3_pairs; and the submitted legacy design, kidney_A 49/147):
  1. design draw   resample TFs (clusters) with replacement; every drawn TF brings its whole block of
                   (target, label) rows; copies of one TF are distinct clusters; targets keep their identity
                   (so target re-use across TFs is preserved).
  2. score draw    binormal: s = d*y + u[target] (+ v[TF] for pooled designs) + e, total null variance 1 ->
                   population AUROC = Phi(d / sqrt 2) for a target value A (d = sqrt2 * Phi^-1(A)). The target
                   intraclass correlation rho_t (and rho_tf) is the empirical one of the observed scGPT scores
                   (within-TF rank-normal scores, one-way ANOVA estimator; median over the 12 layers).
  3. analysis      the benchmark's own statistic (within-TF stratified AUROC; pooled AUROC for the legacy
                   design as submitted) with a two-way TF x target cluster bootstrap SE (B_sim replicates);
                   detected if |A_hat - 0.5| > 1.96 SE (two-sided alpha = 0.05). [A Wald test with bootstrap
                   SE is used inside the simulation because a percentile CI needs far more replicates.]
  4. MDE           smallest A with power 0.8: probit fit of the per-simulation rejection indicator on the
                   effect over a grid around the bootstrap-Wald guess; the grid's empirical power is saved too.
Paired delta (scGPT minus |co-expression|): two scorers with the empirical within-TF correlation rho between
the scGPT score and |co-expression| (median over layers), co-expression at its observed AUROC, scGPT at
AUROC_cx + delta; paired cluster bootstrap SE of the difference on the same resamples.
Hanley-McNeil analytic MDE (stats3.mde_auroc; independence, no clustering) is reported for reference.
Usage: python e4a_power.py <ctx|legacy> ...
"""
from __future__ import annotations
import json, sys, time, zlib
from pathlib import Path
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent; R3 = HERE.parent
sys.path.insert(0, str(R3)); sys.path.insert(0, str(HERE))
import stats3 as S   # noqa: E402
import e4a_lib as L  # noqa: E402

POWER = 0.80; Z = norm.isf(0.025)


def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def normal_scores_within(codes, s):
    """Within-group rank -> normal scores."""
    z = np.empty(s.size)
    for c in np.unique(codes):
        m = np.flatnonzero(codes == c); r = s[m].argsort().argsort() + 0.5
        z[m] = norm.ppf(r / m.size)
    return z


def icc_oneway(groups, x):
    g, inv, n = np.unique(groups, return_inverse=True, return_counts=True)
    k = g.size; N = x.size
    if k < 2 or N - k < 1: return 0.0
    mg = np.bincount(inv, weights=x) / n; m = x.mean()
    msb = (n * (mg - m) ** 2).sum() / (k - 1); msw = ((x - mg[inv]) ** 2).sum() / (N - k)
    n0 = (N - (n ** 2).sum() / N) / (k - 1)
    sb = max(0.0, (msb - msw) / n0)
    return float(sb / (sb + msw)) if sb + msw > 0 else 0.0


def probit_mde(x, rej, target=POWER):
    """Fit P(reject) = Phi(a + b x) by ML; return x where it equals target."""
    x = np.asarray(x, float); rej = np.asarray(rej, float)
    def nll(p):
        q = np.clip(norm.cdf(p[0] + p[1] * x), 1e-9, 1 - 1e-9)
        return -(rej * np.log(q) + (1 - rej) * np.log(1 - q)).sum()
    sx = x.std() if x.std() > 0 else 1.0
    res = minimize(nll, np.array([-1.5, 1.0 / sx]), method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-8, fatol=1e-8))
    a, b = res.x
    return float((norm.ppf(target) - a) / b) if b > 0 else float("nan"), float(a), float(b)


class Design:
    def __init__(self, tf_code, tg_code, y):
        self.tf = np.asarray(tf_code); self.tg = np.asarray(tg_code); self.y = np.asarray(y).astype(bool)
        self.ntf = int(self.tf.max()) + 1; self.ntg = int(self.tg.max()) + 1
        self.blocks = [np.flatnonzero(self.tf == c) for c in range(self.ntf)]

    def draw(self, rng):
        pick = rng.integers(0, self.ntf, self.ntf)
        rows = np.concatenate([self.blocks[c] for c in pick])
        newtf = np.concatenate([np.full(self.blocks[c].size, j) for j, c in enumerate(pick)])
        return newtf, self.tg[rows], self.y[rows], pick


def perm_sd(sp, pooled):
    """Exact SD of the (stratified or pooled) AUROC under random labels within strata (no ties): per stratum
    Var(U) = n1 n0 (n + 1) / 12."""
    if pooled:
        n1 = int(sp.y.sum()); n0 = sp.y.size - n1
        return float(np.sqrt(n1 * n0 * (n1 + n0 + 1) / 12) / (n1 * n0))
    n1 = np.bincount(sp.tf_code, weights=sp.y.astype(float)); n = np.bincount(sp.tf_code); n0 = n - n1
    return float(np.sqrt((n1 * n0 * (n + 1) / 12).sum()) / (n1 * n0).sum())


def tau_hat(tf_code, y, s, min_n=3):
    """Between-TF SD of the binormal effect d_c = sqrt2 * Phi^-1(AUROC_c) beyond Hanley-McNeil sampling noise
    (method of moments, TFs with >= min_n positives and negatives). Returns (tau, d_c array)."""
    d, v = [], []
    for c in np.unique(tf_code):
        m = tf_code == c; n1 = int((m & y).sum()); n0 = int((m & ~y).sum())
        if n1 < min_n or n0 < min_n: continue
        a = float(np.clip(S.auroc(s[m & y], s[m & ~y]), 0.02, 0.98)); q = norm.ppf(a)
        d.append(np.sqrt(2) * q); v.append(S.hanley_mcneil_se(a, n1, n0) ** 2 * (np.sqrt(2) / norm.pdf(q)) ** 2)
    d = np.array(d); v = np.array(v)
    return float(np.sqrt(max(0.0, d.var(ddof=1) - v.mean()))) if d.size > 2 else 0.0, d


def simulate(design, effects, kind, rho_t, rng, nsim, bsim, rho=0.0, a_cx=0.5, pooled=False, rho_tf=0.0,
             with_iid=False, tau1=0.0, tau2=0.0, tau_corr=0.0):
    """Rows: effect, sim, est, se (two-way pigeonhole bootstrap = the pre-registered interval procedure),
    se_tf (TF-only cluster bootstrap), sd_perm (exact within-stratum permutation SD, AUROC only), and the
    corresponding two-sided alpha = 0.05 rejections; reject_iid = DeLong (independence; what the submitted
    analysis assumed) when with_iid. TF-level effect heterogeneity: each drawn TF copy gets its own effect
    d + tau * g (g ~ N(0,1); scGPT and co-expression slopes correlated tau_corr); with it the population mean
    AUROC is Phi(d / sqrt(2 + tau^2)), so d is set from the target value by that formula."""
    out = []
    for eff in effects:
        for i in range(nsim):
            tfc, tgc, y, pick = design.draw(rng)
            sp = S.StratifiedPairs(tfc, tgc, y)
            tfc, tgc, y = sp.tf, sp.tg, sp.y
            def noise():
                u = rng.normal(0, np.sqrt(rho_t), design.ntg)[tgc]
                v = rng.normal(0, np.sqrt(rho_tf), design.ntf)[pick[tfc]] if pooled else 0.0
                return u + v + rng.normal(0, np.sqrt(max(1 - rho_t - (rho_tf if pooled else 0), 1e-6)), y.size)
            fs = L.FastStrat(sp, "pooled" if pooled else "tf")
            wts = fs.boot_weights(bsim, rng, two_way=True); wts1 = fs.boot_weights(bsim, rng, two_way=False)
            g1 = rng.normal(size=design.ntf); g2 = tau_corr * g1 + np.sqrt(max(1 - tau_corr ** 2, 0)) * rng.normal(size=design.ntf)
            sl1 = tau1 * g1[tfc]; sl2 = tau2 * g2[tfc]
            if kind == "auroc":
                d = np.sqrt(2 + tau1 ** 2) * norm.ppf(min(0.5 + eff, 0.995))
                s = noise() + (d + sl1) * y
                est = fs.auc(s)[0] - 0.5
                se = np.nanstd(fs.auc(s, wts), ddof=1); se1 = np.nanstd(fs.auc(s, wts1), ddof=1); sdp = perm_sd(sp, pooled)
                r = dict(effect=eff, sim=i, est=est, se=se, se_tf=se1, sd_perm=sdp, reject=bool(abs(est) > Z * se),
                         reject_tf=bool(abs(est) > Z * se1), reject_perm=bool(abs(est) > Z * sdp))
                if with_iid:
                    dl = S.delong(s[y], s[~y]); r["reject_iid"] = bool(abs(est) > Z * np.sqrt(dl["var_a"]))
                out.append(r)
            else:
                z1 = noise(); z2 = rho * z1 + np.sqrt(1 - rho ** 2) * noise()
                d2 = np.sqrt(2 + tau2 ** 2) * norm.ppf(a_cx); d1 = np.sqrt(2 + tau1 ** 2) * norm.ppf(min(a_cx + eff, 0.995))
                s1 = z1 + (d1 + sl1) * y; s2 = z2 + (d2 + sl2) * y
                est = fs.auc(s1)[0] - fs.auc(s2)[0]
                se = np.nanstd(fs.auc(s1, wts) - fs.auc(s2, wts), ddof=1)
                se1 = np.nanstd(fs.auc(s1, wts1) - fs.auc(s2, wts1), ddof=1)
                r = dict(effect=eff, sim=i, est=est, se=se, se_tf=se1, reject=bool(abs(est) > Z * se), reject_tf=bool(abs(est) > Z * se1))
                if with_iid:
                    dl = S.delong(s1[y], s1[~y], s2[y], s2[~y]); r["reject_iid"] = bool(abs(est) > Z * dl["se_delta"])
                out.append(r)
    return out


def summarise(rows, key="reject"):
    eff = np.array([r["effect"] for r in rows]); rej = np.array([r[key] for r in rows], float)
    k = eff > 0
    mde, a, b = probit_mde(eff[k], rej[k])
    grid = []
    for e in np.unique(eff):
        rr = [r for r in rows if r["effect"] == e]
        g = dict(effect=float(e), power=float(rej[eff == e].mean()), n=len(rr),
                 mean_est=float(np.mean([r["est"] for r in rr])), sd_est=float(np.std([r["est"] for r in rr], ddof=1)),
                 mean_se_two_way=float(np.mean([r["se"] for r in rr])), mean_se_tf_only=float(np.mean([r["se_tf"] for r in rr])))
        if "sd_perm" in rr[0]: g["mean_sd_perm"] = float(np.mean([r["sd_perm"] for r in rr]))
        grid.append(g)
    return mde, grid


def summarise_all(rows, prefix, out, base=0.0):
    """MDE for every test in the rows; oracle MDE from the simulated sampling SD at the null."""
    keys = [("reject", "two_way_bootstrap"), ("reject_tf", "tf_only_bootstrap"), ("reject_perm", "permutation"), ("reject_iid", "iid_delong")]
    for key, lab in keys:
        if key not in rows[0]: continue
        mde, grid = summarise(rows, key)
        out[f"{prefix}_mde_{lab}"] = base + mde; out[f"{prefix}_size_{lab}"] = grid[0]["power"]; out[f"{prefix}_grid_{lab}"] = grid
    g0 = summarise(rows, "reject")[1][0]
    out[f"{prefix}_sampling_sd_at_null"] = g0["sd_est"]
    out[f"{prefix}_mde_oracle_sd"] = base + (Z + norm.ppf(POWER)) * g0["sd_est"]
    out[f"{prefix}_mean_se_two_way_at_null"] = g0["mean_se_two_way"]; out[f"{prefix}_mean_se_tf_only_at_null"] = g0["mean_se_tf_only"]
    if "mean_sd_perm" in g0: out[f"{prefix}_sd_perm_at_null"] = g0["mean_sd_perm"]


def sizes(n):
    if n < 5000: return 200, 200
    if n < 40000: return 100, 80
    return 80, 50


GF_PERM = [0.5, 0.75, 1.0, 1.25, 1.5]
GF_BOOT = [0.6, 0.8, 1.0, 1.2, 1.5]
ZZ = Z + norm.ppf(POWER)


def grid(sd_small, sd_large):
    g = {0.0} | {round(f * ZZ * sd_small, 6) for f in GF_PERM} | {round(f * ZZ * sd_large, 6) for f in GF_BOOT}
    return sorted(g)


def run_ctx(ctx):
    W = HERE / "work" / ctx; f = W / "power.json"
    if f.exists(): log(ctx, "power exists"); return
    h3 = np.load(R3 / "cache" / ctx / "h3_pairs.npz")
    sp = S.StratifiedPairs(h3["tf"], h3["tg"], h3["y"])
    P = np.column_stack([np.load(W / f"scores_L{l:02d}.npz")["prim"] for l in range(12)])
    cx = np.load(W / "coexp.npz")["cx_ext"]
    rows = [json.loads((W / f"rows_L{l:02d}.json").read_text())[0] for l in range(12)]
    a_cx = json.loads((W / "rows_coexp.json").read_text())[0]["auroc"]
    perm = json.loads((W / "perm.json").read_text())
    zc = normal_scores_within(sp.tf_code, cx)
    rho_t_l, rho_l = [], []
    for l in range(12):
        z = normal_scores_within(sp.tf_code, P[:, l])
        rho_t_l.append(icc_oneway(sp.tg_code, z)); rho_l.append(float(np.corrcoef(z, zc)[0, 1]))
    rho_t = float(np.median(rho_t_l)); rho = float(np.median(rho_l))
    taus = [tau_hat(sp.tf_code, sp.y, P[:, l]) for l in range(12)]
    tau1 = float(np.median([x[0] for x in taus])); tau2, dcx = tau_hat(sp.tf_code, sp.y, cx)
    tau_corr = float(np.median([np.corrcoef(x[1], dcx)[0, 1] for x in taus])) if dcx.size > 2 else 0.0
    se0 = float(np.median([r["boot_sd"] for r in rows])); sed = float(np.median([r["delta_cxext_sd"] for r in rows]))
    sdp = float(np.median([s["sd"] for s in perm["null_summary"]]))
    design = Design(sp.tf_code, sp.tg_code, sp.y)
    nsim, bsim = sizes(sp.tf.size)
    rng = np.random.default_rng([20261026, 77, zlib.crc32(ctx.encode())])
    t0 = time.time()
    res = dict(context=ctx, design="h3_pairs (within-TF stratified AUROC)", n_pos=int(sp.y.sum()), n_neg=int((~sp.y).sum()),
               n_tf=int(sp.tfs.size), n_targets=int(sp.tgs.size), rho_target_icc=rho_t, rho_target_icc_by_layer=rho_t_l,
               rho_scgpt_vs_coexp=rho, rho_by_layer=rho_l, a_cx_ext=a_cx, nsim_per_grid=nsim, B_sim=bsim,
               observed_boot_sd_two_way_median=se0, observed_delta_sd_two_way_median=sed, observed_perm_null_sd_median=sdp,
               observed_boot_sd_tf_only_median=float(np.median([(r["tfonly_ci_hi"] - r["tfonly_ci_lo"]) / (2 * Z) for r in rows])),
               tau_tf_effect_sd=tau1, tau_tf_effect_sd_by_layer=[x[0] for x in taus], tau_tf_effect_sd_coexp=tau2,
               tau_corr_scgpt_coexp=tau_corr)
    ra = simulate(design, grid(sdp, se0), "auroc", rho_t, rng, nsim, bsim, tau1=tau1)
    summarise_all(ra, "auroc", res, base=0.5)
    log(ctx, "AUROC MDEs", {k: round(v, 4) for k, v in res.items() if k.startswith("auroc_mde")}, f"{time.time() - t0:.0f}s")
    rd = simulate(design, grid(sed * sdp / se0, sed), "delta", rho_t, rng, nsim, bsim, rho=rho, a_cx=a_cx, tau1=tau1, tau2=tau2, tau_corr=tau_corr)
    summarise_all(rd, "delta", res, base=0.0)
    log(ctx, "delta MDEs", {k: round(v, 4) for k, v in res.items() if k.startswith("delta_mde")}, f"{time.time() - t0:.0f}s")
    n1, n0 = res["n_pos"], res["n_neg"]
    se_hm = S.hanley_mcneil_se(a_cx, n1, n0)
    res.update(mde_auroc_hanley_mcneil=S.mde_auroc(n1, n0), mde_delta_hanley_mcneil=S.mde_paired_delta(se_hm, se_hm, rho),
               mde_auroc_bootstrap_wald=0.5 + ZZ * se0, mde_delta_bootstrap_wald=ZZ * sed, mde_auroc_perm_wald=0.5 + ZZ * sdp,
               seconds=time.time() - t0)
    tmp = f.with_suffix(".tmp"); tmp.write_text(json.dumps(res, indent=1)); tmp.replace(f)


def run_legacy():
    f = HERE / "work" / "power_legacy.json"
    if f.exists(): log("legacy power exists"); return
    z = np.load(HERE / "work" / "legacy_design.npz")
    tfs, tf_code = np.unique(z["tf"], return_inverse=True); tgs, tg_code = np.unique(z["tg"], return_inverse=True)
    y = z["y"].astype(bool)
    atf = np.unique(z["all_tf"], return_inverse=True)[1]; atg = z["all_tg"]; s = z["all_s_full"]; sc = z["all_s_cx"]
    zs = norm.ppf((s.argsort().argsort() + 0.5) / s.size)
    rho_tf = icc_oneway(atf, zs)                        # TF component of the pooled score
    zw = normal_scores_within(atf, s); rho_t = icc_oneway(atg, zw)
    zcw = normal_scores_within(atf, sc); rho = float(np.corrcoef(zw, zcw)[0, 1])
    a_cx = S.auroc(z["s_cx"][y], z["s_cx"][~y]); a_obs = S.auroc(z["s_full"][y], z["s_full"][~y])
    import contexts as X
    _, trp = X.load_trrust(); trs = set(trp)
    ylab = np.array([(a, b) in trs for a, b in zip(z["all_tf"], z["all_tg"])])
    tau1, d1 = tau_hat(atf, ylab, s, min_n=2); tau2, d2 = tau_hat(atf, ylab, sc, min_n=2)
    tau_corr = float(np.corrcoef(d1, d2)[0, 1]) if d1.size > 2 else 0.0
    design = Design(tf_code, tg_code, y)
    rng = np.random.default_rng([20261026, 78])
    n1, n0 = int(y.sum()), int((~y).sum())
    se_hm0 = S.hanley_mcneil_se(0.5, n1, n0); se_hm = S.hanley_mcneil_se(a_cx, n1, n0)
    out = dict(design="submitted legacy design: kidney_A, run_trees.py split", n_pos=n1, n_neg=n0, n_tf=int(tfs.size),
               n_targets=int(tgs.size), rho_tf_icc=rho_tf, rho_target_icc_within_tf=rho_t, rho_full_vs_coexp=rho,
               observed_auroc_full512=a_obs, observed_auroc_coexp=a_cx, mde_auroc_hanley_mcneil=S.mde_auroc(n1, n0),
               mde_delta_hanley_mcneil=S.mde_paired_delta(se_hm, se_hm, rho), nsim_per_grid=400, B_sim=200,
               tau_tf_effect_sd=tau1, tau_tf_effect_sd_coexp=tau2, tau_corr_scgpt_coexp=tau_corr,
               tau_source="all candidate pairs of the 27 held-out TFs among the 187 genes, TRRUST labels (train+held)",
               n_pos_tau_source=int(ylab.sum()))
    t0 = time.time()
    for lab, pooled in (("pooled", True), ("stratified", False)):
        ra = simulate(design, grid(se_hm0, 1.6 * se_hm0), "auroc", rho_t, rng, 400, 200, pooled=pooled, rho_tf=rho_tf, with_iid=pooled, tau1=tau1)
        summarise_all(ra, f"{lab}_auroc", out, base=0.5)
        rd = simulate(design, grid(se_hm0 * 1.2, 2.0 * se_hm0), "delta", rho_t, rng, 400, 200, rho=rho, a_cx=a_cx, pooled=pooled, rho_tf=rho_tf, with_iid=pooled, tau1=tau1, tau2=tau2, tau_corr=tau_corr)
        summarise_all(rd, f"{lab}_delta", out, base=0.0)
        log("legacy", lab, {k: round(v, 4) for k, v in out.items() if "_mde_" in k}, f"{time.time() - t0:.0f}s")
    out["seconds"] = time.time() - t0
    f.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    for c in sys.argv[1:]:
        run_legacy() if c == "legacy" else run_ctx(c)
