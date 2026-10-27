"""E3 primary H2 (CONVENTION3 §4) for one context, plus node-level baselines, logistic adjustment,
restricted-V sensitivity (n_obs >= 3, 5, 10) and the legacy TRRUST TF vs target-only task.

usage: python e3_primary.py <context>
writes partial/primary_<ctx>.json and h2_primary_<ctx>.csv (baseline rows first)
"""
import json, sys, time
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2
import e3_lib as L
S = L.S

ctx = sys.argv[1]
CI = L.CONTEXTS.index(ctx)
N_PRIMARY = 10_000; N_SENS = 2_000; B_BOOT = 1_000; BATCH = 500
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

D = L.load(ctx); n = len(D["genes"]); y = D["is_tf"]; det = D["det"]
tf_pos = np.flatnonzero(y); k = tf_pos.size
log(ctx, "V", n, "TF", k)


def joint_null(En_list, ypos, strata, N, seed, fixed=()):
    """One within-stratum permutation per draw applied jointly to every layer (En_list) and to every
    fixed score (fixed = list of within-stratum rank vectors). Returns (N x len(En_list)), (N x len(fixed))."""
    fs = L.FastStrata(strata, ypos); kk = int(ypos.sum())
    rng = np.random.default_rng(seed)
    out = np.empty((N, len(En_list))); outf = np.empty((N, len(fixed))); done = 0
    while done < N:
        b = min(BATCH, N - done)
        Y = L.perm_batch(ypos, strata, rng, b).astype(np.float64)
        for li, En in enumerate(En_list):
            out[done:done + b, li] = fs.cols(L.loo_scores_cols(En, Y, kk), Y)
        for fi, r in enumerate(fixed):
            outf[done:done + b, fi] = L.strat_cols_fixed(r, Y, fs)
        done += b
    return out, outf


def self_check(En, ypos, strata, obs):
    fs = L.FastStrata(strata, ypos); Y0 = ypos[:, None].astype(float)
    a = fs.cols(L.loo_scores_cols(En, Y0, int(ypos.sum())), Y0)[0]
    a2 = L.strat_cols_varying(L.loo_scores_cols(En, Y0, int(ypos.sum())), Y0, fs)[0]
    assert abs(a - obs) < 1e-9 and abs(a2 - obs) < 1e-9, (a, a2, obs)


def layer_block(X_layers, ypos, strata, N, seed, boot=True, fixed_scores=None):
    """Observed per-layer stats + joint null + max-T. X_layers: list of centred matrices."""
    pos = np.flatnonzero(ypos)
    rows = []; En_list = []
    for l, X in enumerate(X_layers):
        s = S.tf_scores_loo(X, pos)
        a = S.stratified_auroc(s, ypos, strata); u = S.auroc(s[ypos], s[~ypos])
        En = S.l2norm_rows(X); En_list.append(En)
        if l == 0: self_check(En, ypos, strata, a)
        r = dict(layer=l, strat_auroc=a, auroc=u)
        if boot:
            bs = L.boot_strat_fast(s, ypos, strata, L.boot_weights(s.size, B_BOOT, np.random.default_rng(seed + 100 + l)))
            r.update(strat_lo=float(np.quantile(bs[:, 0], .025)), strat_hi=float(np.quantile(bs[:, 0], .975)),
                     auroc_lo=float(np.quantile(bs[:, 1], .025)), auroc_hi=float(np.quantile(bs[:, 1], .975)),
                     strat_boot_sd=float(bs[:, 0].std(ddof=1)))
        r["_score"] = s
        rows.append(r)
    fixed = []
    fixed_obs = []
    if fixed_scores:
        fs = L.FastStrata(strata, ypos)
        for sc in fixed_scores:
            fixed.append(L.strat_ranks_fixed(sc, fs, ypos.size)); fixed_obs.append(S.stratified_auroc(sc, ypos, strata))
    tn, tf_ = joint_null(En_list, ypos, strata, N, seed, fixed)
    t_obs = np.array([r["strat_auroc"] for r in rows])
    z, padj, praw = S.max_t_adjust(t_obs, tn)
    for l, r in enumerate(rows):
        r.update(z=float(z[l]), p_raw=float(praw[l]), p_maxT=float(padj[l]), null=L.null_summary(tn[:, l], t_obs[l]))
    fixed_res = []
    for fi in range(len(fixed)):
        p, b, NN = S.emp_p(fixed_obs[fi], tf_[:, fi])
        fixed_res.append(dict(strat_auroc=float(fixed_obs[fi]), p_raw=p, null=L.null_summary(tf_[:, fi], fixed_obs[fi])))
    return rows, fixed_res, tn


out = dict(context=ctx, conv_hash=L.conv_hash(), n_V=n, n_tf=int(k), n_strata=int(np.unique(det).size))
Xs = [L.centred(D["emb"], l, D["obs_idx"]) for l in range(12)]

# ---------------------------------------------------------------- baselines (scores, no embedding)
lmx, lme = L.log_mean(D["mean_extr"]), L.log_mean(D["mean_expr"])
base_scores = {"log_nobs": D["log_nobs"], "log_mean_expr_extraction": lmx, "log_mean_expr_population": lme}
log("primary null N=10000 (with baselines in the same draws)")
rows, fixed_res, tnull = layer_block(Xs, y, det, N_PRIMARY, 20261026 + 1000 * CI,
                                     fixed_scores=list(base_scores.values()))
S.assert_resolution(N_PRIMARY, 0.05)
np.save(L.HERE / "partial" / f"primary_null_{ctx}.npy", tnull.astype(np.float32))
baselines = []
for (nm, sc), fr in zip(base_scores.items(), fixed_res):
    bs = L.boot_strat_fast(sc, y, det, L.boot_weights(n, B_BOOT, np.random.default_rng(7 + CI)))
    baselines.append(dict(baseline=nm, n_genes=n, n_tf=int(k), auroc=S.auroc(sc[y], sc[~y]),
                          auroc_lo=float(np.quantile(bs[:, 1], .025)), auroc_hi=float(np.quantile(bs[:, 1], .975)),
                          strat_auroc=fr["strat_auroc"], strat_lo=float(np.quantile(bs[:, 0], .025)),
                          strat_hi=float(np.quantile(bs[:, 0], .975)), strat_p_raw=fr["p_raw"], N_perm=N_PRIMARY))
log("baselines (detection / expression) done")

# ---------------------------------------------------------------- static scGPT token embedding / ESM2
paired = []
for name in ("scgpt", "esm2"):
    z = np.load(L.HERE / "partial" / f"static_{name}_{ctx}.npz"); iv = z["idx"]
    E = z["emb"].astype(np.float64); nz = np.linalg.norm(E, axis=1) > 0; iv = iv[nz]; E = E[nz]
    ys, ds = y[iv], det[iv]
    for centre in (True, False):
        Ec = E - E.mean(0, keepdims=True) if centre else E
        s = S.tf_scores_loo(Ec, np.flatnonzero(ys))
        bs = L.boot_strat_fast(s, ys, ds, L.boot_weights(s.size, B_BOOT, np.random.default_rng(11 + CI)))
        rec = dict(baseline=f"static_{name}" + ("" if centre else "_uncentred"), n_genes=int(iv.size), n_tf=int(ys.sum()),
                   auroc=S.auroc(s[ys], s[~ys]), auroc_lo=float(np.quantile(bs[:, 1], .025)),
                   auroc_hi=float(np.quantile(bs[:, 1], .975)), strat_auroc=S.stratified_auroc(s, ys, ds),
                   strat_lo=float(np.quantile(bs[:, 0], .025)), strat_hi=float(np.quantile(bs[:, 0], .975)))
        if centre:
            En = S.l2norm_rows(Ec)
            tn, _ = joint_null([En], ys, ds, N_SENS, 20261026 + 1000 * CI + 77 + (name == "esm2"))
            rec.update(strat_p_raw=S.emp_p(rec["strat_auroc"], tn[:, 0])[0], N_perm=N_SENS)
            # paired difference vs each scGPT layer on the SAME genes and labels (layer V-centred, TF set = matched TFs)
            Wp = L.boot_weights(iv.size, B_BOOT, np.random.default_rng(21 + CI))
            bstat = L.boot_strat_fast(s, ys, ds, Wp)[:, 0]
            for l in range(12):
                sl = S.tf_scores_loo(Xs[l][iv], np.flatnonzero(ys))
                d_obs = rec["strat_auroc"] - S.stratified_auroc(sl, ys, ds)
                dd = bstat - L.boot_strat_fast(sl, ys, ds, Wp)[:, 0]
                paired.append(dict(context=ctx, static=name, layer=l, n_genes=int(iv.size), n_tf=int(ys.sum()),
                                   static_strat_auroc=rec["strat_auroc"], layer_strat_auroc=rec["strat_auroc"] - d_obs,
                                   delta=d_obs, delta_lo=float(np.quantile(dd, .025)), delta_hi=float(np.quantile(dd, .975))))
        baselines.append(rec)
    log("static", name, "done")
out["baselines"] = baselines
pd.DataFrame(paired).to_csv(L.HERE / "partial" / f"baselines_paired_{ctx}.csv", index=False)

# ---------------------------------------------------------------- logistic adjustment per layer
def zs(x): return (x - x.mean()) / x.std()
logit_rows = []
for l, r in enumerate(rows):
    for exname, ex in (("extraction", lmx), ("population", lme)):
        Xr = sm.add_constant(np.column_stack([zs(D["log_nobs"]), zs(ex)]))
        Xf = sm.add_constant(np.column_stack([zs(r["_score"]), zs(D["log_nobs"]), zs(ex)]))
        mr = sm.Logit(y.astype(float), Xr).fit(disp=0); mf = sm.Logit(y.astype(float), Xf).fit(disp=0)
        lr = 2 * (mf.llf - mr.llf); ci = mf.conf_int()[1]
        logit_rows.append(dict(context=ctx, layer=l, expr_covariate=exname, beta_score_per_sd=float(mf.params[1]),
                               beta_lo=float(ci[0]), beta_hi=float(ci[1]), OR_per_sd=float(np.exp(mf.params[1])),
                               LR=float(lr), LR_p=float(chi2.sf(lr, 1)), beta_lognobs=float(mf.params[2]),
                               beta_logexpr=float(mf.params[3]),
                               auroc_full=S.auroc(mf.fittedvalues[y], mf.fittedvalues[~y]),
                               auroc_covariates_only=S.auroc(mr.fittedvalues[y], mr.fittedvalues[~y])))
out["logistic"] = logit_rows
log("logistic done")

# ---------------------------------------------------------------- primary table (baselines first)
tab = [dict(row_type="baseline", **{k_: v for k_, v in b.items()}) for b in baselines]
for r in rows:
    tab.append(dict(row_type="H2_layer", baseline="", layer=r["layer"], n_genes=n, n_tf=int(k), strat_auroc=r["strat_auroc"],
                    strat_lo=r["strat_lo"], strat_hi=r["strat_hi"], auroc=r["auroc"], auroc_lo=r["auroc_lo"],
                    auroc_hi=r["auroc_hi"], z=r["z"], p_raw=r["p_raw"], p_maxT=r["p_maxT"], N_perm=N_PRIMARY,
                    null_mean=r["null"]["mean"], null_sd=r["null"]["sd"]))
cols = ["row_type", "baseline", "layer", "n_genes", "n_tf", "strat_auroc", "strat_lo", "strat_hi", "auroc", "auroc_lo",
        "auroc_hi", "z", "p_raw", "strat_p_raw", "p_maxT", "N_perm", "null_mean", "null_sd"]
pd.DataFrame(tab).reindex(columns=cols).to_csv(L.HERE / f"h2_primary_{ctx}.csv", index=False)
out["layers"] = [{k_: v for k_, v in r.items() if k_ != "_score"} for r in rows]
best = int(np.argmin([r["p_maxT"] for r in rows]))
out["min_p_maxT"] = float(min(r["p_maxT"] for r in rows)); out["best_layer_by_z"] = int(np.argmax([r["z"] for r in rows]))

# ---------------------------------------------------------------- (7) V restricted to n_obs >= 3, 5, 10
restr = []
for thr in (3, 5, 10):
    sub = D["nobs"] >= thr; ys = y[sub]; ds = S.quantile_strata(D["log_nobs"][sub], 10)
    Xr = [L.centred(D["emb"], l, D["obs_idx"], sub=np.flatnonzero(sub)) for l in range(12)]
    rr, _, _ = layer_block(Xr, ys, ds, N_SENS, 20261026 + 1000 * CI + thr, boot=True)
    for r in rr:
        restr.append(dict(context=ctx, min_nobs=thr, n_genes=int(sub.sum()), n_tf=int(ys.sum()), layer=r["layer"],
                          strat_auroc=r["strat_auroc"], strat_lo=r["strat_lo"], strat_hi=r["strat_hi"], auroc=r["auroc"],
                          z=r["z"], p_raw=r["p_raw"], p_maxT=r["p_maxT"], N_perm=N_SENS))
    log("restricted", thr, "done")
out["restricted"] = restr

# ---------------------------------------------------------------- (5) legacy task: TRRUST TF vs TRRUST target-only
U = np.flatnonzero(D["is_trrust_tf"] | D["is_trrust_tgt"]); yt = D["is_trrust_tf"][U]; dt = det[U]
assert not (D["is_trrust_tf"] & D["is_trrust_tgt"]).any()
Xt = [X[U] for X in Xs]                       # V-centred embedding, rows restricted to the task universe
rt, frt, _ = layer_block(Xt, yt, dt, N_PRIMARY, 20261026 + 1000 * CI + 500, boot=True,
                         fixed_scores=[D["log_nobs"][U], lmx[U]])
tr = []
for nm, fr, sc in zip(("log_nobs", "log_mean_expr_extraction"), frt, (D["log_nobs"][U], lmx[U])):
    tr.append(dict(context=ctx, row_type="baseline", baseline=nm, n_genes=int(U.size), n_tf=int(yt.sum()),
                   auroc=S.auroc(sc[yt], sc[~yt]), strat_auroc=fr["strat_auroc"], p_raw=fr["p_raw"], N_perm=N_PRIMARY))
for r in rt:
    tr.append(dict(context=ctx, row_type="layer", baseline="", layer=r["layer"], n_genes=int(U.size), n_tf=int(yt.sum()),
                   strat_auroc=r["strat_auroc"], strat_lo=r["strat_lo"], strat_hi=r["strat_hi"], auroc=r["auroc"],
                   auroc_lo=r["auroc_lo"], auroc_hi=r["auroc_hi"], z=r["z"], p_raw=r["p_raw"], p_maxT=r["p_maxT"],
                   N_perm=N_PRIMARY))
out["trrust_task"] = tr
log("trrust task done")
(L.HERE / "partial" / f"primary_{ctx}.json").write_text(json.dumps(out, indent=1, default=float))
log("DONE")
