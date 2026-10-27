"""E5 main runner: per (context, layer) conditional logit (stratum = TF), two-way cluster bootstrap (B = 500,
refit each time), joint within-TF label permutation draws (N = 500, the same permutations for every layer), and
the within-TF residualised AUROC with the same bootstrap draws. One checkpoint per (context, layer) in ckpt/.

Model (adjusted):   y ~ z(signed cosine, layer L) + z(log(|r_coexp| + 1e-3)) + z(log mean expr target)
                        + z(log n_obs target) + z(STRING-physical genome-wide degree target) + z(GO count target)
TF-level covariates (log mean expression, log n_obs, STRING degree, GO count of the TF) are constant inside every
TF stratum and cancel from the conditional likelihood exactly, i.e. they are adjusted for by the stratification
itself and cannot be (and are not) estimated. TRRUST target in-degree is NOT used (label leakage; non-edges are
matched on it by construction).
Unadjusted:          y ~ z(signed cosine)
Reduced (LRT null):  y ~ covariates only (layer-independent; fitted once per context)
Sensitivity (point estimate + LRT only): degree and GO count entered as log1p instead of raw.

usage: python e5_run.py <context> [<context> ...]
"""
import sys, time, json
from pathlib import Path
import numpy as np
from scipy.stats import chi2
sys.path.insert(0, str(Path(__file__).parent))
import e5_lib as L
S = L.S

CK = Path(__file__).parent / "ckpt"; CK.mkdir(exist_ok=True)
# Refit engine for bootstrap / permutation refits. Layers checkpointed before 15:35 on 2026-10-26 (immune, kidney_A,
# lung L0-L1, lung_external L0-L1, kidney_B L0-L3) used full Newton (clogit_fit); later layers use the chord-Newton
# clogit_fit_fast, which reaches the same maximiser (|beta diff| <= 2.5e-8, check_fast.json) at ~2x lower cost.
# Point estimates and SEs always use clogit_fit. The engine of each layer is stored in its checkpoint.
REFIT = L.clogit_fit_fast


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def run_context(ctx, reverse=False):
    C = L.load_context(ctx); emb = L.load_emb(ctx)
    tf, tg, y = C["tf"], C["tg"], C["y"].astype(np.int64)
    tfs, tfc, tgs, tgc = L.tf_target_codes(tf, tg)
    sp = S.StratifiedPairs(tf, tg, y)
    assert sp.tf.size == y.size and np.array_equal(sp.tf_code, tfc)
    covz = L.zscore(C["cov_raw"])
    cov_s = C["cov_raw"].copy(); cov_s[:, 3] = np.log1p(cov_s[:, 3]); cov_s[:, 4] = np.log1p(cov_s[:, 4])
    covz_s = L.zscore(cov_s)
    # reduced models (layer-independent)
    red_path = CK / f"{ctx}_reduced.json"
    if not red_path.exists():
        fr = L.clogit_fit(L.CLogitData(covz, y, tfc)); frs = L.clogit_fit(L.CLogitData(covz_s, y, tfc))
        D0 = L.CLogitData(np.zeros((y.size, 1)), y, tfc)
        red_path.write_text(json.dumps(dict(ll_reduced=fr["ll"], beta_reduced=fr["beta"].tolist(),
                                            se_reduced=fr["se"].tolist(), conv=fr["converged"],
                                            ll_reduced_sens=frs["ll"], ll_null=D0.ll0(),
                                            n_rows=int(y.size), n_strata=int(D0.groups.size),
                                            n_pos=int(y.sum()), n_targets=int(tgs.size))))
    red = json.loads(red_path.read_text())
    # bootstrap and permutation draws (identical for all layers)
    boots = [L.boot_weights(ctx, b, tfs.size, tgs.size) for b in range(L.B_BOOT)]
    perms = [L.perm_labels(ctx, r, sp).astype(np.int64) for r in range(L.N_PERM)]
    for layer in (range(L.N_LAYERS - 1, -1, -1) if reverse else range(L.N_LAYERS)):
        out = CK / f"{ctx}_L{layer:02d}.npz"
        lock = CK / f"{ctx}_L{layer:02d}.lock"
        if reverse and lock.exists():
            log(ctx, layer, "locked by another worker, skipped"); continue
        if reverse: lock.write_text("reverse worker")
        if out.exists():
            log(ctx, layer, "checkpoint exists, skipped"); continue
        t0 = time.time()
        raw = L.layer_scores(emb, C["obs_idx"], layer, tf, tg)
        sc = L.zscore(raw)
        Xa = np.column_stack([sc, covz]); Xu = sc[:, None]; Xs = np.column_stack([sc, covz_s])
        fa = L.clogit_fit(L.CLogitData(Xa, y, tfc))
        fu = L.clogit_fit(L.CLogitData(Xu, y, tfc))
        fs = L.clogit_fit(L.CLogitData(Xs, y, tfc))
        lrt_a = 2 * (fa["ll"] - red["ll_reduced"]); lrt_u = 2 * (fu["ll"] - red["ll_null"])
        lrt_s = 2 * (fs["ll"] - red["ll_reduced_sens"])
        # residualised score (within-TF OLS) and AUROC
        res, gam = L.residualise_within_tf(raw, C["cov_raw"], tfc)
        auc_raw = sp.auroc(raw); auc_res = sp.auroc(res)
        log(ctx, f"L{layer}", f"adj beta={fa['beta'][0]:+.4f} se={fa['se'][0]:.4f} | unadj {fu['beta'][0]:+.4f} | "
            f"AUROC raw {auc_raw:.4f} resid {auc_res:.4f}")
        bb_a = np.full((L.B_BOOT, Xa.shape[1]), np.nan); bb_u = np.full(L.B_BOOT, np.nan)
        bconv = np.zeros((L.B_BOOT, 2), bool)
        ba_raw = np.empty(L.B_BOOT); ba_res = np.empty(L.B_BOOT)
        for b, (wtf, wtg) in enumerate(boots):
            Xe, ye, ge = L.expand(Xa, y, tfc, wtf, wtg, tgc)
            f1 = REFIT(L.CLogitData(Xe, ye, ge, gw_by_group=wtf), beta0=fa["beta"])
            f2 = REFIT(L.CLogitData(Xe[:, :1], ye, ge, gw_by_group=wtf), beta0=fu["beta"])
            bb_a[b] = f1["beta"]; bb_u[b] = f2["beta"][0]; bconv[b] = (f1["converged"], f2["converged"])
            w = wtf[tfc] * wtg[tgc]
            r_b, _ = L.residualise_within_tf(raw, C["cov_raw"], tfc, w)
            ba_raw[b] = sp.auroc(raw, w); ba_res[b] = sp.auroc(r_b, w)
            if (b + 1) % 100 == 0:
                log(ctx, f"L{layer} boot {b+1}/{L.B_BOOT} ({time.time()-t0:.0f}s)")
        rng_tf = np.random.default_rng([L.SEED, 6, L.CTX_ID[ctx], layer])
        bt_res = sp.bootstrap(res, L.B_BOOT, rng_tf, two_way=False)      # TF-only check (fixed residual)
        pb = np.full((L.N_PERM, Xa.shape[1]), np.nan); pconv = np.zeros(L.N_PERM, bool)
        for r, yp in enumerate(perms):
            f1 = REFIT(L.CLogitData(Xa, yp, tfc))
            pb[r] = f1["beta"]; pconv[r] = f1["converged"]
            if (r + 1) % 100 == 0:
                log(ctx, f"L{layer} perm {r+1}/{L.N_PERM} ({time.time()-t0:.0f}s)")
        np.savez(out, beta_adj=fa["beta"], se_adj=fa["se"], ll_adj=fa["ll"], conv_adj=fa["converged"],
                 beta_unadj=fu["beta"], se_unadj=fu["se"], ll_unadj=fu["ll"], conv_unadj=fu["converged"],
                 beta_sens=fs["beta"], se_sens=fs["se"], ll_sens=fs["ll"],
                 lrt_adj=lrt_a, lrt_unadj=lrt_u, lrt_sens=lrt_s,
                 p_lrt_adj=chi2.sf(lrt_a, 1), p_lrt_unadj=chi2.sf(lrt_u, 1), p_lrt_sens=chi2.sf(lrt_s, 1),
                 boot_adj=bb_a, boot_unadj=bb_u, boot_conv=bconv,
                 perm_adj=pb, perm_conv=pconv,
                 auc_raw=auc_raw, auc_res=auc_res, boot_auc_raw=ba_raw, boot_auc_res=ba_res, boot_auc_res_tfonly=bt_res,
                 resid_ols_coef=gam, score_mean=raw.mean(), score_sd=raw.std(),
                 corr_score_covs=np.array([np.corrcoef(raw, C["cov_raw"][:, j])[0, 1] for j in range(C["cov_raw"].shape[1])]),
                 refit_engine=REFIT.__name__, seconds=time.time() - t0)
        log(ctx, f"L{layer} done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    rev = "--reverse" in sys.argv          # second worker: layers 11 -> 0 (results are deterministic per layer)
    for c in [a for a in sys.argv[1:] if not a.startswith("--")]:
        run_context(c, reverse=rev)
    log("ALL DONE", sys.argv[1:])
