"""Post-merge diagnostics (no new model fits): inference scale of the score coefficient (model SE vs permutation-null
SD vs two-way bootstrap SD), the H3-style decision quantities on the residual AUROC, and stratum-weight concentration
(pooled within-TF AUROC weights each TF by n_pos*n_neg; the conditional logit roughly by k(n-k)/n). Writes
inference_scale.csv, stratum_weights.csv and adds sections 'inference_scale', 'stratum_weights', 'residual_decision'
and 'provenance' to results.json. Run after e5_merge.py."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import e5_lib as L
HERE = Path(__file__).parent; CK = HERE / "ckpt"

rows = []
for ctx in L.CTXS:
    d = pd.read_csv(HERE / f"clogit_{ctx}.csv"); d = d[d.layer != "joint_mean_over_layers"]
    for _, r in d.iterrows():
        z = np.load(CK / f"{ctx}_L{int(r.layer):02d}.npz")
        bb = z["boot_adj"][:, 0]
        rows.append(dict(context=ctx, layer=int(r.layer), beta_adj=r.beta_adj, se_model=r.se_adj,
                         perm_null_sd=r.perm_null_sd_beta, boot_sd=r.boot_sd_beta_adj,
                         ratio_boot_sd_to_model_se=r.boot_sd_beta_adj / r.se_adj,
                         ratio_boot_sd_to_perm_sd=r.boot_sd_beta_adj / r.perm_null_sd_beta,
                         z_model=r.beta_adj / r.se_adj, z_boot=r.beta_adj / r.boot_sd_beta_adj,
                         boot_frac_beta_le0=float(np.mean(bb <= 0)),
                         auroc_boot_sd=float(np.std(z["boot_auc_raw"], ddof=1)),
                         z_boot_auroc=float((z["auc_raw"] - 0.5) / np.std(z["boot_auc_raw"], ddof=1)),
                         boot_frac_auroc_le05=float(np.mean(z["boot_auc_raw"] <= 0.5)),
                         boot_corr_beta_auroc=float(np.corrcoef(bb, z["boot_auc_raw"])[0, 1])))
inf = pd.DataFrame(rows); inf.to_csv(HERE / "inference_scale.csv", index=False)

sw = []
for ctx in L.CTXS:
    P = np.load(L.R3 / "cache" / ctx / "h3_pairs.npz"); tf, y = P["tf"], P["y"].astype(int)
    g = pd.DataFrame(dict(tf=tf, y=y)).groupby("tf").y.agg(["sum", "count"])
    g = g[(g["sum"] > 0) & (g["sum"] < g["count"])]
    k, n = g["sum"].to_numpy(float), g["count"].to_numpy(float)
    w_auc = k * (n - k); w_cl = k * (n - k) / n
    def top(w, m): return float(np.sort(w)[::-1][:m].sum() / w.sum())
    sw.append(dict(context=ctx, n_tf=int(g.shape[0]), max_pos_per_tf=int(k.max()), median_pos_per_tf=float(np.median(k)),
                   auroc_weight_top1=top(w_auc, 1), auroc_weight_top5=top(w_auc, 5), auroc_weight_top10=top(w_auc, 10),
                   clogit_info_top1=top(w_cl, 1), clogit_info_top5=top(w_cl, 5), clogit_info_top10=top(w_cl, 10),
                   eff_n_tf_auroc=float(w_auc.sum() ** 2 / (w_auc ** 2).sum()),
                   eff_n_tf_clogit=float(w_cl.sum() ** 2 / (w_cl ** 2).sum())))
swd = pd.DataFrame(sw); swd.to_csv(HERE / "stratum_weights.csv", index=False)

ra = pd.read_csv(HERE / "residual_auroc.csv")
dec = {}
for ctx, g in ra.groupby("context"):
    dec[ctx] = dict(all_layers_resid_upper95_onesided_lt_060=bool((g.auroc_resid_upper95_onesided < 0.60).all()),
                    max_resid_upper95_onesided=float(g.auroc_resid_upper95_onesided.max()),
                    n_layers_resid_ci_lo_gt_05=int((g.auroc_resid_lo > 0.5).sum()),
                    n_layers_raw_ci_lo_gt_05=int((g.auroc_raw_lo > 0.5).sum()),
                    mean_delta_resid_minus_raw=float(g.delta_resid_minus_raw.mean()),
                    n_layers_delta_ci_excludes_0=int(((g.delta_hi < 0) | (g.delta_lo > 0)).sum()))

res = json.loads((HERE / "results.json").read_text())
res["inference_scale"] = dict(
    file="inference_scale.csv",
    note="model SE = inverse-Hessian SE of the conditional logit (pairs independent given TF); perm_null_sd = SD of the "
         "adjusted coefficient under within-TF label permutation; boot_sd = SD under the two-way (TF x target) pigeonhole "
         "bootstrap. Pigeonhole variance is known to be somewhat upward-biased; pair-level SEs ignore target clustering.",
    by_context={c: dict(ratio_boot_sd_to_model_se_range=[float(g.ratio_boot_sd_to_model_se.min()), float(g.ratio_boot_sd_to_model_se.max())],
                        z_boot_range=[float(g.z_boot.min()), float(g.z_boot.max())],
                        z_model_range=[float(g.z_model.min()), float(g.z_model.max())],
                        max_boot_frac_beta_le0=float(g.boot_frac_beta_le0.max()),
                        z_boot_auroc_range=[float(g.z_boot_auroc.min()), float(g.z_boot_auroc.max())],
                        boot_corr_beta_auroc_range=[float(g.boot_corr_beta_auroc.min()), float(g.boot_corr_beta_auroc.max())])
                for c, g in inf.groupby("context")})
res["stratum_weights"] = dict(file="stratum_weights.csv", rows=swd.to_dict(orient="records"))
res["residual_decision"] = dict(note="H3-style quantities on the residualised score (B = 500 two-way bootstrap, not the "
                                     "pre-registered B = 2,000; descriptive, no decision is re-taken here)", by_context=dec)
res["provenance"] = dict(
    run="checkpointed per (context, layer) in ckpt/; the first agent's background jobs completed all 60 (context, layer) "
        "cells before this resumption; nothing was recomputed, merge and diagnostics re-run on 2026-10-26",
    refit_engines="bootstrap/permutation refits: full Newton for immune, kidney_A, lung L0-L1, lung_external L0-L1, kidney_B L0-L3; "
                  "chord-Newton (clogit_fit_fast, |beta diff| <= 2.5e-8 vs full Newton, check_fast.json) for the other layers; "
                  "point estimates and SEs always full Newton",
    convergence="all 60 x (500 bootstrap x 2 models + 500 permutation) refits converged (boot_conv / perm_conv in ckpt)",
    errata_item3="not applicable: E5 uses no vertex subsampling (edge statistics use the two-way TF x target bootstrap)",
    reductions="none relative to the E5 specification (B = 500, N = 500 as specified)")
(HERE / "results.json").write_text(json.dumps(res, indent=1, default=float))
print(inf.groupby("context")[["ratio_boot_sd_to_model_se", "ratio_boot_sd_to_perm_sd", "z_model", "z_boot", "boot_frac_beta_le0",
                              "z_boot_auroc", "boot_frac_auroc_le05", "boot_corr_beta_auroc"]].agg(["min", "max"]).round(3).T.to_string())
print(swd.round(3).to_string()); print(json.dumps(dec, indent=1))
