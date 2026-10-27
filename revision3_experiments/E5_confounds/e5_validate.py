"""Validate the E5 conditional-logit engine against statsmodels ConditionalLogit.

Where statsmodels can run (all strata of immune and kidney_A; strata with n <= 120 in the large contexts), fit
the adjusted and unadjusted models with both and compare log-likelihoods and coefficients. Also evaluates the
statsmodels log-likelihood and score at the engine's optimum (must agree to ~1e-8). Writes
validation_statsmodels.csv and a timing line for the largest context.
"""
import sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import e5_lib as L
from statsmodels.discrete.conditional_models import ConditionalLogit

warnings.filterwarnings("ignore")
rows = []
for ctx, layers, nmax in [("kidney_A", [0, 5, 11], None), ("immune", [0, 11], None),
                          ("lung", [0], 120), ("kidney_B", [6], 120)]:
    C = L.load_context(ctx); emb = L.load_emb(ctx)
    tfs, tfc, tgs, tgc = L.tf_target_codes(C["tf"], C["tg"])
    sel = np.ones(C["y"].size, bool)
    if nmax is not None:
        sizes = np.bincount(tfc); sel = sizes[tfc] <= nmax
    covz = L.zscore(C["cov_raw"])
    for layer in layers:
        sc = L.zscore(L.layer_scores(emb, C["obs_idx"], layer, C["tf"], C["tg"]))
        for model, Xm in [("adjusted", np.column_stack([sc, covz])), ("unadjusted", sc[:, None])]:
            Xs, ys, gs = Xm[sel], C["y"][sel].astype(int), tfc[sel]
            D = L.CLogitData(Xs, ys, gs)
            t0 = time.time(); fit = L.clogit_fit(D); t_eng = time.time() - t0
            t0 = time.time()
            sm = ConditionalLogit(ys, Xs, groups=gs).fit(method="bfgs", maxiter=500, disp=False, gtol=1e-8)
            t_sm = time.time() - t0
            m = ConditionalLogit(ys, Xs, groups=gs)
            ll_sm_at_eng = m.loglike(fit["beta"]); sc_sm_at_eng = np.max(np.abs(m.score(fit["beta"])))
            rows.append(dict(context=ctx, layer=layer, model=model, n_rows=int(ys.size), n_strata=int(D.groups.size),
                             ll_engine=fit["ll"], ll_statsmodels_fit=float(sm.llf),
                             ll_statsmodels_at_engine_beta=float(ll_sm_at_eng),
                             abs_ll_diff_same_beta=abs(fit["ll"] - ll_sm_at_eng),
                             max_abs_score_statsmodels_at_engine_beta=float(sc_sm_at_eng),
                             beta_score_engine=fit["beta"][0], beta_score_statsmodels=float(sm.params[0]),
                             max_abs_beta_diff=float(np.max(np.abs(fit["beta"] - sm.params))),
                             se_score_engine=fit["se"][0], se_score_statsmodels=float(sm.bse[0]),
                             t_engine_s=t_eng, t_statsmodels_s=t_sm))
            print(rows[-1], flush=True)
df = pd.DataFrame(rows)
df.to_csv(Path(__file__).parent / "validation_statsmodels.csv", index=False)
print("max |ll diff| at same beta:", df.abs_ll_diff_same_beta.max(), " max |beta diff|:", df.max_abs_beta_diff.max())

# timing on the largest context (full strata, adjusted model) -- statsmodels cannot run here
C = L.load_context("lung"); emb = L.load_emb("lung")
tfs, tfc, tgs, tgc = L.tf_target_codes(C["tf"], C["tg"])
sc = L.zscore(L.layer_scores(emb, C["obs_idx"], 0, C["tf"], C["tg"]))
D = L.CLogitData(np.column_stack([sc, L.zscore(C["cov_raw"])]), C["y"], tfc)
t0 = time.time(); fit = L.clogit_fit(D); print("lung L0 adjusted full fit", time.time() - t0, "s", fit["iters"], "iters",
                                                fit["beta"], fit["converged"])
