"""E3 sensitivity (6): the preprint's 0.744 (ledger C28: joint 6-D SV2-SV7, immune, observed SVD).

Label set exactly as recorded in E9_claims_ledger.csv C28: the immune edge dataset
(cycle4_immune_main/cycle1_edge_dataset.tsv); its named genes that are observed; TF = any edge-dataset
source, target-only = the rest (expected 73 vs 222). Procedures:
 (a) the screen's: features (U*S)[:, 1:7] of stats3.svd_observed, StandardScaler + LogisticRegression(C=1),
     5-fold StratifiedKFold(shuffle, random_state=42), fold-mean AUROC per layer, then mean over layers;
     split-seed variability (seeds 0..49) and a 1,000-draw label-permutation null of the whole procedure
     (statistics: mean over layers, max over layers);
 (b) LOO mean cosine to the TF set within the 6 band coordinates (stats3.tf_scores_loo), unstratified and
     detection-stratified (node_labels det_decile), N = 10,000 joint label permutations (plain / within strata),
     max-T over 12 layers and the mean-over-layers statistic;
 (c) the same with the full 512-d V-centred embedding;
 (d) the screen's out-of-fold logistic predictions scored with the detection-stratified AUROC.
"""
import csv, json, time
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
import e3_lib as L
S = L.S
import contexts as X

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)
D = L.load("immune"); genes = D["genes"]; gpos = {g: i for i, g in enumerate(genes)}
rows = list(csv.DictReader(open(X.S38 / "cycle4_immune_main/cycle1_edge_dataset.tsv"), delimiter="\t"))
named = {}
for r in rows:
    named.setdefault(r["source"], int(r["source_idx"])); named.setdefault(r["target"], int(r["target_idx"]))
sources = {r["source"] for r in rows}
obs_set = set(D["obs_idx"].tolist()); rowpos = {r: i for i, r in enumerate(D["obs_idx"])}
obs_named = sorted(g for g in named if named[g] in obs_set)
tf_g = [g for g in obs_named if g in sources]; tg_g = [g for g in obs_named if g not in sources]
assert all(genes[rowpos[named[g]]] == g for g in obs_named)
lab = np.array([rowpos[named[g]] for g in tf_g + tg_g]); y = np.r_[np.ones(len(tf_g)), np.zeros(len(tg_g))].astype(bool)
det = D["det"][lab]
log("labels", len(named), "named;", len(obs_named), "observed;", y.sum(), "TF vs", (~y).sum(), "target-only")

def cv_auc(F, yy, seed=42):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    m = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=1000))
    fold = []
    for tr, te in skf.split(F, yy):
        m.fit(F[tr], yy[tr]); fold.append(roc_auc_score(yy[te], m.predict_proba(F[te])[:, 1]))
    return float(np.mean(fold))

US6, X512 = [], []
for l in range(12):
    US, Sv, Vt = S.svd_observed(D["emb"], l, D["obs_idx"]); US6.append(US[lab, 1:7].copy())
    X512.append(L.centred(D["emb"], l, D["obs_idx"])[lab])
out = dict(conv_hash=L.conv_hash(), n_named=len(named), n_named_observed=len(obs_named), n_tf=int(y.sum()),
           n_target_only=int((~y).sum()), preprint_value=0.744, e9_recomputed_mean=0.747)
yy = y.astype(int)
# (a) the screen's procedure
per = []
for l in range(12):
    F = US6[l]; skf = StratifiedKFold(5, shuffle=True, random_state=42)
    m = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=1000))
    oof = cross_val_predict(m, F, yy, cv=skf, method="predict_proba")[:, 1]
    per.append(dict(layer=l, cv_foldmean=cv_auc(F, yy), oof_auroc=S.auroc(oof[y], oof[~y]),
                    oof_strat_auroc=S.stratified_auroc(oof, y, det)))
seedvar = np.array([[cv_auc(US6[l], yy, seed=s) for l in range(12)] for s in range(50)])
log("screen procedure observed; seeds done")
rng = np.random.default_rng(20261026 + 744); NP = 1000; nullm = np.empty((NP, 12))
for b in range(NP):
    yp = rng.permutation(yy)
    nullm[b] = [cv_auc(US6[l], yp) for l in range(12)]
    if b % 200 == 0: log("screen null", b)
obs_mean = float(np.mean([p["cv_foldmean"] for p in per])); obs_max = float(np.max([p["cv_foldmean"] for p in per]))
out["screen_procedure"] = dict(per_layer=per, mean_over_layers=obs_mean, max_over_layers=obs_max,
    argmax_layer=int(np.argmax([p["cv_foldmean"] for p in per])),
    oof_strat_mean_over_layers=float(np.mean([p["oof_strat_auroc"] for p in per])),
    seed_variability=dict(n_seeds=50, mean_over_layers_by_seed_median=float(np.median(seedvar.mean(1))),
                          lo=float(np.quantile(seedvar.mean(1), .025)), hi=float(np.quantile(seedvar.mean(1), .975)),
                          per_layer_sd_median=float(np.median(seedvar.std(0, ddof=1)))),
    perm_null=dict(N=NP, mean_stat=L.null_summary(nullm.mean(1), obs_mean), max_stat=L.null_summary(nullm.max(1), obs_max),
                   p_mean=S.emp_p(obs_mean, nullm.mean(1))[0], p_max=S.emp_p(obs_max, nullm.max(1))[0],
                   p_floor=1 / (NP + 1)))

# (b), (c) LOO cosine, unstratified and stratified, joint null
def loo_block(Xl, strata, seed, N=10_000):
    k = int(y.sum()); pos = np.flatnonzero(y)
    obs = np.array([S.stratified_auroc(S.tf_scores_loo(Xm, pos), y, strata) for Xm in Xl])
    unst = np.array([S.auroc(*(lambda s: (s[y], s[~y]))(S.tf_scores_loo(Xm, pos))) for Xm in Xl])
    En = [S.l2norm_rows(Xm) for Xm in Xl]; fs = L.FastStrata(strata, y); r = np.random.default_rng(seed)
    nul = np.empty((N, 12)); done = 0
    while done < N:
        b = min(500, N - done); Yb = L.perm_batch(y, strata, r, b).astype(float)
        for l in range(12): nul[done:done + b, l] = fs.cols(L.loo_scores_cols(En[l], Yb, k), Yb)
        done += b
    z, padj, praw = S.max_t_adjust(obs, nul)
    W = L.boot_weights(y.size, 1000, np.random.default_rng(seed + 1))
    bs = [L.boot_strat_fast(S.tf_scores_loo(Xm, pos), y, strata, W)[:, 0] for Xm in Xl]
    return dict(per_layer=[dict(layer=l, auroc=float(obs[l]), lo=float(np.quantile(bs[l], .025)), hi=float(np.quantile(bs[l], .975)),
                                z=float(z[l]), p_raw=float(praw[l]), p_maxT=float(padj[l])) for l in range(12)],
                unstratified_auroc_per_layer=unst.tolist(),
                mean_over_layers=float(obs.mean()), p_mean_over_layers=S.emp_p(obs.mean(), nul.mean(1))[0],
                null_mean_over_layers=L.null_summary(nul.mean(1), obs.mean()), min_p_maxT=float(padj.min()), N=N)
S.assert_resolution(10_000, 0.05)
flat = np.zeros(y.size, int)
out["loo_cos_SV2_7"] = dict(unstratified=loo_block(US6, flat, 20261026 + 7441), stratified=loo_block(US6, det, 20261026 + 7442))
log("LOO 6-D done")
out["loo_cos_full512"] = dict(unstratified=loo_block(X512, flat, 20261026 + 7443), stratified=loo_block(X512, det, 20261026 + 7444))
log("LOO 512 done")
(L.HERE / "ledger_0744.json").write_text(json.dumps(out, indent=1, default=float))
log("DONE", obs_mean, out["loo_cos_SV2_7"]["unstratified"]["mean_over_layers"], out["loo_cos_SV2_7"]["stratified"]["mean_over_layers"])
