"""Verifier: sensitivity of the TRRUST TF/target-only AUROC (estimator i) to the centring set, Tier A, all static models.
(a) centred over all 17,693 Tier A genes (E6 tierA_results convention); (b) centred over the 2,753 TRRUST-labelled genes only
(the convention of the legacy set and of the selection-bias null, where every gene is labelled). Writes verify/centring_sensitivity.csv."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_e6 as V
voc = {k: V.load_model(k) for k in V.MODELS}
A = sorted(set.intersection(*[set(v[0]) for v in voc.values()]))
T = V.trrust(); TF = set(T.TF); TGO = set(T.target) - TF
ttf = np.array([g in TF for g in A]); lab = ttf | np.array([g in TGO for g in A])
LAM = V.lambert(); lam = np.array([g in LAM for g in A])
e6 = pd.read_csv(V.E6 / "tierA_results.csv"); e6 = e6[e6.estimator == "full"]
rows = []
for k, (pos, E) in voc.items():
    X = E[[pos[g] for g in A]].astype(np.float64)
    sa = V.loo_scores(X, ttf)                     # centred over Tier A, TF set = all TRRUST TFs in Tier A
    sb = V.loo_scores(X[lab], ttf[lab])           # centred over labelled genes only
    la = V.loo_scores(X[lab], lam[lab]) if lam[lab].sum() > 2 else None
    r = dict(model=k, trrust_centre_tierA=roc_auc_score(ttf[lab], sa[lab]), trrust_centre_labelled=roc_auc_score(ttf[lab], sb),
             e6_tierA_value=float(e6[(e6.embedding == k) & (e6.task == "node_trrust")].value.iloc[0]))
    r["diff_labelled_minus_tierA"] = r["trrust_centre_labelled"] - r["trrust_centre_tierA"]
    rows.append(r); print(r, flush=True)
df = pd.DataFrame(rows); df.to_csv(V.HERE / "centring_sensitivity.csv", index=False)
from scipy.stats import kendalltau
sb_ = pd.read_csv(V.E6 / "selection_bias.csv"); x = sb_[(sb_.row_type == "subset_null") & (sb_.variant == "legacy178")].set_index("model")
df = df.set_index("model")
print("tau legacy178 vs tierA(centre tierA):", kendalltau(x.legacy_value, df.loc[x.index, "trrust_centre_tierA"])[0])
print("tau legacy178 vs tierA(centre labelled):", kendalltau(x.legacy_value, df.loc[x.index, "trrust_centre_labelled"])[0])
print("tau subset-null mean vs tierA(centre labelled):", kendalltau(x.null_mean, df.loc[x.index, "trrust_centre_labelled"])[0])
print(pd.concat([df, x[["legacy_value", "null_mean"]]], axis=1).round(3).to_string())
