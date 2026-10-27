"""E4c verifier part 7: compare v5 independent refits with the implementer's fold-level values (bookkeeping)."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import t as tdist
H = Path(__file__).resolve().parent; E = H.parent
PF = pd.read_csv(E / "probe_folds.csv")
def nb(d, ratio):
    d = np.asarray(d, float); J = d.size; m = d.mean(); se = np.sqrt((1 / J + ratio) * d.var(ddof=1)); q = tdist.ppf(0.975, J - 1)
    return m, m - q * se, m + q * se
out = []
def imp(ctx, src, feat, layer):
    return PF[(PF.context == ctx) & (PF.model == "M1") & (PF.source == src) & (PF.feature == feat) & (PF.layer == layer) & (PF.regime == "C3")].set_index(["repeat", "fold"]).auc
for u in ("A", "B", "C", "Cn", "D"):
    f = H / f"v5_{u}.csv"
    if not f.exists(): continue
    v = pd.read_csv(f); v = v[v.label_col == 0].set_index(["repeat", "fold"])
    r = v.iloc[0]; i = imp(r.context, r.source, r.feature, int(r.layer)).reindex(v.index)
    out.append(dict(unit=u, context=r.context, source=r.source, feature=r.feature, layer=int(r.layer), n_folds=len(v),
                    verifier_mean=v.auc.mean(), implementer_mean_same_folds=i.mean(), max_abs_fold_diff=(v.auc - i).abs().max(),
                    corr_folds=np.corrcoef(v.auc, i)[0, 1]))
for a, b in (("A", "B"), ("C", "Cn")):
    if (H / f"v5_{a}.csv").exists() and (H / f"v5_{b}.csv").exists():
        va = pd.read_csv(H / f"v5_{a}.csv").set_index(["repeat", "fold"]); vb = pd.read_csv(H / f"v5_{b}.csv").set_index(["repeat", "fold"])
        d = va.auc - np.fmax(vb.auc.reindex(va.index), 0.5); ratio = va.n_test.mean() / va.n_train.mean()
        m, lo, hi = nb(d.values, ratio); m2, lo2, hi2 = nb(d.values, 0.25)
        out.append(dict(unit=f"{a}-{b}", context=va.context.iloc[0], source=va.source.iloc[0], feature="bilinear minus target-only",
                        layer=int(va.layer.iloc[0]), n_folds=len(d), verifier_mean=m, verifier_lo_pairs=lo, verifier_hi_pairs=hi,
                        verifier_lo_genes025=lo2, verifier_hi_genes025=hi2))
# null subset
if (H / "v5_N.csv").exists():
    v = pd.read_csv(H / "v5_N.csv"); ck = pd.read_csv(E / "partial/m1/lung/scgpt_bilinear_L2_C3_r0_null.csv")
    j = v.merge(ck[["fold", "label_col", "auc"]], on=["fold", "label_col"], suffixes=("_verifier", "_implementer"))
    per = j.groupby("label_col")[["auc_verifier", "auc_implementer"]].mean()
    per.to_csv(H / "v7_null_subset_per_draw.csv")
    out.append(dict(unit="N", context="lung", source="scgpt", feature="bilinear rewired draws 1-10 (r0)", layer=2, n_folds=int(j.fold.nunique()),
                    verifier_mean=per.loc[1:, "auc_verifier"].mean(), implementer_mean_same_folds=per.loc[1:, "auc_implementer"].mean(),
                    max_abs_fold_diff=(j.auc_verifier - j.auc_implementer).abs().max(),
                    true_label_verifier=per.loc[0, "auc_verifier"], true_label_implementer=per.loc[0, "auc_implementer"],
                    verifier_draws_ge_true=int((per.loc[1:, "auc_verifier"] >= per.loc[0, "auc_verifier"]).sum())))
O = pd.DataFrame(out); O.to_csv(H / "v7_independent_vs_implementer.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); print(O.round(4).to_string(index=False))
# direction (unit D)
if (H / "v5_D_direction.csv").exists():
    Dv = pd.read_csv(H / "v5_D_direction.csv"); pair = Dv.groupby(["i", "j"])["diff"].mean()
    Di = pd.concat([pd.read_csv(f) for f in sorted((E / "partial/m1/lung_external").glob("direction_L11_r*.csv"))])
    Di = Di[Di["repeat"].isin(Dv["repeat"].unique())]; pi = Di.groupby(["i", "j"])["diff"].mean()
    jj = pd.concat([pair.rename("v"), pi.rename("i")], axis=1, join="inner")
    res = dict(repeats=sorted(Dv["repeat"].unique().tolist()), n_pairs=len(pair), frac_forward_verifier=float((pair > 0).mean()),
               frac_forward_implementer_same_repeats=float((pi > 0).mean()), sign_agreement=float((np.sign(jj.v) == np.sign(jj.i)).mean()))
    pd.DataFrame([res]).to_csv(H / "v7_direction_L11_lung_external.csv", index=False); print(res)
