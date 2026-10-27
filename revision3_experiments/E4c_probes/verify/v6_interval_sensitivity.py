"""E4c verifier part 6 (bookkeeping on the implementer's fold-level CSVs; no fitting).
(a) Nadeau-Bengio ratio sensitivity. The implementer used rho-ratio = n_test/n_train counted in PAIRS. In C3 test pairs need both
    genes in the held-out fifth, so this is ~(1/5)^2/(4/5)^2 = 0.0625 instead of the 1/(K-1) = 0.25 of ordinary 5-fold CV. The
    exchangeable unit in C3 is the gene (gene fifth held out), whose test/train ratio is 0.25. Recompute every C3 M1 interval with
    ratio 0.25 (Bouckaert-Frank corrected repeated k-fold) and count layers whose interval excludes 0, under both ratios.
(b) Direction test: the implementer's exact binomial CI treats the 435-453 pairs as independent, but the pairs share genes (a TF
    appears in many pairs). Gene-cluster bootstrap: multinomial weights on genes, pair weight = w_i * w_j (pigeonhole product),
    B = 2000, percentile 95% CI of the weighted fraction forward > reverse.
Outputs: v6_nb_sensitivity.csv, v6_nb_counts.csv, v6_direction_cluster.csv
"""
from pathlib import Path
import glob
import numpy as np, pandas as pd
from scipy.stats import t as tdist
H = Path(__file__).resolve().parent; E = H.parent
PF = pd.read_csv(E / "probe_folds.csv"); NF = pd.read_csv(E / "null_folds.csv")
CT = ["lung", "lung_external", "kidney_B"]


def nb(d, ratio):
    d = np.asarray(d, float); d = d[~np.isnan(d)]; J = d.size; m = d.mean(); se = np.sqrt((1 / J + ratio) * d.var(ddof=1))
    q = tdist.ppf(0.975, J - 1); return m, m - q * se, m + q * se


def get(ctx, model, src, feat, layer, reg="C3"):
    g = PF[(PF.context == ctx) & (PF.model == model) & (PF.source == src) & (PF.feature == feat) & (PF.layer == layer) & (PF.regime == reg)]
    return g.sort_values(["repeat", "fold"]).set_index(["repeat", "fold"])


rows = []
for ctx in CT:
    for src, layers in (("scgpt", range(12)), ("static", [-1]), ("esm2", [-1])):
        for l in layers:
            tgo = get(ctx, "M1", src, "tg_only", l)
            cos = None
            for feat in ("concat", "hadamard", "absdiff", "bilinear", "tg_bilinear", "tg_hadamard"):
                p = get(ctx, "M1", src, feat, l)
                if not len(p): continue
                r_pair = p.n_test.mean() / p.n_train.mean()
                dn = p.auc - np.fmax(tgo.auc.reindex(p.index), 0.5)
                comps = {"delta_node": dn}
                if src == "scgpt":
                    cs = get(ctx, "ref", "scgpt", "cosine", l).auc.reindex(p.index); comps["delta_cosine"] = p.auc - cs
                cx = get(ctx, "ref", "coexpr", "abs_coexpr", -1).auc.reindex(p.index); comps["delta_abs_coexpr"] = p.auc - cx
                nf = NF[(NF.context == ctx) & (NF.model == "M1") & (NF.source == src) & (NF.feature == feat) & (NF.layer == l) & (NF.regime == "C3")]
                if len(nf):
                    nf = nf.set_index(["repeat", "fold"]).sort_index()
                    comps["delta_rewired"] = nf.obs_auc - nf.null_mean
                    comps["relational_signal"] = nf.obs_auc - np.fmax(nf.null_mean, np.fmax(tgo.auc.reindex(nf.index), 0.5))
                for name, d in comps.items():
                    J = int(d.notna().sum())
                    for lab, ratio in (("pairs", r_pair), ("genes_0.25", 0.25)):
                        m, lo, hi = nb(d.values, ratio)
                        rows.append(dict(context=ctx, source=src, feature=feat, layer=l, comparison=name, ratio_kind=lab, ratio=ratio, J=J,
                                         mean=m, lo=lo, hi=hi))
R = pd.DataFrame(rows); R.to_csv(H / "v6_nb_sensitivity.csv", index=False)
C = R.groupby(["context", "source", "feature", "comparison", "ratio_kind"]).agg(n_layers=("layer", "size"), mean_min=("mean", "min"),
        mean_max=("mean", "max"), n_lo_gt0=("lo", lambda x: int((x > 0).sum())), n_hi_lt0=("hi", lambda x: int((x < 0).sum()))).reset_index()
C.to_csv(H / "v6_nb_counts.csv", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 500)
w = C.pivot_table(index=["context", "source", "feature", "comparison"], columns="ratio_kind", values=["n_lo_gt0", "n_hi_lt0"]).astype(int)
print(w.to_string())
# headline single rows
for ctx, l in (("lung", 2), ("lung", 6)):
    print(R[(R.context == ctx) & (R.source == "scgpt") & (R.feature == "bilinear") & (R.layer == l)].round(4).to_string(index=False))
print(R[(R.context == "lung") & (R.source == "static") & (R.feature == "bilinear")].round(4).to_string(index=False))
# probe AUROC upper bounds vs 0.60 under ratio 0.25
ub = []
for (ctx, model, src, feat, l), g in PF[(PF.regime == "C3") & (PF.model == "M1") & PF.context.isin(CT)].groupby(["context", "model", "source", "feature", "layer"]):
    m, lo, hi = nb(g.auc.values - 0.5, 0.25); ub.append(dict(context=ctx, source=src, feature=feat, layer=l, auc=0.5 + m, hi_025=0.5 + hi))
UB = pd.DataFrame(ub); print("max upper bound (ratio 0.25) per context:", UB.groupby("context").hi_025.max().round(4).to_dict())
UB.to_csv(H / "v6_auc_upper_ratio025.csv", index=False)

# ---------------------------------------------------------------- (b) direction cluster bootstrap
rng = np.random.default_rng(20261027); out = []
for ctx in CT:
    fs = glob.glob(str(E / "partial" / "m1" / ctx / "direction_L*_r*.csv"))
    DD = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    for l, g in DD.groupby("layer"):
        pair = g.groupby(["i", "j"])["diff"].mean().reset_index()
        pair = pair[pair["diff"] != 0]
        genes = np.unique(np.r_[pair.i, pair.j]); gi = np.searchsorted(genes, pair.i); gj = np.searchsorted(genes, pair.j)
        f = (pair["diff"].values > 0).astype(float); obs = f.mean()
        Wg = rng.multinomial(genes.size, np.full(genes.size, 1 / genes.size), size=2000).astype(float)
        Wp = Wg[:, gi] * Wg[:, gj]; bs = (Wp @ f) / Wp.sum(1)
        out.append(dict(context=ctx, layer=int(l), n_pairs=len(pair), n_genes=genes.size, frac_forward=obs,
                        ci_lo_cluster=np.quantile(bs, .025), ci_hi_cluster=np.quantile(bs, .975),
                        n_repeats=int(g["repeat"].nunique())))
DC = pd.DataFrame(out); DT = pd.read_csv(E / "direction_test.csv")[["context", "layer", "frac_forward_higher", "ci_lo", "ci_hi"]]
DC = DC.merge(DT, on=["context", "layer"]); DC.to_csv(H / "v6_direction_cluster.csv", index=False)
print(DC.round(3).to_string(index=False))
print(DC.groupby("context").apply(lambda x: pd.Series(dict(n_lo_gt05_binom=int((x.ci_lo > .5).sum()), n_lo_gt05_cluster=int((x.ci_lo_cluster > .5).sum()),
      fmin=x.frac_forward.min(), fmax=x.frac_forward.max()))))
