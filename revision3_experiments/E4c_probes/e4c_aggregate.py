"""Aggregate E4c checkpoints into probe_results.csv, null_results.csv, relational_signal.csv, direction_test.csv,
probe_folds.csv, null_folds.csv and results.json. Pure bookkeeping: no model fitting."""
import json, glob
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats as sst
import e4c_lib as L
S = L.S
H = L.HERE
KEYS = ["context", "model", "source", "feature", "layer", "regime"]


def load():
    tr, nu = [], []
    for f in glob.glob(str(H / "partial/m1/*/*.csv")):
        if Path(f).name.startswith("direction_"): continue
        (nu if f.endswith("_null.csv") else tr).append(pd.read_csv(f))
    T = pd.concat(tr, ignore_index=True); N = pd.concat(nu, ignore_index=True) if nu else pd.DataFrame()
    nl = []
    for m in ("hgb", "mlp", "hgb_tuned", "mlp_tuned"):
        rows = [json.loads(Path(f).read_text()) for f in glob.glob(str(H / f"partial/{m}/*/*.json"))]
        if rows: nl.append(pd.DataFrame(rows))
    if nl:
        X = pd.concat(nl, ignore_index=True)
        T = pd.concat([T, X[X.label_col == 0]], ignore_index=True)
        N = pd.concat([N, X], ignore_index=True)          # nonlin null rows (incl. col 0 of the same folds)
    ref = [pd.read_csv(f) for f in glob.glob(str(H / "partial/ref/*.csv"))]   # unsupervised scorers, same folds
    if ref: T = pd.concat([T, pd.concat(ref, ignore_index=True)], ignore_index=True)
    return T, N


def nb(d, ratio):
    """Nadeau-Bengio corrected resampled t: var = (1/J + n_test/n_train) s^2, df = J-1."""
    d = np.asarray(d, float); d = d[~np.isnan(d)]; J = d.size
    if J < 2: return dict(mean=np.nan, se=np.nan, lo=np.nan, hi=np.nan, t=np.nan, p=np.nan, J=J)
    m = d.mean(); s2 = d.var(ddof=1); se = np.sqrt((1 / J + ratio) * s2)
    tq = sst.t.ppf(0.975, J - 1)
    t = m / se if se > 0 else np.inf * np.sign(m)
    return dict(mean=m, se=se, lo=m - tq * se, hi=m + tq * se, t=t, p=float(2 * sst.t.sf(abs(t), J - 1)), J=J)


def main():
    T, N = load()
    T["layer"] = T["layer"].astype(int)
    T.sort_values(KEYS + ["repeat", "fold"]).to_csv(H / "probe_folds.csv", index=False)
    # ---------------------------------------------------------------- probe_results.csv
    out = []
    for k, g in T.groupby(KEYS):
        ratio = g.n_test.mean() / g.n_train.mean(); r = nb(g.auc - 0.5, ratio)
        rep_means = g.groupby("repeat").auc.mean()
        out.append(dict(zip(KEYS, k), n_folds=len(g), n_repeats=g["repeat"].nunique(), auc_mean=g.auc.mean(),
                        auc_sd_repeats=(rep_means.std(ddof=1) if rep_means.size > 1 else np.nan),
                        auc_sd_folds=g.auc.std(ddof=1), auc_ci_lo=0.5 + r["lo"], auc_ci_hi=0.5 + r["hi"],
                        nb_t_vs_0p5=r["t"], nb_p_vs_0p5=r["p"], n_within_tf_folds=int((g.metric_type == "within_tf").sum()),
                        n_pooled_folds=int((g.metric_type == "pooled").sum()), mean_n_tf_used=g.n_tf_used.mean(),
                        mean_n_pos=g.n_pos.mean(), mean_n_neg=g.n_neg.mean(), mean_n_test=g.n_test.mean(),
                        mean_n_train=g.n_train.mean(),
                        alpha_mode=(g.alpha.mode().iloc[0] if "alpha" in g and g.alpha.notna().any() else np.nan)))
    P = pd.DataFrame(out); P.to_csv(H / "probe_results.csv", index=False)
    # ---------------------------------------------------------------- null summaries
    nrows, nfold, null_stat = [], [], {}
    if len(N):
        N["layer"] = N["layer"].astype(int)
        for k, g in N.groupby(KEYS):
            piv = g.pivot_table(index=["repeat", "fold"], columns="label_col", values="auc")
            if 0 not in piv.columns or piv.shape[1] < 2: continue
            obs_f = piv[0].values; nul = piv.drop(columns=0)
            obs = float(np.mean(obs_f)); draws = nul.mean(0).values          # statistic: mean AUROC over folds
            Nn = draws.size; p, b, _ = S.emp_p(obs, draws, "greater")
            try: floor = S.assert_resolution(Nn, 0.05); claim = True
            except ValueError: floor = 1 / (Nn + 1); claim = False
            null_stat[k] = (obs, draws)
            # consistency with the true-label run on the same folds (M1: separate K=1 run)
            t_same = T.set_index(KEYS + ["repeat", "fold"]).auc
            try:
                tv = np.array([t_same.loc[k + idx] for idx in piv.index]); dmax = float(np.max(np.abs(tv - obs_f)))
            except KeyError: dmax = np.nan
            nrows.append(dict(zip(KEYS, k), n_folds=len(piv), obs_auc_mean=obs, null_N=Nn, null_mean=draws.mean(),
                              null_sd=draws.std(ddof=1), null_q025=np.quantile(draws, .025), null_q50=np.median(draws),
                              null_q95=np.quantile(draws, .95), null_q975=np.quantile(draws, .975), b=b, p_emp=p,
                              z=(obs - draws.mean()) / draws.std(ddof=1), resolution_floor=floor,
                              significance_claim_allowed=claim, max_abs_diff_vs_true_run=dmax))
            for (rep, fold), row in piv.iterrows():
                v = row.drop(0).values
                nfold.append(dict(zip(KEYS, k), repeat=rep, fold=fold, obs_auc=row[0], null_mean=np.nanmean(v),
                                  null_sd=np.nanstd(v, ddof=1), null_N=int(np.isfinite(v).sum())))
    NR = pd.DataFrame(nrows); NF = pd.DataFrame(nfold)
    # max-T over the 12 layers (M1; same rewiring draws for every layer) per context/source/feature/regime
    if len(NR):
        NR["p_maxT_layers"] = np.nan
        for k, g in NR[NR.model == "M1"].groupby(["context", "model", "source", "feature", "regime"]):
            if g.layer.nunique() < 2: continue
            ks = [tuple(r[c] for c in KEYS) for _, r in g.iterrows()]
            Nmin = min(null_stat[x][1].size for x in ks)
            t_obs = np.array([null_stat[x][0] for x in ks]); t_null = np.column_stack([null_stat[x][1][:Nmin] for x in ks])
            S.assert_resolution(Nmin, 0.05)
            _, padj, _ = S.max_t_adjust(t_obs, t_null, "greater")
            NR.loc[g.index, "p_maxT_layers"] = padj
    NR.to_csv(H / "null_results.csv", index=False); NF.to_csv(H / "null_folds.csv", index=False)
    # ---------------------------------------------------------------- relational signal
    TI = T.set_index(KEYS + ["repeat", "fold"]).auc
    rel = []
    probes = T[~T.feature.isin(["tf_only", "tg_only", "coexpr_marg"]) & (T.model != "ref")]
    for k, g in probes.groupby(KEYS):
        ctx, model, src, feat, layer, reg = k
        def base(f, mdl=model, s=src):
            out = []
            for _, r in g.iterrows():
                key = (ctx, mdl, s, f, layer, reg, r["repeat"], r["fold"])
                out.append(TI.get(key, np.nan))
            return np.array(out, float)
        flags = []
        if src == "coexpr":
            nb_ = base("coexpr_marg", "M1"); tfb = nb_; tgb = nb_
            if model != "M1": flags.append("node baseline = M1 co-expression marginal-only model")
        else:
            tgb = base("tg_only")
            if model == "M1": tfb = base("tf_only")
            else:
                tfb = np.where(g.metric_type.values == "within_tf", 0.5, base("tf_only", "M1"))
                flags.append("TF-only = 0.5 by construction (within-TF metric)")
            if np.all(np.isnan(tgb)) and model != "M1":
                tgb = base("tg_only", "M1"); flags.append("target-only from M1")
        node = np.fmax(tfb, tgb)
        ratio = g.n_test.mean() / g.n_train.mean()
        dn = nb(g.auc.values - node, ratio)
        row = dict(zip(KEYS, k), n_folds=len(g), probe_auc=g.auc.mean(), tf_only_auc=np.nanmean(tfb),
                   tg_only_auc=np.nanmean(tgb), node_max_auc=np.nanmean(node),
                   delta_node=dn["mean"], delta_node_lo=dn["lo"], delta_node_hi=dn["hi"], delta_node_p=dn["p"],
                   notes="; ".join(flags))
        # same-fold comparison with the unsupervised scorers (no training): |co-expression| and the H3 primary
        # scorer (signed cosine, full 512-d) at the probe's layer (layer-free sources: at the best-mean layer)
        def ref_vals(src_, feat_, lay_):
            return np.array([TI.get((ctx, "ref", src_, feat_, lay_, reg, r["repeat"], r["fold"]), np.nan)
                             for _, r in g.iterrows()], float)
        cx = ref_vals("coexpr", "abs_coexpr", -1)
        if layer >= 0: cs = ref_vals("scgpt", "cosine", layer); cs_l = layer
        else:
            cands = {l_: ref_vals("scgpt", "cosine", l_) for l_ in range(12)}
            cs_l = max(cands, key=lambda l_: np.nanmean(cands[l_])); cs = cands[cs_l]
        if np.isfinite(cx).any():
            dc = nb(g.auc.values - cx, ratio); dcs = nb(g.auc.values - cs, ratio)
            row.update(abs_coexpr_auc=np.nanmean(cx), delta_abs_coexpr=dc["mean"], delta_abs_coexpr_lo=dc["lo"],
                       delta_abs_coexpr_hi=dc["hi"], cosine_ref_layer=cs_l, cosine_auc=np.nanmean(cs),
                       delta_cosine=dcs["mean"], delta_cosine_lo=dcs["lo"], delta_cosine_hi=dcs["hi"])
        if len(NF):
            nf = NF[(NF.context == ctx) & (NF.model == model) & (NF.source == src) & (NF.feature == feat) &
                    (NF.layer == layer) & (NF.regime == reg)]
            if len(nf):
                nf = nf.sort_values(["repeat", "fold"])
                gi = g.set_index(["repeat", "fold"])
                idx = [(a, b) for a, b in zip(nf["repeat"], nf["fold"])]
                nodef = pd.Series(node, index=pd.MultiIndex.from_frame(g[["repeat", "fold"]])).reindex(idx).values
                ratio_n = gi.loc[idx].n_test.mean() / gi.loc[idx].n_train.mean()
                dr = nb(nf.obs_auc.values - nf.null_mean.values, ratio_n)
                comb = nb(nf.obs_auc.values - np.fmax(nf.null_mean.values, nodef), ratio_n)
                nr = NR[(NR.context == ctx) & (NR.model == model) & (NR.source == src) & (NR.feature == feat) &
                        (NR.layer == layer) & (NR.regime == reg)].iloc[0]
                row.update(null_folds=len(nf), null_N=int(nr.null_N), rewired_auc=nf.null_mean.mean(),
                           delta_rewired=dr["mean"], delta_rewired_lo=dr["lo"], delta_rewired_hi=dr["hi"],
                           delta_rewired_p_nb=dr["p"], p_emp_rewired=nr.p_emp, p_maxT_layers=nr.get("p_maxT_layers", np.nan),
                           rewired_claim_allowed=bool(nr.significance_claim_allowed),
                           relational_signal=comb["mean"], relational_signal_lo=comb["lo"],
                           relational_signal_hi=comb["hi"], relational_signal_p_nb=comb["p"])
        rel.append(row)
    R = pd.DataFrame(rel)
    for col in ("delta_node_p", "relational_signal_p_nb"):
        if col in R:
            R[col.replace("_p", "_q_bh").replace("_nb", "")] = np.nan
            for _, g in R.groupby(["context", "model", "regime"]):
                ok = g[col].notna()
                if ok.any(): R.loc[g.index[ok], col.replace("_p", "_q_bh").replace("_nb", "")] = S.bh_fdr(g.loc[ok, col].values)
    R.to_csv(H / "relational_signal.csv", index=False)
    # ---------------------------------------------------------------- direction test
    drs = []
    for f in glob.glob(str(H / "partial/m1/*/direction_L*_r*.csv")):
        drs.append(pd.read_csv(f))
    D = []
    if drs:
        DD = pd.concat(drs, ignore_index=True)
        node_info = {}
        for ctx in DD.context.unique():
            e = np.load(L.CACHE / ctx / "edges_trrust.npy"); G = len(json.loads((L.CACHE / ctx / "obs_genes.json").read_text()))
            outdeg = np.bincount(e[:, 0], minlength=G); nobs = np.load(L.CACHE / ctx / "ncells.npy")
            indeg = np.bincount(e[:, 1], minlength=G)
            anyr = np.load(L.CACHE / ctx / "any_regulatory.npy")
            node_info[ctx] = (outdeg, nobs, indeg, {(int(a), int(b)) for a, b in anyr})
        for (ctx, layer), g in DD.groupby(["context", "layer"]):
            pair = g.groupby(["i", "j"])["diff"].mean()
            outdeg, nobs, indeg, anyset = node_info[ctx]
            ii = pair.index.get_level_values(0).values; jj = pair.index.get_level_values(1).values
            # node-level direction baselines (interpretive): TRRUST out-degree within V excluding the pair itself
            # (i -> j counted for i), and detection n_obs
            od_i = outdeg[ii] - 1; od_j = outdeg[jj]
            deg_dec = od_i != od_j; deg_fwd = od_i > od_j
            probe_fwd = pair.values > 0
            no_dec = nobs[ii] != nobs[jj]
            res_mask = deg_dec & ~deg_fwd                        # pairs where the degree baseline points j -> i
            in_i = indeg[ii]; in_j = indeg[jj] - 1; in_dec = in_i != in_j
            extra = dict(frac_indeg_j_gt_i=float((in_j > in_i)[in_dec].mean()), n_indeg_decisive=int(in_dec.sum()),
                         frac_outdeg_i_gt_j=float(deg_fwd[deg_dec].mean()), n_outdeg_decisive=int(deg_dec.sum()),
                         frac_nobs_i_gt_j=float((nobs[ii] > nobs[jj])[no_dec].mean()),
                         probe_frac_when_outdeg_favours_j=float(probe_fwd[res_mask].mean()) if res_mask.any() else np.nan,
                         n_outdeg_favours_j=int(res_mask.sum()),
                         probe_frac_when_outdeg_tied=float(probe_fwd[~deg_dec].mean()) if (~deg_dec).any() else np.nan,
                         n_outdeg_tied=int((~deg_dec).sum()),
                         spearman_probe_diff_vs_outdeg_diff=float(sst.spearmanr(pair.values, od_i - od_j)[0]))
            n = int(pair.size); k = int((pair > 0).sum()); nz = int((pair != 0).sum())
            bt = sst.binomtest(k, nz, 0.5) if nz else None
            ci = bt.proportion_ci(0.95, "exact") if bt else None
            per_rep = g.groupby(["repeat", "i", "j"])["diff"].mean().gt(0).groupby("repeat").mean()
            # strict subset: j -> i absent from every regulatory set (TRRUST and DoRothEA, any level, both files)
            strict = np.array([(int(b_), int(a_)) not in anyset for a_, b_ in zip(ii, jj)])
            ps = pair.values[strict]; ks_ = int((ps > 0).sum()); nzs = int((ps != 0).sum())
            bts = sst.binomtest(ks_, nzs, 0.5) if nzs else None; cis = bts.proportion_ci(0.95, "exact") if bts else None
            extra.update(n_strict=int(strict.sum()), frac_forward_strict=ks_ / nzs if nzs else np.nan,
                         strict_ci_lo=cis.low if cis else np.nan, strict_ci_hi=cis.high if cis else np.nan,
                         strict_binom_p=bts.pvalue if bts else np.nan)
            D.append(dict(context=ctx, model="M1", feature="concat", regime="C3", layer=int(layer), n_pairs=n,
                          n_nonzero=nz, n_forward_higher=k, frac_forward_higher=k / nz if nz else np.nan,
                          ci_lo=ci.low if ci else np.nan, ci_hi=ci.high if ci else np.nan,
                          binom_p_two_sided=bt.pvalue if bt else np.nan, mean_diff=pair.mean(),
                          median_diff=pair.median(), frac_by_repeat_min=per_rep.min(), frac_by_repeat_max=per_rep.max(),
                          n_repeats=int(g["repeat"].nunique()), **extra))
    pd.DataFrame(D).to_csv(H / "direction_test.csv", index=False)
    print("probe rows", len(P), "null rows", len(NR), "relational rows", len(R), "direction rows", len(D))


if __name__ == "__main__":
    main()
