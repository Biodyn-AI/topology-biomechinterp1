"""E10 step 5 (task 4): lineage-restricted TFs and where they are recognised.

Selection (data-driven, fixed before any embedding statistic is computed for these sets):
  For each of the four datasets (TS immune, TS kidney [kidney_B file], TS lung, Krasnow lung) 3,000 cells
  (seeded) are read from raw counts (all 60,606 / 53,514 genes), normalised to counts per 10k, and for every
  Lambert TF (matched by Ensembl ID) the mean CP10k and the fraction of cells with count > 0 are computed.
  Tissue-level expression: immune = TS immune, kidney = TS kidney, lung = min(TS lung, Krasnow) (a lung TF must
  be expressed in BOTH lung studies). Specificity_T = log2((m_T + 0.1) / (max over the other datasets + 0.1)),
  where for lung the 'other' datasets are immune and kidney and for immune / kidney the others include both lung
  datasets. Eligible: fraction expressing >= 5 % in T (in both lung datasets for lung) and Specificity_T >= 1.
  The top 30 per tissue by Specificity_T form the lineage set.
Evaluation in every context and layer: positives = lineage-set TFs observed in V (>= 5 required), negatives =
non-Lambert genes in V, score = LOO mean cosine to the observed lineage set (stats3.tf_scores_loo, full dim,
V-centred), detection-stratified AUROC (cache det_decile); Poisson gene bootstrap (B = 500) for CIs. Reference:
50 random Lambert-TF subsets of the same observed size in the same context (layer-mean).
Output: lineage_tfs.csv (selection table + evaluation rows, column 'part').
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, h5py
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e10_common as K


def cci(v, bt):
    """Centred bootstrap interval: v - (quantile - mean); removes the small downward bias of the
    leave-gene-out bootstrap replicates (see SUMMARY)."""
    bt = np.asarray(bt, float); m = np.nanmean(bt)
    return float(v - (np.nanquantile(bt, .975) - m)), float(v - (np.nanquantile(bt, .025) - m))
S = K.S; X = K.X
N_CELLS, TOP, B, R_RAND = 3000, 30, 500, 50
DATASETS = {"TS_immune": X.H5_IMMUNE, "TS_kidney": X.H5_KIDNEY_B, "TS_lung": X.H5_LUNG, "Krasnow_lung": X.H5_EXT_LUNG}


def raw_profile(path, ens_keep, rng):
    with h5py.File(path, "r") as f:
        r = f["raw"]; ens = np.array([x.decode() for x in r["var"]["_index"][:]])
        col = {e: i for i, e in enumerate(ens)}
        cols = np.array([col.get(e, -1) for e in ens_keep])
        Xg = r["X"]; indptr = Xg["indptr"][:]; nc = indptr.size - 1
        cells = np.sort(rng.choice(nc, min(N_CELLS, nc), replace=False))
        pos = np.full(ens.size, -1); pos[cols[cols >= 0]] = np.flatnonzero(cols >= 0)
        M = np.zeros((cells.size, len(ens_keep)), np.float64); lib = np.zeros(cells.size)
        data, ind = Xg["data"], Xg["indices"]
        for k, cidx in enumerate(cells):
            a, b = indptr[cidx], indptr[cidx + 1]; v = data[a:b].astype(float); ii = ind[a:b]
            lib[k] = v.sum(); p = pos[ii]; m = p >= 0; M[k, p[m]] = v[m]
    cp = M / np.maximum(lib, 1)[:, None] * 1e4
    present = cols >= 0
    return cp.mean(0), (M > 0).mean(0), present


def select():
    lam = K.lambert_all().drop_duplicates("Ensembl ID")
    ens = lam["Ensembl ID"].tolist(); sym = lam["HGNC symbol"].tolist()
    prof = {}
    for i, (nm, p) in enumerate(DATASETS.items()):
        m, fr, pres = raw_profile(p, ens, np.random.default_rng(K.SEED + 50 + i)); prof[nm] = (m, fr, pres)
        print(nm, "TFs present", int(pres.sum()), flush=True)
    df = pd.DataFrame({"ensembl": ens, "symbol": sym})
    for nm, (m, fr, pres) in prof.items():
        df[f"cp10k_{nm}"] = np.where(pres, m, 0.0); df[f"frac_{nm}"] = np.where(pres, fr, 0.0); df[f"present_{nm}"] = pres
    e = 0.1
    m_imm, m_kid = df.cp10k_TS_immune, df.cp10k_TS_kidney
    m_lung = np.minimum(df.cp10k_TS_lung, df.cp10k_Krasnow_lung)
    df["spec_immune"] = np.log2((m_imm + e) / (np.maximum.reduce([m_kid, df.cp10k_TS_lung, df.cp10k_Krasnow_lung]) + e))
    df["spec_kidney"] = np.log2((m_kid + e) / (np.maximum.reduce([m_imm, df.cp10k_TS_lung, df.cp10k_Krasnow_lung]) + e))
    df["spec_lung"] = np.log2((m_lung + e) / (np.maximum(m_imm, m_kid) + e))
    df["elig_immune"] = df.frac_TS_immune >= 0.05
    df["elig_kidney"] = df.frac_TS_kidney >= 0.05
    df["elig_lung"] = (df.frac_TS_lung >= 0.05) & (df.frac_Krasnow_lung >= 0.05)
    sets = {}
    for t in ("immune", "kidney", "lung"):
        d = df[df[f"elig_{t}"] & (df[f"spec_{t}"] >= 1)].sort_values(f"spec_{t}", ascending=False).head(TOP)
        sets[t] = d.symbol.tolist(); df[f"in_set_{t}"] = df.symbol.isin(sets[t])
        df.loc[d.index, f"rank_{t}"] = np.arange(1, len(d) + 1)
        print(t, len(sets[t]), sets[t], flush=True)
    return df, sets


def evaluate(sets):
    allrows = []
    for ci, ctx in enumerate(K.CTXS):
        ck = K.HERE / "partial" / f"lineage_eval_{ctx}.csv"
        if ck.exists():
            allrows += pd.read_csv(ck).to_dict("records"); continue
        rows = []
        c = K.load_cache(ctx); emb = K.load_emb(ctx); genes = c["genes"]; pos = {g: i for i, g in enumerate(genes)}
        nl = c["nl"]; is_lam = nl["is_tf_lambert"].astype(bool); strata = nl["det_decile"]; n = len(genes)
        rng = np.random.default_rng(K.SEED + 60 + ci)
        W = rng.poisson(1.0, (n, B)).astype(float)
        Ens = [K.unit(K.centred(emb, L, c["obs_idx"])) for L in range(K.N_LAYERS)]
        lam_idx = np.flatnonzero(is_lam)
        for t, gl in sets.items():
            P = np.array(sorted(pos[g] for g in gl if g in pos))
            base = dict(part="eval", context=ctx, tissue_set=t, n_set=len(gl), n_observed=int(P.size),
                        median_nobs_set=float(np.median(c["nobs"][P])) if P.size else np.nan,
                        median_nobs_nonTF=float(np.median(c["nobs"][~is_lam])),
                        home=bool(X.CONTEXTS[ctx]["tissue"] == t))
            if P.size < 5:
                rows.append(dict(base, layer="mean", auroc=np.nan, note="fewer than 5 observed")); continue
            rowsel = np.zeros(n, bool); rowsel[P] = True; rowsel |= ~is_lam
            isp = np.zeros(n, bool); isp[P] = True
            vals, boots = [], []
            for L, En in enumerate(Ens):
                sc = S.tf_scores_loo(En, P)
                v = S.stratified_auroc(sc[rowsel], isp[rowsel], strata[rowsel])
                bt = K.wstrat_auroc_var(K.loo_scores_weighted(En, P, W)[rowsel], isp[rowsel], strata[rowsel], W[rowsel])
                vals.append(v); boots.append(bt)
                lo_, hi_ = cci(v, bt)
                rows.append(dict(base, layer=L, auroc=v, lo=lo_, hi=hi_))
            bm = np.mean(boots, 0)
            # reference: random Lambert subsets of the same observed size
            rv = []
            for r in range(R_RAND):
                Q = np.sort(rng.choice(lam_idx, P.size, replace=False))
                sel = np.zeros(n, bool); sel[Q] = True; rs = sel | ~is_lam
                rv.append(np.mean([S.stratified_auroc(S.tf_scores_loo(En, Q)[rs], sel[rs], strata[rs]) for En in Ens]))
            lo_, hi_ = cci(float(np.mean(vals)), bm)
            rows.append(dict(base, layer="mean", auroc=float(np.mean(vals)), lo=lo_, hi=hi_,
                             random_subset_mean=float(np.mean(rv)), random_subset_q025=float(np.quantile(rv, .025)),
                             random_subset_q975=float(np.quantile(rv, .975)),
                             frac_random_ge_obs=float(np.mean(np.array(rv) >= np.mean(vals)))))
            print(ctx, t, P.size, round(float(np.mean(vals)), 3), round(float(np.mean(rv)), 3), flush=True)
        pd.DataFrame(rows).to_csv(ck, index=False); allrows += rows
    return allrows


def main():
    fsel = K.HERE / "partial" / "lineage_selection_full.csv"
    if fsel.exists():
        df = pd.read_csv(fsel); sets = {t: df[df[f"in_set_{t}"]].sort_values(f"rank_{t}").symbol.tolist() for t in ("immune", "kidney", "lung")}
    else:
        df, sets = select(); df.to_csv(fsel, index=False)
    rows = evaluate(sets)
    sel = df[df[[f"in_set_{t}" for t in sets]].any(axis=1)].copy(); sel.insert(0, "part", "selection")
    ev = pd.DataFrame(rows)
    # S2 (all Lambert TFs, from e10_stats) vs the TFs' median n_obs across the 16 contexts
    s2 = []
    for ctx in K.CTXS:
        z = np.load(K.HERE / "partial" / f"stats_{ctx}.npz"); c = K.load_cache(ctx)
        is_lam = c["nl"]["is_tf_lambert"].astype(bool)
        s2.append(dict(part="S2_vs_nobs", context=ctx, layer="mean", auroc=float(z["point"][1].mean()),
                       median_nobs_set=float(np.median(c["nobs"][is_lam])), n_observed=int(is_lam.sum()),
                       median_nobs_nonTF=float(np.median(c["nobs"][~is_lam]))))
    s2 = pd.DataFrame(s2)
    rho = spearmanr(s2.auroc, s2.median_nobs_set)
    m = ev[(ev.layer.astype(str) == "mean") & ev.auroc.notna()]
    rho2 = spearmanr(m.auroc, m.median_nobs_set)
    summ = pd.DataFrame([dict(part="summary", test="spearman_S2_layermean_vs_median_nobs_lambertTF_16ctx", auroc=rho[0], p=rho[1], n_observed=len(s2)),
                         dict(part="summary", test="spearman_lineageset_auroc_vs_median_nobs_set", auroc=rho2[0], p=rho2[1], n_observed=len(m))])
    pd.concat([sel, ev, s2, summ], ignore_index=True).to_csv(K.HERE / "lineage_tfs.csv", index=False)
    print(summ.to_string())


if __name__ == "__main__":
    main()
