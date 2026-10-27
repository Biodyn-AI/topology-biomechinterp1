"""E10 step 3 (task 2): the submitted 'cross-tissue' contrast kidney_A vs kidney_B re-evaluated as WITHIN-tissue
reproducibility (same TS kidney donor TSP2, different HVG axis: 4,803 vs 8,143 genes; different cell draws).

Matched set M = genes observed in both; labels recomputed on M and identical in both contexts (Lambert TF;
TRRUST TF = has >= 1 TRRUST target in M, target-only = TRRUST target in M and not a TRRUST TF).
Per layer and for each context:
  S2b_band_ownV   : SV2-SV7 of the observed-gene SVD of the context's own V (as in the submitted analysis,
                    but zero rows excluded), rows restricted to M
  S2b_band_commonM: SV2-SV7 of an SVD fitted on the M rows only (same genes in both contexts)
  S2b_full        : full-dimension V-centred embedding restricted to M (basis-free)
  S2_ownstrata    : detection-stratified Lambert AUROC (strata = deciles of the context's own log n_obs on M)
  S2_sharedstrata : strata = deciles of the mean of the two contexts' log n_obs (identical strata)
Paired gene bootstrap on M (B = 2,000, the same resample in both contexts): delta = A - B, and the basis-
instability test |delta_band| - |delta_full|. Direct basis diagnostics: principal-angle cosines between the
two contexts' SV2-SV7 directions in the shared 512-d model space (Vt rows) and canonical correlations of their
gene coordinates on M; within-context band SD with the SVD refitted on each gene resample (B = 200) vs with a
fixed basis. Legacy bridge: the submitted computation (curated gene set, SVD on the FULL matrix including the
never-observed zero rows - the documented artefact) is re-run so the submitted numbers are traceable.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e10_common as K
S = K.S
X = K.X

B = 2000; B_REFIT = 200
A, Bc = "kidney_A", "kidney_B"


def trrust_labels(names):
    _, pairs = X.load_trrust(); V = set(names)
    pv = [(a, b) for a, b in pairs if a in V and b in V]
    tfs = {a for a, _ in pv}; tg = {b for _, b in pv} - tfs
    return np.array([g in tfs for g in names]), np.array([g in tg for g in names])


def stat_set(EnA, EnB, tf, tg):
    return S.tf_auroc(EnA, tf, tg), S.tf_auroc(EnB, tf, tg)


def boot_tfauroc(En, tf, tg, W):
    lab = np.r_[tf, tg]; pos = np.r_[np.ones(tf.size, bool), np.zeros(tg.size, bool)]
    Sc = K.loo_scores_weighted(En, tf, W)
    return K.wauroc_var(Sc[lab], pos, W[lab])


def legacy_bridge():
    """Re-run of the submitted M3 computation (revision2_experiments/run_remaining.py, read-only logic):
    curated gene set (cycle1_edge_dataset.tsv), SVD of the FULL matrix including never-observed zero rows
    (documented zero-row artefact), TRRUST TF (genome-wide) vs TRRUST target-only (genome-wide), SV2-SV7, LOO.
    Then the same gene/label sets under CONVENTION3 (observed-gene SVD; full dimension)."""
    tsv = X.REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
    e = pd.read_csv(tsv, sep="\t"); named = set(e["source"].astype(str)) | set(e["target"].astype(str))
    t, _ = X.load_trrust(); tfs_all = set(t.TF.astype(str)); tgts_all = set(t.target.astype(str)) - tfs_all
    rows = []
    ea, eb = K.load_emb(A), K.load_emb(Bc)
    gA = X.h5ad_genes(X.CONTEXTS[A]["h5ad"]); gB = X.h5ad_genes(X.CONTEXTS[Bc]["h5ad"])
    nzA = np.linalg.norm(np.asarray(ea[0], float), axis=1) > 0; nzB = np.linalg.norm(np.asarray(eb[0], float), axis=1) > 0
    subA = np.array([i for i, g in enumerate(gA) if g in named and nzA[i]]); namesA = [gA[i] for i in subA]
    imm_tf = {g for g in namesA if g in tfs_all}
    subB = np.array([i for i, g in enumerate(gB) if g in set(namesA) and nzB[i]]); namesB = [gB[i] for i in subB]
    subA167 = np.array([gA.index(g) for g in namesB])
    sets = {"A_curated187": (ea, subA, namesA, "all_trrust_tf"), "B_curated167": (eb, subB, namesB, "imm_tf"),
            "A_matched167": (ea, subA167, namesB, "imm_tf")}
    cA, cB = K.load_cache(A), K.load_cache(Bc)
    for nm, (e_, sub, names, tfdef) in sets.items():
        tfp = np.array([k for k, g in enumerate(names) if (g in tfs_all if tfdef == "all_trrust_tf" else g in imm_tf)])
        tgp = np.array([k for k, g in enumerate(names) if g in tgts_all])
        c_ = cA if nm.startswith("A") else cB; pos = {int(j): k for k, j in enumerate(c_["obs_idx"])}
        rows_obs = np.array([pos[int(j)] for j in sub])
        for L in range(K.N_LAYERS):
            Xf = np.asarray(e_[L], float); Xf = Xf - Xf.mean(0); U, Sv, _ = np.linalg.svd(Xf, full_matrices=False)
            leg = S.tf_auroc(K.unit((U * Sv)[sub][:, K.BAND]), tfp, tgp)
            US, _, _ = S.svd_observed(e_, L, c_["obs_idx"])
            conv_band = S.tf_auroc(K.unit(US[rows_obs][:, K.BAND]), tfp, tgp)
            conv_full = S.tf_auroc(K.unit(K.centred(e_, L, c_["obs_idx"])[rows_obs]), tfp, tgp)
            rows += [dict(layer=L, statistic=f"legacy_bridge_{nm}_{k}", n_M=len(names), A=v, n_tf=int(tfp.size), n_tg=int(tgp.size))
                     for k, v in (("submitted_zero_row_SVD_band", leg), ("observed_SVD_band", conv_band), ("full_dim", conv_full))]
    return rows


def main():
    out = []; diag = []
    ca, cb = K.load_cache(A), K.load_cache(Bc); ea, eb = K.load_emb(A), K.load_emb(Bc)
    M = sorted(set(ca["genes"]) & set(cb["genes"]))
    pa = {g: i for i, g in enumerate(ca["genes"])}; pb = {g: i for i, g in enumerate(cb["genes"])}
    ia = np.array([pa[g] for g in M]); ib = np.array([pb[g] for g in M])
    lam = set(K.lambert_all()["HGNC symbol"])
    is_lam = np.array([g in lam for g in M]); lam_pos = np.flatnonzero(is_lam)
    is_trtf, is_trtg = trrust_labels(M); trtf = np.flatnonzero(is_trtf); trtg = np.flatnonzero(is_trtg)
    la = np.log(np.maximum(ca["nobs"][ia], 1)); lb = np.log(np.maximum(cb["nobs"][ib], 1))
    st_a, st_b = S.quantile_strata(la, 10), S.quantile_strata(lb, 10); st_s = S.quantile_strata((la + lb) / 2, 10)
    # ---------------- mechanism descriptors
    axA, axB = set(X.h5ad_genes(X.CONTEXTS[A]["h5ad"])), set(X.h5ad_genes(X.CONTEXTS[Bc]["h5ad"]))
    a2 = json.loads((K.R3 / "A2_axis_verification/results.json").read_text())["contexts"]
    from scipy.stats import spearmanr
    exB = np.load(cb["d"] / "expr_extraction.npy", mmap_mode="r"); exA = np.load(ca["d"] / "expr_extraction.npy", mmap_mode="r")
    gB_notA = np.array([g not in axA for g in cb["genes"]]); gA_notB = np.array([g not in axB for g in ca["genes"]])
    fracB = np.array([(r > 0)[gB_notA].sum() / max((r > 0).sum(), 1) for r in np.asarray(exB)])
    fracA = np.array([(r > 0)[gA_notB].sum() / max((r > 0).sum(), 1) for r in np.asarray(exA)])
    mech = dict(axis_A=len(axA), axis_B=len(axB), axis_shared=len(axA & axB), V_A=len(ca["genes"]), V_B=len(cb["genes"]),
                M=len(M), lambert_tf_M=int(is_lam.sum()), trrust_tf_M=int(trtf.size), trrust_target_only_M=int(trtg.size),
                lambert_tf_V_A=int(ca["nl"]["is_tf_lambert"].sum()), lambert_tf_V_B=int(cb["nl"]["is_tf_lambert"].sum()),
                cells_A=a2[A]["n_cells"], cells_B=a2[Bc]["n_cells"],
                cells_truncated_A=a2[A]["cells_with_nnz_gt_max_genes"], cells_truncated_B=a2[Bc]["cells_with_nnz_gt_max_genes"],
                nnz_median_A=a2[A]["nnz_median"], nnz_median_B=a2[Bc]["nnz_median"], max_genes=512,
                spearman_nobs_on_M=float(spearmanr(la, lb)[0]),
                median_nobs_lambert_M_A=float(np.median(ca["nobs"][ia][is_lam])), median_nobs_lambert_M_B=float(np.median(cb["nobs"][ib][is_lam])),
                median_nobs_trrust_tf_M_A=float(np.median(ca["nobs"][ia][is_trtf])), median_nobs_trrust_tf_M_B=float(np.median(cb["nobs"][ib][is_trtf])),
                median_nobs_nonTF_M_A=float(np.median(ca["nobs"][ia][~is_lam])), median_nobs_nonTF_M_B=float(np.median(cb["nobs"][ib][~is_lam])),
                mean_frac_cell_inputs_outside_other_axis_B=float(fracB.mean()), mean_frac_cell_inputs_outside_other_axis_A=float(fracA.mean()),
                drawn_cells_shared=12)
    print(json.dumps(mech, indent=1))
    rng = np.random.default_rng(K.SEED + 30)
    W = rng.multinomial(len(M), np.full(len(M), 1 / len(M)), size=B).T.astype(float)     # (|M|, B), same for A and B
    keep = {}
    for L in range(K.N_LAYERS):
        EA = K.unit(K.centred(ea, L, ca["obs_idx"])[ia]); EB = K.unit(K.centred(eb, L, cb["obs_idx"])[ib])
        USa, _, Vta = S.svd_observed(ea, L, ca["obs_idx"]); USb, _, Vtb = S.svd_observed(eb, L, cb["obs_idx"])
        BA, BB = K.unit(USa[ia][:, K.BAND]), K.unit(USb[ib][:, K.BAND])
        USma, _, Vtma = S.svd_observed(ea, L, ca["obs_idx"][ia]); USmb, _, Vtmb = S.svd_observed(eb, L, cb["obs_idx"][ib])
        CA, CB = K.unit(USma[:, K.BAND]), K.unit(USmb[:, K.BAND])
        specs = {
            "S2b_band_ownV": (BA, BB, "trrust", None), "S2b_band_commonM": (CA, CB, "trrust", None),
            "S2b_full": (EA, EB, "trrust", None), "S2_ownstrata": (EA, EB, "lambert", (st_a, st_b)),
            "S2_sharedstrata": (EA, EB, "lambert", (st_s, st_s)),
        }
        bootd = {}
        for nm, (Ea, Eb, lab, strata) in specs.items():
            if lab == "trrust":
                va, vb = stat_set(Ea, Eb, trtf, trtg)
                ba, bb = boot_tfauroc(Ea, trtf, trtg, W), boot_tfauroc(Eb, trtf, trtg, W)
            else:
                sa, sb = S.tf_scores_loo(Ea, lam_pos), S.tf_scores_loo(Eb, lam_pos)
                va, vb = S.stratified_auroc(sa, is_lam, strata[0]), S.stratified_auroc(sb, is_lam, strata[1])
                ba = K.wstrat_auroc_var(K.loo_scores_weighted(Ea, lam_pos, W), is_lam, strata[0], W)
                bb = K.wstrat_auroc_var(K.loo_scores_weighted(Eb, lam_pos, W), is_lam, strata[1], W)
            d = ba - bb; bootd[nm] = (va, vb, ba, bb); keep.setdefault(nm, []).append((va, vb, ba, bb))
            out.append(dict(layer=L, statistic=nm, n_M=len(M), A=va, B=vb, delta=va - vb,
                            delta_lo=float(np.quantile(d, 0.025)), delta_hi=float(np.quantile(d, 0.975)),
                            se_A=float(ba.std(ddof=1)), se_B=float(bb.std(ddof=1)),
                            A_lo=float(np.quantile(ba, .025)), A_hi=float(np.quantile(ba, .975)),
                            B_lo=float(np.quantile(bb, .025)), B_hi=float(np.quantile(bb, .975)),
                            corr_boot_A_B=float(np.corrcoef(ba, bb)[0, 1])))
        # basis-instability test: |delta_band| - |delta_full|
        for band in ["S2b_band_ownV", "S2b_band_commonM"]:
            va, vb, ba, bb = bootd[band]; fa, fb, fba, fbb = bootd["S2b_full"]
            obs = abs(va - vb) - abs(fa - fb); bd = np.abs(ba - bb) - np.abs(fba - fbb)
            dev = bd - bd.mean()
            out.append(dict(layer=L, statistic=f"absdelta_{band}_minus_absdelta_full", n_M=len(M), delta=obs,
                            delta_lo=float(obs - np.quantile(dev, .975)), delta_hi=float(obs - np.quantile(dev, .025))))
        # direct basis diagnostics
        cos_model = np.linalg.svd(Vta[K.BAND] @ Vtb[K.BAND].T, compute_uv=False)
        cos_modelM = np.linalg.svd(Vtma[K.BAND] @ Vtmb[K.BAND].T, compute_uv=False)
        Ga = USa[ia][:, K.BAND] - USa[ia][:, K.BAND].mean(0); Gb = USb[ib][:, K.BAND] - USb[ib][:, K.BAND].mean(0)
        qa, _ = np.linalg.qr(Ga); qb, _ = np.linalg.qr(Gb); cc = np.linalg.svd(qa.T @ qb, compute_uv=False)
        best = np.abs(Vta[K.BAND] @ Vtb[:20].T).max(1)
        # within-context band variability: SVD refitted on gene resamples vs fixed basis
        rr = np.random.default_rng(K.SEED + 31 + L); ref = {}
        for nm_c, e_, c_ in ((A, ea, ca), (Bc, eb, cb)):
            n = len(c_["genes"]); trt = np.flatnonzero(c_["nl"]["is_trrust_tf"]); trg = np.flatnonzero(c_["nl"]["is_trrust_target_only"])
            Xl = np.asarray(e_[L], dtype=np.float64)[c_["obs_idx"]]
            refit, fixed = [], []
            US0, _, Vt0 = S.svd_observed(e_, L, c_["obs_idx"]); Ef = K.unit(US0[:, K.BAND])
            for r in range(B_REFIT):
                idx = rr.integers(0, n, n)
                Xr = Xl[idx] - Xl[idx].mean(0); _, _, Vt = np.linalg.svd(Xr, full_matrices=False)
                Eb_ = K.unit((Xl - Xl.mean(0)) @ Vt[K.BAND].T)        # resample-fitted basis, all genes projected
                w = np.bincount(idx, minlength=n).astype(float)[:, None]
                refit.append(boot_tfauroc(Eb_, trt, trg, w)[0]); fixed.append(boot_tfauroc(Ef, trt, trg, w)[0])
            ref[nm_c] = (float(np.std(refit, ddof=1)), float(np.std(fixed, ddof=1)))
        diag.append(dict(layer=L, band_principal_cos_ownV=";".join(f"{x:.3f}" for x in cos_model),
                         band_mean_principal_cos_ownV=float(cos_model.mean()),
                         band_mean_principal_cos_commonM=float(cos_modelM.mean()),
                         band_gene_canonical_corr_mean=float(cc.mean()), band_gene_canonical_corr_min=float(cc.min()),
                         band_best_match_cos_in_top20=";".join(f"{x:.3f}" for x in best),
                         sd_band_refit_A=ref[A][0], sd_band_fixedbasis_A=ref[A][1],
                         sd_band_refit_B=ref[Bc][0], sd_band_fixedbasis_B=ref[Bc][1]))
        print(L, {k: round(v[0] - v[1], 4) for k, v in bootd.items()}, diag[-1]["band_mean_principal_cos_ownV"], flush=True)
    # layer-mean rows: mean over layers of the point values and of the (shared) bootstrap replicates
    lm = {}
    for nm, lst in keep.items():
        va = np.mean([x[0] for x in lst]); vb = np.mean([x[1] for x in lst])
        ba = np.mean([x[2] for x in lst], 0); bb = np.mean([x[3] for x in lst], 0); d = ba - bb; lm[nm] = (va, vb, ba, bb)
        out.append(dict(layer="mean", statistic=nm, n_M=len(M), A=va, B=vb, delta=va - vb,
                        delta_lo=float(np.quantile(d, .025)), delta_hi=float(np.quantile(d, .975)),
                        se_A=float(ba.std(ddof=1)), se_B=float(bb.std(ddof=1)),
                        A_lo=float(np.quantile(ba, .025)), A_hi=float(np.quantile(ba, .975)),
                        B_lo=float(np.quantile(bb, .025)), B_hi=float(np.quantile(bb, .975)),
                        corr_boot_A_B=float(np.corrcoef(ba, bb)[0, 1])))
    for band in ["S2b_band_ownV", "S2b_band_commonM"]:
        va, vb, ba, bb = lm[band]; fa, fb, fba, fbb = lm["S2b_full"]
        obs = abs(va - vb) - abs(fa - fb); bd = np.abs(ba - bb) - np.abs(fba - fbb); dev = bd - bd.mean()
        out.append(dict(layer="mean", statistic=f"absdelta_{band}_minus_absdelta_full", n_M=len(M), delta=obs,
                        delta_lo=float(obs - np.quantile(dev, .975)), delta_hi=float(obs - np.quantile(dev, .025))))
        # layer-profile agreement: Pearson of the 12-layer profiles A vs B
        pa_ = np.array([x[0] for x in keep[band]]); pb_ = np.array([x[1] for x in keep[band]])
        fa_ = np.array([x[0] for x in keep["S2b_full"]]); fb_ = np.array([x[1] for x in keep["S2b_full"]])
        out.append(dict(layer="profile", statistic=f"mean_absdelta_over_layers_{band}", n_M=len(M),
                        delta=float(np.mean(np.abs(pa_ - pb_)))))
    fa_ = np.array([x[0] for x in keep["S2b_full"]]); fb_ = np.array([x[1] for x in keep["S2b_full"]])
    out.append(dict(layer="profile", statistic="mean_absdelta_over_layers_S2b_full", n_M=len(M),
                    delta=float(np.mean(np.abs(fa_ - fb_)))))
    out += legacy_bridge()
    df = pd.DataFrame(out)
    df.to_csv(K.HERE / "kidneyA_vs_B.csv", index=False)
    df.to_csv(K.HERE / "partial" / "kidney_layers.csv", index=False)
    pd.DataFrame(diag).to_csv(K.HERE / "partial" / "kidney_basis_diag.csv", index=False)
    (K.HERE / "partial" / "kidney_mech.json").write_text(json.dumps(mech, indent=1))


if __name__ == "__main__":
    main()
