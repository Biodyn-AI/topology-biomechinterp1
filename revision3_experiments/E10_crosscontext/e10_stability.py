"""E10 step 7 (task 6): basis-free representational stability and power-matched subsampling.

(a) For the 20 context pairs that make up the variance-decomposition groups (resample, donor, cell type,
    HVG axis, checkpoint, tissue), per layer, on matched genes M (observed in both): k = 20 nearest neighbours
    of every gene by cosine of the V-centred full-dimension embedding, computed within M in each context;
    per-gene Jaccard of the two neighbour sets. Reported: mean Jaccard of Lambert TFs and non-TFs, and the
    detection-adjusted TF - nonTF difference (within deciles of the pair-mean log n_obs, weighted by TF count).
    RSA: Spearman correlation of the upper-triangle cosine matrices for (i) TF x TF, (ii) size-matched random
    non-TF sets (mean of 10 draws), (iii) a random 1,500-gene subset of M.
(b) Power-matched subsampling: lung, lung_external and kidney_B are subsampled (R = 100, seeded) to immune's
    numbers: |V| genes with exactly immune's Lambert-TF count; STRING edges inside the subsample further
    subsampled to immune's edge count; the H3 benchmark subsampled to immune's number of H3 TFs (all their
    pairs, full-V embedding). Within each subsample, centring, detection deciles, observed-gene SVD and one
    Maslov-Sneppen rewiring are recomputed. S1, S2, S2b_full, S2b_band, S3 (layer means) are compared with
    immune's own values.
Checkpoints: partial/stab_<pair>.csv, partial/powermatch_<ctx>.npz. Output: stability.csv
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e10_common as K
S = K.S
KNN = 20; R_PM = 100
GROUPS = {
    "resample_immune": ["immune", "immune_seed43", "immune_seed44"], "resample_kidney": ["kidney_B", "kidney_ckpt_whole_human"],
    "donor": ["immune_TSP14", "immune_TSP21", "immune_TSP25"], "celltype": ["immune_Bcell", "immune_CD4T", "immune_CD8T"],
    "hvg_axis": ["kidney_A", "kidney_B"], "checkpoint": ["kidney_ckpt_whole_human", "kidney_ckpt_kidney", "kidney_ckpt_brain"],
    "tissue": ["immune", "kidney_B", "lung", "lung_external"],
}


def pairs():
    out = []
    for g, cs in GROUPS.items():
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)): out.append((g, cs[i], cs[j]))
    return out


def knn(En, k):
    Cm = En @ En.T; np.fill_diagonal(Cm, -np.inf)
    return np.argpartition(-Cm, k, axis=1)[:, :k]


def spearman_vec(a, b):
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def stab_pair(g, a, b):
    f = K.HERE / "partial" / f"stab_{a}__{b}.csv"
    if f.exists(): return pd.read_csv(f)
    ca, cb = K.load_cache(a), K.load_cache(b); ea, eb = K.load_emb(a), K.load_emb(b)
    M = sorted(set(ca["genes"]) & set(cb["genes"]))
    pa = {x: i for i, x in enumerate(ca["genes"])}; pb = {x: i for i, x in enumerate(cb["genes"])}
    ia = np.array([pa[x] for x in M]); ib = np.array([pb[x] for x in M])
    lam = ca["nl"]["is_tf_lambert"][ia].astype(bool)
    ln = (np.log(np.maximum(ca["nobs"][ia], 1)) + np.log(np.maximum(cb["nobs"][ib], 1))) / 2
    dec = S.quantile_strata(ln, 10)
    rng = np.random.default_rng(K.SEED + 70)
    tfi = np.flatnonzero(lam); ntf = np.flatnonzero(~lam)
    sub1500 = np.sort(rng.choice(len(M), min(1500, len(M)), replace=False))
    rand_sets = [np.sort(rng.choice(ntf, tfi.size, replace=False)) for _ in range(10)]
    iu_tf = np.triu_indices(tfi.size, 1)
    rows = []
    for L in range(K.N_LAYERS):
        A = K.unit(K.centred(ea, L, ca["obs_idx"])[ia]).astype(np.float32)
        Bm = K.unit(K.centred(eb, L, cb["obs_idx"])[ib]).astype(np.float32)
        na, nb = knn(A, KNN), knn(Bm, KNN)
        inter = np.array([np.intersect1d(x, y, assume_unique=True).size for x, y in zip(na, nb)])
        jac = inter / (2 * KNN - inter)
        adj = []; wts = []
        for d in np.unique(dec):
            m = dec == d
            if (m & lam).sum() and (m & ~lam).sum():
                adj.append(jac[m & lam].mean() - jac[m & ~lam].mean()); wts.append((m & lam).sum())
        def rsa(idx):
            iu = np.triu_indices(idx.size, 1)
            return spearman_vec((A[idx] @ A[idx].T)[iu], (Bm[idx] @ Bm[idx].T)[iu])
        rows.append(dict(part="knn_rsa", group=g, ctx_a=a, ctx_b=b, layer=L, n_M=len(M), n_TF=int(tfi.size),
                         jaccard_TF=float(jac[lam].mean()), jaccard_nonTF=float(jac[~lam].mean()),
                         jaccard_all=float(jac.mean()), jaccard_chance=KNN / (len(M) - 1),
                         jaccard_TF_minus_nonTF_detadj=float(np.average(adj, weights=wts)),
                         rsa_TF=rsa(tfi), rsa_nonTF_sizematched=float(np.mean([rsa(r) for r in rand_sets])),
                         rsa_random1500=rsa(sub1500)))
    df = pd.DataFrame(rows); df.to_csv(f, index=False)
    print(g, a, b, df[["jaccard_TF", "jaccard_nonTF", "rsa_TF", "rsa_nonTF_sizematched"]].mean().round(3).to_dict(), flush=True)
    return df


def stats_on(emb, obs_idx, nl, nobs, edges, sp, Wsel=None, seed=0):
    """Layer-mean S1, S2, S2b_full, S2b_band, S3 for a (sub)universe obs_idx (rows already restricted)."""
    is_tf = nl["is_tf_lambert"].astype(bool); tf = np.flatnonzero(is_tf)
    strata = S.quantile_strata(np.log(np.maximum(nobs, 1)), 10)
    trtf = np.flatnonzero(nl["is_trrust_tf"]); trtg = np.flatnonzero(nl["is_trrust_target_only"])
    epos, eneg, _ = K.rewired_nonpairs(edges, obs_idx.size, seed=seed)
    pa = np.r_[epos[:, 0], eneg[:, 0]]; pb = np.r_[epos[:, 1], eneg[:, 1]]
    yp = np.r_[np.ones(len(epos), bool), np.zeros(len(eneg), bool)]
    out = np.zeros((5, K.N_LAYERS))
    for L in range(K.N_LAYERS):
        En = K.unit(K.centred(emb, L, obs_idx))
        s1 = (En[pa] * En[pb]).sum(1); out[0, L] = S.auroc(s1[yp], s1[~yp])
        out[1, L] = S.stratified_auroc(S.tf_scores_loo(En, tf), is_tf, strata)
        out[2, L] = S.tf_auroc(En, trtf, trtg)
        US, _, _ = S.svd_observed(emb, L, obs_idx); out[3, L] = S.tf_auroc(K.unit(US[:, K.BAND]), trtf, trtg)
        if sp is not None:
            spo, s3_layers, m = sp
            out[4, L] = spo.auroc(s3_layers[L][m])
    return out


def powermatch(ctx, ref):
    f = K.HERE / "partial" / f"powermatch_{ctx}.npz"
    if f.exists(): return dict(np.load(f))
    c = K.load_cache(ctx); emb = K.load_emb(ctx); nl = c["nl"]; n = len(c["genes"])
    is_tf = nl["is_tf_lambert"].astype(bool); tfi = np.flatnonzero(is_tf); ntf = np.flatnonzero(~is_tf)
    edges = np.load(c["d"] / "edges_string.npy")
    z = np.load(c["d"] / "h3_pairs.npz"); sp_full = S.StratifiedPairs(z["tf"], z["tg"], z["y"])
    s3_layers = []
    for L in range(K.N_LAYERS):
        Ef = K.unit(K.centred(emb, L, c["obs_idx"])); s3_layers.append((Ef[sp_full.tf] * Ef[sp_full.tg]).sum(1))
    rng = np.random.default_rng(K.SEED + 80 + K.CTXS.index(ctx))
    res = np.zeros((R_PM, 5, K.N_LAYERS)); meta = []
    for r in range(R_PM):
        sel = np.sort(np.r_[rng.choice(tfi, ref["n_tf"], replace=False), rng.choice(ntf, ref["n_V"] - ref["n_tf"], replace=False)])
        loc = np.full(n, -1); loc[sel] = np.arange(sel.size)
        e = edges[(loc[edges[:, 0]] >= 0) & (loc[edges[:, 1]] >= 0)]; e = np.c_[loc[e[:, 0]], loc[e[:, 1]]]
        if len(e) > ref["n_edges"]: e = e[np.sort(rng.choice(len(e), ref["n_edges"], replace=False))]
        nl_sub = {k: v[sel] for k, v in nl.items()}
        tf_keep = rng.choice(sp_full.tfs, min(ref["h3_tfs"], sp_full.tfs.size), replace=False)
        m = np.isin(sp_full.tf, tf_keep); spo = S.StratifiedPairs(sp_full.tf[m], sp_full.tg[m], sp_full.y[m])
        assert spo.tf.size == m.sum()
        res[r] = stats_on(emb, c["obs_idx"][sel], nl_sub, c["nobs"][sel], e, (spo, s3_layers, m), seed=int(rng.integers(1 << 30)))
        meta.append((len(e), int(spo.y.sum())))
        if r % 10 == 0: print(ctx, r, res[r].mean(1).round(3), flush=True)
    meta = np.array(meta)
    np.savez(f, res=res, n_edges=meta[:, 0], h3_pos=meta[:, 1])
    return dict(res=res, n_edges=meta[:, 0], h3_pos=meta[:, 1])


def main():
    allrows = []
    for g, a, b in pairs():
        allrows.append(stab_pair(g, a, b))
    st = pd.concat(allrows, ignore_index=True)
    lm = st.groupby(["part", "group", "ctx_a", "ctx_b"], as_index=False).mean(numeric_only=True); lm["layer"] = "mean"
    st = pd.concat([st, lm], ignore_index=True)
    # power matching
    ri = json.loads(str(np.load(K.HERE / "partial" / "stats_immune.npz")["meta"]))
    zi = np.load(K.HERE / "partial" / "stats_immune.npz")["point"]
    ref = dict(n_V=ri["n_V"], n_tf=ri["n_lambert_tf"], n_edges=ri["n_string_edges"], h3_tfs=ri["h3_tfs"])
    names = ["S1", "S2", "S2b_full", "S2b_band", "S3"]
    pm = []
    for ctx in ["lung", "lung_external", "kidney_B"]:
        d = powermatch(ctx, ref); full = np.load(K.HERE / "partial" / f"stats_{ctx}.npz")["point"]
        lmv = d["res"].mean(2)                                    # (R, 5) layer means
        for si, nm in enumerate(names):
            v = lmv[:, si]; im = zi[si].mean(); fu = full[si].mean()
            pm.append(dict(part="powermatch", group="powermatch", ctx_a=ctx, ctx_b="immune", layer="mean", statistic=nm,
                           immune_value=im, full_value=fu, matched_mean=float(v.mean()), matched_q025=float(np.quantile(v, .025)),
                           matched_q975=float(np.quantile(v, .975)), frac_matched_gt_immune=float(np.mean(v > im)),
                           shift_full_minus_immune=fu - im, shift_matched_minus_immune=float(v.mean() - im),
                           R=R_PM, n_V=ref["n_V"], n_tf=ref["n_tf"], n_edges_mean=float(d["n_edges"].mean()),
                           h3_tfs=ref["h3_tfs"], h3_pos_mean=float(d["h3_pos"].mean())))
    out = pd.concat([st, pd.DataFrame(pm)], ignore_index=True)
    out.to_csv(K.HERE / "stability.csv", index=False)
    print(pd.DataFrame(pm)[["ctx_a", "statistic", "immune_value", "full_value", "matched_mean", "matched_q025", "matched_q975", "frac_matched_gt_immune"]].round(3).to_string())


if __name__ == "__main__":
    main()
