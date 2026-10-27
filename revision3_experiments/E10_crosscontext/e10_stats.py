"""E10 step 1: per (context, layer) point estimates of S1, S2, S2b (full / SV2-7 band), S3 and a JOINT
Poisson gene bootstrap shared across all 16 contexts (same multiplicity for a gene symbol everywhere, so
the covariance between contexts that share genes is retained).

S2   = H2 primary statistic: detection-stratified (deciles of log n_obs, cache det_decile) AUROC of the
       leave-one-out mean cosine to the Lambert TF set, Lambert TF vs every other gene in V, full dimension,
       V-centred embedding.                                    -> stats3.tf_scores_loo + stats3.stratified_auroc
S2b_full = legacy TRRUST TF vs TRRUST target-only AUROC of the LOO mean cosine to the TRRUST-TF set, full dim.
S2b_band = same on SV2-SV7 coordinates (U_k S_k, columns 1..6) of the observed-gene SVD (stats3.svd_observed).
S3   = H3 primary statistic: within-TF stratified AUROC of signed cosine (full dim) on cache h3_pairs.
S1   = basis-free PPI pair-cosine AUROC: STRING-physical>=700 edges within V vs ONE Maslov-Sneppen rewired
       non-pair set (burn-in 10m, +m swaps, seed fixed per context; rewired pairs equal to a real edge dropped).

Bootstrap (B = 1000, Poisson(1) multiplicities per gene symbol; S3 uses independent TF-role and target-role
multiplicities = two-way pigeonhole; S1 pair weight = w_i w_j). Checkpoint: partial/stats_<ctx>.npz.
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import e10_common as K
S = K.S

B = 1000
OUT = K.HERE / "partial"
STATS = ["S1", "S2", "S2b_full", "S2b_band", "S3"]


def universe():
    U = set()
    for c in K.CTXS: U |= set(K.load_cache(c)["genes"])
    return sorted(U)


def poisson_weights(U, B=B):
    rng = np.random.default_rng(K.SEED + 10)
    Wg = rng.poisson(1.0, (len(U), B)).astype(np.uint8)
    Wtf = rng.poisson(1.0, (len(U), B)).astype(np.uint8)
    Wtg = rng.poisson(1.0, (len(U), B)).astype(np.uint8)
    return Wg, Wtf, Wtg


def run_ctx(ctx, U, Wg_all, Wtf_all, Wtg_all, log):
    f = OUT / f"stats_{ctx}.npz"
    if f.exists():
        log(f"{ctx}: checkpoint exists, skip"); return
    t0 = time.time()
    c = K.load_cache(ctx); emb = K.load_emb(ctx); genes = c["genes"]; n = len(genes)
    uidx = np.searchsorted(U, genes); assert all(U[i] == g for i, g in zip(uidx, genes))
    Wg = Wg_all[uidx].astype(float)
    nl = c["nl"]
    is_tf = nl["is_tf_lambert"].astype(bool); tf = np.flatnonzero(is_tf); strata = nl["det_decile"]
    trtf = np.flatnonzero(nl["is_trrust_tf"]); trtg = np.flatnonzero(nl["is_trrust_target_only"])
    lab2 = np.r_[trtf, trtg]; pos2 = np.r_[np.ones(trtf.size, bool), np.zeros(trtg.size, bool)]
    sp = K.h3_obj(c)
    Wtf = Wtf_all[uidx[sp.tfs]].astype(float); Wtg = Wtg_all[uidx[sp.tgs]].astype(float)
    edges = np.load(c["d"] / "edges_string.npy")
    epos, eneg, acc = K.rewired_nonpairs(edges, n, seed=K.SEED + K.CTXS.index(ctx))
    pa = np.r_[epos[:, 0], eneg[:, 0]]; pb = np.r_[epos[:, 1], eneg[:, 1]]
    ypair = np.r_[np.ones(len(epos), bool), np.zeros(len(eneg), bool)]
    point = np.full((len(STATS), K.N_LAYERS), np.nan); boot = np.full((len(STATS), K.N_LAYERS, B), np.nan)
    for L in range(K.N_LAYERS):
        E = K.centred(emb, L, c["obs_idx"]); En = K.unit(E)
        # S2
        sc = S.tf_scores_loo(En, tf)
        point[1, L] = S.stratified_auroc(sc, is_tf, strata)
        boot[1, L] = K.wstrat_auroc_var(K.loo_scores_weighted(En, tf, Wg), is_tf, strata, Wg)
        # S2b full
        point[2, L] = S.tf_auroc(En, trtf, trtg)
        Sc = K.loo_scores_weighted(En, trtf, Wg)
        boot[2, L] = K.wauroc_var(Sc[lab2], pos2, Wg[lab2])
        # S2b band
        US, _, _ = S.svd_observed(emb, L, c["obs_idx"])
        Eb = K.unit(US[:, K.BAND])
        point[3, L] = S.tf_auroc(Eb, trtf, trtg)
        Sc = K.loo_scores_weighted(Eb, trtf, Wg)
        boot[3, L] = K.wauroc_var(Sc[lab2], pos2, Wg[lab2])
        # S3
        s3 = (En[sp.tf] * En[sp.tg]).sum(1)
        point[4, L] = sp.auroc(s3)
        boot[4, L] = K.w_strat_pairs(sp, s3, Wtf, Wtg)
        # S1
        s1 = (En[pa] * En[pb]).sum(1)
        point[0, L] = S.auroc(s1[ypair], s1[~ypair])
        for b0 in range(0, B, 250):
            Wp = Wg[pa, b0:b0 + 250] * Wg[pb, b0:b0 + 250]
            boot[0, L, b0:b0 + 250] = K.wauroc_fixed(s1, ypair, Wp)
        log(f"{ctx} L{L}: " + " ".join(f"{s}={point[i, L]:.4f}" for i, s in enumerate(STATS)) +
            f"  ({time.time() - t0:.0f}s)")
    meta = dict(context=ctx, n_V=n, n_lambert_tf=int(tf.size), n_trrust_tf=int(trtf.size),
                n_trrust_target_only=int(trtg.size), h3_tfs=int(sp.tfs.size), h3_pos=int(sp.y.sum()),
                h3_neg=int((~sp.y).sum()), n_string_edges=int(len(epos)), n_rewired_nonpairs=int(len(eneg)),
                ms_acceptance=acc)
    np.savez(f, point=point, boot=boot, stats=np.array(STATS), meta=json.dumps(meta))
    log(f"{ctx}: saved {json.dumps(meta)}")


def main():
    OUT.mkdir(exist_ok=True)
    logf = open(K.HERE / "e10_stats.log", "a")
    def log(m):
        print(m, flush=True); logf.write(m + "\n"); logf.flush()
    U = np.array(universe())
    log(f"union universe {U.size} genes; B={B}")
    Wg, Wtf, Wtg = poisson_weights(list(U))
    order = sys.argv[1:] or K.CTXS
    for ctx in order:
        run_ctx(ctx, U, Wg, Wtf, Wtg, log)
    log("done")


if __name__ == "__main__":
    main()
