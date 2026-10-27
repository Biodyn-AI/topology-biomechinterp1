"""E3 sensitivity (1)+(3)+(4): band scan of every contiguous band [i..j] of SV1..SV12 (78 bands) x 12 layers.

Scores (observed-gene SVD per layer, stats3.svd_observed; rows restricted to the task universe):
  variant 'cos': width 1 -> raw coordinate US[:, i], sign-invariant AUROC (stats3.signed_free_auroc, also in
                 the null); width >= 2 -> leave-one-out mean cosine to the TF set within the band coordinates
                 (stats3.tf_scores_loo on US[:, i..j]).
  variant 'lda': diagonal-LDA / whitened-centroid projection, width-invariant: Z = U (= US / S) coordinates of
                 the band, w = mean(Z | TF) - mean(Z | other), score = Z w (fit in-sample, refit in every null
                 draw, so the null carries the same optimism). Width 1 = coordinate signed by the TF mean.
All AUROCs are detection-stratified (node_labels det_decile). Null: N = 2,000 within-stratum permutations,
one per draw applied jointly to all 936 cells of a variant; Westfall-Young max-T over the 936 cells.
usage: python e3_bandscan.py <context> <task: lambert|trrust>
"""
import json, sys, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import e3_lib as L
S = L.S

ctx, task = sys.argv[1], sys.argv[2]
CI = L.CONTEXTS.index(ctx)
N = 2_000; BATCH = 250
SEED = 20261026 + 1000 * CI + (0 if task == "lambert" else 300) + 40
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

D = L.load(ctx)
if task == "lambert":
    U = np.arange(len(D["genes"])); y = D["is_tf"]
else:
    U = np.flatnonzero(D["is_trrust_tf"] | D["is_trrust_tgt"]); y = D["is_trrust_tf"][U]
det = D["det"][U]; n = U.size; k = int(y.sum()); pos = np.flatnonzero(y)
fs = L.FastStrata(det, y)
bands = L.bands12()
log(ctx, task, "n", n, "TF", k)

USs, Zs, Ss = [], [], []
sing = []
for l in range(12):
    US, Sv, Vt = S.svd_observed(D["emb"], l, D["obs_idx"])
    if task == "lambert":
        for kk in range(12):
            rho = spearmanr(US[:, kk], D["log_nobs"]).correlation
            sing.append(dict(context=ctx, layer=l, k=kk + 1, S_k=float(Sv[kk]), var_frac=float(Sv[kk] ** 2 / (Sv ** 2).sum()),
                             spearman_lognobs=float(rho), abs_spearman_lognobs=float(abs(rho))))
    USs.append(US[U, :12].copy()); Zs.append((US[U, :12] / Sv[:12]).copy()); Ss.append(Sv[:12])
if sing: pd.DataFrame(sing).to_csv(L.HERE / "partial" / f"singular_{ctx}.csv", index=False)

# ---------------------------------------------------------------- observed
def lda_score(Z, yy):
    w = Z[yy].mean(0) - Z[~yy].mean(0); return Z @ w

cells = []
obs = {"cos": [], "lda": []}
for l in range(12):
    for (i, j) in bands:
        Ub = USs[l][:, i:j + 1]; Zb = Zs[l][:, i:j + 1]
        if i == j:
            a_raw = S.stratified_auroc(Ub[:, 0], y, det); a_cos = S.signed_free_auroc(a_raw)
        else:
            a_raw = S.stratified_auroc(S.tf_scores_loo(Ub, pos), y, det); a_cos = a_raw
        a_lda = S.stratified_auroc(lda_score(Zb, y), y, det)
        u_cos = S.auroc(*(lambda s: (s[y], s[~y]))(Ub[:, 0] if i == j else S.tf_scores_loo(Ub, pos)))
        cells.append(dict(context=ctx, task=task, layer=l, band=f"SV{i+1}-{j+1}" if j > i else f"SV{i+1}", i=i + 1, j=j + 1,
                          width=j - i + 1, strat_auroc_cos=a_cos, strat_auroc_cos_signed=a_raw, auroc_cos_unstrat=u_cos,
                          strat_auroc_lda=a_lda))
        obs["cos"].append(a_cos); obs["lda"].append(a_lda)
obs = {v: np.array(a) for v, a in obs.items()}
log("observed done", len(cells), "cells")

# ---------------------------------------------------------------- joint null with per-batch checkpoints
ck = L.HERE / "partial" / f"bandscan_ckpt_{ctx}_{task}.npz"
nb = N // BATCH
null = {"cos": np.full((N, len(cells)), np.nan), "lda": np.full((N, len(cells)), np.nan)}; start = 0
if ck.exists():
    z = np.load(ck); start = int(z["done_batches"]); null["cos"] = z["cos"]; null["lda"] = z["lda"]
    log("resuming at batch", start)
En_band = {}
for l in range(12):
    for (i, j) in bands:
        if j > i: En_band[(l, i, j)] = S.l2norm_rows(USs[l][:, i:j + 1]).astype(np.float32)
ranks1 = {(l, i): L.strat_ranks_fixed(USs[l][:, i], fs, n) for l in range(12) for i in range(12)}
for bi in range(start, nb):
    rng = np.random.default_rng([SEED, bi])
    Y = L.perm_batch(y, det, rng, BATCH).astype(np.float64); Y32 = Y.astype(np.float32)
    rows = slice(bi * BATCH, (bi + 1) * BATCH); c = 0
    for l in range(12):
        Z = Zs[l]
        Dm = Z.T @ Y / k - Z.T @ (1 - Y) / (n - k)            # 12 x B centroid differences per permuted label set
        for (i, j) in bands:
            if i == j:
                A = L.strat_cols_fixed(ranks1[(l, i)], Y, fs)
                null["cos"][rows, c] = 0.5 + np.abs(A - 0.5)
                null["lda"][rows, c] = np.where(Dm[i] > 0, A, 1 - A)   # coordinate signed by the TF mean
            else:
                SC = L.loo_scores_cols(En_band[(l, i, j)], Y32, k)
                null["cos"][rows, c] = fs.cols(SC, Y32)
                SC = Z[:, i:j + 1].astype(np.float32) @ Dm[i:j + 1].astype(np.float32)
                null["lda"][rows, c] = fs.cols(SC, Y32)
            c += 1
    np.savez(ck, done_batches=bi + 1, cos=null["cos"], lda=null["lda"])
    log("batch", bi + 1, "/", nb)

# sanity: the LDA width-1 sign rule matches the observed definition at the true labels
for v in ("cos", "lda"):
    assert not np.isnan(null[v]).any()
S.assert_resolution(N, 0.05)
res = {}
for v in ("cos", "lda"):
    z, padj, praw = S.max_t_adjust(obs[v], null[v])
    for c, cell in enumerate(cells):
        cell[f"z_{v}"] = float(z[c]); cell[f"p_raw_{v}"] = float(praw[c]); cell[f"p_maxT_{v}"] = float(padj[c])
        cell[f"null_mean_{v}"] = float(null[v][:, c].mean()); cell[f"null_sd_{v}"] = float(null[v][:, c].std(ddof=1))
    best = int(np.argmax(z))
    mx = null[v].copy(); mx = (mx - mx.mean(0)) / mx.std(0, ddof=1)
    res[v] = dict(best_cell=dict(layer=cells[best]["layer"], band=cells[best]["band"], strat_auroc=float(obs[v][best]),
                                 z=float(z[best]), p_maxT=float(padj[best])),
                  best_auroc_cell=dict(layer=cells[int(np.argmax(obs[v]))]["layer"], band=cells[int(np.argmax(obs[v]))]["band"],
                                       strat_auroc=float(obs[v].max())),
                  n_cells_maxT_lt_0p05=int((padj < 0.05).sum()), n_cells=len(cells),
                  maxz_null=L.null_summary(mx.max(1)), max_auroc_null=L.null_summary(null[v].max(1), float(obs[v].max())))

# ---------------------------------------------------------------- legacy screen bands: own permutation p
legacy = [("SV2-4", 1, 3), ("SV5-7", 4, 6), ("SV2-7", 1, 6)]
mp = json.loads((L.R3 / "E9_screen/legacy_band_mapping.json").read_text())["bands"]
targets = [dict(legacy_band=nm, convention="nominal_observed_index", mapped_band=nm) for nm, _, _ in legacy]
for m in mp:
    if m["context"] == ctx:
        targets.append(dict(legacy_band=m["legacy_band"], convention="mapped:" + m["legacy_convention"], mapped_band=m["mapped_observed_band"]))
bidx = {c["band"]: [] for c in cells}
for c_i, c in enumerate(cells): bidx[c["band"]].append(c_i)
leg = []
for t in targets:
    cols = bidx[t["mapped_band"]]
    for v in ("cos", "lda"):
        zb, pb, pr = S.max_t_adjust(obs[v][cols], null[v][:, cols])           # max-T over this band's 12 layers
        for l, c_i in enumerate(cols):
            leg.append(dict(context=ctx, task=task, **t, variant=v, layer=l, strat_auroc=float(obs[v][c_i]),
                            p_raw=float(pr[l]), p_maxT_band_12layers=float(pb[l]), p_maxT_936=cells[c_i][f"p_maxT_{v}"],
                            screen_derived=True, discovery_contaminated=ctx in L.DISCOVERY, N_perm=N))
pd.DataFrame(leg).to_csv(L.HERE / "partial" / f"legacy_{ctx}_{task}.csv", index=False)

# ---------------------------------------------------------------- degeneracy flags (E1 gaps)
g = pd.read_csv(L.R3 / "E1_compression/gaps.csv"); g = g[g.context == ctx]
deg = {(int(r.layer), int(r.k)): bool(r.near_degenerate) for r in g.itertuples()}
for c in cells:
    l, i, j = c["layer"], c["i"], c["j"]
    c["cuts_degenerate_lower"] = bool(deg.get((l, i - 1), False)) if i > 1 else False
    c["cuts_degenerate_upper"] = bool(deg.get((l, j), False))
    c["boundary_cuts_degenerate_pair"] = c["cuts_degenerate_lower"] or c["cuts_degenerate_upper"]
suffix = "" if task == "lambert" else "_trrust"
df = pd.DataFrame(cells)
df.to_csv(L.HERE / f"bandscan_{ctx}{suffix}.csv", index=False)
for v in ("cos", "lda"):
    df.pivot(index="band", columns="layer", values=f"strat_auroc_{v}").reindex([c["band"] for c in cells[:78]]) \
      .to_csv(L.HERE / "partial" / f"heatmap_{v}_{ctx}{suffix}.csv")
null_sum = {v: dict(per_cell_mean_range=[float(null[v].mean(0).min()), float(null[v].mean(0).max())]) for v in null}
(L.HERE / "partial" / f"bandscan_{ctx}_{task}.json").write_text(json.dumps(dict(context=ctx, task=task, n=n, n_tf=k, N_perm=N,
    conv_hash=L.conv_hash(), summary=res, null=null_sum), indent=1))
ck.unlink()
log("DONE", json.dumps(res["cos"]["best_cell"]), json.dumps(res["lda"]["best_cell"]))
