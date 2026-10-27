"""E1 - compression of the observed-gene representation (CONVENTION3 section 5, sensitivity; answers D2).

For every context in contexts.CONTEXTS and every layer L0..L11, on the mean-centred V x 512 matrix
(V = observed genes, cache/<ctx>/obs_idx.npy; SVD via stats3.svd_observed):
  1. effective rank (stats3.effective_rank), participation ratio (stats3.participation_ratio), variance
     fraction on SV1 and cumulative on SV1..10                                    -> per_context_layer.csv
  2. gap ratios S_k/S_{k+1}, k = 1..11, gene bootstrap B = 200 (resample genes with replacement, recentre,
     recompute the spectrum). Because sorted singular values make every bootstrap ratio >= 1, the percentile
     CI can never contain 1; the near-degeneracy flag therefore uses the basic (reverse-percentile) CI
     [2r - q97.5, 2r - q2.5], which corrects the upward sorting bias. Eigenvector stability of SVk / SVk+1
     and of their 2-D span is reported alongside                                   -> gaps.csv
     (the same bootstrap gives 95% CIs for ER, PR and SV1 variance in per_context_layer.csv)
  3. Spearman rho (signed and |.|) of each SVk coordinate (k = 1..10) with log n_obs -> detection_corr.csv
  4. feature-shuffle control at L0 and L11: each of the 512 columns permuted independently across genes,
     20 draws                                                                       -> shuffle.csv
  5. depth trend: Spearman rho of ER and PR against layer index (descriptive)     -> results.json / SUMMARY
  6. detection-matched comparability, five core contexts, n = 1,300 genes, 200 draws, ER/PR/SV1 at L0, L11:
       (a) 'pooled'  : common absolute detection distribution - bins = deciles of pooled log n_obs over the
                       five contexts; per-bin quota proportional to min over contexts of genes available
                       in that bin, scaled to 1,300 (sampling without replacement)
       (b) 'within'  : proportional allocation over the context's own log n_obs deciles (size equalised,
                       the context's detection distribution kept)
                                                                                    -> matched.csv
  7. zero-row artefact (kidney_A, immune; L0, L11): the same statistics on the FULL G x 512 matrix
     including never-observed all-zero rows, and corr(SV1, observed indicator)     -> zero_row_artefact.csv
  extra (CONVENTION3 section 2 detection sensitivity): ER/PR/SV1 at L0, L11 on V restricted to
     n_obs >= 3, 5, 10                                                              -> detection_threshold.csv
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import rankdata, spearmanr, pearsonr

HERE = Path(__file__).resolve().parent
R3 = HERE.parent
sys.path.insert(0, str(R3))
import contexts as X  # noqa: E402
import stats3 as S  # noqa: E402

EMB = Path.home() / "biomechinterp-local/subproject41_r3/emb"
CACHE = R3 / "cache"
PART = HERE / "partial"
SEED = 20261026
B_GAP = 200
N_SHUF = 20
N_MATCH = 200
N_SUB = 1300
CORE = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
CTX = list(X.CONTEXTS)
LAYERS = list(range(12))
KMAX = 12          # gaps for k = 1..11 need S_1..S_12
KDET = 10


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def load(ctx):
    emb = np.load(EMB / f"{ctx}.npy", mmap_mode="r")
    obs_idx = np.load(CACHE / ctx / "obs_idx.npy")
    ncell = np.load(CACHE / ctx / "ncells.npy").astype(float)
    assert emb.shape[0] == 12 and emb.shape[2] == 512 and ncell.size == obs_idx.size
    assert ncell.min() >= 1, f"{ctx}: observed gene with n_obs < 1"
    return emb, obs_idx, ncell


def spec_stats(Sv):
    v = Sv ** 2
    return dict(er=S.effective_rank(Sv), pr=S.participation_ratio(Sv), sv1_var=float(v[0] / v.sum()),
                cum10_var=float(v[:10].sum() / v.sum()))


def gram_spectrum(Xb, vectors=False):
    """Singular values (descending) of the column-centred Xb via the 512 x 512 Gram matrix
    (agrees with the thin SVD to ~1e-15 relative on the top components; checked in run())."""
    Xc = Xb - Xb.mean(axis=0, keepdims=True)
    if vectors:
        w, v = np.linalg.eigh(Xc.T @ Xc)
        return np.sqrt(np.clip(w[::-1], 0, None)), v[:, ::-1].T
    w = np.linalg.eigvalsh(Xc.T @ Xc)
    return np.sqrt(np.clip(w[::-1], 0, None))


def qsum(a):
    a = np.asarray(a, float)
    return dict(mean=float(a.mean()), sd=float(a.std(ddof=1)), q025=float(np.quantile(a, 0.025)),
                q50=float(np.quantile(a, 0.5)), q975=float(np.quantile(a, 0.975)))


# --------------------------------------------------------------------------------------------- per context
def run_context(ci, ctx):
    out_p = PART / f"{ctx}.json"
    if out_p.exists():
        return json.loads(out_p.read_text())
    emb, obs_idx, ncell = load(ctx)
    lognobs = np.log(ncell)
    r_det = rankdata(lognobs)
    V = obs_idx.size
    rows, gaps, dets, shuf, thr = [], [], [], [], []
    for L in LAYERS:
        US, Sv, Vt = S.svd_observed(emb, L, obs_idx)
        st = spec_stats(Sv)
        # 2. gene bootstrap of the spectrum
        Xr = np.asarray(emb[L], dtype=np.float64)[obs_idx]
        rng = np.random.default_rng([SEED, ci, L, 2])
        bs = np.empty((B_GAP, KMAX)); bst = {k: np.empty(B_GAP) for k in ("er", "pr", "sv1_var")}
        cos_ind = np.empty((B_GAP, KMAX)); cos_pair = np.empty((B_GAP, KMAX - 1))
        for b in range(B_GAP):
            sb, vb = gram_spectrum(Xr[rng.integers(0, V, V)], vectors=True)
            bs[b] = sb[:KMAX]
            sstat = spec_stats(sb)
            for k in bst: bst[k][b] = sstat[k]
            M = Vt[:KMAX] @ vb[:KMAX].T                  # |cos| between original and bootstrap directions
            cos_ind[b] = np.abs(np.diag(M))
            for k in range(KMAX - 1):                    # smallest principal-angle cosine of span(k, k+1)
                cos_pair[b, k] = np.linalg.svd(M[k:k + 2, k:k + 2], compute_uv=False).min()
        row = dict(context=ctx, layer=L, n_V=V, **st)
        for k in bst:
            row[f"{k}_boot_lo"] = float(np.quantile(bst[k], 0.025))
            row[f"{k}_boot_hi"] = float(np.quantile(bst[k], 0.975))
        row["sv2_var"] = float(Sv[1] ** 2 / (Sv ** 2).sum())
        row["sv3_var"] = float(Sv[2] ** 2 / (Sv ** 2).sum())
        row["n_sv_below_1e-3_rel"] = int((Sv / Sv[0] < 1e-3).sum())
        rows.append(row)
        rb = bs[:, :-1] / bs[:, 1:]
        for k in range(KMAX - 1):
            r = float(Sv[k] / Sv[k + 1]); qlo, qhi = np.quantile(rb[:, k], [0.025, 0.975])
            gaps.append(dict(context=ctx, layer=L, k=k + 1, S_k=float(Sv[k]), S_k1=float(Sv[k + 1]), ratio=r,
                             boot_mean=float(rb[:, k].mean()), boot_sd=float(rb[:, k].std(ddof=1)),
                             pct_lo=float(qlo), pct_hi=float(qhi), basic_lo=float(2 * r - qhi),
                             basic_hi=float(2 * r - qlo), near_degenerate=bool(2 * r - qhi <= 1.0),
                             med_abscos_svk=float(np.median(cos_ind[:, k])),
                             med_abscos_svk1=float(np.median(cos_ind[:, k + 1])),
                             med_pair_span_mincos=float(np.median(cos_pair[:, k]))))
        # 3. detection tracking of SV1..10
        for k in range(KDET):
            rho = float(np.corrcoef(rankdata(US[:, k]), r_det)[0, 1])
            dets.append(dict(context=ctx, layer=L, k=k + 1, spearman=rho, abs_spearman=abs(rho),
                             var_frac=float(Sv[k] ** 2 / (Sv ** 2).sum())))
        # 4. feature shuffle at L0, L11
        if L in (0, 11):
            rng_s = np.random.default_rng([SEED, ci, L, 4])
            ss = [spec_stats(gram_spectrum(rng_s.permuted(Xr, axis=0))) for _ in range(N_SHUF)]
            d = dict(context=ctx, layer=L, n_V=V, n_draws=N_SHUF, er_obs=st["er"], pr_obs=st["pr"],
                     sv1_var_obs=st["sv1_var"])
            for k in ("er", "pr", "sv1_var"):
                a = np.array([s[k] for s in ss])
                d[f"{k}_shuf_mean"] = float(a.mean()); d[f"{k}_shuf_sd"] = float(a.std(ddof=1))
                d[f"{k}_shuf_min"] = float(a.min()); d[f"{k}_shuf_max"] = float(a.max())
            d["er_shuf_over_obs"] = d["er_shuf_mean"] / st["er"]
            shuf.append(d)
            # extra: detection threshold sensitivity
            for t in (1, 3, 5, 10):
                keep = ncell >= t
                Xt = Xr[keep]
                thr.append(dict(context=ctx, layer=L, min_nobs=t, n_genes=int(keep.sum()),
                                **spec_stats(np.linalg.svd(Xt - Xt.mean(0), compute_uv=False))))
    res = dict(rows=rows, gaps=gaps, dets=dets, shuf=shuf, thr=thr)
    out_p.write_text(json.dumps(res))
    log(f"done {ctx} (V={V})")
    return res


# ------------------------------------------------------------------------------------------ matched (core)
def run_matched():
    out_p = PART / "matched.json"
    if out_p.exists():
        return json.loads(out_p.read_text())
    data = {c: load(c) for c in CORE}
    pooled = np.concatenate([np.log(data[c][2]) for c in CORE])
    edges = np.unique(np.quantile(pooled, np.linspace(0, 1, 11)))
    nb = len(edges) - 1

    def binof(x):
        return np.clip(np.searchsorted(edges, x, side="right") - 1, 0, nb - 1)

    avail = np.array([[int((binof(np.log(data[c][2])) == b).sum()) for b in range(nb)] for c in CORE])
    mn = avail.min(0)
    assert mn.sum() >= N_SUB, f"pooled matching infeasible: {mn.sum()} < {N_SUB}"
    quota = np.floor(mn / mn.sum() * N_SUB).astype(int)
    rem = N_SUB - quota.sum()
    frac = mn / mn.sum() * N_SUB - quota
    for j in np.argsort(-frac)[:rem]: quota[j] += 1
    assert (quota <= mn).all() and quota.sum() == N_SUB
    info = dict(pooled_bin_edges_log_nobs=edges.tolist(), pooled_bin_edges_nobs=np.exp(edges).tolist(),
                avail_per_context_bin={c: avail[i].tolist() for i, c in enumerate(CORE)},
                pooled_quota=quota.tolist())
    rows = []
    for ci, c in enumerate(CORE):
        emb, obs_idx, ncell = data[c]
        lg = np.log(ncell); pb = binof(lg); wb = S.quantile_strata(lg, 10)
        nwb = int(wb.max()) + 1
        # within-context deciles: proportional allocation of 1300 over the realised strata (keeps the
        # context's own detection distribution; only |V| is equalised)
        wavail = np.bincount(wb, minlength=nwb)
        fr = wavail / wavail.sum() * N_SUB; wq = np.floor(fr).astype(int)
        for j in np.argsort(-(fr - wq))[: N_SUB - wq.sum()]: wq[j] += 1
        info.setdefault("within_decile_counts", {})[c] = wavail.tolist()
        for L in (0, 11):
            Xr = np.asarray(emb[L], dtype=np.float64)[obs_idx]
            full = spec_stats(S.svd_observed(emb, L, obs_idx)[1])
            for scheme, strata, q in (("pooled", pb, quota), ("within", wb, wq)):
                rng = np.random.default_rng([SEED, ci, L, 6, 0 if scheme == "pooled" else 1])
                draws = {k: [] for k in ("er", "pr", "sv1_var")}
                with_repl = False
                for _ in range(N_MATCH):
                    idx = []
                    for b in range(len(q)):
                        pool = np.flatnonzero(strata == b)
                        rep = q[b] > pool.size
                        with_repl |= rep
                        idx.append(rng.choice(pool, q[b], replace=rep))
                    st = spec_stats(gram_spectrum(Xr[np.concatenate(idx)]))
                    for k in draws: draws[k].append(st[k])
                d = dict(context=c, layer=L, scheme=scheme, n_sub=N_SUB, n_draws=N_MATCH, n_V=obs_idx.size,
                         any_bin_with_replacement=bool(with_repl), er_full_V=full["er"], pr_full_V=full["pr"],
                         sv1_var_full_V=full["sv1_var"],
                         median_nobs_sub=float(np.median(ncell[np.concatenate(idx)])))
                for k, a in draws.items():
                    s = qsum(a)
                    for kk, vv in s.items(): d[f"{k}_{kk}"] = vv
                rows.append(d)
        log(f"matched {c}")
    res = dict(rows=rows, info=info)
    out_p.write_text(json.dumps(res))
    return res


# ------------------------------------------------------------------------------------------ zero-row artefact
def run_zero_rows():
    rows = []
    for ci, c in enumerate(["kidney_A", "immune"]):
        emb = np.load(EMB / f"{c}.npy", mmap_mode="r")
        obs_idx = np.load(CACHE / c / "obs_idx.npy")
        G = emb.shape[1]
        ind = np.zeros(G); ind[obs_idx] = 1.0
        for L in (0, 11):
            for matrix, idx in (("full_incl_zero_rows", np.arange(G)), ("observed_only", obs_idx)):
                US, Sv, Vt = S.svd_observed(emb, L, idx)
                st = spec_stats(Sv)
                d = dict(context=c, layer=L, matrix=matrix, n_rows=int(idx.size), n_zero_rows=int(G - obs_idx.size)
                         if matrix.startswith("full") else 0, **st)
                if matrix.startswith("full"):
                    d["pearson_sv1_observed_indicator"] = float(pearsonr(US[:, 0], ind)[0])
                    d["spearman_sv1_observed_indicator"] = float(spearmanr(US[:, 0], ind)[0])
                    d["max_abs_corr_sv2to10_observed_indicator"] = float(max(abs(pearsonr(US[:, k], ind)[0])
                                                                           for k in range(1, 10)))
                    Xr = np.asarray(emb[L], dtype=np.float64)[idx]
                    rng = np.random.default_rng([SEED, 100 + ci, L, 7])
                    sh = [spec_stats(gram_spectrum(rng.permuted(Xr, axis=0)))["er"] for _ in range(N_SHUF)]
                    d["er_feature_shuffle_mean"] = float(np.mean(sh)); d["er_feature_shuffle_sd"] = float(np.std(sh, ddof=1))
                rows.append(d)
    return rows


def run():
    PART.mkdir(exist_ok=True)
    # numerical check: Gram spectrum == thin SVD (lung, L11)
    emb, obs_idx, _ = load("lung")
    _, Sv, _ = S.svd_observed(emb, 11, obs_idx)
    sg = gram_spectrum(np.asarray(emb[11], np.float64)[obs_idx])
    chk = dict(max_rel_diff_top12=float(np.abs(sg[:12] / Sv[:12] - 1).max()),
               er_svd=S.effective_rank(Sv), er_gram=S.effective_rank(sg))
    log(f"gram check {chk}")
    t0 = time.time()
    res = Parallel(n_jobs=3, verbose=0)(delayed(run_context)(i, c) for i, c in enumerate(CTX))
    log(f"contexts done in {time.time() - t0:.0f}s")
    per = pd.DataFrame([r for x in res for r in x["rows"]])
    gaps = pd.DataFrame([r for x in res for r in x["gaps"]])
    dets = pd.DataFrame([r for x in res for r in x["dets"]])
    shuf = pd.DataFrame([r for x in res for r in x["shuf"]])
    thr = pd.DataFrame([r for x in res for r in x["thr"]])
    per.to_csv(HERE / "per_context_layer.csv", index=False)
    gaps.to_csv(HERE / "gaps.csv", index=False)
    dets.to_csv(HERE / "detection_corr.csv", index=False)
    shuf.to_csv(HERE / "shuffle.csv", index=False)
    thr.to_csv(HERE / "detection_threshold.csv", index=False)
    m = run_matched()
    pd.DataFrame(m["rows"]).to_csv(HERE / "matched.csv", index=False)
    z = pd.DataFrame(run_zero_rows())
    z.to_csv(HERE / "zero_row_artefact.csv", index=False)
    log("zero-row done")

    # 5. depth trend + monotonicity
    trend = {}
    for c, g in per.groupby("context", sort=False):
        g = g.sort_values("layer")
        er, pr = g.er.to_numpy(), g.pr.to_numpy()
        der, dpr = np.diff(er), np.diff(pr)
        trend[c] = dict(
            spearman_er_layer=float(spearmanr(g.layer, er)[0]), spearman_pr_layer=float(spearmanr(g.layer, pr)[0]),
            spearman_sv1var_layer=float(spearmanr(g.layer, g.sv1_var)[0]),
            er_L0=float(er[0]), er_L11=float(er[-1]), er_fold_L0_over_L11=float(er[0] / er[-1]),
            pr_L0=float(pr[0]), pr_L11=float(pr[-1]), pr_fold_L0_over_L11=float(pr[0] / pr[-1]),
            sv1_var_L0=float(g.sv1_var.iloc[0]), sv1_var_L11=float(g.sv1_var.iloc[-1]),
            er_n_increases=int((der > 0).sum()), er_max_increase=float(der.max()),
            er_strictly_decreasing=bool((der < 0).all()), er_argmin_layer=int(g.layer.iloc[er.argmin()]),
            pr_n_increases=int((dpr > 0).sum()), pr_strictly_decreasing=bool((dpr < 0).all()),
            pr_argmin_layer=int(g.layer.iloc[pr.argmin()]),
            er_half_reached_layer=(int(g.layer.iloc[np.flatnonzero(er <= er[0] / 2)[0]])
                                   if (er <= er[0] / 2).any() else None))
    tr = pd.DataFrame(trend).T.reset_index().rename(columns={"index": "context"})
    tr.to_csv(HERE / "depth_trend.csv", index=False)

    # headline numbers
    fh = (R3 / "FROZEN.txt").read_text().split()
    frozen_hash = fh[fh.index("sha256") + 1]
    actual = hashlib.sha256((R3 / "CONVENTION3.md").read_bytes()).hexdigest()
    dg = gaps[gaps.near_degenerate]
    det_max = dets.loc[dets.groupby(["context", "layer"]).abs_spearman.idxmax()]
    results = dict(
        experiment="E1_compression", convention3_sha256_frozen=frozen_hash, convention3_sha256_actual=actual,
        hash_matches=bool(frozen_hash == actual), seed=SEED, numpy=np.__version__,
        specs=dict(B_gap=B_GAP, n_shuffle=N_SHUF, n_matched_draws=N_MATCH, n_sub=N_SUB, core=CORE,
                   degeneracy_rule="basic bootstrap CI lower bound 2r - q97.5 <= 1 (percentile CI is bounded "
                                   "below by 1 because bootstrap singular values are sorted)",
                   embeddings=str(EMB), gram_check=chk),
        depth_trend=trend,
        core_table={c: {k: trend[c][k] for k in ("er_L0", "er_L11", "er_fold_L0_over_L11", "pr_L0", "pr_L11",
                                                 "sv1_var_L0", "sv1_var_L11", "spearman_er_layer",
                                                 "spearman_pr_layer", "er_strictly_decreasing",
                                                 "er_n_increases", "pr_strictly_decreasing", "pr_n_increases")}
                    for c in CORE},
        near_degenerate_pairs=dict(n_total=int(len(dg)), n_tested=int(len(gaps)),
                                   per_core_context={c: dg[dg.context == c][["layer", "k"]].values.tolist()
                                                     for c in CORE}),
        detection=dict(
            max_abs_spearman_by_context_layer=det_max[["context", "layer", "k", "abs_spearman"]].to_dict("records"),
            sv1_abs_spearman_core={c: dets[(dets.context == c) & (dets.k == 1)].set_index("layer").abs_spearman
                                   .round(4).to_dict() for c in CORE}),
        shuffle=shuf[["context", "layer", "er_obs", "er_shuf_mean", "er_shuf_sd", "er_shuf_over_obs"]]
        .to_dict("records"),
        matched=m["info"], zero_row=z.to_dict("records"))
    (HERE / "results.json").write_text(json.dumps(results, indent=1, default=float))
    log("all done")


if __name__ == "__main__":
    run()
