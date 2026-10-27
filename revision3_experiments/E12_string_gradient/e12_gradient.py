"""E12 - STRING confidence gradient under the round-3 convention.

Question: among physical STRING links (combined score >= 150) between observed genes, do higher-confidence
links join genes whose centred, full-dimensional scGPT vectors are more similar, beyond what the
endpoints' network degree, detection and co-expression explain?

The screen reported a perfect rank correlation (rho = 1 over five score quintiles, p = 1.4e-24, exact 0.0167)
for a co-pole z on one direction of the kidney_A 209-gene subset. Here, per context and layer:
  rho      Spearman(score, cosine) over all links
  rho_adj  partial Spearman controlling for rank(log deg_i + log deg_j), rank(log nobs_i + log nobs_j),
           rank(|co-expression|) (deg = degree in the >= 150 physical network within V)
  bins     STRING bands 150-399 / 400-699 / 700-899 / 900-999: mean cosine and AUROC vs degree-matched
           non-linked pairs (one per link, same endpoint degree deciles, no physical link at any score)
  sv_rho   Spearman(score, -|x_i - x_j|) for SV1..SV3 coordinates (descriptive; screen-band mapping)
Intervals: vertex subsampling without replacement, m = n/2, delete-d scaling sqrt(m/(n-m)) = 1
(CONVENTION3_ERRATA item 3), B = 500, percentile of (stat_b - stat) added to stat. In each subsample the statistic is
computed on the full-sample ranks restricted to the retained links (Pearson correlation of ranks; the partial version by
residualisation on the restricted covariate ranks), which avoids re-ranking 500 times per layer.
"""
import sys, json, time
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import contexts as X
import stats3 as S

OUT = Path(__file__).parent
CACHE = OUT.parent / "cache"
CTXS = sys.argv[1:] or ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
BANDS = [(150, 400), (400, 700), (700, 900), (900, 1001)]
B = 500
SEED = 20261026


def spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra -= ra.mean(); rb -= rb.mean()
    return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))


def partial_spearman(a, b, Z):
    """Pearson correlation of rank(a), rank(b) after OLS residualisation on ranked covariates Z."""
    R = np.column_stack([np.ones(len(a))] + [rankdata(z) for z in Z.T])
    ra, rb = rankdata(a), rankdata(b)
    ea = ra - R @ np.linalg.lstsq(R, ra, rcond=None)[0]
    eb = rb - R @ np.linalg.lstsq(R, rb, rcond=None)[0]
    return float(ea @ eb / np.sqrt((ea @ ea) * (eb @ eb)))


def rank_stats_masked(ra, rb, R, k):
    """Subsample versions on the full-sample ranks restricted to mask k (Pearson of ranks; partial via normal equations)."""
    a, b, M = ra[k], rb[k], R[k]
    ac, bc = a - a.mean(), b - b.mean()
    r = float(ac @ bc / np.sqrt((ac @ ac) * (bc @ bc)))
    G = M.T @ M
    ea = a - M @ np.linalg.solve(G, M.T @ a); eb = b - M @ np.linalg.solve(G, M.T @ b)
    return r, float(ea @ eb / np.sqrt((ea @ ea) * (eb @ eb)))


print("loading STRING physical >= 150", flush=True)
STG = X.load_string_genomewide(150)
rows = []
for ctx in CTXS:
    t0 = time.time()
    d = CACHE / ctx
    og = json.load(open(d / "obs_genes.json")); obs_idx = np.load(d / "obs_idx.npy")
    pos = {g: k for k, g in enumerate(og)}; n = len(og)
    st = [((pos[a], pos[b]), s) for (a, b), s in STG.items() if a in pos and b in pos]
    E_ = np.array([p for p, _ in st], dtype=np.int64); sc = np.array([s for _, s in st], dtype=float)
    pi, pj = E_[:, 0], E_[:, 1]
    deg = np.bincount(np.r_[pi, pj], minlength=n).astype(float)
    nobs = np.load(d / "nobs_reproduced.npy").astype(float)
    Xe = np.load(d / "expr.npy").astype(np.float64)
    Xe = (Xe - Xe.mean(0)) / (Xe.std(0) + 1e-12)
    coexp = np.abs(np.einsum("ci,ci->i", Xe[:, pi], Xe[:, pj]) / Xe.shape[0]); del Xe
    Z = np.column_stack([np.log(deg[pi]) + np.log(deg[pj]), np.log1p(nobs[pi]) + np.log1p(nobs[pj]), coexp])

    # degree-matched non-linked pairs (one per link), drawn once
    rng = np.random.default_rng(SEED)
    dec = np.minimum((rankdata(deg, method="ordinal") - 1) * 10 // n, 9).astype(int)
    linked = set(map(tuple, np.sort(E_, axis=1)))
    pool = [np.flatnonzero(dec == k) for k in range(10)]
    ni, nj = np.empty(len(pi), int), np.empty(len(pi), int)
    for e in range(len(pi)):
        while True:
            a = rng.choice(pool[dec[pi[e]]]); b = rng.choice(pool[dec[pj[e]]])
            if a != b and (min(a, b), max(a, b)) not in linked: break
        ni[e], nj[e] = a, b

    emb, _ = X.load_context(ctx)
    masks = []
    for b in range(B):
        keep = np.zeros(n, bool); keep[np.random.default_rng(SEED + 1 + b).choice(n, n // 2, replace=False)] = True
        masks.append(keep[pi] & keep[pj])
    rsc = rankdata(sc); RZ = np.column_stack([np.ones(len(sc))] + [rankdata(z) for z in Z.T])
    print(f"{ctx}: {n} genes, {len(sc)} links, mean links per subsample {np.mean([m.sum() for m in masks]):.0f}", flush=True)
    for L in range(12):
        Xl = np.asarray(emb[L], dtype=np.float64)[obs_idx]; Xl -= Xl.mean(0)
        En = S.l2norm_rows(Xl)
        cos = np.einsum("ij,ij->i", En[pi], En[pj]); cosn = np.einsum("ij,ij->i", En[ni], En[nj])
        US, _, _ = S.svd_observed(emb, L, obs_idx)
        rho = spearman(sc, cos); rho_adj = partial_spearman(sc, cos, Z)
        rb, rab = np.empty(B), np.empty(B)
        rcos = rankdata(cos)
        for b, k in enumerate(masks):
            rb[b], rab[b] = rank_stats_masked(rsc, rcos, RZ, k)
        row = dict(ctx=ctx, layer=L, n_genes=n, n_links=len(sc), rho=rho,
                   rho_lo=rho + np.quantile(rb - rho, 0.025), rho_hi=rho + np.quantile(rb - rho, 0.975),
                   rho_adj=rho_adj, rho_adj_lo=rho_adj + np.quantile(rab - rho_adj, 0.025),
                   rho_adj_hi=rho_adj + np.quantile(rab - rho_adj, 0.975),
                   rho_score_coexp=spearman(sc, coexp), rho_score_degprod=spearman(sc, Z[:, 0]))
        for lo, hi in BANDS:
            m = (sc >= lo) & (sc < hi)
            row[f"n_{lo}"] = int(m.sum()); row[f"meancos_{lo}"] = float(cos[m].mean())
            row[f"auroc_{lo}"] = S.auroc(cos[m], cosn[m])
        for k in range(3):
            row[f"sv{k+1}_rho"] = spearman(sc, -np.abs(US[pi, k] - US[pj, k]))
        rows.append(row)
        print(f"{ctx} L{L} rho {rho:+.3f} [{row['rho_lo']:+.3f},{row['rho_hi']:+.3f}] adj {rho_adj:+.3f} "
              f"[{row['rho_adj_lo']:+.3f},{row['rho_adj_hi']:+.3f}] auroc " +
              " ".join(f"{row[f'auroc_{lo}']:.3f}" for lo, _ in BANDS), flush=True)
    print(f"{ctx} done {time.time()-t0:.0f}s, {len(sc)} links", flush=True)
    import pandas as pd
    pd.DataFrame(rows).to_csv(OUT / ("gradient.csv" if len(sys.argv) == 1 else f"gradient_{'_'.join(CTXS)}.csv"), index=False)
