"""E7 SENSITIVITY (not part of the H4 primary; effect sizes with intervals, no pass/fail).

Paired CIs, with the SAME resamples as the primary Delta (same seeds -> same gene bootstrap weights, same two-way
TF x target weights, same vertex subsets), for:
  family R   remove-top-j (drop SV1..SVj) minus full (T512), j in {1, 2, 3, 5, 10}
  family N   screen band minus same-width top-k: nominal SV2-4, SV5-7 (vs T3), SV2-7 (vs T6); kidney_A-rule SV4-6
             (vs T3). The kidney_A-rule SV1-3 / SV1-6 bands ARE T3 / T6 (difference 0 by construction).
Contexts: the four H4 evaluation contexts + kidney_A (where the legacy mapping is defined).
Output: sens_R_N_ci.csv (per layer and layer-mean).

  python e7_sens.py --ctx lung
"""
import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import rankdata

import e7lib as L

S = L.S
HERE = L.HERE
CONTRASTS = [("R1", "T512"), ("R2", "T512"), ("R3", "T512"), ("R5", "T512"), ("R10", "T512"),
             ("Nnom_SV2-4", "T3"), ("Nnom_SV5-7", "T3"), ("Nnom_SV2-7", "T6"), ("Nmap_SV5-7", "T3")]


def ppi_vertex_draws_fast(D, Sc, B, seed, frac=0.5):
    """Same node subsets as stats3.vertex_subsample_ci / e7lib.ppi_vertex_draws; AUROC of every column of Sc on each
    subset through mid-ranks (identical to stats3.auroc incl. ties)."""
    n = D["n"]; m = int(frac * n); rng = np.random.default_rng(seed)
    pa, pb, y = D["pa"], D["pb"], D["py"]; out = np.full((Sc.shape[1], B), np.nan)
    for r in range(B):
        nodes = np.sort(rng.choice(n, m, replace=False)); keep = np.zeros(n, bool); keep[nodes] = True
        e = keep[pa] & keep[pb]
        if e.sum() <= 2: continue
        yy = y[e]; n1 = yy.sum(); n0 = yy.size - n1
        if n1 == 0 or n0 == 0: continue
        R = rankdata(Sc[e], axis=0)
        out[:, r] = (R[yy].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)
    return out, m


def run(ctx):
    out = L.WORK / f"sens_{ctx}.csv"
    if out.exists():
        print(ctx, "exists"); return
    t0 = time.time()
    D = L.load_ctx(ctx); emb = L.load_emb(ctx); svds = L.layer_svds(emb, D["obs_idx"])
    reps = {r[0]: r for r in L.representations(ctx)}
    need = sorted({a for a, b in CONTRASTS} | {b for a, b in CONTRASTS})
    FS, _ = L.load_faststrat()
    rows = []
    for t in L.TASKS:
        pts = {r: np.empty(L.N_LAYERS) for r in need}; dr = {r: np.empty((L.N_LAYERS, L.B_BOOT)) for r in need}
        if t == "NODE":
            W = L.gene_boot_weights(D["n"], L.B_BOOT, L.SEED + 101)
        elif t == "EDGE":
            fs = FS.FastStrat(D["sp"]); wts = fs.boot_weights(L.B_BOOT, np.random.default_rng(L.SEED + 102), True)
        cols = []
        for l in range(L.N_LAYERS):
            for r in need:
                En = L.unit_rep(svds[l][0], reps[r][3])
                v, s = L.METRIC[t](En, D, True); pts[r][l] = v
                if t == "NODE": dr[r][l] = L.node_boot(En, D, W)
                elif t == "EDGE": dr[r][l] = fs.auc(s, wts)
                else: cols.append(s)
            print(f"[{time.strftime('%H:%M:%S')}] {ctx} {t} layer {l} ({time.time()-t0:.0f}s)", flush=True)
        if t == "PPI":
            Sc = np.column_stack(cols); del cols
            dd, m = ppi_vertex_draws_fast(D, Sc, L.B_BOOT, L.SEED + 103)
            # verify the fast path against e7lib.ppi_vertex_draws (stats3-verified) on 5 draws, 2 columns
            chk, _ = L.ppi_vertex_draws(D, [Sc[:, 0], Sc[:, -1]], 5, L.SEED + 103)
            assert np.nanmax(np.abs(chk - dd[[0, -1], :5])) < 1e-12, "fast ppi draws mismatch"
            k = 0
            for l in range(L.N_LAYERS):
                for r in need: dr[r][l] = dd[k]; k += 1
        n = D["n"]; m = int(0.5 * n)
        for a, b in CONTRASTS:
            d = dr[a] - dr[b]
            for lay in list(range(L.N_LAYERS)) + ["mean"]:
                if lay == "mean": est = pts[a].mean() - pts[b].mean(); dd_ = d.mean(0)
                else: est = pts[a][lay] - pts[b][lay]; dd_ = d[lay]
                if t == "PPI": lo, hi = L.vertex_ci_from_draws(est, dd_, n, m)
                else: lo, hi = float(np.nanquantile(dd_, 0.025)), float(np.nanquantile(dd_, 0.975))
                rows.append(dict(ctx=ctx, task=t, rep=a, reference=b, family="R" if a.startswith("R") else "N",
                                 layer=lay, value=float(pts[a].mean() if lay == "mean" else pts[a][lay]),
                                 reference_value=float(pts[b].mean() if lay == "mean" else pts[b][lay]),
                                 diff=float(est), ci_lo=lo, ci_hi=hi, B=int(np.isfinite(dd_).sum())))
    pd.DataFrame(rows).to_csv(out, index=False)
    print(ctx, "sens done", f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--ctx"); a = ap.parse_args()
    if a.ctx == "merge":
        fr = [pd.read_csv(L.WORK / f"sens_{c}.csv") for c in L.EVAL_CTX + L.EXTRA_CTX]
        pd.concat(fr).to_csv(HERE / "sens_R_N_ci.csv", index=False); print("merged")
    else:
        run(a.ctx)
