"""E2 pass 2 (needs pass-1 _ckpt/<ctx>_ms_<N>.npz and _label_<N>.npz). For the top-3 cells of the primary
family (largest |z| under Maslov-Sneppen):
  * replay of the IDENTICAL Maslov-Sneppen chain (same seed) and label permutations (same seed), now also
    evaluating the co-pole rate at K/n in {0.10,0.20,0.27,0.33,0.50} (poles within V) for the top cells and
    for two fixed diagnostic cells (L0 SV1, L11 SV1); replayed r must equal pass 1 exactly;
  * vertex-subsampling CI (stats3.vertex_subsample_ci, B=500) for r and for the co-pole ratios
    (observed / analytic null mean: stats3.copole_expected_configuration for MS; 2K(K-1)/(n(n-1)) for label);
  * hub exclusion (top 5% / 10% by within-V degree removed) with a fresh MS chain, N=2000;
  * leave-one-gene-out jackknife of r for the top cell.
Writes _ckpt/<ctx>_pass2.json.
"""
import argparse, json, time
import numpy as np
import e2_common as E

S = E.S
DIAG = [0 * 10 + 0, 11 * 10 + 0]          # L0 SV1, L11 SV1


def top_cells(ctx, N, k=3):
    z = np.load(E.CKPT / f"{ctx}_ms_{N}.npz")
    obs, null = z["obs"], z["null"]
    zz = (obs - null.mean(0)) / null.std(0, ddof=1)
    order = np.argsort(-np.abs(zz), kind="mergesort")
    return [int(c) for c in order[:k]], obs, null


def label_mean(n, K):
    return 2.0 * K * (K - 1) / (n * (n - 1))


def run(ctx, N, B=500, N_hub=2000):
    log = E.logger(E.HERE / f"e2_pass2_{ctx}.log")
    out_p = E.CKPT / f"{ctx}_pass2.json"
    if out_p.exists():
        log(f"[{ctx}] pass2 exists"); return
    t0 = time.time()
    c = E.load_ctx(ctx); n = c["n"]; pi, pj = c["edges"][:, 0].copy(), c["edges"][:, 1].copy()
    tops, obs_all, null_ms = top_cells(ctx, N)
    lab_ck = np.load(E.CKPT / f"{ctx}_label_{N}.npz")["null"]
    Xc = E.coords(ctx, c["obs_idx"])
    cells = tops + [d for d in DIAG if d not in tops]
    Xs = Xc[:, cells]
    Lm = E.pole_matrix(Xs)                                   # column = fi*len(cells) + j
    nc = len(cells); F = len(E.POLE_FRACS)
    E.check_vectorised(Xs, Lm, pi, pj, cols=tuple(range(nc)))
    cp_obs = E.copole_cols(Lm, pi, pj)
    cp_exp_cfg = np.array([S.copole_expected_configuration(Lm[:, i], pi, pj) for i in range(Lm.shape[1])])
    cp_exp_lab = np.array([label_mean(n, max(1, int(f * n))) for f in E.POLE_FRACS for _ in cells])
    log(f"[{ctx}] top cells {[E.cell_name(x) for x in tops]}; replaying MS N={N}")
    # --- replay MS
    R, C, acc = E.run_ms(Xs, Lm, pi, pj, n, N, seed=E.seed_for(ctx, "ms"), log=log, every=2000)
    assert np.allclose(R, null_ms[:, cells], rtol=0, atol=1e-10), "MS replay differs from pass 1"
    rng = np.random.default_rng(E.seed_for(ctx, "label"))
    RL, CL = E.run_perm(Xs, Lm, pi, pj, n, N, rng)
    assert np.allclose(RL, lab_ck[:, cells], rtol=0, atol=1e-10), "label replay differs from pass 1"
    log(f"[{ctx}] replays identical to pass 1 ({time.time()-t0:.0f}s)")
    res = dict(context=ctx, N=N, top_cells=tops, cells=cells, fracs=list(E.POLE_FRACS))
    cp = []
    for fi, f in enumerate(E.POLE_FRACS):
        for j, cell in enumerate(cells):
            col = fi * nc + j
            for null_name, nul, ana in (("ms", C[:, col], cp_exp_cfg[col]), ("label", CL[:, col], cp_exp_lab[col])):
                p, b, NN = S.emp_p(cp_obs[col], nul, "greater")
                cp.append(dict(cell=cell, layer=E.cell_name(cell)[0], sv=E.cell_name(cell)[1], pole_frac=f,
                               K=max(1, int(f * n)), null=null_name, copole_obs=float(cp_obs[col]),
                               null_mean=float(nul.mean()), null_sd=float(nul.std(ddof=1)),
                               analytic_null_mean=float(ana), ratio_obs_null=float(cp_obs[col] / nul.mean()),
                               z=float((cp_obs[col] - nul.mean()) / nul.std(ddof=1)), p_greater=float(p), b=b, N=NN))
    res["copole"] = cp
    res["ms_acc_replay_mean"] = float(acc.mean())
    # --- vertex subsampling CIs
    rng = np.random.default_rng(E.seed_for(ctx, "vsub"))
    vci = []
    for j, cell in enumerate(tops):
        x = Xc[:, cell]
        est, lo, hi, dr = S.vertex_subsample_ci(lambda nodes, a, b: S.assortativity(x, a, b), n, pi, pj, B, rng)
        vci.append(dict(cell=cell, layer=E.cell_name(cell)[0], sv=E.cell_name(cell)[1], stat="r",
                        pole_frac=None, est=est, lo=lo, hi=hi, B=B, n_nan=int(np.isnan(dr).sum())))
        for f in E.POLE_FRACS:
            for kind in ("ms", "label"):
                def fn(nodes, a, b, x=x, f=f, kind=kind):
                    pos = np.full(n, -1); pos[nodes] = np.arange(nodes.size)
                    lab = S.pole_labels(x[nodes], f); aa, bb = pos[a], pos[b]
                    rate = S.copole_rate(lab, aa, bb)
                    den = S.copole_expected_configuration(lab, aa, bb) if kind == "ms" else \
                        label_mean(nodes.size, max(1, int(f * nodes.size)))
                    return rate / den
                est, lo, hi, dr = S.vertex_subsample_ci(fn, n, pi, pj, B, rng)
                vci.append(dict(cell=cell, layer=E.cell_name(cell)[0], sv=E.cell_name(cell)[1],
                                stat=f"copole_ratio_{kind}_analytic", pole_frac=f, est=est, lo=lo, hi=hi, B=B,
                                n_nan=int(np.isnan(dr).sum())))
    res["vertex_ci"] = vci
    log(f"[{ctx}] vertex CIs done ({time.time()-t0:.0f}s)")
    # --- hub exclusion
    deg = np.bincount(np.r_[pi, pj], minlength=n)
    hub = []
    for frac in (0.05, 0.10):
        nh = int(np.ceil(frac * n))
        hubs = np.argsort(-deg, kind="mergesort")[:nh]
        keep = np.ones(n, bool); keep[hubs] = False
        e = keep[pi] & keep[pj]; a, b = pi[e], pj[e]
        o = E.assort_cols(Xc[:, tops], a, b)
        Rh, _, acch = E.run_ms(Xc[:, tops], None, a, b, n, N_hub, seed=E.seed_for(ctx, "hub", frac))
        sm = E.null_summary(o, Rh)
        S.assert_resolution(N_hub, 0.05, 1)
        for j, cell in enumerate(tops):
            hub.append(dict(cell=cell, layer=E.cell_name(cell)[0], sv=E.cell_name(cell)[1], hub_frac=frac,
                            n_hubs_removed=nh, min_hub_degree=int(deg[hubs].min()), m_remaining=int(e.sum()),
                            edge_frac_remaining=float(e.mean()), r_full=float(obs_all[cell]), r_nohub=float(o[j]),
                            null_mean=float(sm["null_mean"][j]), null_sd=float(sm["null_sd"][j]), z=float(sm["z"][j]),
                            p_raw=float(sm["p_raw"][j]), p_adj_top3=float(sm["p_adj"][j]), N=N_hub,
                            acc_mean=float(acch.mean())))
    res["hub"] = hub
    log(f"[{ctx}] hub exclusion done ({time.time()-t0:.0f}s)")
    # --- leave-one-gene-out jackknife of r (top cell)
    x = Xc[:, tops[0]]; r0 = S.assortativity(x, pi, pj)
    q95, q90 = np.quantile(deg, 0.95), np.quantile(deg, 0.90)
    rank = np.empty(n, int); rank[np.argsort(-deg, kind="mergesort")] = np.arange(n)
    jk = np.full(n, np.nan)
    for g in np.flatnonzero(deg > 0):
        e = (pi != g) & (pj != g)
        jk[g] = S.assortativity(x, pi[e], pj[e])
    infl = r0 - jk
    order = np.argsort(-np.nan_to_num(infl, nan=-np.inf), kind="mergesort")[:10]
    res["jackknife"] = [dict(rank=i + 1, gene=c["genes"][g], degree_withinV=int(deg[g]), degree_rank=int(rank[g]) + 1,
                             hub_top5pct=bool(rank[g] < int(np.ceil(0.05 * n))),
                             hub_top10pct=bool(rank[g] < int(np.ceil(0.10 * n))),
                             coord=float(x[g]), r_without_gene=float(jk[g]), influence=float(infl[g]))
                        for i, g in enumerate(order)]
    finite = infl[np.isfinite(infl)]
    res["jackknife_summary"] = dict(cell=tops[0], r_full=float(r0), max_influence=float(finite.max()),
                                    sum_top10_influence=float(infl[order].sum()),
                                    jackknife_se=float(np.sqrt((finite.size - 1) / finite.size *
                                                               ((jk[np.isfinite(jk)] - np.nanmean(jk)) ** 2).sum())),
                                    n_genes_with_edges=int(finite.size), deg_q95=float(q95), deg_q90=float(q90))
    out_p.write_text(json.dumps(res, indent=1, default=float))
    log(f"[{ctx}] pass2 written ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--ctx", required=True); a = ap.parse_args()
    for ctx in a.ctx.split(","):
        run(ctx, E.N_PRIMARY[ctx])
