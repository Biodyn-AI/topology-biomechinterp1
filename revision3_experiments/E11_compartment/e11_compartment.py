"""E11 part 1: compartment axis under the paper convention (CONVENTION3 sec. 2).

For each context in {lung, lung_external, kidney_B, immune, kidney_A}, layer L0..L11 and direction SV1..SV12
of the observed-gene SVD (stats3.svd_observed; centred over V):
  * poles at K = floor(0.27 n) within the universe; for each GO term and each pole a two-sided Fisher exact test
    (stats3.fisher_2x2) of "term-annotated in pole vs outside pole"; OR (Haldane) and p recorded;
  * K-free statistic: sign-invariant AUROC (stats3.signed_free_auroc) of the SV coordinate, term vs non-term;
  * detection-controlled sensitivity: AUROC pooled within detection deciles (cache node_labels det_decile),
    null = label permutation within deciles.
Family (primary) = 12 directions x 12 layers x 2 terms (extracellular space GO:0005615, cytosol GO:0005829);
secondary family = same for mitochondrion GO:0005739 and ER lumen GO:0005788.
Correction: Westfall-Young max-T (stats3.max_t_adjust) over N = 2,000 joint permutations of the annotation
labels across genes: one permutation of the gene axis per draw, applied to all four term labels (term sizes and
term overlaps preserved) and to every direction and layer.
Universes: 'V' (all observed genes; un-annotated genes count as not carrying the term) = primary;
'annotated' (V restricted to genes with >= 1 GO annotation in gene2go_all.pkl; poles recomputed) = sensitivity.
Outputs: compartment_tests.csv, maxT.csv, part1_meta.json.
"""
import json, time
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from e11_common import (S, HERE, CTXS, TERMS, PRIMARY, SECONDARY, POLE_FRAC, N_LAYERS, N_DIRS,
                        load_emb, load_cache, load_gene2go, sv_coords, convention_hash)

N_PERM = 2000
BATCH = 100
SEED = 20261026
TERM_IDS = PRIMARY + SECONDARY


def within_strata_ranks(C, strata):
    R = np.empty_like(C)
    for s in np.unique(strata):
        m = strata == s
        R[m] = rankdata(C[m], axis=0)
    return R


def fisher_table(K, n, n1):
    """(OR, p) of the pole test for every possible in-pole count a = 0..min(K, n1)."""
    amax = min(K, n1)
    ors = np.full(amax + 1, np.nan); ps = np.full(amax + 1, np.nan)
    for a in range(max(0, K - (n - n1)), amax + 1):
        o, p = S.fisher_2x2(a, K - a, n1 - a, (n - K) - (n1 - a))
        ors[a], ps[a] = o, p
    return ors, ps


def run_universe(ctx, uname, C, Y, det, log_nobs, rng_seed):
    n, m = C.shape; T = Y.shape[1]
    K = int(POLE_FRAC * n)
    n1 = Y.sum(0).astype(int); n0 = n - n1
    R = rankdata(C, axis=0)
    Rs = within_strata_ranks(C, det)
    strata = np.unique(det)
    n1s = np.array([[Y[det == s, t].sum() for t in range(T)] for s in strata])         # strata x T
    ns = np.array([(det == s).sum() for s in strata])
    const_s = (n1s * (n1s + 1) / 2).sum(0); den_s = (n1s * (ns[:, None] - n1s)).sum(0)
    lab = np.stack([S.pole_labels(C[:, c], POLE_FRAC) for c in range(m)], 1)
    Top = (lab == 1).astype(float); Bot = (lab == -1).astype(float)
    tabs = [fisher_table(K, n, int(n1[t])) for t in range(T)]

    def stats_from(Yf):                       # Yf: n x (B*T) float, column = b*T + t
        B = Yf.shape[1] // T
        rs = (R.T @ Yf).reshape(m, B, T)
        auc = (rs - (n1 * (n1 + 1) / 2)[None, None, :]) / (n1 * n0)[None, None, :]
        aT = np.rint(Top.T @ Yf).astype(int).reshape(m, B, T)
        aB = np.rint(Bot.T @ Yf).astype(int).reshape(m, B, T)
        pT = np.stack([tabs[t][1][aT[:, :, t]] for t in range(T)], 2)
        pB = np.stack([tabs[t][1][aB[:, :, t]] for t in range(T)], 2)
        fis = -np.log10(np.minimum(pT, pB))
        return auc, fis, aT, aB, pT, pB

    def strat_from(Yf):
        B = Yf.shape[1] // T
        rs = (Rs.T @ Yf).reshape(m, B, T)
        return (rs - const_s[None, None, :]) / den_s[None, None, :]

    Yf = Y.astype(float)
    auc_o, fis_o, aT_o, aB_o, pT_o, pB_o = [x[:, 0, :] for x in stats_from(Yf)]
    sauc_o = strat_from(Yf)[:, 0, :]
    # sanity checks against the shared estimators
    for c, t in [(0, 0), (m - 1, T - 1), (13, 1)]:
        a_ref = S.auroc(C[Y[:, t], c], C[~Y[:, t], c]); assert abs(a_ref - auc_o[c, t]) < 1e-9
        s_ref = S.stratified_auroc(C[:, c], Y[:, t], det); assert abs(s_ref - sauc_o[c, t]) < 1e-9
        a = int(((lab[:, c] == 1) & Y[:, t]).sum())
        assert abs(S.fisher_2x2(a, K - a, n1[t] - a, n - K - (n1[t] - a))[1] - pT_o[c, t]) < 1e-12

    rng = np.random.default_rng(rng_seed); rng_s = np.random.default_rng(rng_seed + 1)
    auc_n = np.empty((N_PERM, m, T)); fis_n = np.empty((N_PERM, m, T)); sauc_n = np.empty((N_PERM, m, T))
    groups = [np.flatnonzero(det == s) for s in strata]
    done = 0
    while done < N_PERM:
        b = min(BATCH, N_PERM - done)
        P = np.stack([rng.permutation(n) for _ in range(b)])                     # b x n
        Yp = Y[P].transpose(1, 0, 2).reshape(n, b * T).astype(float)
        a, f, *_ = stats_from(Yp)
        auc_n[done:done + b] = a.transpose(1, 0, 2); fis_n[done:done + b] = f.transpose(1, 0, 2)
        Ps = np.broadcast_to(np.arange(n), (b, n)).copy()
        for g in groups:
            Ps[:, g] = rng_s.permuted(np.broadcast_to(g, (b, g.size)), axis=1)
        Yps = Y[Ps].transpose(1, 0, 2).reshape(n, b * T).astype(float)
        sauc_n[done:done + b] = strat_from(Yps).transpose(1, 0, 2)
        done += b

    sf = lambda a: 0.5 + np.abs(a - 0.5)                                         # = stats3.signed_free_auroc
    stat = {"auroc_signfree": (sf(auc_o), sf(auc_n)), "fisher_neglog10p": (fis_o, fis_n),
            "strat_auroc_signfree": (sf(sauc_o), sf(sauc_n))}
    rho = np.array([spearmanr(C[:, c], log_nobs).correlation for c in range(m)])
    fam_cols = {"primary": [0, 1], "secondary": [2, 3]}
    padj = {}; praw = {}; maxrows = []
    for fam, tcols in fam_cols.items():
        for sname, (o, nl) in stat.items():
            t_obs = o[:, tcols].T.reshape(-1)                      # term-major: idx = j*m + c
            t_null = nl[:, :, tcols].transpose(0, 2, 1).reshape(N_PERM, -1)
            z, pa, pr = S.max_t_adjust(t_obs, t_null, tail="greater")
            for j, t in enumerate(tcols):
                padj[(sname, t)] = pa[j * m:(j + 1) * m]; praw[(sname, t)] = pr[j * m:(j + 1) * m]
            best = int(np.argmax(z)); j, c = divmod(best, m); L, k = divmod(c, N_DIRS)
            mu = t_null.mean(0); sd = t_null.std(0, ddof=1); sd[sd == 0] = np.inf
            zmax = ((t_null - mu) / sd).max(1)
            q = np.quantile(zmax, [0.5, 0.95, 0.99])
            maxrows.append(dict(context=ctx, universe=uname, family=fam, statistic=sname, N=N_PERM,
                                n_genes=n, K=K, family_size=t_obs.size,
                                best_term=TERM_IDS[tcols[j]], best_term_name=TERMS[TERM_IDS[tcols[j]]],
                                best_layer=L, best_direction=f"SV{k + 1}", best_observed=float(t_obs[best]),
                                best_z=float(z[best]), best_p_raw=float(pr[best]), best_p_adj=float(pa[best]),
                                n_cells_padj_lt_0p05=int((pa < 0.05).sum()),
                                maxnull_z_q50=float(q[0]), maxnull_z_q95=float(q[1]), maxnull_z_q99=float(q[2]),
                                p_floor=1.0 / (N_PERM + 1)))
    rows = []
    for c in range(m):
        L, k = divmod(c, N_DIRS)
        for t, tid in enumerate(TERM_IDS):
            oT, pT = tabs[t][0][aT_o[c, t]], pT_o[c, t]; oB, pB = tabs[t][0][aB_o[c, t]], pB_o[c, t]
            rows.append(dict(
                context=ctx, universe=uname, layer=L, direction=f"SV{k + 1}", term=tid, term_name=TERMS[tid],
                family="primary" if tid in PRIMARY else "secondary", n_genes=n, n_term=int(n1[t]), K=K,
                a_top=int(aT_o[c, t]), OR_top=float(oT), p_top=float(pT),
                a_bottom=int(aB_o[c, t]), OR_bottom=float(oB), p_bottom=float(pB),
                best_pole="top" if pT <= pB else "bottom", p_best_pole=float(min(pT, pB)),
                OR_best_pole=float(oT if pT <= pB else oB),
                auroc_signed=float(auc_o[c, t]), auroc_signfree=float(sf(auc_o[c, t])),
                auroc_null_mean=float(sf(auc_n[:, c, t]).mean()), auroc_null_sd=float(sf(auc_n[:, c, t]).std(ddof=1)),
                p_raw_auroc=float(praw[("auroc_signfree", t)][c]), p_adj_auroc=float(padj[("auroc_signfree", t)][c]),
                p_raw_fisher_perm=float(praw[("fisher_neglog10p", t)][c]),
                p_adj_fisher=float(padj[("fisher_neglog10p", t)][c]),
                strat_auroc_signed=float(sauc_o[c, t]), strat_auroc_signfree=float(sf(sauc_o[c, t])),
                p_raw_strat=float(praw[("strat_auroc_signfree", t)][c]),
                p_adj_strat=float(padj[("strat_auroc_signfree", t)][c]),
                spearman_coord_vs_log_nobs=float(rho[c])))
    return rows, maxrows


def main():
    t0 = time.time()
    g2go = load_gene2go()
    frozen, actual, ok = convention_hash(); assert ok, "CONVENTION3 hash mismatch"
    S.assert_resolution(N_PERM, 0.05, 1)
    all_rows, all_max, meta = [], [], {}
    for ci, ctx in enumerate(CTXS):
        genes, obs_idx, nl, _ = load_cache(ctx)
        emb = load_emb(ctx)
        assert emb.shape[1] >= obs_idx.max() + 1
        C, Ss = sv_coords(emb, obs_idx)
        Y = np.array([[tid in g2go.get(g, ()) for tid in TERM_IDS] for g in genes], bool)
        det = nl["det_decile"]; log_nobs = nl["log_nobs"]
        ann = np.array([len(g2go.get(g, ())) > 0 for g in genes])
        meta[ctx] = dict(n_V=len(genes), n_annotated=int(ann.sum()),
                         term_counts_V={TERMS[t]: int(Y[:, j].sum()) for j, t in enumerate(TERM_IDS)},
                         sv1_var_frac_by_layer=[float(s[0] ** 2 / (s ** 2).sum()) for s in Ss])
        for ui, (uname, mask) in enumerate([("V", np.ones(len(genes), bool)), ("annotated", ann)]):
            rows, mx = run_universe(ctx, uname, C[mask], Y[mask], det[mask], log_nobs[mask],
                                    SEED + 1000 * ci + 100 * ui)
            all_rows += rows; all_max += mx
            best = [r for r in mx if r["family"] == "primary"]
            print(f"[{time.time() - t0:6.0f}s] {ctx:14s} {uname:9s} n={mask.sum()} " +
                  " | ".join(f"{r['statistic']}: {r['best_term_name']} L{r['best_layer']} {r['best_direction']} "
                             f"obs={r['best_observed']:.3f} p_adj={r['best_p_adj']:.4f}" for r in best), flush=True)
        pd.DataFrame(all_rows).to_csv(HERE / "compartment_tests.csv", index=False)   # checkpoint
        pd.DataFrame(all_max).to_csv(HERE / "maxT.csv", index=False)
    meta.update(convention_sha256=frozen, n_perm=N_PERM, seed=SEED, pole_frac=POLE_FRAC,
                runtime_s=time.time() - t0)
    (HERE / "part1_meta.json").write_text(json.dumps(meta, indent=1))
    print("done", time.time() - t0, flush=True)


if __name__ == "__main__":
    main()
