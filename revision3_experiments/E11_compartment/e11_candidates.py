"""E11 part 3: re-scoping of four screen candidates that the submitted paper lists as failed.

(i)   iter_0011 H03 "591 GO terms at SV2 poles, 0 significant" (matrix: cycle1_main = kidney_A).
      (a) own convention re-run with a valid test: SVD of the 209-named-gene submatrix (zero rows included, as in
          the screen), SV2, K = 52 top / bottom, terms with 3..50 named genes (all GO aspects), one-sided Fisher per
          pole (best pole), BH over the term family per layer; plus a check that the screen's shuffle null is
          degenerate (it applies one bijection to pole sets AND the term set, so every overlap is invariant).
      (b) paper convention in lung, immune and kidney_A: observed-gene SVD, SV2 coordinate, poles K = floor(0.27 n)
          within V, GO Biological Process terms only (namespace from go-basic.obo) with 10..500 genes in V,
          two-sided Fisher (stats3.fisher_2x2) at each pole; BH within the (context, layer) family of
          terms x 2 poles, and within the context over all 12 layers.
(iv)  iter_0058 H03 / iter_0059 H03 effective rank vs AUROC across layers, confounded by depth (matrix: true immune
      cycle4_immune main/seed43/seed44, non-zero genes). Paper convention in lung and immune (+ immune_seed43/44 for
      the pooled 36-observation design): ER = stats3.effective_rank of the observed-gene SVD (also the screen's
      s^2-entropy variant); AUROC = (A) H2 statistic: detection-stratified AUROC of the LOO mean cosine to the
      Lambert TF set, full dimension, centred; (B) legacy: LOO mean-cosine AUROC in the SV2-SV7 band, TRRUST TF vs
      target-only. Raw Spearman, partial Spearman on layer (rank residuals), and a depth-free within-layer test on
      the immune triplet (layer fixed effects; permutation of ER across replicates within layer, N = 10,000).
(ii)  iter_0061 H03 FFL intermediacy (matrix: true immune cycle4_immune_main; 264 FFLs). Paper convention in lung (replication: lung_external)
      and immune: TRRUST FFLs A->B, A->C, B->C within V; t = <B-A, C-A>/|C-A|^2 in SV2-4, SV5-7, SV2-7 of each layer;
      statistic frac(0 <= t <= 1) (primary, one-sided greater) and mean t; null = gene-label permutation over V
      (one permutation per draw applied to all 36 cells), N = 2,000, max-T over the 36 cells.
(iii) subproject_40 iter_0006-0008 H07/H09/H12 PH rewiring null (matrix: true immune cycle4 main/seed43/seed44,
      180 genes drawn from the full gene axis incl. zero rows): reported as a screen result only; this script
      reproduces the screen's gene draws to document the fraction of never-observed (zero) rows they contained.
Outputs: further_candidates.csv, cand_i_gobp_sv2.csv, cand_iv_effrank.csv, cand_ii_ffl.csv, cand_iii_ph_draws.csv,
         part3_meta.json.
"""
import csv, json, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, spearmanr, pearsonr
from e11_common import S, X, HERE, load_emb, load_cache, load_gene2go, go_namespaces, convention_hash

SEED = 20261026
EDGE = X.REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
t0 = time.time()
summary = []          # rows for further_candidates.csv
meta = {}


def log(*a):
    print(f"[{time.time() - t0:6.0f}s]", *a, flush=True)


def fisher_lookup(K, n, n1):
    lo = max(0, K - (n - n1)); hi = min(K, n1)
    ors = np.full(hi + 1, np.nan); ps = np.full(hi + 1, np.nan)
    for a in range(lo, hi + 1):
        ors[a], ps[a] = S.fisher_2x2(a, K - a, n1 - a, (n - K) - (n1 - a))
    return ors, ps


# ====================================================================================== (i)
def cand_i():
    g2go = load_gene2go(); ns = go_namespaces()
    rows = []
    # ---- (a) own convention (kidney_A cycle1_main, 209 named genes, submatrix SVD incl. zero rows) ----------
    g2i = {}
    with open(EDGE) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            s_, t_ = r["source"].strip().upper(), r["target"].strip().upper()
            g2i.setdefault(s_, int(r["source_idx"])); g2i.setdefault(t_, int(r["target_idx"]))
    names = list(g2i.keys()); idx = [g2i[g] for g in names]
    embA = load_emb("kidney_A")
    t2g = defaultdict(set)
    for g in names:
        for t in g2go.get(g, ()): t2g[t].add(g)
    valid = {t: gs for t, gs in t2g.items() if 3 <= len(gs) <= 50}
    n_bp_valid = sum(1 for t in valid if ns.get(t) == "biological_process")
    own = {}
    for L in (7, 8, 11):
        Xs = np.asarray(embA[L])[idx]; Xs = Xs - Xs.mean(0)
        U, Sv, _ = np.linalg.svd(Xs, full_matrices=False)
        proj = U[:, 1] * Sv[1]
        order = np.argsort(-proj); top = {names[i] for i in order[:52]}; bot = {names[i] for i in order[-52:]}
        N = len(names); res = []
        for t, gs in valid.items():
            best = (0.0, 1.0)
            for pole in (top, bot):
                a = len(gs & pole); tab = [[a, 52 - a], [len(gs) - a, N - len(gs) - 52 + a]]
                o, p = fisher_exact(tab, alternative="greater")
                if p < best[1] or (p == best[1] and o > best[0]): best = (o, p)
            res.append((t, best[0], best[1]))
        # degenerate-null check: a bijection applied to the pole sets and the term set preserves |term & pole|
        rng = np.random.default_rng(SEED)
        t_chk = max(res, key=lambda r: r[1])[0]
        mp = dict(zip(names, rng.permutation(names)))
        inv = len({mp[g] for g in valid[t_chk]} & {mp[g] for g in top}) == len(valid[t_chk] & top)
        q = S.bh_fdr([r[2] for r in res])
        own[L] = dict(n_terms=len(res), min_p=float(min(r[2] for r in res)), n_p_lt_0p05=int(sum(r[2] < 0.05 for r in res)),
                      n_bh_lt_0p05=int((q < 0.05).sum()), min_q=float(q.min()), shuffle_null_degenerate=bool(inv),
                      top_term=min(res, key=lambda r: r[2])[0])
        log(f"(i-own) L{L}: {own[L]}")
    meta["i_own"] = dict(n_named=len(names), n_terms_3_50=len(valid), n_terms_3_50_BP_only=n_bp_valid, by_layer=own,
                         screen_record="iter_0011/h03_gobp_sv2_poles.json: emp_p = 1.0 and null_std_or = 0.0 for every "
                                       "reported hit at L7, L8, L11")
    for L, v in own.items():
        summary.append(dict(candidate="i_GO_terms_SV2_poles", screen_iteration="iter_0011 H03",
                            screen_matrix="cycle1_main (kidney_A; labelled immune)",
                            screen_convention="SVD of 209-named-gene submatrix incl. 22 zero rows; K=52; terms 3-50 named genes, all GO aspects; shuffle null",
                            rerun_context="kidney_A", rerun_convention="own convention, Fisher one-sided best pole, BH over terms",
                            layer=L, statistic="n terms BH q<0.05", value=v["n_bh_lt_0p05"], n_tests=v["n_terms"],
                            min_p=v["min_p"], p_adj=v["min_q"], note=f"screen null degenerate={v['shuffle_null_degenerate']}"))
    # ---- (b) paper convention -----------------------------------------------------------------------
    for ctx in ("lung", "immune", "kidney_A"):
        genes, obs_idx, nl, _ = load_cache(ctx); emb = load_emb(ctx); n = len(genes)
        tg = defaultdict(list)
        for k, g in enumerate(genes):
            for t in g2go.get(g, ()):
                if ns.get(t) == "biological_process": tg[t].append(k)
        terms = sorted(t for t, ks in tg.items() if 10 <= len(ks) <= 500)
        M = np.zeros((n, len(terms)), bool)
        for j, t in enumerate(terms): M[tg[t], j] = True
        n1 = M.sum(0); K = int(0.27 * n)
        tabs = {int(v): fisher_lookup(K, n, int(v)) for v in np.unique(n1)}
        allp = []
        for L in range(12):
            US, _, _ = S.svd_observed(emb, L, obs_idx)
            lab = S.pole_labels(US[:, 1], 0.27)
            aT = M[lab == 1].sum(0); aB = M[lab == -1].sum(0)
            for j, t in enumerate(terms):
                oT, pT = tabs[int(n1[j])][0][aT[j]], tabs[int(n1[j])][1][aT[j]]
                oB, pB = tabs[int(n1[j])][0][aB[j]], tabs[int(n1[j])][1][aB[j]]
                rows.append(dict(context=ctx, layer=L, term=t, n_term=int(n1[j]), n_genes=n, K=K,
                                 a_top=int(aT[j]), OR_top=oT, p_top=pT, a_bottom=int(aB[j]), OR_bottom=oB, p_bottom=pB))
                allp += [pT, pB]
        df = pd.DataFrame([r for r in rows if r["context"] == ctx])
        q_ctx = S.bh_fdr(np.array(allp)).reshape(-1, 2)
        df["q_top_ctx"] = q_ctx[:, 0]; df["q_bottom_ctx"] = q_ctx[:, 1]
        qt, qb = np.empty(len(df)), np.empty(len(df))
        for L in range(12):
            m = (df.layer == L).values
            qq = S.bh_fdr(np.concatenate([df.p_top[m].values, df.p_bottom[m].values]))
            qt[m] = qq[:m.sum()]; qb[m] = qq[m.sum():]
        df["q_top_layer"] = qt; df["q_bottom_layer"] = qb
        rows = [r for r in rows if r["context"] != ctx] + df.to_dict("records")
        for L in range(12):
            d = df[df.layer == L]
            nsig = int(((d.q_top_layer < 0.05) | (d.q_bottom_layer < 0.05)).sum())
            enr = d[((d.q_top_layer < 0.05) & (d.OR_top > 1)) | ((d.q_bottom_layer < 0.05) & (d.OR_bottom > 1))]
            best = d.assign(pm=np.minimum(d.p_top, d.p_bottom)).sort_values("pm").iloc[0]
            summary.append(dict(candidate="i_GO_terms_SV2_poles", screen_iteration="iter_0011 H03",
                                screen_matrix="cycle1_main (kidney_A; labelled immune)",
                                screen_convention="SVD of 209-named-gene submatrix incl. 22 zero rows; K=52; terms 3-50 named genes, all GO aspects; shuffle null",
                                rerun_context=ctx, rerun_convention="paper: observed-gene SVD SV2, K=0.27n in V, GO BP 10-500 genes, two-sided Fisher, BH within layer (terms x 2 poles)",
                                layer=L, statistic="n terms BH q<0.05 (either pole)", value=nsig, n_tests=2 * len(d),
                                min_p=float(best.pm), p_adj=float(min(d.q_top_layer.min(), d.q_bottom_layer.min())),
                                note=f"n enriched (OR>1) q<0.05: {len(enr)}; top term {best.term}"))
        log(f"(i-paper) {ctx}: {len(terms)} BP terms; per-layer BH hits: "
            + str([s["value"] for s in summary if s["rerun_context"] == ctx and s["candidate"].startswith("i_")]))
        meta.setdefault("i_paper", {})[ctx] = dict(n_terms=len(terms), n_V=n, K=K,
                                                   n_cells_bh_ctx_lt_0p05=int((q_ctx < 0.05).sum()))
    pd.DataFrame(rows).to_csv(HERE / "cand_i_gobp_sv2.csv", index=False)


# ====================================================================================== (iv)
def er_s2(Sv):
    p = Sv ** 2 / (Sv ** 2).sum(); p = p[p > 0]
    return float(np.exp(-(p * np.log(p)).sum()))


def partial_spearman(x, y, z):
    from scipy.stats import rankdata
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    A = np.c_[np.ones_like(rz), rz]
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]; ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    # [verifier fix] if x (or y) is rank-collinear with z (e.g. ER strictly monotone in layer, rho = -1), the
    # residuals are floating-point noise (~1e-15) and their Pearson r is meaningless: the partial correlation is
    # not estimable. Return NaN instead of a noise correlation.
    tol = 1e-8 * np.sqrt(len(rx))
    if np.linalg.norm(ex) < tol or np.linalg.norm(ey) < tol:
        return (float("nan"), float("nan"))
    return pearsonr(ex, ey)


def cand_iv():
    rows = []
    for ctx in ("lung", "immune", "immune_seed43", "immune_seed44"):
        genes, obs_idx, nl, _ = load_cache(ctx); emb = load_emb(ctx)
        tf = np.flatnonzero(nl["is_tf_lambert"]); tr_tf = np.flatnonzero(nl["is_trrust_tf"])
        tr_tg = np.flatnonzero(nl["is_trrust_target_only"]); det = nl["det_decile"]
        is_tf = nl["is_tf_lambert"].astype(bool)
        for L in range(12):
            US, Sv, _ = S.svd_observed(emb, L, obs_idx)
            Xc = np.asarray(emb[L], dtype=np.float64)[obs_idx]; Xc = Xc - Xc.mean(0)
            sc = S.tf_scores_loo(Xc, tf)
            h2 = S.stratified_auroc(sc, is_tf, det)
            leg = S.tf_auroc(US[:, 1:7], tr_tf, tr_tg)
            rows.append(dict(context=ctx, layer=L, eff_rank=S.effective_rank(Sv), eff_rank_s2=er_s2(Sv),
                             h2_strat_auroc_full=h2, legacy_auroc_sv2_7=leg, n_V=len(genes), n_tf_lambert=tf.size,
                             n_trrust_tf=tr_tf.size, n_trrust_target_only=tr_tg.size))
        log(f"(iv) {ctx} done")
    df = pd.DataFrame(rows); df.to_csv(HERE / "cand_iv_effrank.csv", index=False)
    res = {}
    for ctx in ("lung", "immune"):
        d = df[df.context == ctx]
        for er in ("eff_rank", "eff_rank_s2"):
            for au in ("h2_strat_auroc_full", "legacy_auroc_sv2_7"):
                r_raw = spearmanr(d[er], d[au]); r_le = spearmanr(d.layer, d[er]); r_la = spearmanr(d.layer, d[au])
                pr = partial_spearman(d[er].values, d[au].values, d.layer.values.astype(float))
                res[(ctx, er, au)] = dict(rho_raw=r_raw.correlation, rho_layer_er=r_le.correlation,
                                          rho_layer_auroc=r_la.correlation, partial_rho=pr[0], partial_p=pr[1], n=12,
                                          partial_estimable=bool(np.isfinite(pr[0])))
    # pooled immune triplet (the screen's 36-observation design) + depth-free within-layer test
    d = df[df.context.isin(["immune", "immune_seed43", "immune_seed44"])].copy()
    rng = np.random.default_rng(SEED + 4)
    for er in ("eff_rank", "eff_rank_s2"):
        for au in ("h2_strat_auroc_full", "legacy_auroc_sv2_7"):
            r_raw = spearmanr(d[er], d[au]).correlation
            pr = partial_spearman(d[er].values, d[au].values, d.layer.values.astype(float))
            e = d[er] - d.groupby("layer")[er].transform("mean"); a = d[au] - d.groupby("layer")[au].transform("mean")
            r_w = spearmanr(e, a).correlation
            ev = e.values.reshape(-1); lay = d.layer.values
            grp = [np.flatnonzero(lay == L) for L in range(12)]
            null = np.empty(10000)
            for b in range(null.size):
                ep = ev.copy()
                for g in grp: ep[g] = ev[rng.permutation(g)]
                null[b] = spearmanr(ep, a.values).correlation
            p_w, bcount, N = S.emp_p(r_w, null, "two-sided")
            res[("immune_triplet", er, au)] = dict(rho_raw=r_raw, partial_rho=pr[0], partial_p=pr[1], n=36,
                                                   partial_estimable=bool(np.isfinite(pr[0])),
                                                   within_layer_rho=r_w, within_layer_p_perm=p_w, within_b=bcount,
                                                   within_N=N, within_null_sd=float(null.std(ddof=1)))
    for (ctx, er, au), v in res.items():
        summary.append(dict(candidate="iv_effrank_vs_auroc", screen_iteration="iter_0058 H03 / iter_0059 H03",
                            screen_matrix="cycle4_immune main/seed43/seed44 (true immune), non-zero genes [2039]",
                            screen_convention="ER on s^2 entropy; AUROC = 5-fold LR on SV2-7, TRRUST TF vs target; 36 (seed,layer) obs",
                            rerun_context=ctx, rerun_convention=f"paper: observed-gene SVD; {er}; {au}",
                            layer="all", statistic="partial Spearman (layer)", value=v["partial_rho"], n_tests=v["n"],
                            min_p=v["partial_p"], p_adj=v.get("within_layer_p_perm", np.nan),
                            note=json.dumps({k: (round(float(x), 4) if isinstance(x, (float, np.floating)) else x)
                                             for k, x in v.items()})))
    meta["iv"] = {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in res.items()}
    log("(iv) summary written")


# ====================================================================================== (ii)
def cand_ii():
    """Two nulls, each N = 2,000 joint draws over the 36 cells: 'label' = gene-label permutation over V (the screen's
    null); 'role' = permutation within TRRUST-TF genes and within all other genes separately, so that A and B stay
    TFs and TF-identity clustering (H2) cannot by itself produce intermediacy."""
    N_PERM = 2000
    rows = []
    for ci, ctx in enumerate(("lung", "immune", "lung_external")):
        genes, obs_idx, nl, d = load_cache(ctx); emb = load_emb(ctx); n = len(genes)
        e = np.load(d / "edges_trrust.npy"); out = defaultdict(set)
        for a, b in e: out[int(a)].add(int(b))
        F = np.array(sorted((a, b, c) for a in out for b in out[a] if b in out and b != a
                            for c in (out[a] & out[b]) if c not in (a, b)), dtype=np.int64)
        subs = {"SV2_4": slice(1, 4), "SV5_7": slice(4, 7), "SV2_7": slice(1, 7)}
        coords = []
        for L in range(12):
            US, _, _ = S.svd_observed(emb, L, obs_idx)
            for nm, sl in subs.items(): coords.append((L, nm, US[:, sl].copy()))

        def stat(Z, A, B, C):
            ca = Z[C] - Z[A]; t = ((Z[B] - Z[A]) * ca).sum(-1) / (ca * ca).sum(-1)
            return np.mean((t >= 0) & (t <= 1), -1), np.mean(t, -1), np.median(t, -1)

        A, B, C = F[:, 0], F[:, 1], F[:, 2]
        obs = np.array([stat(Z, A, B, C) for _, _, Z in coords])            # cells x 3
        is_tf = nl["is_trrust_tf"].astype(bool); g_tf = np.flatnonzero(is_tf); g_ot = np.flatnonzero(~is_tf)
        res = {}
        for nk, null_name in enumerate(("label", "role")):
            rng = np.random.default_rng(SEED + 10 * ci + nk)
            nullf = np.empty((N_PERM, len(coords))); nullt = np.empty((N_PERM, len(coords)))
            for b in range(N_PERM):
                if null_name == "label":
                    p = rng.permutation(n)
                else:
                    p = np.arange(n); p[g_tf] = g_tf[rng.permutation(g_tf.size)]; p[g_ot] = g_ot[rng.permutation(g_ot.size)]
                Ap, Bp, Cp = p[A], p[B], p[C]
                for j, (_, _, Z) in enumerate(coords):
                    f_, m_, _ = stat(Z, Ap, Bp, Cp); nullf[b, j] = f_; nullt[b, j] = m_
            zf, paf, prf = S.max_t_adjust(obs[:, 0], nullf, "greater")
            zt, pat, prt = S.max_t_adjust(obs[:, 1], nullt, "two-sided")
            res[null_name] = (nullf, nullt, paf, prf, pat, prt)
        for j, (L, nm, _) in enumerate(coords):
            r = dict(context=ctx, layer=L, subspace=nm, n_ffl=len(F), frac_between=obs[j, 0], mean_t=obs[j, 1],
                     median_t=obs[j, 2], N=N_PERM)
            for null_name, (nullf, nullt, paf, prf, pat, prt) in res.items():
                r.update({f"frac_null_mean_{null_name}": nullf[:, j].mean(), f"frac_null_sd_{null_name}": nullf[:, j].std(ddof=1),
                          f"p_raw_frac_{null_name}": prf[j], f"p_adj_frac_{null_name}": paf[j],
                          f"mean_t_null_mean_{null_name}": nullt[:, j].mean(), f"mean_t_null_sd_{null_name}": nullt[:, j].std(ddof=1),
                          f"p_raw_mean_t_2s_{null_name}": prt[j], f"p_adj_mean_t_2s_{null_name}": pat[j]})
            rows.append(r)
        dd = pd.DataFrame([r for r in rows if r["context"] == ctx])
        for null_name in ("label", "role"):
            bf = dd.sort_values(f"p_adj_frac_{null_name}").iloc[0]; bt = dd.sort_values(f"p_adj_mean_t_2s_{null_name}").iloc[0]
            summary.append(dict(candidate="ii_FFL_intermediacy", screen_iteration="iter_0061 H03 (lineage iter_0060 H03)",
                                screen_matrix="cycle4_immune_main (true immune), non-zero genes; 264 FFLs (226 valid)",
                                screen_convention="t-interpolation in SV2-4/SV5-7/SV2-7; gene-position permutation N=2000 at peak layers only; statistic mean t",
                                rerun_context=ctx,
                                rerun_convention=f"paper: observed-gene SVD; 12 layers x 3 subspaces; {null_name} permutation N=2000; max-T over 36 cells",
                                layer=int(bf.layer), statistic=f"frac(0<=t<=1) best cell {bf.subspace}", value=float(bf.frac_between),
                                n_tests=len(dd), min_p=float(bf[f"p_raw_frac_{null_name}"]), p_adj=float(bf[f"p_adj_frac_{null_name}"]),
                                note=f"null mean {bf[f'frac_null_mean_{null_name}']:.4f}; n_ffl={len(F)}; best mean-t cell L{int(bt.layer)} "
                                     f"{bt.subspace} mean_t={bt.mean_t:.3f} (null {bt[f'mean_t_null_mean_{null_name}']:.3f}) "
                                     f"p_adj_2s={bt[f'p_adj_mean_t_2s_{null_name}']:.4f}"))
            log(f"(ii) {ctx} {null_name}: {len(F)} FFLs; best frac cell L{int(bf.layer)} {bf.subspace} {bf.frac_between:.3f} "
                f"vs {bf[f'frac_null_mean_{null_name}']:.3f} p_adj={bf[f'p_adj_frac_{null_name}']:.4f}")
    pd.DataFrame(rows).to_csv(HERE / "cand_ii_ffl.csv", index=False)


# ====================================================================================== (iii)
def cand_iii():
    """Reproduce subproject_40 iter_0006 gene draws: run_seed = 620000 + seed_index*10000 + split_index*1000 + layer*10,
    candidate pool = first / second half of the 4,941-gene axis, 180 genes, rng.choice(pool, 180, replace=False)."""
    rows = []
    for si, ctx in enumerate(("immune", "immune_seed43", "immune_seed44")):
        genes, obs_idx, nl, _ = load_cache(ctx); G = 4941
        obs = np.zeros(G, bool); obs[obs_idx] = True
        for spi, pool in enumerate((np.arange(0, G // 2), np.arange(G // 2, G))):
            for L in range(12):
                rng = np.random.default_rng(620000 + 0 * 100000 + si * 10000 + spi * 1000 + L * 10)
                idx = rng.choice(pool, size=180, replace=False)
                rows.append(dict(context=ctx, split=("source_disjoint", "target_disjoint")[spi], layer=L,
                                 n_drawn=180, n_zero_rows=int((~obs[idx]).sum()),
                                 frac_zero_rows=float((~obs[idx]).mean()), pool_frac_zero=float((~obs[pool]).mean())))
    df = pd.DataFrame(rows); df.to_csv(HERE / "cand_iii_ph_draws.csv", index=False)
    summary.append(dict(candidate="iii_PH_rewiring", screen_iteration="subproject_40 iter_0006 H07 (+iter_0007 H09, iter_0008 H12)",
                        screen_matrix="cycle4_immune main/seed43/seed44 (true immune); 180 genes drawn from the FULL 4,941-gene axis",
                        screen_convention="PCA(14) of 180 drawn genes; H1 lifetime sum (Euclidean) vs degree-preserving kNN-rewired geodesic H1; 8 null draws; Fisher across seeds; 24 = 12 layers x 2 splits",
                        rerun_context="none (screen result only)", rerun_convention="not re-run; draws reproduced to count zero rows",
                        layer="all", statistic="mean fraction of never-observed (zero) rows among drawn genes",
                        value=float(df.frac_zero_rows.mean()), n_tests=len(df), min_p=np.nan, p_adj=np.nan,
                        note=f"range {df.frac_zero_rows.min():.3f}-{df.frac_zero_rows.max():.3f}; zero rows are identical points "
                             f"(coincident after centring), which collapse into one vertex cluster before PH"))
    meta["iii"] = dict(mean_frac_zero=float(df.frac_zero_rows.mean()), min=float(df.frac_zero_rows.min()),
                       max=float(df.frac_zero_rows.max()))
    log(f"(iii) zero-row fraction in screen draws: {meta['iii']}")


if __name__ == "__main__":
    frozen, actual, ok = convention_hash(); assert ok
    meta["convention_sha256"] = frozen
    which = sys.argv[1:] or ["iii", "i", "iv", "ii"]
    prefix = {"i": "i_", "ii": "ii_", "iii": "iii_", "iv": "iv_"}
    fc, mp = HERE / "further_candidates.csv", HERE / "part3_meta.json"
    keep = pd.read_csv(fc) if fc.exists() else pd.DataFrame()          # rows of parts not re-run are kept
    if len(keep): keep = keep[~keep.candidate.str.split("_").str[0].isin(which)]
    old_meta = json.loads(mp.read_text()) if mp.exists() else {}
    for w in which:
        {"i": cand_i, "ii": cand_ii, "iii": cand_iii, "iv": cand_iv}[w]()
        pd.concat([keep, pd.DataFrame(summary)], ignore_index=True).to_csv(fc, index=False)   # checkpoint
        mp.write_text(json.dumps({**old_meta, **meta}, indent=1, default=float))
    log("done")
