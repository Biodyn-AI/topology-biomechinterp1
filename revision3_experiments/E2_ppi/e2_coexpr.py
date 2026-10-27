"""E2 sensitivity (f): does the STRING-physical enrichment survive holding co-expression fixed?

Positives: STRING-physical >= 700 edges within V. Candidate non-edges: random gene pairs in V with NO physical
STRING link at any score (pool of 30 m pairs, fixed seed). Co-expression = Pearson correlation of log-normalised
expression on (i) the exact extraction cells (cache expr_extraction.npy) and (ii) the 3,000-cell tier (expr.npy).
Matching: deciles of |co-expression| over positives + pool; for every positive 5 non-edges drawn from the same
decile (without replacement). Scorers: top-cell coordinate similarity -|x_i - x_j| (top cell of the primary
family under Maslov-Sneppen), full-dimension cosine of the V-centred embedding at every layer, and |co-expression|
itself. Statistics: unmatched AUROC (vs 5 m random non-edges), matched pooled AUROC, and decile-stratified
AUROC (stats3.stratified_auroc; positives compared only with non-edges of the same |co-expression| decile).
Vertex-subsampling CIs (stats3.vertex_subsample_ci, B=500) for the stratified AUROC of the top-cell score,
the cosine at the top-cell layer and |co-expression|. Writes _ckpt/<ctx>_coexpr.json.
"""
import json, time
import numpy as np
import e2_common as E
from e2_pass2 import top_cells

S = E.S
RATIO = 5


def zcols(M):
    M = np.asarray(M, np.float64); M = M - M.mean(0); sd = M.std(0)
    Z = np.where(sd > 0, M / np.where(sd > 0, sd, 1), 0.0)
    return np.ascontiguousarray(Z.T), int((sd == 0).sum())     # genes x cells


def pair_corr(ZT, a, b):
    """Pearson correlation of gene pairs (a, b) from the standardised genes x cells matrix."""
    nc = ZT.shape[1]; chunk = max(500, int(4e6 // nc))
    out = np.empty(a.size)
    for s in range(0, a.size, chunk):
        out[s:s + chunk] = np.einsum("ij,ij->i", ZT[a[s:s + chunk]], ZT[b[s:s + chunk]]) / nc
    return out


def run(ctx, B=500):
    out_p = E.CKPT / f"{ctx}_coexpr.json"
    log = E.logger(E.HERE / f"e2_coexpr_{ctx}.log")
    if out_p.exists():
        log(f"[{ctx}] exists"); return
    t0 = time.time()
    c = E.load_ctx(ctx); n = c["n"]; pi, pj = c["edges"][:, 0].copy(), c["edges"][:, 1].copy(); m = pi.size
    tops, _, _ = top_cells(ctx, E.N_PRIMARY[ctx]); top = tops[0]; Ltop, ktop = E.cell_name(top)
    rng = np.random.default_rng(E.seed_for(ctx, "coexpr"))
    ph = np.load(E.CKPT / f"{ctx}_phys_any.npy")
    linked = set((ph[:, 0] * n + ph[:, 1]).tolist()) | set((ph[:, 1] * n + ph[:, 0]).tolist())
    a = rng.integers(0, n, 40 * m); b = rng.integers(0, n, 40 * m)
    k = a != b; a, b = np.minimum(a[k], b[k]), np.maximum(a[k], b[k])
    key = a * n + b; _, first = np.unique(key, return_index=True); first = np.sort(first)
    a, b, key = a[first], b[first], key[first]
    ok = np.array([q not in linked for q in key.tolist()])
    a, b = a[ok][:30 * m], b[ok][:30 * m]
    log(f"[{ctx}] pool {a.size} non-edges (no physical link at any score)")
    Xc, Vt, _ = E.coords(ctx, c["obs_idx"], want_vt=True)
    xs = Xc[:, top]
    emb = E.load_emb(ctx)
    res = dict(context=ctx, top_cell=dict(cell=top, layer=Ltop, sv=ktop), m=m, pool=int(a.size), ratio=RATIO, B=B,
               sources={})
    for src in ("extraction", "tier3000"):
        M = np.load(E.CACHE / ctx / ("expr_extraction.npy" if src == "extraction" else "expr.npy"))
        Z, nzero = zcols(M)
        ce_pos = np.abs(pair_corr(Z, pi, pj)); ce_pool = np.abs(pair_corr(Z, a, b))
        dec_all = S.quantile_strata(np.r_[ce_pos, ce_pool], 10)
        dpos, dpool = dec_all[:m], dec_all[m:]
        sel, short = [], 0
        for d in np.unique(dpos):
            need = RATIO * int((dpos == d).sum()); cand = np.flatnonzero(dpool == d)
            take = rng.choice(cand, min(need, cand.size), replace=False); short += need - take.size
            sel.append(take)
        sel = np.concatenate(sel)
        unm = rng.choice(a.size, min(RATIO * m, a.size), replace=False)
        # all pairs: positives then matched negatives
        PA = np.r_[pi, a[sel]]; PB = np.r_[pj, b[sel]]; y = np.r_[np.ones(m, bool), np.zeros(sel.size, bool)]
        strata = np.r_[dpos, dpool[sel]]
        UA = np.r_[pi, a[unm]]; UB = np.r_[pj, b[unm]]; yu = np.r_[np.ones(m, bool), np.zeros(unm.size, bool)]
        scorers = {"top_cell_neg_absdiff": lambda A, Bq: -np.abs(xs[A] - xs[Bq]),
                   "abs_coexpr": lambda A, Bq: np.abs(pair_corr(Z, A, Bq))}
        for L in range(E.N_LAYERS):
            Ec = np.asarray(emb[L], np.float64)[c["obs_idx"]]; Ec = S.l2norm_rows(Ec - Ec.mean(0))
            scorers[f"cosine_L{L}"] = (lambda A, Bq, Ec=Ec: np.einsum("ij,ij->i", Ec[A], Ec[Bq]))
        rows = {}
        for name, fn in scorers.items():
            s_m = fn(PA, PB); s_u = fn(UA, UB)
            rows[name] = dict(auroc_unmatched=S.auroc(s_u[yu], s_u[~yu]), auroc_matched_pooled=S.auroc(s_m[y], s_m[~y]),
                              auroc_matched_stratified=S.stratified_auroc(s_m, y, strata))
            if name in ("top_cell_neg_absdiff", "abs_coexpr", f"cosine_L{Ltop}"):
                def stat(nodes, _a, _b, s_m=s_m):
                    inS = np.zeros(n, bool); inS[nodes] = True; kk = inS[PA] & inS[PB]
                    return S.stratified_auroc(s_m[kk], y[kk], strata[kk])
                est, lo, hi, dr = S.vertex_subsample_ci(stat, n, PA, PB, B, np.random.default_rng(E.seed_for(ctx, src, name)))
                rows[name].update(strat_lo=lo, strat_hi=hi, strat_B=B)
            log(f"  [{ctx}/{src}] {name}: {rows[name]}")
        res["sources"][src] = dict(n_cells=int(M.shape[0]), zero_var_genes=nzero, n_matched_neg=int(sel.size),
                                   shortfall=int(short), mean_abs_coexpr_pos=float(ce_pos.mean()),
                                   mean_abs_coexpr_pool=float(ce_pool.mean()),
                                   mean_abs_coexpr_matched=float(ce_pool[sel].mean()), scorers=rows)
    out_p.write_text(json.dumps(res, indent=1, default=float))
    log(f"[{ctx}] done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    import sys
    for ctx in sys.argv[1].split(","):
        run(ctx)
