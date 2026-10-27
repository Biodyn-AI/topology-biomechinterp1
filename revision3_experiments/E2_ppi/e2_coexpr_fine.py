"""E2 sensitivity (f, finer): as e2_coexpr.py but matching / stratifying on 50 quantile bins of |co-expression|
instead of deciles (point estimates only). Added after the decile analysis showed residual within-decile
co-expression signal (|co-expression| stratified AUROC 0.53-0.66); labelled sensitivity, no pass/fail.
Same pool construction and seeds as e2_coexpr.py. Writes sens_coexpr_fine.csv.
"""
import numpy as np, pandas as pd
import e2_common as E
from e2_pass2 import top_cells
from e2_coexpr import zcols, pair_corr

S = E.S
RATIO, NBINS = 5, 50


def run(ctx):
    c = E.load_ctx(ctx); n = c["n"]; pi, pj = c["edges"][:, 0].copy(), c["edges"][:, 1].copy(); m = pi.size
    top = top_cells(ctx, E.N_PRIMARY[ctx])[0][0]; Ltop, ktop = E.cell_name(top)
    rng = np.random.default_rng(E.seed_for(ctx, "coexpr_fine"))
    ph = np.load(E.CKPT / f"{ctx}_phys_any.npy")
    linked = set((ph[:, 0] * n + ph[:, 1]).tolist()) | set((ph[:, 1] * n + ph[:, 0]).tolist())
    a = rng.integers(0, n, 40 * m); b = rng.integers(0, n, 40 * m)
    k = a != b; a, b = np.minimum(a[k], b[k]), np.maximum(a[k], b[k])
    key = a * n + b; _, first = np.unique(key, return_index=True); first = np.sort(first)
    a, b, key = a[first], b[first], key[first]
    ok = np.array([q not in linked for q in key.tolist()]); a, b = a[ok][:30 * m], b[ok][:30 * m]
    Xc = E.coords(ctx, c["obs_idx"]); xs = Xc[:, top]
    emb = E.load_emb(ctx)
    Ec = np.asarray(emb[Ltop], np.float64)[c["obs_idx"]]; Ec = S.l2norm_rows(Ec - Ec.mean(0))
    rows = []
    for src in ("extraction", "tier3000"):
        Z, _ = zcols(np.load(E.CACHE / ctx / ("expr_extraction.npy" if src == "extraction" else "expr.npy")))
        ce_pos = np.abs(pair_corr(Z, pi, pj)); ce_pool = np.abs(pair_corr(Z, a, b))
        bins = S.quantile_strata(np.r_[ce_pos, ce_pool], NBINS); bpos, bpool = bins[:m], bins[m:]
        sel, short = [], 0
        for d in np.unique(bpos):
            need = RATIO * int((bpos == d).sum()); cand = np.flatnonzero(bpool == d)
            take = rng.choice(cand, min(need, cand.size), replace=False); short += need - take.size; sel.append(take)
        sel = np.concatenate(sel)
        PA = np.r_[pi, a[sel]]; PB = np.r_[pj, b[sel]]; y = np.r_[np.ones(m, bool), np.zeros(sel.size, bool)]
        strata = np.r_[bpos, bpool[sel]]
        for name, s in (("top_cell_neg_absdiff", -np.abs(xs[PA] - xs[PB])),
                        (f"cosine_L{Ltop}", np.einsum("ij,ij->i", Ec[PA], Ec[PB])),
                        ("abs_coexpr", np.abs(pair_corr(Z, PA, PB)))):
            rows.append(dict(context=ctx, coexpr_source=src, n_bins=NBINS, scorer=name, top_layer=Ltop, top_sv=ktop, m=m,
                             n_matched_neg=int(sel.size), shortfall=int(short),
                             mean_abs_coexpr_pos=float(ce_pos.mean()), mean_abs_coexpr_matched=float(ce_pool[sel].mean()),
                             auroc_matched_pooled=S.auroc(s[y], s[~y]), auroc_matched_stratified=S.stratified_auroc(s, y, strata)))
            print(rows[-1], flush=True)
    return rows


if __name__ == "__main__":
    out = []
    for ctx in E.PRIMARY_CTX:
        out += run(ctx)
    pd.DataFrame(out).to_csv(E.HERE / "sens_coexpr_fine.csv", index=False)
