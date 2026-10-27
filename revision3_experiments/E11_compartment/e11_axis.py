"""E11 part 1b: is it one compartment AXIS? Descriptive follow-up on compartment_tests.csv (universe V).

For each context: the best extracellular-space cell (max sign-invariant AUROC over 12 layers x SV1..SV12), the
AUROC of the other three terms on the same coordinate oriented so that extracellular > 0.5 (cytosol < 0.5 means
cytosolic genes sit at the opposite pole), a gene-bootstrap 95% CI of the extracellular AUROC (B = 2,000, SVD
held fixed: conditional on the estimated direction), and the |Pearson r| between this coordinate and the lung best
cell on matched genes (directions compared only on matched genes, never by index), together with the extracellular
orientation agreement. Also the SV1 extracellular AUROC / p_adj per layer.
Outputs: axis_summary.csv, sv1_by_layer.csv.
"""
import numpy as np, pandas as pd
from e11_common import S, HERE, CTXS, load_emb, load_cache, load_gene2go

B = 2000
SEED = 20261026 + 7


def main():
    ct = pd.read_csv(HERE / "compartment_tests.csv"); v = ct[ct.universe == "V"]
    g2go = load_gene2go(); rng = np.random.default_rng(SEED)
    coord = {}; rows = []
    for ctx in CTXS:
        x = v[(v.context == ctx) & (v.term == "GO:0005615")].sort_values("auroc_signfree").iloc[-1]
        L, k = int(x.layer), int(x.direction[2:]) - 1
        genes, obs_idx, nl, _ = load_cache(ctx); emb = load_emb(ctx)
        US, _, _ = S.svd_observed(emb, L, obs_idx); z = US[:, k]
        sgn = 1.0 if x.auroc_signed >= 0.5 else -1.0; z = z * sgn           # extracellular-high orientation
        coord[ctx] = dict(zip(genes, z))
        y = np.array(["GO:0005615" in g2go.get(g, ()) for g in genes])
        boot = np.empty(B); n = len(genes)
        for b in range(B):
            i = rng.integers(0, n, n); yy = y[i]
            boot[b] = S.auroc(z[i][yy], z[i][~yy]) if 0 < yy.sum() < n else np.nan
        r = dict(context=ctx, layer=L, direction=x.direction, extracellular_auroc=float(S.auroc(z[y], z[~y])),
                 ci_lo=float(np.nanquantile(boot, 0.025)), ci_hi=float(np.nanquantile(boot, 0.975)),
                 p_adj_auroc=float(x.p_adj_auroc), strat_auroc_oriented=float(0.5 + sgn * (x.strat_auroc_signed - 0.5)),
                 spearman_vs_log_nobs_oriented=float(sgn * x.spearman_coord_vs_log_nobs))
        for t, nm in (("GO:0005829", "cytosol"), ("GO:0005739", "mitochondrion"), ("GO:0005788", "ER_lumen")):
            yt = np.array([t in g2go.get(g, ()) for g in genes])
            q = v[(v.context == ctx) & (v.term == t) & (v.layer == L) & (v.direction == x.direction)].iloc[0]
            r[f"{nm}_auroc_oriented"] = float(S.auroc(z[yt], z[~yt])); r[f"{nm}_p_adj_auroc"] = float(q.p_adj_auroc)
        rows.append(r)
    for r in rows:
        a, b = coord["lung"], coord[r["context"]]
        common = sorted(set(a) & set(b))
        rr = np.corrcoef([a[g] for g in common], [b[g] for g in common])[0, 1]
        r.update(n_matched_with_lung=len(common), r_with_lung_best_cell=float(rr))
    # lung-matched cell: in every context, the (layer, SVk) whose coordinate best matches the lung axis on matched genes
    from e11_common import sv_coords
    for r in rows:
        ctx = r["context"]; genes, obs_idx, nl, _ = load_cache(ctx); C, _ = sv_coords(load_emb(ctx), obs_idx)
        a = coord["lung"]; pos = [k for k, g in enumerate(genes) if g in a]; ref = np.array([a[genes[k]] for k in pos])
        cors = np.array([np.corrcoef(ref, C[pos, c])[0, 1] for c in range(C.shape[1])])
        c = int(np.argmax(np.abs(cors))); Lm, km = divmod(c, 12)
        q = v[(v.context == ctx) & (v.term == "GO:0005615") & (v.layer == Lm) & (v.direction == f"SV{km + 1}")].iloc[0]
        qc = v[(v.context == ctx) & (v.term == "GO:0005829") & (v.layer == Lm) & (v.direction == f"SV{km + 1}")].iloc[0]
        sg = np.sign(cors[c])
        r.update(lung_matched_layer=Lm, lung_matched_direction=f"SV{km + 1}", lung_matched_abs_r=float(abs(cors[c])),
                 lung_matched_extracellular_auroc_oriented=float(0.5 + sg * (q.auroc_signed - 0.5)),
                 lung_matched_extracellular_p_adj=float(q.p_adj_auroc),
                 lung_matched_cytosol_auroc_oriented=float(0.5 + sg * (qc.auroc_signed - 0.5)),
                 lung_matched_cytosol_p_adj=float(qc.p_adj_auroc),
                 lung_matched_strat_extracellular_oriented=float(0.5 + sg * (q.strat_auroc_signed - 0.5)))
    pd.DataFrame(rows).to_csv(HERE / "axis_summary.csv", index=False)
    s1 = v[(v.direction == "SV1") & (v.term.isin(["GO:0005615", "GO:0005829"]))][
        ["context", "layer", "term_name", "auroc_signed", "auroc_signfree", "p_adj_auroc", "strat_auroc_signed",
         "p_adj_strat", "OR_top", "OR_bottom", "p_best_pole", "p_adj_fisher", "spearman_coord_vs_log_nobs"]]
    s1.to_csv(HERE / "sv1_by_layer.csv", index=False)
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
