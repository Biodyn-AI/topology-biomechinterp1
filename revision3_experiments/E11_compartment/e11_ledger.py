"""E11 part 2: ledger reproduction of "SV1 extracellular OR = 6.37" (screen iter_0006 H02) and of the submitted
retraction (round-2 rebuild_headline.py: OR 0.46 / p 0.13 at L11; best L8 OR 2.41 / p 0.031), then the same
tests under the paper convention (observed-gene SVD, CONVENTION3 sec. 2) in kidney_A and in every E11 context.

Conventions documented here (zero rows are used ONLY to reproduce / document the artefact, as instructed):
  OWN  (iter_0006/run_iter0006_screen.py::run_h02_svd_biology): layer-11 matrix of cycle1_main (= kidney_A,
       float32, all 4,803 rows incl. never-observed zero rows), centred over all rows, np.linalg.svd, SV1 = U[:,0]
       (no S scaling, numpy's arbitrary sign); the 209 genes named in cycle1_edge_dataset.tsv; signed-loading
       quartiles (n25 = 209 // 4 = 52) top vs bottom; GO terms with 5..200 named genes (all GO aspects in
       gene2go_all.pkl); scipy fisher_exact one-sided 'greater' (sample OR).
  R2   (revision2_experiments/rebuild_headline.py, the submitted retraction): SVD of the full 4,803-row matrix
       incl. zero rows (float64), restricted to the curated genes non-zero at L0 (n = 187), K = int(0.27 n),
       each pole vs the rest, two-sided Fisher (Haldane OR), the smaller p of the two poles; 12 layers.
  PAPER: stats3.svd_observed on V; the same two test statistics applied to the named genes in V; plus the
       full-V pole test from part 1 (compartment_tests.csv).
Output: ledger_iter0006.json, ledger_tests.csv.
"""
import csv, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from e11_common import S, X, HERE, CTXS, load_emb, load_cache, load_gene2go, convention_hash

EDGE = X.REPO / "subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
EXT, CYT = "GO:0005615", "GO:0005829"


def named_gene_map():
    """iter_0006 load_gene_index_map(): dict insertion order, later rows overwrite the index."""
    g2i = {}
    with open(EDGE) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g2i[row["source"]] = int(row["source_idx"]); g2i[row["target"]] = int(row["target_idx"])
    return g2i


def own_test(proj, names, gene2go, term, flip=False):
    """iter_0006 quartile test on a signed projection (list order = names)."""
    rows = [dict(gene=g, sv1_loading=float(p), sv1_abs=float(abs(p))) for g, p in zip(names, proj)]
    rows.sort(key=lambda x: -x["sv1_abs"])
    if flip:
        for r in rows: r["sv1_loading"] = -r["sv1_loading"]
    srt = sorted(rows, key=lambda x: -x["sv1_loading"])
    n25 = max(5, len(srt) // 4)
    top = set(r["gene"] for r in srt[:n25]); bot = set(r["gene"] for r in srt[-n25:])
    tg = {g for g in names if g in gene2go and term in gene2go[g]}
    a = len(tg & top); c = len(tg & bot)
    OR, p = fisher_exact([[a, n25 - a], [c, n25 - c]], alternative="greater")
    return dict(n=len(names), n25=n25, a_top=a, c_bottom=c, OR=float(OR), p=float(p), top=top, bot=bot)


def own_term_screen(proj, names, gene2go, flip=False):
    """Full iter_0006 term screen (to confirm GO:0005615 is the top hit, 284 terms)."""
    from collections import defaultdict
    t2g = defaultdict(set)
    for g in names:
        for t in gene2go.get(g, ()): t2g[t].add(g)
    valid = [t for t, gs in t2g.items() if 5 <= len(gs) <= 200]
    r = own_test(proj, names, gene2go, EXT, flip); top, bot, n25 = r["top"], r["bot"], r["n25"]
    out = []
    for t in valid:
        a = len(t2g[t] & top); c = len(t2g[t] & bot)
        if a + c == 0: continue
        OR, p = fisher_exact([[a, n25 - a], [c, n25 - c]], alternative="greater")
        out.append((p, t, OR))
    out.sort()
    return len(out), out[:3]


def pole_test(proj, has_term):
    """R2 / paper pole test: K = int(0.27 n), each pole vs rest, two-sided Fisher; returns both poles."""
    n = proj.size; lab = S.pole_labels(proj, 0.27); K = int(0.27 * n); n1 = int(has_term.sum())
    res = {}
    for name, v in (("top", 1), ("bottom", -1)):
        a = int((has_term & (lab == v)).sum())
        res[name] = S.fisher_2x2(a, K - a, n1 - a, (n - K) - (n1 - a)) + (a,)
    best = min(res, key=lambda k: res[k][1])
    return dict(n=n, K=K, n_term=n1, best_pole=best, OR=res[best][0], p=res[best][1],
                OR_top=res["top"][0], p_top=res["top"][1], a_top=res["top"][2],
                OR_bottom=res["bottom"][0], p_bottom=res["bottom"][1], a_bottom=res["bottom"][2])


def main():
    frozen, actual, ok = convention_hash(); assert ok
    gene2go = load_gene2go()
    g2i = named_gene_map(); names = list(g2i.keys()); nidx = np.array([g2i[g] for g in names])
    embA = load_emb("kidney_A"); genesA_all = X.h5ad_genes(X.CONTEXTS["kidney_A"]["h5ad"])
    assert all(genesA_all[i] == g for g, i in g2i.items()), "edge-file indices disagree with kidney_A gene axis"
    nzA = S.observed_mask(embA)
    out = dict(convention_sha256=frozen, source=dict(
        script="iterations/iter_0006/run_iter0006_screen.py::run_h02_svd_biology",
        record="iterations/iter_0006/executor_hypothesis_screen.json H02; iter0006_results.json H02_svd_biology",
        matrix="subproject_38 .../cycle1_main/layer_gene_embeddings.npy (12 x 4803 x 512) = context kidney_A "
               "(Tabula Sapiens KIDNEY h5ad; the screen labelled it immune)",
        named_genes=str(EDGE), n_named=len(names), n_named_zero_rows=int((~nzA[nidx]).sum())))
    rows = []

    # ---------------- 1. OWN convention, exact reproduction -------------------------------------------------
    X11 = np.asarray(embA[11]); Xc = X11 - X11.mean(axis=0, keepdims=True)          # float32, as in iter_0006
    U, Sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    projA = U[nidx, 0]
    rep = {}
    for flip in (False, True):
        r = own_test(projA, names, gene2go, EXT, flip)
        nt, top3 = own_term_screen(projA, names, gene2go, flip)
        rep["flipped" if flip else "numpy_sign"] = dict(OR=r["OR"], p=r["p"], a_top=r["a_top"], c_bottom=r["c_bottom"],
                                                        n25=r["n25"], n_terms_tested=nt,
                                                        top3=[dict(term=t, OR=o, p=p) for p, t, o in top3],
                                                        zero_rows_in_top=int(sum(not nzA[g2i[g]] for g in r["top"])),
                                                        zero_rows_in_bottom=int(sum(not nzA[g2i[g]] for g in r["bot"])))
    match = [k for k, v in rep.items() if abs(v["OR"] - 6.367741935483871) < 1e-9 and abs(v["p"] - 0.00025746446015579904) < 1e-12]
    obs_ind = nzA.astype(float)
    r_ind = float(np.corrcoef(U[:, 0], obs_ind)[0, 1])
    zero_val = np.unique(np.round(U[~nzA, 0], 7))
    ext_named = np.array([EXT in gene2go.get(g, ()) for g in names]); zr = ~nzA[nidx]
    out["own_convention"] = dict(reproduced=bool(match), matching_orientation=match, orientations=rep,
                                 sv1_var_frac_L11=float(Sv[0] ** 2 / (Sv ** 2).sum()),
                                 r_SV1_vs_observed_indicator_L11=r_ind,
                                 n_distinct_SV1_values_on_zero_rows=int(zero_val.size),
                                 extracellular_rate_zero_row_named=float(ext_named[zr].mean()),
                                 extracellular_rate_observed_named=float(ext_named[~zr].mean()),
                                 n_zero_row_named=int(zr.sum()))
    # the same OWN statistic with the zero-row named genes dropped (SVD unchanged)
    keep = [g for g in names if nzA[g2i[g]]]
    pk = U[[g2i[g] for g in keep], 0]
    for flip in (False, True):
        r = own_test(pk, keep, gene2go, EXT, flip)
        rows.append(dict(test="own", convention="own_zero_row_SVD_named_nonzero_only", context="kidney_A", layer=11,
                         term=EXT, orientation="flipped" if flip else "numpy_sign", n=r["n"], OR=r["OR"], p=r["p"],
                         a_top=r["a_top"], c_bottom=r["c_bottom"]))
    for k, v in rep.items():
        rows.append(dict(test="own", convention="own_zero_row_SVD_209_named", context="kidney_A", layer=11, term=EXT,
                         orientation=k, n=len(names), OR=v["OR"], p=v["p"], a_top=v["a_top"], c_bottom=v["c_bottom"]))

    # ---------------- 2. R2 retraction reproduction (zero-row SVD) -------------------------------------------
    cur = sorted(i for i, g in enumerate(genesA_all) if g in set(names))
    nz0 = np.linalg.norm(np.asarray(embA[0], dtype=np.float64), axis=1) > 0
    sub = np.array([i for i in cur if nz0[i]]); subn = [genesA_all[i] for i in sub]
    r2 = {}
    for term in (EXT, CYT):
        ht = np.array([term in gene2go.get(g, ()) for g in subn])
        per = []
        for L in range(12):
            XL = np.asarray(embA[L], dtype=np.float64); XL = XL - XL.mean(0, keepdims=True)
            UL, SL, _ = np.linalg.svd(XL, full_matrices=False)
            US = (UL * SL)[sub]
            t = pole_test(US[:, 0], ht); per.append(t)
            rows.append(dict(test="pole_vs_rest", convention="R2_zero_row_SVD_curated_nonzero", context="kidney_A",
                             layer=L, term=term, orientation="sign_invariant", n=t["n"], OR=t["OR"], p=t["p"],
                             best_pole=t["best_pole"], K=t["K"]))
            if term == EXT and L == 11:
                r_ind_L = float(np.corrcoef(UL[:, 0], nzA.astype(float))[0, 1])
        b = int(np.argmin([t["p"] for t in per]))
        r2[term] = dict(L11=dict(OR=per[11]["OR"], p=per[11]["p"], pole=per[11]["best_pole"]),
                        best_layer=b, best=dict(OR=per[b]["OR"], p=per[b]["p"], pole=per[b]["best_pole"]), n=len(sub),
                        K=per[0]["K"])
    out["R2_retraction_convention"] = dict(
        source="revision2_experiments/rebuild_headline.py [B]; headline_rebuilt.csv rows SV1_extracellular_L11 / _best_layer",
        submitted_values=dict(L11_OR=0.4644303065355697, L11_p=0.12909553549175648, best_layer=8,
                              best_OR=2.4141201264488936, best_p=0.031046878529548395),
        reproduced=r2, r_SV1_vs_observed_indicator_L11=r_ind_L,
        reproduced_exactly=bool(abs(r2[EXT]["L11"]["OR"] - 0.4644303065355697) < 1e-6 and
                                abs(r2[EXT]["L11"]["p"] - 0.12909553549175648) < 1e-9 and r2[EXT]["best_layer"] == 8 and
                                abs(r2[EXT]["best"]["OR"] - 2.4141201264488936) < 1e-6))

    # ---------------- 3. PAPER convention: kidney_A and every context ----------------------------------------
    paper = {}
    for ctx in CTXS:
        genes, obs_idx, nl, _ = load_cache(ctx); emb = load_emb(ctx)
        pos = {g: k for k, g in enumerate(genes)}
        nm = [g for g in names if g in pos]; ni = np.array([pos[g] for g in nm])
        res = {}
        for L in range(12):
            US, _, _ = S.svd_observed(emb, L, obs_idx)
            for term in (EXT, CYT):
                ht = np.array([term in gene2go.get(g, ()) for g in nm])
                t = pole_test(US[ni, 0], ht)
                rows.append(dict(test="pole_vs_rest", convention="paper_observed_SVD_named_in_V", context=ctx, layer=L,
                                 term=term, orientation="sign_invariant", n=t["n"], OR=t["OR"], p=t["p"],
                                 best_pole=t["best_pole"], K=t["K"]))
                ow = [own_test(US[ni, 0], nm, gene2go, term, f) for f in (False, True)]
                for f, r in zip(("svd_observed_sign", "flipped"), ow):
                    rows.append(dict(test="own", convention="paper_observed_SVD_named_in_V", context=ctx, layer=L,
                                     term=term, orientation=f, n=r["n"], OR=r["OR"], p=r["p"],
                                     a_top=r["a_top"], c_bottom=r["c_bottom"]))
                res[(L, term)] = dict(pole=t, own=min(ow, key=lambda r: r["p"]))
        paper[ctx] = dict(
            n_named_in_V=len(nm),
            L11_SV1_extracellular_pole_vs_rest=dict(OR=res[(11, EXT)]["pole"]["OR"], p=res[(11, EXT)]["pole"]["p"],
                                                     pole=res[(11, EXT)]["pole"]["best_pole"]),
            L11_SV1_extracellular_own_best_orientation=dict(OR=res[(11, EXT)]["own"]["OR"], p=res[(11, EXT)]["own"]["p"]),
            best_layer_extracellular_pole_vs_rest=min(range(12), key=lambda L: res[(L, EXT)]["pole"]["p"]),
            min_p_extracellular_pole_vs_rest=min(res[(L, EXT)]["pole"]["p"] for L in range(12)),
            L11_SV1_cytosol_pole_vs_rest=dict(OR=res[(11, CYT)]["pole"]["OR"], p=res[(11, CYT)]["pole"]["p"]))
    out["paper_convention_named_genes"] = paper
    # full-V paper-convention test at SV1 (from part 1)
    ct = pd.read_csv(HERE / "compartment_tests.csv")
    q = ct[(ct.universe == "V") & (ct.direction == "SV1") & (ct.term.isin([EXT, CYT]))]
    out["paper_convention_full_V_SV1"] = {
        ctx: {t: dict(L11=q[(q.context == ctx) & (q.term == t) & (q.layer == 11)][
            ["OR_best_pole", "p_best_pole", "best_pole", "auroc_signfree", "p_adj_auroc"]].iloc[0].to_dict(),
                      best_layer_by_auroc=int(q[(q.context == ctx) & (q.term == t)].sort_values("auroc_signfree").iloc[-1].layer))
              for t in (EXT, CYT)} for ctx in CTXS}
    pd.DataFrame(rows).to_csv(HERE / "ledger_tests.csv", index=False)
    (HERE / "ledger_iter0006.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: out[k] for k in ("own_convention", "R2_retraction_convention")}, indent=1, default=float)[:4000])
    print(json.dumps(out["paper_convention_named_genes"], indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
