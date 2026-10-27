"""Build Supplementary Table S2 (compartment axis) from revision3_experiments/E11_compartment.

Reads only CSV outputs written by the E11 scripts. Quantities computed here are trivial derived ones:
  * z of a test cell = (sign-free AUROC - mean of its permutation null) / SD of its null, from the null mean and SD
    stored in compartment_tests.csv (this is the standardisation used by stats3.max_t_adjust; checked below against
    the best_z recorded in maxT.csv);
  * detection-stratified AUROC oriented like the raw AUROC, 0.5 + sign(A - 0.5) * (A_strat - 0.5);
  * the enriched / opposite Fisher pole (the larger / smaller of OR_top and OR_bottom), ranges and maxima.
Writes  sagmb_revision/supplement/S2.tex                         (LaTeX fragment, one longtable whose rows are the panels)
        sagmb_revision/supplementary_tables/Table_S2_compartment.csv     (section-stacked values)
Run:    /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python sagmb_revision/supplement/build_S2.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REV = HERE.parent
ROOT = REV.parent
E11 = ROOT / "revision3_experiments" / "E11_compartment"
OUT_TEX = HERE / "S2.tex"
OUT_CSV = REV / "supplementary_tables" / "Table_S2_compartment.csv"

CORE = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
PROBE_CTX = ["lung", "lung_external", "kidney_B", "immune"]
LABEL = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune",
         "kidney_A": "kidney~A"}
PLAIN = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney B", "immune": "immune",
         "kidney_A": "kidney A"}
TERMS = [("GO:0005615", "extracellular_space", "extracellular"), ("GO:0005829", "cytosol", "cytosol")]
COMP_DIR = {"lung": "SV1", "lung_external": "SV1", "kidney_B": "SV1", "immune": "SV2", "kidney_A": "SV1"}
REPS = [("scGPT_L0", r"scGPT $L_0$"), ("scGPT_L3", r"scGPT $L_3$"), ("scGPT_L11", r"scGPT $L_{11}$"),
        ("scGPT_static_token", "scGPT input-token embedding"), ("ESM2_15B", "ESM2-15B per-protein embedding")]
RT = dict(float_precision="round_trip")

ct = pd.read_csv(E11 / "compartment_tests.csv", **RT)
mx = pd.read_csv(E11 / "maxT.csv", **RT)
ax = pd.read_csv(E11 / "axis_summary.csv", **RT)
sb = pd.read_csv(E11 / "static_baselines.csv", **RT)
sp = pd.read_csv(E11 / "static_baselines_paired.csv", **RT)
FLOOR = float(mx.p_floor.iloc[0])
assert abs(FLOOR - 1 / 2001) < 1e-12

ct["z_maxT"] = (ct.auroc_signfree - ct.auroc_null_mean) / ct.auroc_null_sd
ct["strat_auroc_oriented"] = 0.5 + np.sign(ct.auroc_signed - 0.5) * (ct.strat_auroc_signed - 0.5)
# z check: the best cell of each (context, universe, primary family, sign-free AUROC) must reproduce maxT.csv best_z
for r in mx[(mx.family == "primary") & (mx.statistic == "auroc_signfree")].itertuples():
    q = ct[(ct.context == r.context) & (ct.universe == r.universe) & (ct.term == r.best_term)
           & (ct.layer == r.best_layer) & (ct.direction == r.best_direction)].iloc[0]
    assert abs(q.z_maxT - r.best_z) < 1e-9 and abs(q.auroc_signfree - r.best_observed) < 1e-12, r.context


def one(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).to_numpy()
    out = df[m]
    assert len(out) == 1, (kw, len(out))
    return out.iloc[0]


# ------------------------------------------------------------------------------------------ Panel A: best cells
V = ct[(ct.universe == "V") & (ct.family == "primary")]
best = []
for c in CORE:
    for tid, tname, tlab in TERMS:
        g = V[(V.context == c) & (V.term == tid)]
        assert len(g) == 144
        b = g.loc[g.auroc_signfree.idxmax()]
        an = one(ct, context=c, universe="annotated", term=tid, layer=b.layer, direction=b.direction)
        enr_top = b.OR_top >= b.OR_bottom
        d = dict(context=c, universe="V", term=tid, term_name=tname, layer=int(b.layer), direction=b.direction,
                 n_genes=int(b.n_genes), n_term=int(b.n_term), auroc_signed=b.auroc_signed,
                 auroc_signfree=b.auroc_signfree, auroc_null_mean=b.auroc_null_mean, auroc_null_sd=b.auroc_null_sd,
                 z_maxT=b.z_maxT, p_adj_auroc=b.p_adj_auroc, strat_auroc_signed=b.strat_auroc_signed,
                 strat_auroc_oriented=b.strat_auroc_oriented, p_adj_strat=b.p_adj_strat,
                 spearman_coord_vs_log_nobs=b.spearman_coord_vs_log_nobs,
                 OR_top=b.OR_top, p_top=b.p_top, OR_bottom=b.OR_bottom, p_bottom=b.p_bottom,
                 OR_enriched_pole=b.OR_top if enr_top else b.OR_bottom,
                 p_enriched_pole=b.p_top if enr_top else b.p_bottom,
                 OR_opposite_pole=b.OR_bottom if enr_top else b.OR_top,
                 annotated_n_genes=int(an.n_genes), annotated_auroc_signfree=an.auroc_signfree,
                 annotated_p_adj_auroc=an.p_adj_auroc,
                 maxnull_z_q95=one(mx, context=c, universe="V", family="primary",
                                   statistic="auroc_signfree").maxnull_z_q95)
        if tid == "GO:0005615":
            a = one(ax, context=c)
            assert a.layer == b.layer and a.direction == b.direction
            assert abs(a.extracellular_auroc - b.auroc_signfree) < 1e-9
            d.update(ci_lo=a.ci_lo, ci_hi=a.ci_hi,
                     interval="95% gene-bootstrap percentile interval, B=2000, direction held fixed")
        best.append(d)
best = pd.DataFrame(best)

# ------------------------------------------------------------------------------------------ Panel B: by layer
bylayer = ct[(ct.universe == "V") & (ct.term.isin(["GO:0005615", "GO:0005829"]))
             & (((ct.context != "immune") & (ct.direction == "SV1"))
                | ((ct.context == "immune") & (ct.direction.isin(["SV1", "SV2"]))))].copy()
bylayer = bylayer.sort_values(["context", "direction", "term", "layer"])
# cross-check against sv1_by_layer.csv
s1 = pd.read_csv(E11 / "sv1_by_layer.csv", **RT)
chk = bylayer[bylayer.direction == "SV1"].merge(s1, on=["context", "layer", "term_name"], suffixes=("", "_s1"))
assert len(chk) == 120 and (chk.auroc_signfree - chk.auroc_signfree_s1).abs().max() < 1e-12

# ------------------------------------------------------------------------------------------ CSV
sec = []
f1 = best.copy(); f1.insert(0, "section", "A_best_cell_per_term"); sec.append(f1)
f2 = bylayer.copy(); f2.insert(0, "section", "B_compartment_direction_by_layer"); sec.append(f2)
GOID = {n: t for t, n, _ in TERMS}
f3 = sb.rename(columns={"term": "term_name"}); f3.insert(1, "term", f3.term_name.map(GOID))
f3.insert(0, "section", "C_probe_and_own_svd")
f3["interval"] = "probe_lo/probe_hi: 95% gene-bootstrap percentile interval of the out-of-fold AUROC, B=2000"
sec.append(f3)
f4 = sp.rename(columns={"lo": "delta_lo", "hi": "delta_hi", "term": "term_name"})
f4.insert(1, "term", f4.term_name.map(GOID))
f4.insert(0, "section", "C_probe_paired_difference")
f4["interval"] = "delta_lo/delta_hi: 95% paired gene-bootstrap percentile interval, B=2000 (same resamples)"
sec.append(f4)
f5 = ax.copy(); f5.insert(0, "section", "D_axis_cell_other_terms_and_replication")
f5["interval"] = "ci_lo/ci_hi: 95% gene-bootstrap percentile interval of extracellular AUROC, B=2000"
sec.append(f5)
f6 = mx.copy(); f6.insert(0, "section", "maxT_family_best_cells"); sec.append(f6)
f7 = ct.copy(); f7.insert(0, "section", "all_tests_both_universes"); sec.append(f7)
out = pd.concat(sec, ignore_index=True, sort=False)
out.insert(2, "extraction", out.context.map(PLAIN))
out["source"] = out.section.map({
    "A_best_cell_per_term": "E11_compartment/compartment_tests.csv; axis_summary.csv (CI); maxT.csv (null q95)",
    "B_compartment_direction_by_layer": "E11_compartment/compartment_tests.csv (= sv1_by_layer.csv for SV1)",
    "C_probe_and_own_svd": "E11_compartment/static_baselines.csv",
    "C_probe_paired_difference": "E11_compartment/static_baselines_paired.csv",
    "D_axis_cell_other_terms_and_replication": "E11_compartment/axis_summary.csv",
    "maxT_family_best_cells": "E11_compartment/maxT.csv",
    "all_tests_both_universes": "E11_compartment/compartment_tests.csv"})
for col in ("layer", "n_genes", "n_term", "K", "a_top", "a_bottom", "n_pos", "best_sv", "n_matched_with_lung",
            "lung_matched_layer", "N", "family_size", "best_layer", "n_cells_padj_lt_0p05", "annotated_n_genes"):
    v = pd.to_numeric(out[col], errors="raise")
    assert ((v.dropna() % 1) == 0).all(), col
    out[col] = v.astype("Int64")
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_CSV, index=False)


# ------------------------------------------------------------------------------------------ LaTeX helpers
def f(x, d=3):
    return f"{x:,.{d}f}".replace(",", "{,}")


def pa(p):
    return f"{p:.4f}"


def sci(p):
    m, e = f"{p:.1e}".split("e")
    return rf"${m}\times10^{{{int(e)}}}$"


def cell(L, d):
    return rf"$L_{{{int(L)}}}$ $\SV{{{d[2:]}}}$"




q95 = {c: one(mx, context=c, universe="V", family="primary", statistic="auroc_signfree").maxnull_z_q95 for c in CORE}
ann_n = {c: int(best[(best.context == c)].annotated_n_genes.iloc[0]) for c in CORE}
nV = {c: int(best[(best.context == c)].n_genes.iloc[0]) for c in CORE}
strat_not_floor = best[best.p_adj_strat > FLOOR + 1e-12]
assert len(strat_not_floor) == 1
snf = strat_not_floor.iloc[0]
null_means = ct[(ct.universe == "V") & (ct.term == "GO:0005615")
                & ct.apply(lambda r: r.direction == COMP_DIR[r.context], axis=1)].auroc_null_mean
ext = sb[sb.term == "extracellular_space"]
max_strat_diff = (ext.probe_auroc - ext.probe_strat_auroc).abs().max()
ngenes = {c: (int(one(ext, context=c, representation="scGPT_L3").n_genes),
              int(one(ext, context=c, representation="scGPT_L3").n_pos)) for c in PROBE_CTX}
assert (best.annotated_p_adj_auroc <= FLOOR + 1e-12).all()

cap = (
    r"\caption[The early compartment axis]{\textbf{The early compartment axis.} Each of the first twelve directions "
    r"($\SV{1}$--$\SV{12}$) of the observed-gene decomposition at each layer, tested for separating genes directly annotated "
    r"to extracellular space (GO:0005615) or cytosol (GO:0005829) from all other observed genes. AUROC, sign-free AUROC "
    r"$\max(A, 1 - A)$; $z$, AUROC standardised by the mean and SD of its permutation null; $p_{\mathrm{adj}}$, "
    r"Westfall--Young max-$T$ adjusted $p$ over the 288 tests of one extraction (12 directions $\times$ 12 layers "
    r"$\times$ 2 terms; 2{,}000 joint permutations of the gene labels), 0.0005 $= 1/2{,}001$ being its floor. "
    r"\textbf{A}, the cell (layer, direction) with the largest AUROC for each term; $n$, observed genes annotated to the "
    r"term; under each "
    r"extraction, $|V|$ / genes with any GO annotation. CI, 95\% gene-bootstrap percentile interval ($B = 2{,}000$, "
    r"direction held fixed), extracellular cells only. 95th percentile of the null maximum of $z$: "
    + ", ".join(f(q95[c], 2) for c in CORE) + r" (row order). "
    r"Strat., AUROC within deciles of $\log n_{\mathrm{obs}}$, oriented like the AUROC (own max-$T$ family, labels permuted "
    r"within deciles; $p_{\mathrm{adj}}$ at the floor "
    rf"except kidney~A cytosol, {pa(snf.p_adj_strat)}); $|\rho|$, absolute Spearman correlation of the coordinate with "
    r"$\log n_{\mathrm{obs}}$; Annot., AUROC of the same cell among annotated genes only ($p_{\mathrm{adj}}$ at the floor); "
    r"OR, Fisher odds ratio (Haldane) for the term in the 27\% pole where it is enriched (two-sided $p$) and in the "
    r"opposite pole. "
    r"\textbf{B}, extracellular AUROC of the compartment direction ($\SV{1}$; $\SV{2}$ in immune) at every layer; for "
    r"kidney~A the row shows $\SV{1}$, which is not the lung compartment direction ($r$ and lung-matched cell in D); "
    r"$p_{\mathrm{adj}}$ as in A; each layer is decomposed separately; permutation-null mean "
    rf"{f(null_means.min())}--{f(null_means.max())}. "
    r"\textbf{C}, extracellular space on one common gene set per extraction (observed genes in the scGPT vocabulary and "
    r"the ESM2 file; $n$, with extracellular genes in parentheses). Probe, out-of-fold AUROC of L2 logistic regression on "
    r"the first 256 principal components (5 outer folds, regularisation by inner 5-fold cross-validation), with 95\% "
    r"gene-bootstrap percentile interval ($B = 2{,}000$ resamples of the out-of-fold scores); $\Delta$, paired difference "
    r"from scGPT $L_3$ on the same resamples; detection-stratified probe AUROCs differ from these by "
    rf"$\leq {f(max_strat_diff)}$ (cytosol probes in the CSV). Own-SVD, sign-free AUROC of $\SV{{1}}$ of the "
    r"representation's own decomposition of these genes, and the best of its first twelve directions (index; optimistic by "
    r"construction). Input-token embedding, scGPT's gene-embedding table, before any cell is processed; ESM2-15B, the "
    r"5{,}120-dimensional per-protein embedding of the 15-billion-parameter protein language model ESM2. "
    r"\textbf{D}, on the extracellular cell of A, oriented so that extracellular genes score high: AUROC of cytosol, "
    r"mitochondrion (GO:0005739) and endoplasmic-reticulum lumen (GO:0005788) against all other genes (below 0.5, the "
    r"opposite end; $p_{\mathrm{adj}}$ for cytosol as in A, for the other two over the 288 tests of these two terms, "
    r"same permutations); $r$, Pearson correlation with the lung coordinate on the $n$ shared "
    r"genes; lung-matched cell, the extraction's cell most correlated with the lung coordinate ($|r|$), and its "
    r"extracellular AUROC oriented like the lung coordinate. "
    r"Full values: \texttt{Table\_S2\_compartment.csv}.}\label{tab:S2}\\"
)

L = []
A = L.append
TW = r"\textwidth"


def begin_panel(spec, ncol, title):
    A(r"\scriptsize\setlength{\tabcolsep}{2.4pt}\renewcommand{\arraystretch}{1.05}%")
    A(r"\begin{tabular*}{%s}{@{\extracolsep{\fill}}%s@{}}" % (TW, spec))
    A(r"\multicolumn{%d}{@{}l}{\textbf{%s}}\\" % (ncol, title))
    A(r"\midrule")


def end_panel(last=False):
    A(r"\end{tabular*}\\")
    if not last:
        A(r"\midrule")


A(r"\setcounter{table}{1}")
A(r"\begingroup")
A(r"\footnotesize")
A(r"\setlength{\LTleft}{0pt}\setlength{\LTright}{0pt}")
# one-column longtable: each panel is one unbreakable row holding its own tabular, so the table may break between panels
A(r"\begin{longtable}{@{}p{\textwidth}@{}}")
A(cap)
A(r"\toprule")
A(r"\endfirsthead")
A(r"\multicolumn{1}{@{}l}{\scriptsize\itshape Table~\thetable\ (continued)}\\")
A(r"\toprule")
A(r"\endhead")
A(r"\multicolumn{1}{r@{}}{\scriptsize\itshape continued on next page}\\")
A(r"\endfoot")
A(r"\bottomrule")
A(r"\endlastfoot")

# ---- Panel A
begin_panel("lllrcrrrrrrr", 12, r"A\quad Best cell for each term (288 tests per extraction)")
A(r"Extraction & Term ($n$) & Cell & AUROC & 95\% CI & $z$ & $p_{\mathrm{adj}}$ & Strat. & $|\rho|$ & Annot. & "
  r"OR, enriched pole ($p$) & OR, opposite \\")
A(r"\midrule")
for i, c in enumerate(CORE):
    for j, (tid, tname, tlab) in enumerate(TERMS):
        b = best[(best.context == c) & (best.term == tid)].iloc[0]
        au = f(b.auroc_signfree) + " & " + (f"{f(b.ci_lo)}--{f(b.ci_hi)}" if tid == "GO:0005615" else "--")
        first = LABEL[c] if j == 0 else f"{f(nV[c], 0)} / {f(ann_n[c], 0)}"
        A(f"{first} & {tlab} ({f(b.n_term, 0)}) & {cell(b.layer, b.direction)} & {au} & "
          f"{f(b.z_maxT, 1)} & {pa(b.p_adj_auroc)} & {f(b.strat_auroc_oriented)} & "
          f"{f(abs(b.spearman_coord_vs_log_nobs), 2)} & {f(b.annotated_auroc_signfree)} & "
          f"{f(b.OR_enriched_pole, 2)} ({sci(b.p_enriched_pole)}) & {f(b.OR_opposite_pole, 2)} \\\\")
    if i < len(CORE) - 1:
        A(r"\addlinespace[1.5pt]")
end_panel()

# ---- Panel B
begin_panel("ll*{12}{r}", 14, r"B\quad Extracellular AUROC of the compartment direction at every layer")
A(r"Extraction & & " + " & ".join(rf"$L_{{{l}}}$" for l in range(12)) + r" \\")
A(r"\midrule")
for i, c in enumerate(CORE):
    g = ct[(ct.context == c) & (ct.universe == "V") & (ct.term == "GO:0005615")
           & (ct.direction == COMP_DIR[c])].sort_values("layer")
    assert list(g.layer) == list(range(12))
    A(f"{LABEL[c]} & AUROC & " + " & ".join(f(x) for x in g.auroc_signfree) + r" \\")
    A(rf"\quad $\SV{{{COMP_DIR[c][2:]}}}$ & $p_{{\mathrm{{adj}}}}$ & " + " & ".join(pa(x) for x in g.p_adj_auroc)
      + r" \\")
    if i < len(CORE) - 1:
        A(r"\addlinespace[1.5pt]")
end_panel()

# ---- Panel C
begin_panel("lrrrr", 5, r"C\quad Linear probe and own-decomposition AUROC for extracellular space, common gene sets")
A(r"Representation & " + " & ".join(LABEL[c] for c in PROBE_CTX) + r" \\")
A(r" & " + " & ".join(f"$n$ = {f(ngenes[c][0], 0)} ({f(ngenes[c][1], 0)})" for c in PROBE_CTX) + r" \\")
A(r"\midrule")
REPS_SHORT = {"scGPT_L0": r"scGPT $L_0$", "scGPT_L3": r"scGPT $L_3$", "scGPT_L11": r"scGPT $L_{11}$",
              "scGPT_static_token": "scGPT input-token embedding", "ESM2_15B": "ESM2-15B"}
A(r"\multicolumn{5}{@{}l}{\itshape Probe, out-of-fold AUROC (95\% CI)}\\")
for rep, _ in REPS:
    A(rf"\quad {REPS_SHORT[rep]} & " + " & ".join(
        f"{f(r.probe_auroc)} ({f(r.probe_lo)}--{f(r.probe_hi)})"
        for r in [one(ext, context=c, representation=rep) for c in PROBE_CTX]) + r" \\")
A(r"\addlinespace[1.5pt]")
A(r"\multicolumn{5}{@{}l}{\itshape $\Delta$ probe AUROC from scGPT $L_3$ (95\% CI)}\\")
spx = sp[sp.term == "extracellular_space"]
sg = lambda x: (f"$+${f(x)}" if x >= 0 else f"$-${f(-x)}")
for rep, _ in REPS:
    if rep == "scGPT_L3":
        continue
    A(rf"\quad {REPS_SHORT[rep]} & " + " & ".join(
        f"{sg(r.delta)} ({sg(r.lo)}, {sg(r.hi)})"
        for r in [one(spx, context=c, representation=rep) for c in PROBE_CTX]) + r" \\")
A(r"\addlinespace[1.5pt]")
A(r"\multicolumn{5}{@{}l}{\itshape Own-SVD: $\SV{1}$ AUROC; best of $\SV{1}$--$\SV{12}$ (index)}\\")
for rep, _ in REPS:
    A(rf"\quad {REPS_SHORT[rep]} & " + " & ".join(
        f"{f(r.sv1_auroc)}; {f(r.best_sv_auroc)} ({int(r.best_sv)})"
        for r in [one(ext, context=c, representation=rep) for c in PROBE_CTX]) + r" \\")
end_panel()

# ---- Panel D
begin_panel("llrrrrlr", 8, r"D\quad Other compartments on the extracellular cell, and replication of the direction")
A(r" & & \multicolumn{3}{c}{Oriented AUROC ($p_{\mathrm{adj}}$)} & & \multicolumn{2}{c}{Lung-matched cell} \\")
A(r"\cmidrule(lr){3-5}\cmidrule(l){7-8}")
A(r"Extraction & Cell & Cytosol & Mitochondrion & ER lumen & $r$ with lung ($n$) & Cell ($|r|$) & "
  r"Extracellular ($p_{\mathrm{adj}}$) \\")
A(r"\midrule")
for c in CORE:
    a = one(ax, context=c)
    rr = "--" if c == "lung" else f"{f(a.r_with_lung_best_cell)} ({f(a.n_matched_with_lung, 0)})"
    lm = "--" if c == "lung" else f"{cell(a.lung_matched_layer, a.lung_matched_direction)} ({f(a.lung_matched_abs_r)})"
    le = "--" if c == "lung" else (f"{f(a.lung_matched_extracellular_auroc_oriented)} "
                                   f"({pa(a.lung_matched_extracellular_p_adj)})")
    A(f"{LABEL[c]} & {cell(a.layer, a.direction)} & {f(a.cytosol_auroc_oriented)} ({pa(a.cytosol_p_adj_auroc)}) & "
      f"{f(a.mitochondrion_auroc_oriented)} ({pa(a.mitochondrion_p_adj_auroc)}) & "
      f"{f(a.ER_lumen_auroc_oriented)} ({pa(a.ER_lumen_p_adj_auroc)}) & {rr} & {lm} & {le} \\\\")
end_panel(last=True)
A(r"\end{longtable}")
A(r"\endgroup")
OUT_TEX.write_text("\n".join(L) + "\n")
print(f"wrote {OUT_TEX} and {OUT_CSV} ({len(out)} rows, {out.shape[1]} columns)")
print(best[["context", "term_name", "layer", "direction", "auroc_signfree", "z_maxT", "strat_auroc_oriented",
            "OR_enriched_pole", "p_enriched_pole"]].round(4).to_string())
print("caption characters:", len(cap))
