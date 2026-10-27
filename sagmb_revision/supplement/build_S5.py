"""Build Supplementary Table S5 (H1 conditioned on compartment and on the detection threshold) and its CSV.

Sources (read only): revision3_experiments/F1_h1_compartment/
  classes.csv, strat_null.csv, residual_assort.csv, resid_controls.csv, resid_controls_summary.csv, within_between.csv,
  nobs_threshold.csv, nobs_threshold_cells.csv, intervals.csv, chain_diagnostics.csv
and revision3_experiments/E2_ppi/h1_summary.csv (full-V reference rows of panel E).
Writes sagmb_revision/supplement/S5.tex and sagmb_revision/supplementary_tables/Table_S5_compartment_conditioning.csv.
All intervals are the delete-d columns (*_deleted) of intervals.csv / nobs_threshold.csv (CONVENTION3_ERRATA item 3).
Derived quantities (own computation): counts of cells with max-T adjusted p < 0.05, medians, minima, maxima and
ranges over cells or seeds.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
F1 = ROOT / "revision3_experiments" / "F1_h1_compartment"
E2 = ROOT / "revision3_experiments" / "E2_ppi"
OUT_TEX = ROOT / "sagmb_revision" / "supplement" / "S5.tex"
OUT_CSV = ROOT / "sagmb_revision" / "supplementary_tables" / "Table_S5_compartment_conditioning.csv"

CTX = ["lung", "lung_external", "kidney_B"]
LAB = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B"}
TOP3 = {"lung": [(3, 1), (6, 2), (8, 3)], "lung_external": [(8, 2), (3, 1), (7, 2)], "kidney_B": [(9, 3), (8, 3), (7, 2)]}


def f(x, d=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    s = f"{float(x):.{d}f}"
    if s.startswith("-"):
        s = s[1:] if float(s) == 0 else "$-$" + s[1:]
    return s


def ci(lo, hi, d=3):
    return f"{f(lo, d)}--{f(hi, d)}"


def rng(lo, hi, d=2):
    """Range over contexts; a single value when both ends round to the same number."""
    a, b = f(lo, d), f(hi, d)
    return a if a == b else f"{a}--{b}"


def wci(est, lo, hi, d=3):
    return f"{f(est, d)} ({ci(lo, hi, d)})"


def n_(x):
    return f"{int(round(float(x))):,}"


def cell(L, k):
    return f"$L_{{{int(L)}}}\\,\\SV{{{int(k)}}}$"


def parse(c):
    a, b = c.split()
    return int(a[1:]), int(b[2:])


def panel_rows(letter, title, spec, head, body, pre=""):
    t = [f"\\textbf{{{letter}}}\\enspace {title}\\\\*",
         "{\\scriptsize" + pre + "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}" + spec + "@{}}"]
    t += head + ["\\midrule"] + body + ["\\end{tabular*}}\\\\"]
    return t


csv_parts = []


def add_csv(df, panel, source):
    d = df.copy()
    d.insert(0, "source_file", source)
    d.insert(0, "panel", panel)
    csv_parts.append(d)


iv = pd.read_csv(F1 / "intervals.csv")


def interval(c, analysis, statistic, L, k):
    t = iv[(iv.context == c) & (iv.analysis == analysis) & (iv.statistic == statistic) & (iv.layer == L) & (iv.sv == k)]
    assert len(t) == 1, (c, analysis, statistic, L, k)
    t = t.iloc[0]
    return t.est, t.ci_lo_deleted, t.ci_hi_deleted


add_csv(iv, "intervals_delete_d_and_prereg", "F1_h1_compartment/intervals.csv")

# ------------------------------------------------------------------ Panel A: classes and same-class edge share
cls = pd.read_csv(F1 / "classes.csv")
add_csv(cls, "A_classes", "F1_h1_compartment/classes.csv")
CLS = [("extracellular", "Extracellular space"), ("plasma_membrane", "Plasma membrane"), ("ER_Golgi", "ER or Golgi"),
       ("mitochondrion", "Mitochondrion"), ("nucleus", "Nucleus"), ("cytosol", "Cytosol"),
       ("other_annotated", "Other annotated"), ("unannotated", "Unannotated")]
A_body = []
for key, lab in CLS:
    row = [lab]
    for c in CTX:
        t = cls[(cls.context == c) & (cls.cls == key)].iloc[0]
        row += [n_(t.n_genes), f(t.frac_edges_within_this_class), f(t.config_expected_frac_edges_within_this_class)]
    A_body.append(" & ".join(row) + " \\\\")
A_body.append("\\midrule")
row_all, row_ms = ["All classes (same-class edges)"], ["MS null: mean (SD); $z$"]
for c in CTX:
    t = cls[(cls.context == c) & (cls.cls == "ALL_same_class")].iloc[0]
    row_all += [n_(t.n_genes), f(t.frac_edges_within_this_class), f(t.config_expected_frac_edges_within_this_class)]
    row_ms += ["", f"{f(t.ms_null_mean)} ({f(t.ms_null_sd, 4)})", f"$z = {t.ms_z:.1f}$"]
A_body += [" & ".join(row_all) + " \\\\", " & ".join(row_ms) + " \\\\"]
A_head = ["& \\multicolumn{3}{c}{lung} & \\multicolumn{3}{c}{lung (ext.)} & \\multicolumn{3}{c}{kidney~B} \\\\",
          "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(l){8-10}",
          "Class & Genes & Edge share & Expected & Genes & Edge share & Expected & Genes & Edge share & Expected \\\\"]
A = panel_rows("A", "Compartment classes and the share of STRING edges whose two genes fall in the same class",
               "lrrrrrrrrr", A_head, A_body)

# ------------------------------------------------------------------ Panel B: class-aware nulls
sn = pd.read_csv(F1 / "strat_null.csv")
add_csv(sn, "B_class_aware_nulls_cells", "F1_h1_compartment/strat_null.csv")
B_body = []
for c in CTX:
    g = sn[sn.context == c]
    cells = TOP3[c] + ([(7, 6)] if c == "lung" else [])
    for i, (L, k) in enumerate(cells):
        t = g[(g.layer == L) & (g.sv == k)].iloc[0]
        row = [LAB[c] if i == 0 else "", cell(L, k) + ("$^{*}$" if (c, L, k) == ("lung", 7, 6) else ""), f(t.r),
               f(t.f1ms_null_mean), f(t.strat_null_mean), f(t.strat_null_sd, 4), f(t.strat_z, 1),
               f(t.strat_excess_over_ms_excess), f(t.cpms_null_mean), f(t.cpms_null_sd, 4), f(t.cpms_z, 1),
               f(t.cpms_excess_over_ms_excess)]
        B_body.append(" & ".join(row) + " \\\\")
    row = ["", "120 cells", "--", "--",
           f"max {f(g.strat_null_mean.max())}", "", f"{int((g.strat_p_adj_maxT < 0.05).sum())}/120",
           f"{f(g.strat_excess_over_ms_excess.median())} ({f(g.strat_excess_over_ms_excess.min())})",
           f"max {f(g.cpms_null_mean.max())}", "", f"{int((g.cpms_p_adj_maxT < 0.05).sum())}/120",
           f"{f(g.cpms_excess_over_ms_excess.median())} ({f(g.cpms_excess_over_ms_excess.min())})"]
    assert (g.strat_p_adj_maxT.min() <= 1 / 2001 + 1e-12) and (g.cpms_p_adj_maxT.min() <= 1 / 2001 + 1e-12)
    assert ((g.strat_p_adj_maxT < 0.05) <= (g.r > 0)).all() and ((g.cpms_p_adj_maxT < 0.05) <= (g.r > 0)).all()
    B_body.append(" & ".join(row) + " \\\\")
    if c != CTX[-1]:
        B_body.append("\\addlinespace[2pt]")
B_head = ["& & & & \\multicolumn{4}{c}{Permutation within classes} & \\multicolumn{4}{c}{Class-pair-preserving MS} \\\\",
          "\\cmidrule(lr){5-8}\\cmidrule(l){9-12}",
          "Extraction & Cell & $r$ & MS null & Null mean & Null SD & $z$ & Retained & Null mean & Null SD & $z$ & Retained \\\\"]
B = panel_rows("B", "Compartment-aware nulls ($N = 2{,}000$) at the three strongest cells of each extraction",
               "llrrrrrrrrrr", B_head, B_body)

# ------------------------------------------------------------------ Panel C: residualisation
ra = pd.read_csv(F1 / "residual_assort.csv")
rc = pd.read_csv(F1 / "resid_controls.csv")
rcs = pd.read_csv(F1 / "resid_controls_summary.csv")
add_csv(ra, "C_residualised_cells", "F1_h1_compartment/residual_assort.csv")
add_csv(rc, "C_residualisation_controls_cells", "F1_h1_compartment/resid_controls.csv")
add_csv(rcs, "C_residualisation_controls_summary", "F1_h1_compartment/resid_controls_summary.csv")
DESIGNS = [("class8_dummies+lnobs", "residclass", "Class factor (8 levels) $+ \\log n_{\\mathrm{obs}}$"),
           ("cc30+lnobs", "resid30cc", "30 GO-CC terms $+ \\log n_{\\mathrm{obs}}$ (pre-specified)"),
           (None, None, "\\itshape Post hoc controls (no test)"),
           ("lnobs", None, "$\\log n_{\\mathrm{obs}}$ alone"),
           ("cc30", None, "30 GO-CC terms, without $\\log n_{\\mathrm{obs}}$"),
           ("compartment8_multilabel+lnobs", None, "8 class-defining terms, multi-label $+ \\log n_{\\mathrm{obs}}$"),
           ("bp30+lnobs", None, "30 GO-BP terms $+ \\log n_{\\mathrm{obs}}$"),
           ("mf30+lnobs", None, "30 GO-MF terms $+ \\log n_{\\mathrm{obs}}$"),
           ("random30+lnobs", None, "30 random gene sets $+ \\log n_{\\mathrm{obs}}$ (5 seeds)")]
C_body = []
for key, rakey, lab in DESIGNS:
    if key is None:
        C_body.append(f"\\multicolumn{{10}}{{@{{}}l}}{{{lab}}}\\\\")
        continue
    row = [lab]
    for c in CTX:
        if key == "random30+lnobs":
            s = rcs[(rcs.context == c) & rcs.design.str.startswith("random30+lnobs_seed")]
            assert len(s) == 5
            row += [ci(s.top1_ratio.min(), s.top1_ratio.max(), 2), ci(s.median_ratio.min(), s.median_ratio.max(), 2), "--"]
            continue
        s = rcs[(rcs.context == c) & (rcs.design == key)].iloc[0]
        if rakey is not None:
            g = ra[(ra.context == c) & (ra.design == rakey)]
            L, k = parse(s.top1_cell)
            est, lo, hi = interval(c, rakey, "ratio_resid_over_raw", L, k)
            assert abs(est - s.top1_ratio) < 1e-9 and abs(g.ratio_resid_over_raw.median() - s.median_ratio) < 1e-9
            assert g.p_adj_maxT.min() <= 1 / 2001 + 1e-12
            row += [wci(est, lo, hi), f(s.median_ratio), f"{int((g.p_adj_maxT < 0.05).sum())}"]
        else:
            row += [f(s.top1_ratio), f(s.median_ratio), "--"]
    C_body.append(" & ".join(row) + " \\\\")
C_head = ["& \\multicolumn{3}{c}{lung, $L_3\\,\\SV{1}$} & \\multicolumn{3}{c}{lung (ext.), $L_8\\,\\SV{2}$} & "
          "\\multicolumn{3}{c}{kidney~B, $L_9\\,\\SV{3}$} \\\\",
          "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(l){8-10}",
          "Design (regressors) & Strongest cell & Median & Sig. & Strongest cell & Median & Sig. & "
          "Strongest cell & Median & Sig. \\\\"]
C = panel_rows("C", "Residualisation: residual $r$ / raw $r$ after regressing each coordinate on gene annotations (MS null, $N = 2{,}000$)",
               ">{\\raggedright\\arraybackslash}p{3.5cm}crrcrrcrr", C_head, C_body, pre="\\setlength{\\tabcolsep}{2.5pt}")
share = {}
for key in ("cc30", "bp30+lnobs", "mf30+lnobs"):
    s = rcs[rcs.design == key]
    share[key] = (s.frac_edges_sharing_a_term.min(), s.frac_edges_sharing_a_term.max(),
                  s.frac_random_pairs_sharing_a_term.min(), s.frac_random_pairs_sharing_a_term.max())
r2 = rcs[rcs.design == "cc30+lnobs"].top1_R2

# ------------------------------------------------------------------ Panel D: within- vs between-class edges
wb = pd.read_csv(F1 / "within_between.csv")
add_csv(wb, "D_within_between_cells", "F1_h1_compartment/within_between.csv")
D_body = []
for c in CTX:
    L, k = TOP3[c][0]
    for i, (es, lab) in enumerate((("within", "Same class"), ("between", "Different classes"))):
        g = wb[(wb.context == c) & (wb.edges == es)]
        graw, gcc = g[g.coords == "raw"], g[g.coords == "class_centred"]
        e1 = interval(c, f"{es}_class_edges_raw", "r", L, k)
        e2 = interval(c, f"{es}_class_edges_class_centred", "r", L, k)
        assert abs(e1[0] - graw[(graw.layer == L) & (graw.sv == k)].r.iloc[0]) < 1e-9
        assert abs(e2[0] - gcc[(gcc.layer == L) & (gcc.sv == k)].r.iloc[0]) < 1e-9
        D_body.append(" & ".join([LAB[c] if i == 0 else "", cell(L, k) if i == 0 else "", lab, n_(graw.m.iloc[0]),
                                  wci(*e1), wci(*e2), f(graw.r.median()), f(gcc.r.median()),
                                  f"{int((graw.p_adj_maxT < 0.05).sum())}", f"{int((gcc.p_adj_maxT < 0.05).sum())}"])
                       + " \\\\")
D_head = ["& & & & \\multicolumn{2}{c}{Strongest cell, $r$ (95\\% CI)} & \\multicolumn{2}{c}{Median $r$ over 120 cells} & "
          "\\multicolumn{2}{c}{Sig.} \\\\",
          "\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(l){9-10}",
          "Extraction & Cell & Edges & $m$ & Raw & Class-centred & Raw & Class-centred & Raw & Class-centred \\\\"]
D = panel_rows("D", "Edges within one class and between classes, each against its own MS null ($N = 1{,}000$)",
               "lllrccrrrr", D_head, D_body)

# ------------------------------------------------------------------ Panel E: detection threshold
nt = pd.read_csv(F1 / "nobs_threshold.csv")
ntc = pd.read_csv(F1 / "nobs_threshold_cells.csv")
add_csv(nt, "E_nobs_threshold_summary", "F1_h1_compartment/nobs_threshold.csv")
add_csv(ntc, "E_nobs_threshold_cells", "F1_h1_compartment/nobs_threshold_cells.csv")
hs = pd.read_csv(E2 / "h1_summary.csv").set_index("context")
E_body = []
for c in CTX:
    L0, k0 = TOP3[c][0]
    s = hs.loc[c]
    est, lo, hi = interval(c, "raw_all_edges", "r", L0, k0)
    E_body.append(" & ".join([LAB[c], "all ($V$)", n_(s.n_V), n_(s.m), f"{int(s.n_cells_p_adj_lt_05)}", cell(L0, k0),
                              wci(est, lo, hi), f(s.median_r), cell(L0, k0), "1", wci(est, lo, hi), f(s.top_z, 1)]) + " \\\\")
    for thr in (3, 5, 10):
        t = nt[(nt.context == c) & (nt.nobs_min == thr)].iloc[0]
        assert abs(t.family_min_p_adj - 1 / (t.N + 1)) < 1e-12
        Lt, kt = parse(t.top_cell)
        Lm, km = parse(t.e2top1_matched)
        em, lom, him = interval(c, f"nobs_ge_{thr}", "r", Lm, km)
        assert abs(em - t.e2top1_r) < 1e-9
        g = ntc[(ntc.context == c) & (ntc.nobs_min == thr)]
        assert int((g.p_adj_maxT < 0.05).sum()) == int(t.n_cells_p_adj_lt_05)
        E_body.append(" & ".join(["", f"$\\geq {thr}$", n_(t.n_V), n_(t.m), f"{int(t.n_cells_p_adj_lt_05)}", cell(Lt, kt),
                                  wci(t.top_r, t.top_r_ci_lo_deleted, t.top_r_ci_hi_deleted), f(t.median_r),
                                  cell(Lm, km), f(t.e2top1_abs_corr), wci(em, lom, him), f(t.e2top1_z, 1)]) + " \\\\")
    if c != CTX[-1]:
        E_body.append("\\addlinespace[2pt]")
E_head = ["& & & & & \\multicolumn{3}{c}{Strongest cell} & \\multicolumn{4}{c}{Strongest direction of the full-$V$ analysis} \\\\",
          "\\cmidrule(lr){6-8}\\cmidrule(l){9-12}",
          "Extraction & $n_{\\mathrm{obs}}$ & $n_V$ & $m$ & Sig. & Cell & $r$ (95\\% CI) & Median $r$ & Cell & $|$corr$|$ & "
          "$r$ (95\\% CI) & $z$ \\\\"]
E = panel_rows("E", "Gene set restricted to genes observed in at least 3, 5 or 10 sampled cells (SVD recomputed; MS null, $N = 2{,}000$)",
               "llrrrlcrlrcr", E_head, E_body)

# ------------------------------------------------------------------ chain diagnostics
chd = pd.read_csv(F1 / "chain_diagnostics.csv")
add_csv(chd, "chain_diagnostics", "F1_h1_compartment/chain_diagnostics.csv")
cp_acc = chd[chd.null == "cpms:null"].acc_mean
cp_lag = chd[chd.null == "cpms:null"].lag1_top_cell

cap = (
    "\\textbf{H1 conditioned on compartment co-membership and on the detection threshold.} "
    "Lung, lung (ext.) and kidney~B; $r$, scalar assortativity of a spectral coordinate over the STRING physical edges "
    "(combined score $\\geq 700$) among observed genes, as in Table~S3; MS, Maslov--Sneppen degree-preserving rewiring; Sig., number "
    "of the 120 cells ($\\SV{1}$--$\\SV{10}$ $\\times$ $L_0$--$L_{11}$) with two-sided max-$T$ adjusted $p < 0.05$ (the family minimum "
    "is at the permutation floor $1/(N+1)$ in every family; every significant cell has $r > 0$). "
    "95\\% CI: vertex subsampling ($n_V/2$ genes, $B = 500$) with delete-$d$ scaling, coordinates and residuals held fixed. "
    "\\textbf{(A)} Each gene is assigned one class from its direct GO cellular-component annotations, in the priority order "
    "listed (extracellular space GO:0005615, plasma membrane GO:0005886, endoplasmic reticulum, ER lumen or Golgi, mitochondrion, "
    "nucleus, cytosol, any other direct GO annotation, none). Edge share, fraction of all STRING edges with both genes in that class; "
    "expected, its stub-matching (configuration-model) expectation; MS null, the same-class share in $N = 2{,}000$ rewired graphs. "
    "\\textbf{(B)} Cells, the three cells with the largest $|z|$ under MS in Table~S3 ($^{*}$lung $L_7\\,\\SV{6}$, a cell with little "
    "compartment content, chosen post hoc). MS null, null mean of $r$ in an MS replica with the same $N$; permutation within "
    "classes, one permutation of coordinates within each class applied to all 120 cells (preserves class, not degree); "
    "class-pair-preserving MS, double-edge swaps accepted only if every edge keeps its pair of classes (preserves every degree "
    "and the multiset of edge class pairs; acceptance "
    f"{f(cp_acc.min())}--{f(cp_acc.max())}, lag-1 autocorrelation of the strongest cell's null statistic "
    f"{f(cp_lag.min(), 2)}--{f(cp_lag.max(), 2)}); "
    "retained, $(r - \\text{null mean})/(r - \\text{MS null mean})$, the share of the "
    "excess over MS that remains (no interval). Rows `120 cells': largest null mean over cells; Sig.; median (minimum) retained share. "
    "\\textbf{(C)} Ordinary least squares of each coordinate on the regressors listed; entries are residual $r$ / raw $r$ at the "
    "strongest cell (with 95\\% CI for the two tested designs) and the median over the 120 cells; Sig. for the residual coordinates. "
    "Class factor, seven class indicators and an intercept; 30 GO-CC/BP/MF terms, indicators of the 30 most frequent direct "
    "cellular-component, biological-process or molecular-function annotations among the genes of $V$; 8 class-defining terms, "
    "the eight class definitions of panel A as overlapping indicators; random gene sets, 30 sets of the same sizes as the CC terms; "
    "for these, ranges over the five seeds (two decimals). The 30-CC design explains "
    f"{100 * r2.min():.0f}--{100 * r2.max():.0f}\\% of the strongest coordinate's variance. Share of STRING edges whose two genes "
    f"share at least one term: CC {rng(*share['cc30'][:2])} (random pairs {rng(*share['cc30'][2:])}), "
    f"BP {rng(*share['bp30+lnobs'][:2])} ({rng(*share['bp30+lnobs'][2:])}), "
    f"MF {rng(*share['mf30+lnobs'][:2])} ({rng(*share['mf30+lnobs'][2:])}; its most frequent "
    "term, `protein binding', is largely interaction-derived). "
    "\\textbf{(D)} STRING edges split into those whose two genes share a class and the rest, each evaluated with its own MS chain "
    "(which does not keep edges within or between classes); class-centred, coordinates centred within each class before computing $r$. "
    "\\textbf{(E)} $V$ restricted to genes with $n_{\\mathrm{obs}}$ (number of sampled cells contributing to the gene) at or above "
    "the threshold, SVD recomputed on the restricted genes and edges among them; rows `all ($V$)', the unrestricted analysis "
    "(Table~S3; MS $N = 10{,}000$ for lung and lung (ext.), interval from the draws of this analysis); strongest direction of the "
    "full-$V$ analysis, the restricted direction at the same layer with the largest $|$correlation$|$ of coordinates on the kept genes. "
    "Full values: \\texttt{Table\\_S5\\_compartment\\_conditioning.csv}."
)

lines = ["\\setcounter{table}{4}",
         "\\begingroup\\footnotesize",
         "\\setlength{\\tabcolsep}{3pt}",
         "\\begin{longtable}{@{}p{\\textwidth}@{}}",
         "\\caption{" + cap + "}\\label{tab:S5}\\\\",
         "\\toprule",
         "\\endfirsthead",
         "\\multicolumn{1}{@{}l}{\\footnotesize\\itshape Table~\\thetable\\ (continued)}\\\\",
         "\\toprule",
         "\\endhead",
         "\\midrule",
         "\\multicolumn{1}{@{}r}{\\footnotesize\\itshape continued on next page}\\\\",
         "\\endfoot",
         "\\bottomrule",
         "\\endlastfoot"]
lines += A + ["\\midrule"] + B + ["\\midrule"] + C + ["\\midrule"] + D + ["\\midrule"] + E
lines += ["\\end{longtable}", "\\endgroup", ""]
OUT_TEX.write_text("\n".join(lines))

out = pd.concat(csv_parts, ignore_index=True, sort=False)
front = ["panel", "source_file", "context"]
out = out[front + [c for c in out.columns if c not in front]]
out.to_csv(OUT_CSV, index=False)
print("wrote", OUT_TEX, OUT_CSV, out.shape)
