"""Build Supplementary Table S3 (H1 interaction sensitivity analyses) and its CSV.

Sources (read only):
  revision3_experiments/E2_ppi/            h1_summary.csv, h1_primary_<ctx>.csv, verify/verify_c.json (delete-d CIs),
                                           sens_nulls_summary.csv, sens_{label,degree_withinV,degree_genomewide,detection}.csv,
                                           sens_robustness_summary.csv, sens_robustness.csv, sens_hub_exclusion.csv,
                                           sens_jackknife_summary.csv, sens_jackknife_top10.csv, sens_coexpr_matched.csv,
                                           sens_coexpr_fine.csv, sens_complexes_summary.csv, sens_complexes.csv,
                                           sens_copole_K.csv, chain_diagnostics.csv
  revision3_experiments/E12_string_gradient/gradient.csv
Writes sagmb_revision/supplement/S3.tex and sagmb_revision/supplementary_tables/Table_S3_interaction.csv.

Delete-d intervals (CONVENTION3_ERRATA item 3): top-cell r intervals are the verifier's delete-d intervals
(verify_c.json); vertex-subsampling intervals stored with the pre-registered sqrt(m/n) scaling (co-expression AUROC,
co-pole ratio) are rescaled to delete-d by multiplying both half-widths by sqrt(2); the E12 gradient intervals are
already delete-d.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
E2 = ROOT / "revision3_experiments" / "E2_ppi"
E12 = ROOT / "revision3_experiments" / "E12_string_gradient"
OUT_TEX = ROOT / "sagmb_revision" / "supplement" / "S3.tex"
OUT_CSV = ROOT / "sagmb_revision" / "supplementary_tables" / "Table_S3_interaction.csv"
SQ2 = np.sqrt(2.0)

CTX = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
LAB = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune", "kidney_A": "kidney~A"}
ROB = [("immune_seed43", "immune, draw (seed 43)"), ("immune_seed44", "immune, draw (seed 44)"),
       ("immune_TSP14", "immune, donor TSP14"), ("immune_TSP21", "immune, donor TSP21"),
       ("immune_TSP25", "immune, donor TSP25"), ("immune_Bcell", "immune, B cells"),
       ("immune_CD4T", "immune, CD4 T cells"), ("immune_CD8T", "immune, CD8 T cells"),
       ("kidney_ckpt_whole_human", "kidney cells, whole-human ckpt."), ("kidney_ckpt_kidney", "kidney cells, kidney ckpt."),
       ("kidney_ckpt_brain", "kidney cells, brain ckpt.")]


# ------------------------------------------------------------------ formatting
def f(x, d=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    s = f"{float(x):.{d}f}"
    if s.startswith("-"):
        s = s[1:] if float(s) == 0 else "$-$" + s[1:]
    return s


def ci(lo, hi, d=3):
    return f"{f(lo, d)}--{f(hi, d)}"


def wci(est, lo, hi, d=3):
    return f"{f(est, d)} ({ci(lo, hi, d)})"


def n_(x):
    return f"{int(round(float(x))):,}"


def cell(L, k):
    return f"$L_{{{int(L)}}}\\,\\SV{{{int(k)}}}$"


def deld(est, lo, hi):
    """Pre-registered sqrt(m/n) interval -> delete-d interval (half-widths x sqrt(2))."""
    return est - SQ2 * (est - lo), est + SQ2 * (hi - est)


def panel_rows(letter, title, spec, head, body):
    """One outer longtable row with the panel title, one with the nested full-width tabular."""
    t = [f"\\textbf{{{letter}}}\\enspace {title}\\\\*",
         "{\\scriptsize\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}" + spec + "@{}}"]
    t += head + ["\\midrule"] + body + ["\\end{tabular*}}\\\\"]
    return t


csv_parts = []


def add_csv(df, panel, source):
    d = df.copy()
    d.insert(0, "source_file", source)
    d.insert(0, "panel", panel)
    csv_parts.append(d)


# ------------------------------------------------------------------ Panel A: H1 family per extraction
hs = pd.read_csv(E2 / "h1_summary.csv").set_index("context")
vc = json.loads((E2 / "verify" / "verify_c.json").read_text())
rs = pd.read_csv(E2 / "sens_robustness_summary.csv")
rs_ms = rs[rs.null == "ms"].set_index("context")

A_body = []
for c in CTX:
    s = hs.loc[c]
    lo, hi = vc[c]["ci_delete_d"]
    assert vc[c]["cell"] == f"L{int(s.top_layer)} SV{int(s.top_sv)}"
    assert abs(vc[c]["r"] - s.top_r) < 1e-9
    A_body.append(" & ".join([LAB[c], n_(s.n_V), n_(s.m), n_(s.N), f"{int(s.n_cells_p_adj_lt_05)}",
                              cell(s.top_layer, s.top_sv), f(s.top_r), ci(lo, hi), f(s.top_null_mean),
                              f(s.top_z, 1), f(s.median_r)]) + " \\\\")
A_body.append("\\addlinespace[2pt]\\multicolumn{11}{@{}l}{\\itshape Robustness extractions (MS, $N = 1{,}000$)}\\\\")
for c, lab in ROB:
    s = rs_ms.loc[c]
    A_body.append(" & ".join([lab, n_(s.n_V), n_(s.m), n_(s.N), f"{int(s.n_cells_p_adj_lt_05)}",
                              cell(s.top_layer, s.top_sv), f(s.top_r), "--", f(s.top_null_mean),
                              f(s.top_z, 1), f(s.median_r)]) + " \\\\")
A_head = ["Extraction & $n_V$ & $m$ & $N$ & Sig. & Strongest cell & $r$ & 95\\% CI & Null mean & $z$ & Median $r$ \\\\"]
A = panel_rows("A", "H1 family against Maslov--Sneppen rewiring (MS), per extraction", "lrrrrlrcrrr", A_head, A_body)

summ = hs.reset_index().copy()
summ["top_r_ci_lo_deleted"] = [vc[c]["ci_delete_d"][0] for c in summ.context]
summ["top_r_ci_hi_deleted"] = [vc[c]["ci_delete_d"][1] for c in summ.context]
summ = summ.rename(columns={"top_r_ci_lo": "top_r_ci_lo_prereg", "top_r_ci_hi": "top_r_ci_hi_prereg",
                            "top_excess_r_lo": "top_excess_r_lo_prereg", "top_excess_r_hi": "top_excess_r_hi_prereg",
                            "second_r_ci": "second_r_ci_prereg", "third_r_ci": "third_r_ci_prereg"})
add_csv(summ, "A_summary_main", "E2_ppi/h1_summary.csv + verify/verify_c.json")
add_csv(rs_ms.reset_index(), "A_summary_robustness", "E2_ppi/sens_robustness_summary.csv")
for c in CTX:
    add_csv(pd.read_csv(E2 / f"h1_primary_{c}.csv").assign(null="maslov_sneppen"), "A_cells_main",
            f"E2_ppi/h1_primary_{c}.csv")
rob_cells = pd.read_csv(E2 / "sens_robustness.csv")
add_csv(rob_cells[rob_cells.null == "ms"], "A_cells_robustness", "E2_ppi/sens_robustness.csv")

# ------------------------------------------------------------------ Panel B: alternative nulls at the strongest cell
ns = pd.read_csv(E2 / "sens_nulls_summary.csv")
NULLS = [("label", "sens_label.csv", "Label shuffle over $V$"),
         ("deg_v", "sens_degree_withinV.csv", "Degree strata ($V$ network)"),
         ("deg_g", "sens_degree_genomewide.csv", "Degree strata (genome-wide STRING)"),
         ("det", "sens_detection.csv", "Detection strata ($\\log n_{\\mathrm{obs}}$)"),
         ("maslov_sneppen", None, "MS rewiring (primary)")]
cells_by_null = {}
for key, fn, _ in NULLS:
    if fn is not None:
        d = pd.read_csv(E2 / fn)
        cells_by_null[key] = d
        add_csv(d, "B_cells_altnull", f"E2_ppi/{fn}")
add_csv(ns, "B_summary_altnull", "E2_ppi/sens_nulls_summary.csv")
B_body = []
for key, fn, lab in NULLS:
    row = [lab]
    for c in ("lung", "lung_external"):
        L, k = int(hs.loc[c, "top_layer"]), int(hs.loc[c, "top_sv"])
        s = ns[(ns.context == c) & (ns.null == key)].iloc[0]
        if key == "maslov_sneppen":
            d = pd.read_csv(E2 / f"h1_primary_{c}.csv")
        else:
            d = cells_by_null[key][cells_by_null[key].context == c]
        t = d[(d.layer == L) & (d.sv == k)].iloc[0]
        if key != "maslov_sneppen":
            assert abs(t.null_mean - s.ms_top_null_mean_this) < 1e-12 and abs(t.z - s.ms_top_z_under_this_null) < 1e-9
        else:
            assert abs(t.null_mean - s.top_null_mean) < 1e-12
        row += [f(t.null_mean), f(t.null_sd), f(t.z, 1), f"{int(s.n_cells_p_adj_lt_05)}"]
    B_body.append(" & ".join(row) + " \\\\")
B_head = ["& \\multicolumn{4}{c}{lung, $L_3\\,\\SV{1}$ ($r = " + f(hs.loc['lung', 'top_r']) + "$)} & "
          "\\multicolumn{4}{c}{lung (ext.), $L_8\\,\\SV{2}$ ($r = " + f(hs.loc['lung_external', 'top_r']) + "$)} \\\\",
          "\\cmidrule(lr){2-5}\\cmidrule(l){6-9}",
          "Null & Null mean & Null SD & $z$ & Sig. & Null mean & Null SD & $z$ & Sig. \\\\"]
B = panel_rows("B", "Null mean of $r$ at the strongest cell under each null ($N = 10{,}000$)", "lrrrrrrrr", B_head, B_body)

# ------------------------------------------------------------------ Panel C: hub removal and jackknife
hub = pd.read_csv(E2 / "sens_hub_exclusion.csv")
jks = pd.read_csv(E2 / "sens_jackknife_summary.csv").set_index("context")
add_csv(hub, "C_hub_exclusion", "E2_ppi/sens_hub_exclusion.csv")
add_csv(jks.reset_index(), "C_jackknife_summary", "E2_ppi/sens_jackknife_summary.csv")
add_csv(pd.read_csv(E2 / "sens_jackknife_top10.csv"), "C_jackknife_top10_genes", "E2_ppi/sens_jackknife_top10.csv")
C_body = []
for c in CTX:
    L, k = int(hs.loc[c, "top_layer"]), int(hs.loc[c, "top_sv"])
    row = [LAB[c], cell(L, k), f(hs.loc[c, "top_r"])]
    for fr in (0.05, 0.10):
        h = hub[(hub.context == c) & (hub.layer == L) & (hub.sv == k) & np.isclose(hub.hub_frac, fr)].iloc[0]
        row += [n_(h.n_hubs_removed), f(h.edge_frac_remaining, 2), f(h.r_nohub), f(h.z, 1)]
    j = jks.loc[c]
    assert int(j.layer) == L and int(j.sv) == k
    row += [f(j.max_influence, 4), f(j.jackknife_se)]
    C_body.append(" & ".join(row) + " \\\\")
C_head = ["& & & \\multicolumn{4}{c}{Top 5\\% by degree removed} & \\multicolumn{4}{c}{Top 10\\% by degree removed} & "
          "\\multicolumn{2}{c}{Leave one gene out} \\\\",
          "\\cmidrule(lr){4-7}\\cmidrule(lr){8-11}\\cmidrule(l){12-13}",
          "Extraction & Cell & $r$ & Genes & Edges kept & $r$ & $z$ & Genes & Edges kept & $r$ & $z$ & Max.\\ change & Jackknife SE \\\\"]
C = panel_rows("C", "Hub removal (MS null on the remaining graph, $N = 2{,}000$) and leave-one-gene-out influence, strongest cell",
               "llrrrrrrrrrrr", C_head, C_body)

# ------------------------------------------------------------------ Panel D: co-expression-matched non-edges, co-complex pairs
ce = pd.read_csv(E2 / "sens_coexpr_matched.csv")
cf = pd.read_csv(E2 / "sens_coexpr_fine.csv")
ce["strat_lo_deleted"] = np.where(ce.strat_lo.notna(), ce.auroc_matched_stratified - SQ2 * (ce.auroc_matched_stratified - ce.strat_lo), np.nan)
ce["strat_hi_deleted"] = np.where(ce.strat_hi.notna(), ce.auroc_matched_stratified + SQ2 * (ce.strat_hi - ce.auroc_matched_stratified), np.nan)
ce = ce.rename(columns={"strat_lo": "strat_lo_prereg", "strat_hi": "strat_hi_prereg"})
add_csv(ce, "D_coexpr_matched_deciles", "E2_ppi/sens_coexpr_matched.csv")
add_csv(cf, "D_coexpr_matched_50bins", "E2_ppi/sens_coexpr_fine.csv")
imm50 = cf[(cf.context == "immune") & (cf.coexpr_source == "extraction")].iloc[0]
D1_body = []
for c in CTX:
    e = ce[(ce.context == c) & (ce.coexpr_source == "extraction")]
    L, k = int(e.top_layer.iloc[0]), int(e.top_sv.iloc[0])
    assert (L, k) == (int(hs.loc[c, "top_layer"]), int(hs.loc[c, "top_sv"]))
    row = [LAB[c], cell(L, k), f"{f(e.mean_abs_coexpr_pos.iloc[0])}/{f(e.mean_abs_coexpr_matched.iloc[0])}"]
    for sc in (f"cosine_L{L}", "top_cell_neg_absdiff", "abs_coexpr"):
        t = e[e.scorer == sc].iloc[0]
        row.append(wci(t.auroc_matched_stratified, t.strat_lo_deleted, t.strat_hi_deleted))
    fb = cf[(cf.context == c) & (cf.coexpr_source == "extraction")]
    t3 = ce[(ce.context == c) & (ce.coexpr_source == "tier3000")]
    row += [f(fb[fb.scorer == f"cosine_L{L}"].auroc_matched_stratified.iloc[0]),
            f(fb[fb.scorer == "abs_coexpr"].auroc_matched_stratified.iloc[0]),
            f(t3[t3.scorer == f"cosine_L{L}"].auroc_matched_stratified.iloc[0]),
            f(t3[t3.scorer == "abs_coexpr"].auroc_matched_stratified.iloc[0])]
    D1_body.append(" & ".join(row) + " \\\\")
D1_head = ["& & & \\multicolumn{3}{c}{Deciles of $|$co-expression$|$, extraction cells} & \\multicolumn{2}{c}{50 bins} & "
           "\\multicolumn{2}{c}{3,000 cells} \\\\",
           "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}\\cmidrule(l){9-10}",
           "Extraction & Cell & Mean $|$coexpr$|$ & Cosine & $-|\\Delta x|$ & $|$coexpr$|$ & Cosine & $|$coexpr$|$ & Cosine & $|$coexpr$|$ \\\\"]
D1 = panel_rows("D", "STRING edges against co-expression-matched non-edges (stratified AUROC), and co-complex pairs as edge sets",
                "llcccccccc", D1_head, D1_body)

cx = pd.read_csv(E2 / "sens_complexes_summary.csv")
add_csv(cx, "D_complexes_summary", "E2_ppi/sens_complexes_summary.csv")
add_csv(pd.read_csv(E2 / "sens_complexes.csv"), "D_complexes_cells", "E2_ppi/sens_complexes.csv")
SETS = [("ComplexPortal", "ComplexPortal"), ("ComplexPortal minus STRING700", "ComplexPortal, not STRING $\\geq 700$"),
        ("hu.MAP 2.0", "hu.MAP 2.0"), ("hu.MAP 2.0 minus STRING700", "hu.MAP 2.0, not STRING $\\geq 700$")]
D2_body = []
for key, lab in SETS:
    row = [lab]
    for c in ("lung", "lung_external"):
        s = cx[(cx.context == c) & (cx.edge_set == key)].iloc[0]
        assert s.string_top_cell == f"L{int(hs.loc[c, 'top_layer'])} SV{int(hs.loc[c, 'top_sv'])}"
        row += [n_(s.m), f(s.r_at_string_top), f(s.z_at_string_top, 1), f"{int(s.n_cells_p_adj_lt_05)}"]
    D2_body.append(" & ".join(row) + " \\\\")
D2_head = ["& \\multicolumn{4}{c}{lung, $L_3\\,\\SV{1}$} & \\multicolumn{4}{c}{lung (ext.), $L_8\\,\\SV{2}$} \\\\",
           "\\cmidrule(lr){2-5}\\cmidrule(l){6-9}",
           "Co-complex edge set (MS, $N = 2{,}000$) & Pairs & $r$ & $z$ & Sig. & Pairs & $r$ & $z$ & Sig. \\\\"]
D2 = ["{\\scriptsize\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrrrr@{}}"] + D2_head + ["\\midrule"] + \
     D2_body + ["\\end{tabular*}}\\\\"]

# ------------------------------------------------------------------ Panel E: co-pole rate vs pole size; label vs MS null
cp = pd.read_csv(E2 / "sens_copole_K.csv")
cp["ratio_ci_lo_deleted"] = cp.ratio_analytic_est - SQ2 * (cp.ratio_analytic_est - cp.ratio_ci_lo)
cp["ratio_ci_hi_deleted"] = cp.ratio_analytic_est + SQ2 * (cp.ratio_ci_hi - cp.ratio_analytic_est)
cp = cp.rename(columns={"ratio_ci_lo": "ratio_ci_lo_prereg", "ratio_ci_hi": "ratio_ci_hi_prereg"})
add_csv(cp, "E_copole_K", "E2_ppi/sens_copole_K.csv")
FR = [0.10, 0.20, 0.27, 0.33, 0.50]
E1_body = []
relgap = []
for c in ("lung", "lung_external", "kidney_B"):
    L, k = int(hs.loc[c, "top_layer"]), int(hs.loc[c, "top_sv"])
    row = [LAB[c], cell(L, k)]
    for fr in FR:
        t = cp[(cp.context == c) & (cp.layer == L) & (cp.sv == k) & (cp.null == "ms") & np.isclose(cp.pole_frac, fr)].iloc[0]
        row.append(wci(t.ratio_obs_null, t.ratio_ci_lo_deleted, t.ratio_ci_hi_deleted, 2))
        relgap.append(abs(t.analytic_null_mean / t.null_mean - 1))
    E1_body.append(" & ".join(row) + " \\\\")
E1_head = ["& & \\multicolumn{5}{c}{Observed / MS-null co-pole rate (95\\% CI) at pole size $K/n$} \\\\",
           "\\cmidrule(l){3-7}",
           "Extraction & Cell & 0.10 & 0.20 & 0.27 & 0.33 & 0.50 \\\\"]
E1 = panel_rows("E", "Co-pole rate as a function of pole size, and label shuffle against rewiring at $K/n = 0.27$",
                "llccccc", E1_head, E1_body)
E2_cells = [("lung", 3, 1), ("lung_external", 8, 2), ("kidney_B", 9, 3), ("kidney_B", 11, 1), ("immune", 4, 6),
            ("kidney_A", 5, 2), ("kidney_A", 2, 2), ("kidney_A", 4, 2), ("kidney_A", 0, 1), ("kidney_A", 11, 1)]
E2_body = []
for c, L, k in E2_cells:
    sub = cp[(cp.context == c) & (cp.layer == L) & (cp.sv == k) & np.isclose(cp.pole_frac, 0.27)]
    tl = sub[sub.null == "label"].iloc[0]; tm = sub[sub.null == "ms"].iloc[0]
    assert abs(tl.copole_obs - tm.copole_obs) < 1e-12
    E2_body.append(" & ".join([LAB[c], cell(L, k), f(tm.copole_obs),
                               f(tl.null_mean), f(tl.null_sd, 4), f(tl.ratio_obs_null, 2), f(tl.z, 1),
                               f(tm.null_mean), f(tm.null_sd, 4), f(tm.ratio_obs_null, 2), f(tm.z, 1)]) + " \\\\")
E2_head = ["& & & \\multicolumn{4}{c}{Label shuffle} & \\multicolumn{4}{c}{MS rewiring} \\\\",
           "\\cmidrule(lr){4-7}\\cmidrule(l){8-11}",
           "Extraction & Cell & Observed & Null mean & Null SD & Ratio & $z$ & Null mean & Null SD & Ratio & $z$ \\\\"]
E2t = ["{\\scriptsize\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llrrrrrrrrr@{}}"] + E2_head + ["\\midrule"] + \
      E2_body + ["\\end{tabular*}}\\\\"]

# ------------------------------------------------------------------ Panel F: STRING confidence gradient
gr = pd.read_csv(E12 / "gradient.csv")
add_csv(gr.rename(columns={"ctx": "context"}), "F_string_gradient", "E12_string_gradient/gradient.csv")
BANDS = ["150", "400", "700", "900"]
F_body = []
mono_counts = {}
for c in CTX:
    g = gr[gr.ctx == c].sort_values("layer")
    mono = int(((g.auroc_150 < g.auroc_400) & (g.auroc_400 < g.auroc_700) & (g.auroc_700 < g.auroc_900)).sum())
    mono_counts[c] = mono
    for i, L in enumerate((0, 11)):
        t = g[g.layer == L].iloc[0]
        F_body.append(" & ".join([LAB[c] if i == 0 else "", n_(t.n_links) if i == 0 else "", f"$L_{{{L}}}$",
                                  wci(t.rho, t.rho_lo, t.rho_hi), wci(t.rho_adj, t.rho_adj_lo, t.rho_adj_hi)]
                                 + [f(t[f"auroc_{b}"]) for b in BANDS] + [""]) + " \\\\")
    F_body.append(" & ".join(["", "", "all 12", ci(g.rho.min(), g.rho.max()),
                              f"{ci(g.rho_adj.min(), g.rho_adj.max())} [{f(g.rho_adj_lo.min())}]"]
                             + [ci(g[f'auroc_{b}'].min(), g[f'auroc_{b}'].max()) for b in BANDS]
                             + [f"{mono}/12"]) + " \\\\")
    if c != CTX[-1]:
        F_body.append("\\addlinespace[2pt]")
F_head = ["& & & & & \\multicolumn{4}{c}{AUROC by STRING combined-score band} & \\\\",
          "\\cmidrule(lr){6-9}",
          "Extraction & Links & Layer & $\\rho$ (95\\% CI) & $\\rho_{\\mathrm{adj}}$ (95\\% CI) & 150--399 & 400--699 & 700--899 & 900--999 & Monotone \\\\"]
F = panel_rows("F", "STRING confidence gradient (physical links, combined score $\\geq 150$, within $V$)",
               "lrlccccccr", F_head, F_body)

# ------------------------------------------------------------------ chain diagnostics (caption numbers, CSV)
ch = pd.read_csv(E2 / "chain_diagnostics.csv")
add_csv(ch, "chain_diagnostics_MS", "E2_ppi/chain_diagnostics.csv")
acc_lo, acc_hi = ch.acc_mean.min(), ch.acc_mean.max()
lag_lo, lag_hi = ch.lag1_top_cell.min(), ch.lag1_top_cell.max()
gap_lo, gap_hi = min(relgap), max(relgap)

# ------------------------------------------------------------------ caption and assembly
cap = (
    "\\textbf{H1: physical-interaction partners share spectral coordinates (significance, nulls, hubs, co-expression, "
    "co-complex pairs, pole size and STRING confidence).} "
    "$r$, scalar assortativity of a spectral coordinate over the STRING v12 physical edges (combined score $\\geq 700$) whose two "
    "ends are observed genes ($V$; $n_V$ genes, $m$ edges), i.e.\\ the Pearson correlation of the coordinate across the two ends "
    "of every edge. Cell, one direction at one layer; the family is $\\SV{1}$--$\\SV{10}$ at $L_0$--$L_{11}$ (120 cells). "
    "MS, Maslov--Sneppen degree-preserving rewiring of the same edges (double-edge swaps; burn-in $10m$ accepted swaps, then $m$ "
    "between draws); $N$, joint null draws applied to all 120 cells; Sig., number of cells (of 120) with two-sided Westfall--Young "
    "max-$T$ adjusted $p < 0.05$ (the family minimum is at the permutation floor $1/(N+1)$ in every family shown, and every "
    "significant cell has $r > 0$); "
    "$z = (r - \\text{null mean})/\\text{null SD}$; strongest cell, the cell with the largest $|z|$ under MS. "
    "95\\% CI: vertex subsampling ($n_V/2$ genes without replacement, $B = 500$) with delete-$d$ scaling. "
    "\\textbf{(A)} Robustness extractions: two further 256-cell immune draws, one draw per donor for three donors, three immune "
    "cell types, and one further 192-cell draw of kidney~B's cells embedded with each of three scGPT checkpoints (whole-human, "
    "kidney, brain); no interval was computed for them. "
    "\\textbf{(B)} Label shuffle, permutation of coordinates across $V$; degree strata, permutation within deciles of degree in the "
    "$V$ network or in the genome-wide STRING network; detection strata, permutation within deciles of $\\log n_{\\mathrm{obs}}$, the "
    "number of sampled cells that contributed to a gene. "
    "\\textbf{(C)} The 5\\% or 10\\% most connected genes of $V$ are removed with their edges; genes, number removed; edges kept, "
    "fraction of the $m$ edges remaining; max.\\ change, largest decrease of $r$ when one gene and its edges are removed; jackknife SE, "
    "delete-1 jackknife standard error of $r$. "
    "\\textbf{(D)} Upper part: every STRING edge is compared with five gene pairs of $V$ that have no physical STRING link at any score, "
    "drawn from the same decile of $|$co-expression$|$ (Pearson correlation of log-normalised expression over the extraction's cells); "
    "AUROC stratified by decile, for three scores: cosine of the two genes' centred 512-dimensional vectors at the strongest cell's "
    "layer, $-|\\Delta x|$ (negative absolute difference of the two genes' coordinates on the strongest cell), and $|$co-expression$|$ "
    "itself (the co-expression left after matching); 95\\% CI as above; mean $|$coexpr$|$, over edges / matched non-edges; 50 bins, "
    f"the same with 50 quantile bins (immune: {n_(imm50.n_matched_neg)} of {n_(5 * imm50.m)} non-edges available); 3,000 cells, matching on co-expression computed on 3,000 cells of the same population. "
    "Lower part: co-complex pairs within $V$ (ComplexPortal; hu.MAP~2.0; each also without the pairs that are STRING physical "
    "$\\geq 700$ edges) as the edge set of the same 120-cell family; $r$ and $z$ at the STRING strongest cell. "
    "\\textbf{(E)} Co-pole rate, fraction of edges whose two ends fall in the same pole (top or bottom $K$ genes of the coordinate, "
    "$K = \\lfloor (K/n)\\,n_V \\rfloor$); ratio, observed rate / null mean ($N$ as in A). Upper part: MS null; the interval is "
    "for the ratio of the observed rate to its configuration-model expectation, which differs from the MS null mean by "
    f"{100 * gap_lo:.1f}--{100 * gap_hi:.1f}\\% in these cells. Lower part: $K/n = 0.27$ under both nulls. "
    "\\textbf{(F)} $\\rho$, Spearman correlation of a link's STRING score with the cosine of its two genes' centred, L2-normalised "
    "512-dimensional vectors; $\\rho_{\\mathrm{adj}}$, partial Spearman correlation controlling for the ranks of the summed log degree "
    "(score $\\geq 150$ network within $V$), of the summed $\\log(1 + n_{\\mathrm{obs}})$ and of $|$co-expression$|$ (3,000 cells); "
    "95\\% CI as above. Band AUROC, links of each combined-score band against one non-linked pair per link (no physical link at any "
    "score) matched on the degree deciles of both ends, drawn once (no interval). Rows `all 12': ranges over layers; for "
    "$\\rho_{\\mathrm{adj}}$ the smallest lower confidence bound over layers in brackets; monotone, layers at which the four band "
    "AUROCs increase strictly. "
    f"MS chains: mean acceptance {f(acc_lo, 2)}--{f(acc_hi, 2)}, lag-1 autocorrelation of the strongest cell's null statistic "
    f"{f(lag_lo, 2)}--{f(lag_hi, 2)}. "
    "Vertex-subsampling intervals with the pre-registered scaling factor $\\sqrt{m/n}$, which omits the finite-population "
    "term and so gives half-widths smaller by about a factor $\\sqrt{2}$, are in the CSV (columns ending "
    "\\texttt{\\_prereg}). "
    "Full values: \\texttt{Table\\_S3\\_interaction.csv}."
)

lines = ["\\setcounter{table}{2}",
         "\\begingroup\\footnotesize",
         "\\setlength{\\tabcolsep}{3pt}",
         "\\begin{longtable}{@{}p{\\textwidth}@{}}",
         "\\caption{" + cap + "}\\label{tab:S3}\\\\",
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
lines += A + ["\\midrule"] + B + ["\\midrule"] + C + ["\\midrule"] + D1 + ["\\addlinespace[4pt]"] + D2 + ["\\midrule"] + \
    E1 + ["\\addlinespace[4pt]"] + E2t + ["\\midrule"] + F
lines += ["\\end{longtable}", "\\endgroup", ""]
OUT_TEX.write_text("\n".join(lines))

out = pd.concat(csv_parts, ignore_index=True, sort=False)
front = ["panel", "source_file", "context"]
out = out[front + [c for c in out.columns if c not in front]]
out.to_csv(OUT_CSV, index=False)
print("wrote", OUT_TEX, OUT_CSV, out.shape, "monotone", mono_counts, "copole gap", gap_lo, gap_hi)
