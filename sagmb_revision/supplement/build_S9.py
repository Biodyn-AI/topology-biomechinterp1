"""Build Supplementary Table S9 (cross-extraction analyses) from revision3_experiments/E10_crosscontext.

Writes
  sagmb_revision/supplement/S9.tex                                (LaTeX fragment, one longtable)
  sagmb_revision/supplementary_tables/Table_S9_cross_extraction.csv

All numbers are read from the E10 outputs (varcomp.csv, transfer.csv, lineage_tfs.csv, detection.csv,
kidneyA_vs_B.csv, stability.csv, results.json). Quantities computed here: square roots of the stored Q-profile
upper bounds of tau^2 (to report them as tau), ranges and counts over layers/extractions.
Run:  /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python build_S9.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
E10 = ROOT / "revision3_experiments" / "E10_crosscontext"
OUT_TEX = HERE / "S9.tex"
OUT_CSV = HERE.parent / "supplementary_tables" / "Table_S9_cross_extraction.csv"

MAIN = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
ROBUST = ["immune_seed43", "immune_seed44", "immune_TSP14", "immune_TSP21", "immune_TSP25", "immune_Bcell",
          "immune_CD4T", "immune_CD8T", "kidney_ckpt_whole_human", "kidney_ckpt_kidney", "kidney_ckpt_brain"]
ALL = MAIN + ROBUST
IMM8 = ROBUST[:8]
CKPT = ROBUST[8:]
NAME = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune",
        "kidney_A": "kidney~A", "immune_seed43": "immune (seed 43)", "immune_seed44": "immune (seed 44)",
        "immune_TSP14": "immune (donor TSP14)", "immune_TSP21": "immune (donor TSP21)",
        "immune_TSP25": "immune (donor TSP25)", "immune_Bcell": "immune (B cells)",
        "immune_CD4T": "immune (CD4 T cells)", "immune_CD8T": "immune (CD8 T cells)",
        "kidney_ckpt_whole_human": "kidney (whole-human ckpt)", "kidney_ckpt_kidney": "kidney (kidney ckpt)",
        "kidney_ckpt_brain": "kidney (brain ckpt)"}
STATS = ["S1", "S2", "S2b_full", "S2b_band", "S3"]
STATNAME = {"S1": "S1 interaction", "S2": "S2 TF identity", "S2b_full": "S2b TRRUST, full",
            "S2b_band": "S2b TRRUST, band", "S3": "S3 edges"}
GROUPS = [("resample_immune", "Cell draw, immune"), ("resample_kidney", "Cell draw, kidney"),
          ("donor", "Donor (immune)"), ("celltype", "Cell type (immune)"),
          ("hvg_axis", "Gene axis (kidney~A vs B)"), ("checkpoint", "Checkpoint (same cells)"),
          ("tissue", "Tissue"), ("tissue3_no_lung_external", "Tissue, without lung (ext.)"),
          ("study_lung_vs_lung_external", "Study (lung vs lung (ext.))")]


def num(x, d=3, sign=False):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "--"
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    if s.startswith("-"):
        return "$-$" + s[1:]
    if s.startswith("+"):
        return "$+$" + s[1:]
    return s


def ci(lo, hi, d=3):
    return f"[{num(lo, d)}, {num(hi, d)}]"


def rng(a, b, d=3, sign=False):
    if a >= 0 and b >= 0:
        return f"{num(a, d)}--{num(b, d)}"
    return f"{num(a, d, sign)} to {num(b, d, sign)}"


def pval(p):
    if p >= 0.001:
        return f"{p:.3f}" if p < 0.1 else f"{p:.2f}"
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def cell2(top, bottom):
    return f"{top}\\newline {bottom}"


def intc(n):
    return f"{int(round(n)):,}".replace(",", "{,}")


def medc(x):
    """Median count with thousands separator; keeps a .5 when the median falls between two counts."""
    s = f"{x:,.0f}" if float(x).is_integer() else f"{x:,.1f}"
    return s.replace(",", "{,}")


# ----------------------------------------------------------------------------------------------- load
vc = pd.read_csv(E10 / "varcomp.csv")
tr = pd.read_csv(E10 / "transfer.csv")
lin = pd.read_csv(E10 / "lineage_tfs.csv")
det = pd.read_csv(E10 / "detection.csv")
kab = pd.read_csv(E10 / "kidneyA_vs_B.csv")
stab = pd.read_csv(E10 / "stability.csv")
res = json.loads((E10 / "results.json").read_text())

NC = 7
rows = []


def panel(title):
    if rows:
        rows.append("\\midrule")
    rows.append(f"\\multicolumn{{{NC}}}{{@{{}}p{{16.4cm}}@{{}}}}{{\\footnotesize\\textbf{{{title}}}}}\\\\*")


def header(cells):  # '\\*': no page break between a panel title and its column header
    rows.append(" & ".join(cells) + "\\\\*")
    rows.append("\\cmidrule(lr){1-7}")


def line(cells):
    rows.append(" & ".join(cells) + "\\\\")


def mc(n, text, spec="l"):
    return f"\\multicolumn{{{n}}}{{>{{\\scriptsize}}{spec}}}{{{text}}}"


def mcp(n, width, text):
    return f"\\multicolumn{{{n}}}{{>{{\\scriptsize\\raggedright\\arraybackslash}}p{{{width}}}}}{{{text}}}"


# ======================================================================================= A: variance components
panel("A. Variance components of the layer-mean statistic across the 16 extractions: between-extraction SD $\\tau$ "
      "within each group, and the tissue decision")
header(["Component", "$k$", STATNAME["S1"], STATNAME["S2"], STATNAME["S2b_full"], STATNAME["S2b_band"],
        STATNAME["S3"]])
VM = vc[vc.layer.astype(str) == "mean"]
for g, lab in GROUPS:
    x = VM[VM.group == g].set_index("stat")
    if g == "tissue":
        rows.append("\\cmidrule(lr){1-7}")
    cells = [lab, str(int(x.loc["S1"].k))]
    for s in STATS:
        r = x.loc[s]
        cells.append(num(r.tau))
    line(cells)
TIS = VM[VM.group == "tissue"].set_index("stat")
line(["Tissue $\\tau$, upper Q-profile 95\\% bound", ""] +
     [num(float(np.sqrt(TIS.loc[s].tau2_qprof_hi))) for s in STATS])
DD = VM[VM.group == "D_tissue_minus_within_pooled"].set_index("stat")
rows.append("\\cmidrule(lr){1-7}")
K3 = 1e3  # variance-scale rows are shown x 10^3
line(["$\\tau^2_{\\mathrm{tissue}}$ ($\\times 10^{3}$)", ""] + [num(K3 * DD.loc[s].tau2_tissue, 2) for s in STATS])
line(["Pooled within-tissue floor $\\tau^2_{\\mathrm{within}}$ ($\\times 10^{3}$)", ""] +
     [num(K3 * DD.loc[s].tau2_within_pooled, 2) for s in STATS])
line(["$D = \\tau^2_{\\mathrm{tissue}} - \\tau^2_{\\mathrm{within}}$ ($\\times 10^{3}$) [95\\% interval]", ""] +
     [cell2(num(K3 * DD.loc[s].D, 2, sign=True), ci(K3 * DD.loc[s].D_param_lo, K3 * DD.loc[s].D_param_hi, 2))
      for s in STATS])
line(["$p$ ($H_0$: tissue $=$ floor); tissue effect claimed", ""] +
     [f"{pval(DD.loc[s].p_param_H0_tissue_eq_floor)}; {'yes' if bool(DD.loc[s].tissue_dependence_claimed) else 'no'}"
      for s in STATS])
for s in STATS:
    assert res["tissue_dependence_decision"][s]["claimed"] is False
nlay = {s: res["tissue_dependence_per_layer_D_lower_gt0"][s] for s in STATS}

# ======================================================================================= B: transfer
panel("B. A TF classifier trained on immune genes, applied unchanged to the other extractions (layer mean)")
header(["Test extraction", "TFs / non-TFs (novel genes)", "Novel genes: AUROC (detection-stratified)",
        "Novel genes: range over layers", "In-context classifier, same genes (detection-stratified)",
        "Transfer minus in-context (layers CI $<0$)", "Shared genes: AUROC"])
TM = tr[(tr.centred == True) & (tr.layer.astype(str) == "mean")]
T0 = tr[(tr.centred == True) & (tr.layer.astype(str) == "0")]
oofm = TM[(TM.test == "immune") & (TM.subset == "oof_5fold")].iloc[0]
o0 = T0[(T0.test == "immune") & (T0.subset == "oof_5fold")].iloc[0]
line(["immune (training; 5-fold out-of-fold)", f"{intc(o0.n_pos)} / {intc(o0.n_neg)}",
      f"{num(oofm.auroc)} ({num(oofm.auroc_detstrat)})", "--", "--",
      "--", "--"])
below = res["transfer_novel_layers_transfer_below_incontext_ci"]
for c in ["kidney_B", "lung", "lung_external", "kidney_A"]:
    nov = TM[(TM.test == c) & (TM.subset == "novel")].iloc[0]
    sh = TM[(TM.test == c) & (TM.subset == "shared")].iloc[0]
    n0 = T0[(T0.test == c) & (T0.subset == "novel")].iloc[0]
    perlayer = tr[(tr.centred == True) & (tr.test == c) & (tr.subset == "novel") & (tr.layer.astype(str) != "mean")]
    assert len(perlayer) == 12
    nb = int((perlayer.diff_hi < 0).sum())
    assert nb == below[c]["n_layers_ci_below0"]
    line([NAME[c], f"{intc(n0.n_pos)} / {intc(n0.n_neg)}", f"{num(nov.auroc)} ({num(nov.auroc_detstrat)})",
          cell2(rng(perlayer.auroc.min(), perlayer.auroc.max()), f"lower $\\geq$ {num(perlayer.lo.min())}"),
          f"{num(nov.incontext_auroc)} ({num(nov.incontext_auroc_detstrat)})",
          f"{num(nov.transfer_minus_incontext, sign=True)} ({nb}/12)", num(sh.auroc)])

# ======================================================================================= C: lineage TFs
panel("C. Lineage-specific TF sets: detection-stratified AUROC of set members against non-TF genes (layer mean)")
header(["Extraction", "Immune set: AUROC [95\\% CI]", "Immune set: observed (random)", "Kidney set: AUROC [95\\% CI]",
        "Kidney set: observed (random)", "Lung set: AUROC [95\\% CI]", "Lung set: observed (random)"])
LE = lin[(lin.part == "eval") & (lin.layer.astype(str) == "mean")]


def lin_cells(c, setname):
    r = LE[(LE.context == c) & (LE.tissue_set == setname)]
    assert len(r) == 1
    r = r.iloc[0]
    nobs = int(r.n_observed)
    if np.isnan(r.auroc):
        return ["--", f"{nobs}"]
    mark = "" if bool(r.ci_reliable) else "$^{d}$"
    home = "$^{h}$" if bool(r.home) else ""
    top = f"{num(r.auroc)}{home}"
    return [cell2(top, ci(r.lo_clip, r.hi_clip) + mark), f"{nobs} ({num(r.random_subset_mean)})"]


for c in MAIN + CKPT:
    if c == CKPT[0]:
        rows.append("\\cmidrule(lr){1-7}")
    cells = [NAME[c]]
    for sname in ["immune", "kidney", "lung"]:
        cells += lin_cells(c, sname)
    line(cells)
x = LE[(LE.context.isin(IMM8)) & (LE.tissue_set == "immune")]
cells = ["immune, 8 robustness extractions", rng(x.auroc.min(), x.auroc.max()),
         f"{int(x.n_observed.min())}--{int(x.n_observed.max())}"]
for sname in ["kidney", "lung"]:
    y = LE[(LE.context.isin(IMM8)) & (LE.tissue_set == sname)]
    ok = y.dropna(subset=["auroc"])
    lo_, hi_ = int(y.n_observed.min()), int(y.n_observed.max())
    cells += [("--" if ok.empty else f"{len(ok)} of 8 evaluable"), (f"{lo_}" if lo_ == hi_ else f"{lo_}--{hi_}")]
line(cells)
sp = res["lineage_summary"]["spearman_lineageset_auroc_vs_median_nobs_set"]

# ======================================================================================= D: detection
panel("D. Detection alone: AUROC of $\\log n_{\\mathrm{obs}}$ for separating TFs, and the extraction regime")
header(["Extraction", "$|V|$; cells (input $>512$ genes)", "Median non-zero genes per cell",
        "Median $n_{\\mathrm{obs}}$, Lambert TFs / non-TFs", "Lambert TFs vs non-TFs: AUROC [95\\% CI]",
        "TRRUST TFs vs target-only: AUROC [95\\% CI]", "Lambert TFs with $n_{\\mathrm{obs}} \\geq 10$"])
DA = det[det.part == "a"].set_index(["context", "test"])
DS = det[det.part == "descriptor"].set_index("context")
for c in ALL:
    if c == ROBUST[0]:
        rows.append("\\cmidrule(lr){1-7}")
    a1 = DA.loc[(c, "lognobs_auroc_lambertTF_vs_nonTF")]
    a2 = DA.loc[(c, "lognobs_auroc_trrustTF_vs_targetonly")]
    s = DS.loc[c]
    line([NAME[c], f"{intc(s.n_V)}; {intc(s.n_cells)} ({intc(s.cells_truncated_at_512)})", medc(s.nnz_median),
          f"{s.median_nobs_lambertTF:g} / {s.median_nobs_nonTF:g}", cell2(num(a1.value), ci(a1.lo, a1.hi)),
          cell2(num(a2.value), ci(a2.lo, a2.hi)), f"{100 * s.frac_TF_nobs_ge10:.0f}\\%"])

# ======================================================================================= E: kidney~A vs B
km = res["kidney_mechanism"]
panel("E. The two kidney extractions (one donor, different gene axes): mechanism, matched re-evaluation and the "
      "curated-gene comparison")
rows.append(" & ".join(["Quantity", mcp(2, "4.2cm", "kidney~A"), mcp(2, "4.2cm", "kidney~B"),
                        mcp(2, "4.2cm", "Shared / difference A $-$ B [95\\% CI]")])
            + "\\\\")
rows.append("\\cmidrule(lr){1-7}")


def krow(label, a, b, d):
    rows.append(" & ".join([label, mcp(2, "4.2cm", a), mcp(2, "4.2cm", b), mcp(2, "4.2cm", d)]) + "\\\\")


krow("Gene axis $G$", intc(km["axis_A"]), intc(km["axis_B"]), f"{intc(km['axis_shared'])} shared")
krow("Observed genes $|V|$", intc(km["V_A"]), intc(km["V_B"]), f"{intc(km['M'])} observed in both ($M$)")
krow("Lambert TFs in $V$", intc(km["lambert_tf_V_A"]), intc(km["lambert_tf_V_B"]),
     f"{km['lambert_tf_M']} in $M$; {km['trrust_tf_M']} TRRUST TFs / {km['trrust_target_only_M']} target-only")
krow("Cells (input $>512$ genes)", f"{km['cells_A']} ({km['cells_truncated_A']})",
     f"{km['cells_B']} ({km['cells_truncated_B']})", f"{km['drawn_cells_shared']} cells drawn in both")
krow("Median non-zero genes per cell", medc(km["nnz_median_A"]), medc(km["nnz_median_B"]), "--")
krow("Mean fraction of a cell's detected genes not on the other axis",
     num(km["mean_frac_cell_inputs_outside_other_axis_A"]), num(km["mean_frac_cell_inputs_outside_other_axis_B"]), "--")
krow("Median $n_{\\mathrm{obs}}$ on $M$: Lambert TFs / TRRUST TFs / non-TFs",
     f"{km['median_nobs_lambert_M_A']:g} / {km['median_nobs_trrust_tf_M_A']:g} / {km['median_nobs_nonTF_M_A']:g}",
     f"{km['median_nobs_lambert_M_B']:g} / {km['median_nobs_trrust_tf_M_B']:g} / {km['median_nobs_nonTF_M_B']:g}",
     f"Spearman $\\rho$ of $n_{{\\mathrm{{obs}}}}$ on $M$ {num(km['spearman_nobs_on_M'])}")
rows.append("\\cmidrule(lr){1-7}")
KM = kab[kab.layer.astype(str) == "mean"].set_index("statistic")
for st, lab in [("S2b_band_ownV", "S2b, band $\\SV{2}$--$\\SV{7}$, each extraction's own SVD"),
                ("S2b_band_commonM", "S2b, band $\\SV{2}$--$\\SV{7}$, SVD fitted on $M$"),
                ("S2b_full", "S2b, full dimension"),
                ("S2_ownstrata", "S2, own detection strata"),
                ("S2_sharedstrata", "S2, shared detection strata")]:
    r = KM.loc[st]
    krow(f"Matched genes $M$: {lab}", num(r.A), num(r.B), f"{num(r.delta, sign=True)} {ci(r.delta_lo, r.delta_hi)}")
r = KM.loc["absdelta_S2b_band_ownV_minus_absdelta_full"]
krow("Basis instability: $|\\Delta_{\\mathrm{band}}| - |\\Delta_{\\mathrm{full}}|$ (own SVD)", "--", "--",
     f"{num(r.delta, sign=True)} {ci(r.delta_lo, r.delta_hi)}")
rows.append("\\cmidrule(lr){1-7}")
rows.append(" & ".join(["TRRUST comparison on the curated genes", mcp(2, "4.2cm", "kidney~A, 187 genes (65 / 122)"),
                        mcp(2, "4.2cm", "kidney~A on kidney~B's 167 genes"),
                        mcp(2, "4.2cm", "kidney~B, 167 genes (59 / 108)")]) + "\\\\")
rows.append("\\cmidrule(lr){1-7}")
LBm = res["legacy_bridge_layer_mean"]


def lb_layers(stat, layers):
    x = kab[kab.statistic == stat].copy()
    x["L"] = x.layer.astype(int)
    return x.set_index("L").loc[layers, "A"].tolist()


bz = lb_layers("legacy_bridge_B_curated167_submitted_zero_row_SVD_band", [0, 1, 2, 3])
bz_late = lb_layers("legacy_bridge_B_curated167_submitted_zero_row_SVD_band", list(range(4, 12)))
bo = lb_layers("legacy_bridge_B_curated167_observed_SVD_band", [0, 1, 2, 3])
krow("Band $\\SV{2}$--$\\SV{7}$, SVD including never-observed zero rows (the screen's convention), layer mean",
     num(LBm["legacy_bridge_A_curated187_submitted_zero_row_SVD_band"]),
     num(LBm["legacy_bridge_A_matched167_submitted_zero_row_SVD_band"]),
     cell2(num(LBm["legacy_bridge_B_curated167_submitted_zero_row_SVD_band"]),
           "$L_0$--$L_3$: " + ", ".join(num(v) for v in bz) + f"; $L_4$--$L_{{11}}$: {rng(min(bz_late), max(bz_late))}"))
krow("Same band, observed-gene SVD, layer mean", num(LBm["legacy_bridge_A_curated187_observed_SVD_band"]),
     num(LBm["legacy_bridge_A_matched167_observed_SVD_band"]),
     cell2(num(LBm["legacy_bridge_B_curated167_observed_SVD_band"]), "$L_0$--$L_3$: " + ", ".join(num(v) for v in bo)))
krow("Full dimension, layer mean", num(LBm["legacy_bridge_A_curated187_full_dim"]),
     num(LBm["legacy_bridge_A_matched167_full_dim"]), num(LBm["legacy_bridge_B_curated167_full_dim"]))
hm = res["legacy_bridge_hanley_mcneil_se"]
krow("Hanley--McNeil SE of one AUROC", f"{num(hm['auc_0p61_n65_122'])} (at 0.61)", "--",
     f"{num(hm['auc_0p55_n59_108'])} (at 0.55)")

# ----------------------------------------------------------------------------------------------- caption
caption = (
    "\\textbf{Differences between extractions: variance components, transfer of a TF read-out, lineage TFs, detection, "
    "and the two kidney extractions.} "
    "Statistics (all AUROCs, V-centred embeddings, layer mean over the 12 layers unless stated): S1, cosine of STRING "
    "physical pairs against one degree-preserving rewired non-pair set; S2, the H2 statistic (detection-stratified "
    "LOO-cosine AUROC, Lambert TFs against all other observed genes); S2b, TRRUST TFs against target-only genes, in "
    "the full representation or in the band $\\SV{2}$--$\\SV{7}$; S3, the H3 statistic (within-TF stratified AUROC of the "
    "TF--target cosine). "
    "\\textbf{A}: moment estimate of the between-extraction SD $\\tau$ (AUROC units; 0 when the estimate of $\\tau^2$ is "
    "negative) among the $k$ extractions of each group: cell draws (immune, seed 43, seed 44; kidney~B and the "
    "whole-human checkpoint run), donors TSP14/21/25, B/CD4 T/CD8 T cells, kidney~A vs B, the three checkpoints on the "
    "same 192 cells, and the four tissue-level extractions (immune, kidney~B, lung, lung (ext.)). Sampling covariances "
    "come from one Poisson gene bootstrap shared by all extractions ($B = 1{,}000$; S3 two-way TF $\\times$ target; S1 "
    "node bootstrap, whose SE is within 5\\% of the delete-$d$ vertex-subsampling SE). A tissue effect is claimed only if $D$, the tissue component minus the pooled within-tissue floor (cell "
    "draws, donors, gene axis), has a 95\\% lower bound $>0$ under a parametric random-effects bootstrap ($R = 10{,}000$); "
    f"the lower bound exceeds 0 neither for the layer mean nor at any single layer ({max(nlay.values())} of 12 layers "
    "for every statistic). "
    "Rows marked $\\times 10^{3}$ are in units of $10^{-3}$ AUROC$^2$. Q-profile bounds treat the extractions as "
    "random draws; with 2--4 extractions per group they are wide. "
    "\\textbf{B}: L2-regularised logistic regression (balanced classes; $C$ chosen by 5-fold CV from "
    "$\\{0.01, 0.1, 1, 10\\}$) on unit-normalised $V$-centred embeddings, fitted "
    "per layer on immune genes; novel genes are not observed in immune, shared genes are; the in-context reference is a "
    "5-fold out-of-fold classifier within the test extraction scored on the same genes; intervals are 95\\% gene-bootstrap "
    "intervals ($B = 2{,}000$); AUROCs without stratification, with the detection-stratified AUROC (within deciles of "
    "the test extraction's $\\log n_{\\mathrm{obs}}$, as for S2) in parentheses; below the range over layers, the smallest lower bound "
    "over layers; the difference, its intervals and, in parentheses, the layers whose interval of the difference lies "
    "below 0 are without stratification; per-layer values in the CSV. "
    "\\textbf{C}: sets of the 30 Lambert TFs most specific to each tissue in raw expression (18 for lung, which must be "
    "expressed in both lung studies), selected from raw expression alone, without any embedding statistic; for example immune BACH2, ETS1, BCL11B, "
    "IKZF1, PAX5, EBF1, IRF4; kidney PAX2, PAX8, HNF1B, HNF4A, EMX2, SALL1; lung NKX2-1, FOXA2, FOXF1, GATA6, EPAS1 "
    "(full lists in the CSV). Score: LOO cosine to the observed set members; at least 5 observed members required; centred "
    "Poisson gene-bootstrap intervals ($B = 500$) clipped to $[0, 1]$; observed, members in $V$; random, mean AUROC of "
    "50 random Lambert-TF subsets of the same size; $^{h}$the set's own tissue; $^{d}$interval unreliable (fewer than 10 observed members). "
    f"Across sets and extractions the AUROC falls with the set's median $n_{{\\mathrm{{obs}}}}$ (Spearman $\\rho = "
    f"{sp['rho']:.2f}$, $n = {int(sp['n'])}$; rows not independent). "
    "\\textbf{D}: AUROC of $\\log n_{\\mathrm{obs}}$ alone (no stratification), gene-bootstrap 95\\% intervals "
    "($B = 2{,}000$); values above 0.5 mean TFs are observed in more cells. Cells whose input exceeded scGPT's 512-gene "
    "cap are truncated. "
    "\\textbf{E}: $M$, genes observed in both extractions, with labels recomputed on $M$; paired gene bootstrap on $M$ "
    "($B = 2{,}000$), $\\Delta = $ kidney~A $-$ kidney~B. Curated-gene comparison: TRRUST TFs / targets of the curated gene "
    "list derived from the screen's kidney~A edge list, with the SVD fitted on the full matrix including genes no sampled cell contained (zero rows); the "
    "Hanley--McNeil SE is that of a single AUROC at the stated value and sample sizes. "
    "Full values: \\texttt{Table\\_S9\\_cross\\_extraction.csv}."
)

colspec = ("@{}>{\\scriptsize\\raggedright\\arraybackslash}p{3.1cm}"
           "*{6}{>{\\scriptsize\\raggedright\\arraybackslash}p{2.1cm}}@{}")
tex = ["\\setcounter{table}{8}",
       "\\begingroup", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}",
       "\\renewcommand{\\arraystretch}{1.05}", f"\\begin{{longtable}}{{{colspec}}}",
       f"\\caption{{{caption}}}\\label{{tab:S9}}\\\\", "\\toprule", "\\endfirsthead",
       f"\\multicolumn{{{NC}}}{{@{{}}l}}{{\\footnotesize Table~\\thetable\\ (continued)}}\\\\", "\\toprule", "\\endhead",
       "\\midrule", f"\\multicolumn{{{NC}}}{{r@{{}}}}{{\\footnotesize continued on next page}}\\\\", "\\endfoot",
       "\\bottomrule", "\\endlastfoot"]
tex += rows
tex += ["\\end{longtable}", "\\endgroup"]
OUT_TEX.write_text("\n".join(tex) + "\n")

# ----------------------------------------------------------------------------------------------- CSV
parts = []


def add(df, panel_id, source, context=None):
    df = df.copy()
    if "context" not in df.columns:
        if context is not None:
            df.insert(0, "context", context)
        elif "contexts" in df.columns:
            df.insert(0, "context", df["contexts"])
        elif "test" in df.columns:
            df.insert(0, "context", df["test"])
        elif "ctx_a" in df.columns:
            df.insert(0, "context", df["ctx_a"].astype(str) + " vs " + df["ctx_b"].astype(str))
    df.insert(0, "source_file", source)
    df.insert(0, "panel", panel_id)
    parts.append(df)


vc2 = vc.copy()
vc2["tau_qprof_hi"] = np.sqrt(vc2.tau2_qprof_hi.clip(lower=0))
add(vc2, "A_variance_components", "E10_crosscontext/varcomp.csv (tau_qprof_hi = sqrt(tau2_qprof_hi), added)")
add(tr, "B_transfer", "E10_crosscontext/transfer.csv")
lin2 = lin.copy()
lin2.loc[lin2.part == "selection", "context"] = "selection (raw expression of 4 datasets)"
lin2.loc[lin2.part == "summary", "context"] = "all extractions"
add(lin2, "C_lineage_tfs", "E10_crosscontext/lineage_tfs.csv (part = selection / eval / S2_vs_nobs / summary)")
add(det, "D_detection", "E10_crosscontext/detection.csv (part a: log n_obs AUROC; b: logistic interaction; "
                        "c: S2 restricted by n_obs; descriptor)")
kab2 = kab.copy()  # neutral name for the screen's zero-row SVD convention (values unchanged)
kab2["statistic"] = (kab2["statistic"].str.replace("submitted_zero_row", "screen_zero_row", regex=False)
                     .str.replace("legacy_bridge_", "curated_bridge_", regex=False))
add(kab2, "E_kidneyA_vs_B", "E10_crosscontext/kidneyA_vs_B.csv (statistic names: 'curated_bridge_' = the curated-gene "
    "comparison; '_screen_zero_row_' = the screen's SVD convention including never-observed zero rows)",
    context="kidney_A vs kidney_B")
kmd = pd.DataFrame([dict(quantity=k, value=v) for k, v in km.items()])
add(kmd, "E_kidney_mechanism", "E10_crosscontext/results.json kidney_mechanism", context="kidney_A vs kidney_B")
hmd = pd.DataFrame([dict(quantity=f"hanley_mcneil_se_{k}", value=v) for k, v in hm.items()])
add(hmd, "E_kidney_mechanism", "E10_crosscontext/results.json (Hanley-McNeil SE, curated-gene comparison)",
    context="kidney_A vs kidney_B")
add(stab, "extra_stability", "E10_crosscontext/stability.csv (neighbour Jaccard / RSA; not shown in the PDF table)")
out = pd.concat(parts, ignore_index=True, sort=False)
_txt = " ".join(out.astype(str).values.ravel()).lower() + " " + " ".join(out.columns).lower()
for bad in ("legacy", "submitted", "revision", "verifier", "round-"):
    assert bad not in _txt, bad
out.to_csv(OUT_CSV, index=False)
print("wrote", OUT_TEX, OUT_CSV, out.shape)
