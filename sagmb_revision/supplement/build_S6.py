"""Build Supplementary Table S6 (TF identity, H2) from revision3_experiments/E3_node.

Writes
  sagmb_revision/supplement/S6.tex                         (LaTeX fragment, one longtable)
  sagmb_revision/supplementary_tables/Table_S6_tf_identity.csv

Every number is read from the E3 output files; the only quantities computed here are ranges, minima, maxima and
counts over layers/extractions, and the layer mean of the stored per-layer out-of-fold AUROCs of the ledger.
Run:  /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python build_S6.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
E3 = ROOT / "revision3_experiments" / "E3_node"
OUT_TEX = HERE / "S6.tex"
OUT_CSV = HERE.parent / "supplementary_tables" / "Table_S6_tf_identity.csv"

MAIN = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
ROBUST = ["immune_seed43", "immune_seed44", "immune_TSP14", "immune_TSP21", "immune_TSP25", "immune_Bcell",
          "immune_CD4T", "immune_CD8T", "kidney_ckpt_whole_human", "kidney_ckpt_kidney", "kidney_ckpt_brain"]
ALL = MAIN + ROBUST
NAME = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune",
        "kidney_A": "kidney~A", "immune_seed43": "immune (seed 43)", "immune_seed44": "immune (seed 44)",
        "immune_TSP14": "immune (donor TSP14)", "immune_TSP21": "immune (donor TSP21)",
        "immune_TSP25": "immune (donor TSP25)", "immune_Bcell": "immune (B cells)",
        "immune_CD4T": "immune (CD4 T cells)", "immune_CD8T": "immune (CD8 T cells)",
        "kidney_ckpt_whole_human": "kidney (whole-human ckpt)", "kidney_ckpt_kidney": "kidney (kidney ckpt)",
        "kidney_ckpt_brain": "kidney (brain ckpt)"}


# ----------------------------------------------------------------------------------------------- formatting
def num(x, d=3, sign=False):
    """Number in math mode with a true minus sign; sign=True prints '+' for positive values."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "--"
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    if s.startswith("-"):
        return "$-$" + s[1:]
    if s.startswith("+"):
        return "$+$" + s[1:]
    return s


def ci(lo, hi, d=3, sign=False):
    """Interval [lower, upper]; only negative bounds carry a sign."""
    return f"[{num(lo, d)}, {num(hi, d)}]"


def rng(a, b, d=3, sign=False):
    """Range over layers/extractions: en-dash when both bounds are non-negative, 'to' otherwise."""
    if a >= 0 and b >= 0:
        return f"{num(a, d)}--{num(b, d)}"
    return f"{num(a, d, sign)} to {num(b, d, sign)}"


def pval(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "--"
    if p >= 0.001:
        return f"{p:.3f}" if p < 0.1 else f"{p:.2f}"
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def cell2(top, bottom):
    return f"{top}\\newline {bottom}"


def intc(n):
    return f"{int(n):,}".replace(",", "{,}")


def band_tex(b):
    """'SV3-12' -> 'SV$_{3}$--SV$_{12}$' style using the preamble macro."""
    b = str(b).replace("SV", "")
    if "-" in b:
        i, j = b.split("-")
        return f"$\\SV{{{i}}}$--$\\SV{{{j}}}$"
    return f"$\\SV{{{b}}}$"


def choice_tex(s):
    """'SV4-12@L2 (57)' -> band at layer with count."""
    first = s.split(";")[0].strip()
    cell, cnt = first.split(" ")
    band, lay = cell.split("@")
    return f"{band_tex(band)} at $L_{{{lay[1:]}}}$ {cnt}"


# ----------------------------------------------------------------------------------------------- load
summ = pd.read_csv(E3 / "h2_summary.csv").set_index("context")
prim = {c: pd.read_csv(E3 / f"h2_primary_{c}.csv") for c in ALL}
base = pd.read_csv(E3 / "baselines.csv")
paired = pd.read_csv(E3 / "baselines_paired.csv")
logit = pd.read_csv(E3 / "logistic_adjustment.csv")
restr = pd.read_csv(E3 / "restricted_V.csv")
bmax = pd.read_csv(E3 / "bandscan_maxT.csv")
bands = {c: pd.read_csv(E3 / f"bandscan_{c}.csv") for c in ALL}
legacy = pd.read_csv(E3 / "screen_bands.csv")   # copy of legacy_bands.csv
cfit = pd.read_csv(E3 / "crossfit.csv")
trr = pd.read_csv(E3 / "trrust_task.csv")
ledger = json.loads((E3 / "ledger_0744.json").read_text())

rows = []          # LaTeX body lines
NC = 7             # columns of the longtable


def panel(title):
    rows.append("\\midrule")
    rows.append(f"\\multicolumn{{{NC}}}{{@{{}}p{{16.4cm}}@{{}}}}{{\\footnotesize\\textbf{{{title}}}}}\\\\")


def header(cells):
    rows.append(" & ".join(cells) + "\\\\")
    rows.append("\\cmidrule(lr){1-7}")


def line(cells):
    assert len(cells) == NC, cells
    rows.append(" & ".join(cells) + "\\\\")


EXT_HDR = ["", "lung", "lung (ext.)", "kidney~B", "immune", "kidney~A", "all 16"]

# ======================================================================================= Panel A: primary H2
rows.append(f"\\multicolumn{{{NC}}}{{@{{}}p{{16.4cm}}@{{}}}}{{\\footnotesize\\textbf{{A. Pre-registered H2 statistic "
            f"(full 512-dimensional representation), per extraction}}}}\\\\")
header(["Extraction", "$|V|$ / TFs", "AUROC, 12 layers", "Best layer: AUROC [95\\% CI]", "$z$, 12 layers",
        "Smallest max-$T$ $p$", "Layers max-$T$ $p<0.05$"])
for c in ALL:
    s = summ.loc[c]
    lay = prim[c][prim[c].row_type == "H2_layer"]
    zmin, zmax = lay.z.min(), lay.z.max()
    b = int(s.best_layer)
    best = lay[lay.layer == b].iloc[0]
    assert abs(best.strat_auroc - s.strat_auroc_max) < 1e-12
    if c == ROBUST[0]:
        rows.append("\\cmidrule(lr){1-7}")
    line([NAME[c], f"{intc(s.n_V)} / {intc(s.n_tf)}", rng(s.strat_auroc_min, s.strat_auroc_max),
          cell2(f"$L_{{{b}}}$: {num(best.strat_auroc)}", ci(s.best_strat_lo, s.best_strat_hi)),
          f"{zmin:.1f}--{zmax:.1f}", pval(s.min_p_maxT), f"{int(s.n_layers_maxT_lt_0p05)}/12"])

# ======================================================================================= Panel B: baselines
panel("B. Baselines: detection-stratified AUROC of single-gene scores and static embeddings on the same genes and strata")
header(EXT_HDR)
B = base.set_index(["context", "baseline"])


def brow(label, key, with_ci=True, show_unstrat=False):
    cells = [label]
    for c in MAIN:
        r = B.loc[(c, key)]
        top = num(r.strat_auroc)
        if show_unstrat:
            top += f" ({num(r.auroc)})"
        cells.append(cell2(top, ci(r.strat_lo, r.strat_hi)) if with_ci else top)
    vals = base[base.baseline == key].strat_auroc
    cells.append(rng(vals.min(), vals.max()))
    line(cells)


brow("Detection alone, $\\log n_{\\mathrm{obs}}$ (unstratified AUROC)", "log_nobs", with_ci=True, show_unstrat=True)
brow("Mean expression, extraction cells", "log_mean_expr_extraction")
brow("Mean expression, population", "log_mean_expr_population")
brow("scGPT input gene embedding, $V$-centred", "static_scgpt")
brow("scGPT input gene embedding, uncentred", "static_scgpt_uncentred", with_ci=False)
brow("ESM2-15B per-protein embedding, $V$-centred", "static_esm2")
brow("ESM2-15B per-protein embedding, uncentred", "static_esm2_uncentred", with_ci=False)
cells = ["Genes with an ESM2-15B vector"]
for c in MAIN:
    cells.append(intc(B.loc[(c, "static_esm2")].n_genes))
v = base[base.baseline == "static_esm2"].n_genes
cells.append(f"{intc(v.min())}--{intc(v.max())}")
line(cells)
for key, label in [("scgpt", "Input embedding minus layer (paired; 12 layers)"),
                   ("esm2", "ESM2-15B minus layer (paired; 12 layers)")]:
    cells = [label]
    for c in MAIN:
        p = paired[(paired.context == c) & (paired.static == key)]
        cells.append(cell2(rng(p.delta.min(), p.delta.max(), sign=True), f"lower $\\geq$ {num(p.delta_lo.min())}"))
    p = paired[paired.static == key]
    n0 = int(((p.delta_lo <= 0) & (p.delta_hi >= 0)).sum())
    cells.append(cell2(rng(p.delta.min(), p.delta.max(), sign=True), f"{n0}/{len(p)} CIs include 0"))
    line(cells)

# ======================================================================================= Panel C: logistic + restricted V
panel("C. Adjustment: logistic regression TF $\\sim$ score + $\\log n_{\\mathrm{obs}}$ + $\\log$ mean expression "
      "(standardised), per layer; and $V$ restricted by detection")
header(EXT_HDR)
for cov, lab in [("extraction", "extraction-cell expression"), ("population", "population expression")]:
    L = logit[logit.expr_covariate == cov]
    cells = [f"$\\beta$ per SD of score, {lab}"]
    for c in MAIN:
        x = L[L.context == c]
        cells.append(rng(x.beta_score_per_sd.min(), x.beta_score_per_sd.max()))
    cells.append(rng(L.beta_score_per_sd.min(), L.beta_score_per_sd.max()))
    line(cells)
    if cov == "extraction":
        cells = ["Smallest lower 95\\% bound of $\\beta$"]
        for c in MAIN:
            cells.append(num(L[L.context == c].beta_lo.min()))
        cells.append(num(L.beta_lo.min()))
        line(cells)
    cells = [f"Largest likelihood-ratio $p$, {lab}"]
    for c in MAIN:
        cells.append(pval(L[L.context == c].LR_p.max()))
    cells.append(pval(L.LR_p.max()))
    line(cells)
RL = restr[restr.row_type == "H2_layer"]
for m in [3, 5, 10]:
    cells = [f"H2 AUROC, genes with $n_{{\\mathrm{{obs}}}} \\geq {m}$ (12 layers)"]
    for c in MAIN:
        x = RL[(RL.context == c) & (RL.min_nobs == m)]
        cells.append(cell2(rng(x.strat_auroc.min(), x.strat_auroc.max()),
                           f"{intc(x.n_genes.iloc[0])} / {intc(x.n_tf.iloc[0])}"))
    x = RL[RL.min_nobs == m]
    cells.append(rng(x.strat_auroc.min(), x.strat_auroc.max()))
    line(cells)

# ======================================================================================= Panel D: single directions
panel("D. Single directions $\\SV{1}$--$\\SV{12}$ (raw coordinate, sign-free AUROC), best of the 12 layers")
header(["Direction", "lung: AUROC (layer)", "lung: $z$", "lung: $p$ (936 cells)", "lung (ext.): AUROC (layer)",
        "lung (ext.): $z$", "lung (ext.): $p$ (936 cells)"])
single = {}
for c in ["lung", "lung_external"]:
    bc = bands[c]
    w1 = bc[(bc.row_type == "band_cell") & (bc.task == "lambert") & (bc.width == 1)]
    best = w1.loc[w1.groupby("i").strat_auroc_cos.idxmax()].set_index("i")
    single[c] = best
for k in range(1, 13):
    cells = [f"$\\SV{{{k}}}$"]
    for c in ["lung", "lung_external"]:
        r = single[c].loc[float(k)]
        cells += [f"{num(r.strat_auroc_cos)} ($L_{{{int(r.layer)}}}$)", f"{r.z_cos:.1f}", pval(r.p_maxT_cos)]
    line(cells)
cells = ["Full 512, best layer$^{a}$"]
for c in ["lung", "lung_external"]:
    s = summ.loc[c]
    cells += [f"{num(s.strat_auroc_max)} ($L_{{{int(s.best_layer)}}}$)", f"{s.best_z:.1f}", pval(s.min_p_maxT)]
line(cells)

# ======================================================================================= Panel E: band scan, cross-fit
panel("E. All 78 contiguous bands of $\\SV{1}$--$\\SV{12}$ at 12 layers (936 cells), screen bands, and cross-fitted band choice")
header(EXT_HDR)
bm = bmax[(bmax.task == "lambert")]
COS = bm[bm.variant == "cos"].set_index("context")
LDA = bm[bm.variant == "lda"].set_index("context")
cells = ["Cells with max-$T$ $p<0.05$ (of 936)"] + [str(int(COS.loc[c].n_cells_maxT_lt_0p05)) for c in MAIN]
cells.append(f"{int(COS.loc[ALL].n_cells_maxT_lt_0p05.min())}--{int(COS.loc[ALL].n_cells_maxT_lt_0p05.max())}")
line(cells)
cells = ["Best cell by AUROC"]
for c in MAIN:
    r = COS.loc[c]
    cells.append(cell2(f"{band_tex(r.best_auroc_band)}, $L_{{{int(r.best_auroc_layer)}}}$", num(r.best_auroc)))
cells.append(rng(COS.loc[ALL].best_auroc.min(), COS.loc[ALL].best_auroc.max()))
line(cells)
cells = ["Cell with the largest $z$: AUROC ($z$)"]
for c in MAIN:
    r = COS.loc[c]
    cells.append(cell2(f"{band_tex(r.best_z_band)}, $L_{{{int(r.best_z_layer)}}}$",
                       f"{num(r.best_z_strat_auroc)} ({r.best_z:.1f})"))
cells.append("--")
line(cells)
cells = ["Null 97.5\\% quantile of the maximum AUROC"] + [num(COS.loc[c].null_max_auroc_q97_5) for c in MAIN]
cells.append(rng(COS.loc[ALL].null_max_auroc_q97_5.min(), COS.loc[ALL].null_max_auroc_q97_5.max()))
line(cells)
cells = ["Centroid (LDA) variant, best cell (in-sample)"]
for c in MAIN:
    r = LDA.loc[c]
    cells.append(cell2(f"{band_tex(r.best_auroc_band)}, $L_{{{int(r.best_auroc_layer)}}}$", num(r.best_auroc)))
cells.append(rng(LDA.loc[ALL].best_auroc.min(), LDA.loc[ALL].best_auroc.max()))
line(cells)
# screen bands: nominal observed-index convention (kidney~A: the bands read in the reference decomposition of the full
# matrix, zero rows included, mapped to the observed basis; the screen's own later kidney~A band is quoted in the caption)
LB = legacy[(legacy.row_type == "legacy_band") & (legacy.task == "lambert") & (legacy.variant == "cos")]
for band in ["SV2-4", "SV5-7", "SV2-7"]:
    cells = [f"Screen band {band_tex(band)}: AUROC, 12 layers (layers $p<0.05$; smallest $p$ of 936)"]
    for c in MAIN:
        conv = "mapped:legacy_full_zero" if c == "kidney_A" else "nominal_observed_index"
        x = LB[(LB.context == c) & (LB.legacy_band == band) & (LB.convention == conv)]
        assert len(x) == 12, (c, band, len(x))
        lab = rng(x.strat_auroc.min(), x.strat_auroc.max())
        if c == "kidney_A":
            lab = band_tex(x.mapped_band.iloc[0]) + "$^{b}$: " + lab
        cells.append(cell2(lab, f"({int((x.p_maxT_band_12layers < 0.05).sum())}; {pval(x.p_maxT_936.min())})"))
    cells.append("--")
    line(cells)
CF = cfit[cfit.task == "lambert"].set_index("context")
for lab, m, lo, hi, sg in [("Cross-fitted band, held-out AUROC", "heldout_median", "heldout_lo", "heldout_hi", False),
                           ("Full 512, held-out AUROC", "full512_heldout_median", "full512_heldout_lo",
                            "full512_heldout_hi", False),
                           ("Band minus full 512, held-out", "band_minus_full512_median", "band_minus_full512_lo",
                            "band_minus_full512_hi", True)]:
    cells = [lab + " (median [2.5--97.5\\% of 200 halvings])"]
    for c in MAIN:
        r = CF.loc[c]
        cells.append(cell2(num(r[m], sign=sg), ci(r[lo], r[hi], sign=sg)))
    cells.append(rng(CF.loc[ALL][m].min(), CF.loc[ALL][m].max(), sign=sg))
    line(cells)
cells = ["Optimism (in-sample minus held-out), median"] + [num(CF.loc[c].optimism_median) for c in MAIN]
cells.append(rng(CF.loc[ALL].optimism_median.min(), CF.loc[ALL].optimism_median.max()))
line(cells)
cells = ["Most frequent choice (of 200)"] + [choice_tex(CF.loc[c].top_choices) for c in MAIN] + ["--"]
line(cells)

# ======================================================================================= Panel F: 0.744 ledger and TRRUST
panel("F. The reported 0.744 (immune; screen labels: 73 TFs against 222 target-only genes) and the TRRUST task")
sp = ledger["screen_procedure"]
oof_pooled = float(np.mean([r["oof_auroc"] for r in sp["per_layer"]]))
lc = ledger["loo_cos_SV2_7"]
lf = ledger["loo_cos_full512"]
strat_layers = [r["auroc"] for r in lc["stratified"]["per_layer"]]
imax, imin = int(np.argmax(strat_layers)), int(np.argmin(strat_layers))
sv = sp["seed_variability"]


def frow(label, value, note):
    rows.append(f"{label} & {value} & \\multicolumn{{5}}{{>{{\\scriptsize\\raggedright\\arraybackslash}}p{{11.2cm}}@{{}}}}{{{note}}}\\\\")


rows.append(f"Step & Layer-mean AUROC & \\multicolumn{{5}}{{>{{\\scriptsize}}l}}{{Procedure and remarks}}\\\\")
rows.append("\\cmidrule(lr){1-7}")
frow("Value reported by the screen", num(ledger["preprint_value"]), "As recorded in the screen's ledger.")
frow("Screen procedure, reproduced", num(sp["mean_over_layers"]),
     "Logistic regression ($C = 1$) on standardised $\\SV{2}$--$\\SV{7}$ coordinates, 5-fold stratified CV (seed 42), "
     f"fold-mean AUROC averaged over 12 layers; over CV seeds 0--49 median {num(sv['mean_over_layers_by_seed_median'])}, "
     f"2.5--97.5\\% {num(sv['lo'])}--{num(sv['hi'])}; permutation $p = 1/1{{,}}001$ ($N = 1{{,}}000$, the floor); "
     f"best layer $L_{{{sp['argmax_layer']}}}$ {num(sp['max_over_layers'])}.")
frow("Same out-of-fold predictions, pooled", num(oof_pooled), "Unstratified AUROC of the pooled out-of-fold predictions.")
frow("Same out-of-fold predictions, stratified", num(sp["oof_strat_mean_over_layers"]),
     "Detection-stratified AUROC (deciles of $\\log n_{\\mathrm{obs}}$).")
frow("LOO cosine in $\\SV{2}$--$\\SV{7}$, unstratified", num(lc["unstratified"]["mean_over_layers"]),
     f"Unsupervised score; layer-mean permutation $p$ {pval(lc['unstratified']['p_mean_over_layers'])} ($N = 10{{,}}000$).")
frow("LOO cosine in $\\SV{2}$--$\\SV{7}$, stratified", num(lc["stratified"]["mean_over_layers"]),
     f"Layer-mean permutation $p$ {pval(lc['stratified']['p_mean_over_layers'])}; per layer highest at $L_{{{imax}}}$ "
     f"({num(strat_layers[imax])}), lowest at $L_{{{imin}}}$ ({num(strat_layers[imin])}).")
frow("LOO cosine, full 512, unstratified", num(lf["unstratified"]["mean_over_layers"]),
     f"Layer-mean permutation $p$ {pval(lf['unstratified']['p_mean_over_layers'])}.")
frow("LOO cosine, full 512, stratified", num(lf["stratified"]["mean_over_layers"]),
     f"Layer-mean permutation $p$ {pval(lf['stratified']['p_mean_over_layers'])}.")
frow("Attribution of the 0.747", "--",
     "Order-dependent: removing supervision first costs 0.10 and then stratifying 0.09 (0.747 $\\to$ 0.644 $\\to$ 0.553); "
     "stratifying first costs 0.06 and then removing supervision 0.13 (0.747 / 0.745 $\\to$ 0.683 $\\to$ 0.553). "
     "Supervision accounts for 0.10--0.13, detection for 0.06--0.09.")
rows.append("\\cmidrule(lr){1-7}")
header(["TRRUST task: TRRUST TFs against target-only genes"] + EXT_HDR[1:])
TL = trr[trr.row_type == "layer"]
TB = trr[trr.row_type == "baseline"].set_index(["context", "baseline"])
cells = ["Genes / TFs"]
for c in MAIN:
    x = TL[TL.context == c].iloc[0]
    cells.append(f"{intc(x.n_genes)} / {intc(x.n_tf)}")
cells.append(f"{intc(TL.n_genes.min())}--{intc(TL.n_genes.max())} / {intc(TL.n_tf.min())}--{intc(TL.n_tf.max())}")
line(cells)
cells = ["Full 512, AUROC over 12 layers (layers max-$T$ $p<0.05$)"]
for c in MAIN:
    x = TL[TL.context == c]
    cells.append(cell2(rng(x.strat_auroc.min(), x.strat_auroc.max()), f"({int((x.p_maxT < 0.05).sum())}/12)"))
cells.append(rng(TL.strat_auroc.min(), TL.strat_auroc.max()))
line(cells)
for key, lab in [("log_mean_expr_extraction", "Mean expression alone"), ("log_nobs", "Detection alone")]:
    cells = [lab] + [num(TB.loc[(c, key)].strat_auroc) for c in MAIN]
    v = trr[(trr.row_type == "baseline") & (trr.baseline == key)].strat_auroc
    cells.append(rng(v.min(), v.max()))
    line(cells)
CT = cfit[cfit.task == "trrust"].set_index("context")
cells = ["Cross-fitted band minus full 512 (median [2.5--97.5\\%])"]
for c in MAIN:
    r = CT.loc[c]
    cells.append(cell2(num(r.band_minus_full512_median, sign=True),
                       ci(r.band_minus_full512_lo, r.band_minus_full512_hi, sign=True)))
cells.append(rng(CT.band_minus_full512_median.min(), CT.band_minus_full512_median.max(), sign=True))
line(cells)

# ----------------------------------------------------------------------------------------------- caption, assembly
# the screen's later kidney~A iterations used SV5-SV7 of the observed genes, which maps to itself (Table S14)
KA57 = legacy[(legacy.row_type == "legacy_band") & (legacy.task == "lambert") & (legacy.variant == "cos")
              & (legacy.context == "kidney_A") & (legacy.legacy_band == "SV5-7") & (legacy.convention == "mapped:screen_nonzero")]
assert len(KA57) == 12 and set(KA57.mapped_band) == {"SV5-7"}
KA57_TXT = (f"AUROC {rng(KA57.strat_auroc.min(), KA57.strat_auroc.max())} over the 12 layers "
            f"({int((KA57.p_maxT_band_12layers < 0.05).sum())} layers with $p<0.05$; smallest $p$ of 936, "
            f"{pval(KA57.p_maxT_936.min())})")
caption = (
    "\\textbf{Transcription-factor identity (H2): primary statistic, baselines, adjustment, spectral bands and the "
    "reported 0.744.} "
    "Score: leave-one-out (LOO) mean cosine of a gene's $V$-centred vector to the census TFs of Lambert et al.\\ (excluding "
    "itself); statistic: AUROC of TFs against all other observed genes $V$, stratified by deciles of "
    "$\\log n_{\\mathrm{obs}}$ ($n_{\\mathrm{obs}}$, cells whose input contained the gene) unless marked unstratified. "
    "$|V|$, observed genes; $L_\\ell$, layer; $z$, standardised distance from the within-stratum permutation null; "
    "max-$T$, Westfall--Young correction over the 12 layers (primary, $N = 10{,}000$ joint permutations, floor "
    "$1/10{,}001 = 1.0\\times10^{-4}$) or over the 936 band cells (Panels D, E; $N = 2{,}000$, floor $5.0\\times10^{-4}$). "
    "Intervals are 95\\% gene-bootstrap percentile intervals ($B = 1{,}000$; LOO scores held fixed), written [lower, upper]; "
    "ranges over layers or extractions use an en-dash (or `to' when a bound is negative). "
    "Because the TF centroid is not re-estimated in the bootstrap, a half-sampling SE that re-estimates it is 0.83--0.84 times "
    "the bootstrap SE in lung and lung (ext.) but 1.15--1.32 times it in the two immune extractions checked (immune, CD4 T "
    "cells), so intervals in the immune extractions are about 15--30\\% too narrow. "
    "\\textbf{A}: all 16 extractions (the 11 robustness extractions below the rule; the three checkpoints are applied to the "
    "same 192 kidney cells); the decision uses lung and lung (ext.). "
    "\\textbf{B}: the same statistic for single-gene scores (higher value scores higher) and for static embeddings "
    "scored by the same LOO cosine; unstratified AUROC in parentheses for detection; ESM2-15B is the 5{,}120-dimensional "
    "per-protein embedding (UCE file), matched for 94--96\\% of $V$ in the lung extractions; paired differences (static "
    "minus residual-stream layer) use identical genes, labels and bootstrap draws, with the range over the 12 layers and "
    "the lowest lower bound; `all 16' gives the range over extractions (and, for the paired rows, how many of the "
    "192 extraction--layer intervals include 0). "
    "\\textbf{C}: standardised coefficient of the score, its smallest lower 95\\% Wald bound and the largest "
    "likelihood-ratio $p$ over the 12 layers; restricted-$V$ rows re-centre and re-stratify on the retained genes "
    "($N = 2{,}000$; genes / TFs below the range). "
    "\\textbf{D}: each direction scored by its raw coordinate with a sign-free AUROC, at the layer where it is highest; "
    "$p$ is the 936-cell max-$T$ $p$ of that cell; $^{a}$primary analysis (Panel A): $z$ and 12-layer max-$T$ $p$ from $N = 10{,}000$ permutations. "
    "\\textbf{E}: cos variant (LOO cosine within the band) unless stated; the centroid (LDA) variant is fitted in-sample "
    "and refitted in each null draw, so its maximum is optimistic. Screen bands are those of the automated screen, "
    "indexed on the observed-gene basis; in parentheses, the number of layers with $p<0.05$ on the band's own 12-layer "
    "max-$T$ and the smallest $p$ in the 936-cell max-$T$; $^{b}$kidney~A: the band read in the reference decomposition of the "
    "full kidney~A matrix, zero rows included, and mapped to the observed-gene basis, in which its $\\SV{k}$ is observed "
    "$\\SV{k-1}$ (Table~S14). The screen assembled $\\SV{2}$--$\\SV{4}$ on a 209-gene decomposition whose band also maps "
    "to observed $\\SV{1}$--$\\SV{3}$, and its later kidney~A iterations on mean embeddings used $\\SV{5}$--$\\SV{7}$ of the "
    "observed genes, which maps to itself: " + KA57_TXT + ". Cross-fit: the best cell chosen on a random half of the genes and scored on the other half, 200 halvings; "
    "the brackets are the 2.5--97.5\\% quantiles over halvings of one data set, not sampling intervals. "
    "Immune and kidney~A are the extractions in which the screen's results were found. "
    "\\textbf{F}: 361 genes named by the screen, 295 observed; per-layer values in the CSV. TRRUST task: TRRUST TFs against "
    "genes that are only ever TRRUST targets, same machinery as Panel A. "
    "Per-layer values of every analysis and all 936 band cells of the five main extractions are in the CSV. "
    "Full values: \\texttt{Table\\_S6\\_tf\\_identity.csv}."
)

colspec = ("@{}>{\\scriptsize\\raggedright\\arraybackslash}p{3.1cm}"
           "*{6}{>{\\scriptsize\\raggedright\\arraybackslash}p{2.1cm}}@{}")
tex = []
tex.append("\\setcounter{table}{5}")
tex.append("\\begingroup")
tex.append("\\footnotesize")
tex.append("\\setlength{\\tabcolsep}{3pt}")
tex.append("\\renewcommand{\\arraystretch}{1.05}")
tex.append(f"\\begin{{longtable}}{{{colspec}}}")
tex.append(f"\\caption{{{caption}}}\\label{{tab:S6}}\\\\")
tex.append("\\toprule")
tex.append("\\endfirsthead")
tex.append(f"\\multicolumn{{{NC}}}{{@{{}}l}}{{\\footnotesize Table~\\thetable\\ (continued)}}\\\\")
tex.append("\\toprule")
tex.append("\\endhead")
tex.append("\\midrule")
tex.append(f"\\multicolumn{{{NC}}}{{r@{{}}}}{{\\footnotesize continued on next page}}\\\\")
tex.append("\\endfoot")
tex.append("\\bottomrule")
tex.append("\\endlastfoot")
tex.extend(rows)
tex.append("\\end{longtable}")
tex.append("\\endgroup")
OUT_TEX.write_text("\n".join(tex) + "\n")

# ----------------------------------------------------------------------------------------------- CSV
parts = []


def add(df, panel_id, source, **extra):
    df = df.copy()
    for k, v in extra.items():
        df.insert(0, k, v)
    df.insert(0, "source_file", source)
    df.insert(0, "panel", panel_id)
    parts.append(df)


for c in ALL:
    add(prim[c], "A_primary_H2_per_layer", f"E3_node/h2_primary_{c}.csv", context=c)
add(summ.reset_index(), "A_summary", "E3_node/h2_summary.csv")
add(base, "B_baselines", "E3_node/baselines.csv")
add(paired, "B_paired_static_minus_layer", "E3_node/baselines_paired.csv")
add(logit, "C_logistic_adjustment", "E3_node/logistic_adjustment.csv")
add(restr, "C_restricted_V", "E3_node/restricted_V.csv")
for c in MAIN:  # all 936 cells of the five main extractions; the 16-extraction summary is bandscan_maxT.csv
    bc = bands[c]
    add(bc[bc.row_type == "band_cell"], "D_E_bandscan_cells", f"E3_node/bandscan_{c}.csv")
add(bmax, "E_bandscan_maxT_summary", "E3_node/bandscan_maxT.csv")
add(legacy[legacy.row_type == "legacy_band"], "E_screen_bands", "E3_node/screen_bands.csv")
add(cfit, "E_F_crossfit", "E3_node/crossfit.csv")
add(trr, "F_trrust_task", "E3_node/trrust_task.csv")
# ledger (flattened)
led = []
for r in sp["per_layer"]:
    led.append(dict(context="immune", ledger_quantity="screen_procedure", layer=r["layer"], cv_foldmean=r["cv_foldmean"],
                    oof_auroc=r["oof_auroc"], oof_strat_auroc=r["oof_strat_auroc"]))
for key, d in [("loo_cos_SV2_7", lc), ("loo_cos_full512", lf)]:
    for strat in ["unstratified", "stratified"]:
        for r in d[strat]["per_layer"]:
            led.append(dict(context="immune", ledger_quantity=f"{key}_{strat}", layer=r["layer"], ledger_auroc=r["auroc"],
                            ledger_lo=r["lo"], ledger_hi=r["hi"], z=r["z"], p_raw=r["p_raw"], p_maxT=r["p_maxT"]))
        led.append(dict(context="immune", ledger_quantity=f"{key}_{strat}", layer="mean",
                        ledger_auroc=d[strat]["mean_over_layers"], p_mean_over_layers=d[strat]["p_mean_over_layers"],
                        N_perm=d[strat]["N"]))
led.append(dict(context="immune", ledger_quantity="screen_procedure", layer="mean", cv_foldmean=sp["mean_over_layers"],
                oof_auroc=oof_pooled, oof_strat_auroc=sp["oof_strat_mean_over_layers"],
                p_mean_over_layers=sp["perm_null"]["p_mean"], N_perm=sp["perm_null"]["N"],
                seed_median=sv["mean_over_layers_by_seed_median"], seed_q025=sv["lo"], seed_q975=sv["hi"],
                n_seeds=sv["n_seeds"]))
led.append(dict(context="immune", ledger_quantity="label_set", n_named=ledger["n_named"],
                n_named_observed=ledger["n_named_observed"], n_tf=ledger["n_tf"], n_target_only=ledger["n_target_only"],
                screen_reported_value=ledger["preprint_value"], screen_recomputed_mean=ledger["e9_recomputed_mean"]))
add(pd.DataFrame(led), "F_ledger_0744", "E3_node/ledger_0744.json (oof_auroc at layer=mean is the mean over layers)")
out = pd.concat(parts, ignore_index=True, sort=False)
# neutral names for the public CSV (values unchanged): the screen's bands are called "screen", as in the fragment
out = out.rename(columns={"legacy_band": "screen_band"})
for col_ in ("row_type", "convention"):
    out[col_] = out[col_].map(lambda x: x.replace("legacy_", "screen_") if isinstance(x, str) else x)
_txt = " ".join(out.astype(str).values.ravel()).lower() + " " + " ".join(out.columns).lower()
for bad in ("legacy", "submitted", "revision", "verifier", "round-"):
    assert bad not in _txt, bad
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_CSV, index=False)
print("wrote", OUT_TEX, "and", OUT_CSV, out.shape)
print("oof pooled layer mean", oof_pooled)
