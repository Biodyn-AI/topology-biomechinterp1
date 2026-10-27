"""Build Supplementary Table S11 (co-occurrence instrument, Word2Vec) from revision3_experiments/E8_word2vec.

Writes
  sagmb_revision/supplement/S11.tex                          (LaTeX fragment, one longtable)
  sagmb_revision/supplementary_tables/Table_S11_cooccurrence.csv

All numbers are read from the E8 outputs (w2v_results.csv, paired_vs_scgpt.csv, h1_assortativity_scgpt.csv,
saturation.csv, calibration.csv, cka.csv, models/*.json). H1b (STRING pair-cosine AUROC) intervals are the
delete-d vertex-subsampling intervals (columns *_deleted; CONVENTION3_ERRATA item 3). Quantities computed here:
ranges, best layers, counts of layers whose interval excludes 0, and maxima over SVs/seeds.
Run:  /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python build_S11.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
E8 = ROOT / "revision3_experiments" / "E8_word2vec"
OUT_TEX = HERE / "S11.tex"
OUT_CSV = HERE.parent / "supplementary_tables" / "Table_S11_cooccurrence.csv"

MAIN = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
NAME = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune",
        "kidney_A": "kidney~A"}


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


def rng(a, b, d=3):
    if a >= 0 and b >= 0:
        return f"{num(a, d)}--{num(b, d)}"
    return f"{num(a, d)} to {num(b, d)}"


def pval(p):
    if p >= 0.001:
        return f"{p:.4f}" if p < 0.001 else (f"{p:.3f}" if p < 0.1 else f"{p:.2f}")
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def sci(x):
    m, e = f"{x:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def cell2(top, bottom):
    return f"{top}\\newline {bottom}"


# ----------------------------------------------------------------------------------------------- load
w = pd.read_csv(E8 / "w2v_results.csv")
pv = pd.read_csv(E8 / "paired_vs_scgpt.csv")
h1s = pd.read_csv(E8 / "h1_assortativity_scgpt.csv")
sat = pd.read_csv(E8 / "saturation.csv").set_index(["context", "model"])
cal = pd.read_csv(E8 / "calibration.csv").set_index(["context", "model"])
cka = pd.read_csv(E8 / "cka.csv")
models = {}
for f in sorted((E8 / "models").glob("*.json")):
    if f.name.startswith("._"):
        continue
    models[f.stem] = json.loads(f.read_text())

NC = 7
rows = []


def panel(title):
    if rows:
        rows.append("\\midrule")
    rows.append(f"\\multicolumn{{{NC}}}{{@{{}}p{{16.4cm}}@{{}}}}{{\\footnotesize\\textbf{{{title}}}}}\\\\")


def lab(n, text):
    width = {1: "3.1cm", 2: "5.41cm", 3: "7.72cm"}[n]
    if n == 1:
        return text
    return f"\\multicolumn{{{n}}}{{@{{}}>{{\\scriptsize\\raggedright\\arraybackslash}}p{{{width}}}}}{{{text}}}"


def sub(text):
    rows.append(f"\\multicolumn{{{NC}}}{{@{{}}l}}{{\\scriptsize\\textit{{{text}}}}}\\\\")


def line(cells):
    rows.append(" & ".join(cells) + "\\\\")


def mean_row(c, corpus, task):
    r = w[(w.context == c) & (w.corpus == corpus) & (w.seed.astype(str) == "mean") & (w.task == task)]
    assert len(r) == 1, (c, corpus, task, len(r))
    return r.iloc[0]


def seed_row(c, corpus, window, epochs, seed, task):
    r = w[(w.context == c) & (w.corpus == corpus) & (w.window == window) & (w.epochs == epochs)
          & (w.seed.astype(str) == str(seed)) & (w.task == task)]
    assert len(r) == 1, (c, corpus, window, epochs, task, len(r))
    return r.iloc[0]


def paired(c, model, task):
    x = pv[(pv.context == c) & (pv.w2v == model) & (pv.task == task)].copy()
    x["layer"] = x.layer.astype(int)
    assert len(x) == 12, (c, model, task, len(x))
    return x.set_index("layer").sort_index()


# ======================================================================================= A: per extraction
panel("A. Word2Vec trained on the full population corpus (five-token window, 10 epochs, mean of 5 seeds) against "
      "scGPT's residual stream, per extraction")
line([lab(2, "Task and quantity")] + [NAME[c] for c in MAIN])
rows.append("\\cmidrule(lr){1-7}")
summ = {}

sub("H1b: STRING pair-cosine AUROC against rewired non-pairs (delete-$d$ vertex-subsampling intervals)")
cells = [lab(2, "Word2Vec [95\\% CI]")]
for c in MAIN:
    r = mean_row(c, "full_w5_e10", "H1b_pair_cos_auroc")
    cells.append(cell2(num(r.value), ci(r.ci_lo_deleted, r.ci_hi_deleted)))
line(cells)
cells = [lab(2, "scGPT, 12 layers (best layer)")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H1b_pair_cos_auroc")
    cells.append(f"{rng(p.scgpt.min(), p.scgpt.max())} ($L_{{{int(p.scgpt.idxmax())}}}$)")
line(cells)
cells = [lab(2, "Word2Vec minus scGPT at its best layer")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H1b_pair_cos_auroc")
    b = p.loc[p.scgpt.idxmax()]
    cells.append(cell2(num(b.delta, sign=True), ci(b.delta_ci_lo_deleted, b.delta_ci_hi_deleted)))
line(cells)
cells = [lab(2, "Word2Vec minus scGPT at $L_{11}$")]
for c in MAIN:
    b = paired(c, "full_w5_e10", "H1b_pair_cos_auroc").loc[11]
    cells.append(cell2(num(b.delta, sign=True), ci(b.delta_ci_lo_deleted, b.delta_ci_hi_deleted)))
line(cells)
cells = [lab(2, "Layers with the difference's CI $<0$ (scGPT higher) / $>0$")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H1b_pair_cos_auroc")
    cells.append(f"{int((p.delta_ci_hi_deleted < 0).sum())} / {int((p.delta_ci_lo_deleted > 0).sum())}")
line(cells)

sub("H2: detection-stratified TF AUROC (LOO cosine; gene-bootstrap intervals)")
cells = [lab(2, "Word2Vec [95\\% CI]")]
for c in MAIN:
    r = mean_row(c, "full_w5_e10", "H2_strat_tf_auroc")
    cells.append(cell2(num(r.value), ci(r.ci_lo, r.ci_hi)))
line(cells)
cells = [lab(2, "Word2Vec permutation $p$ (null mean)")]
for c in MAIN:
    r = mean_row(c, "full_w5_e10", "H2_strat_tf_auroc")
    cells.append(cell2(pval(r.p_raw), f"({num(r.null_mean)})"))
line(cells)
cells = [lab(2, "scGPT, 12 layers (best layer)")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H2_strat_tf_auroc")
    cells.append(f"{rng(p.scgpt.min(), p.scgpt.max())} ($L_{{{int(p.scgpt.idxmax())}}}$)")
line(cells)
cells = [lab(2, "Word2Vec minus scGPT at its best layer, raw")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H2_strat_tf_auroc")
    b = p.loc[p.scgpt.idxmax()]
    cells.append(cell2(num(b.delta, sign=True), ci(b.delta_ci_lo, b.delta_ci_hi)))
line(cells)
cells = [lab(2, "Same, null-centred (each instrument's permutation-null mean subtracted)")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H2_strat_tf_auroc")
    b = p.loc[p.scgpt.idxmax()]
    cells.append(cell2(num(b.delta_nullcentred, sign=True), ci(b.delta_nullcentred_ci_lo, b.delta_nullcentred_ci_hi)))
line(cells)
cells = [lab(2, "Layers with the null-centred CI $<0$ (scGPT higher) / $>0$")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H2_strat_tf_auroc")
    cells.append(f"{int((p.delta_nullcentred_ci_hi < 0).sum())} / {int((p.delta_nullcentred_ci_lo > 0).sum())}")
line(cells)

sub("H3: within-TF stratified edge AUROC (two-way TF $\\times$ target bootstrap intervals)")
cells = [lab(2, "Word2Vec [95\\% CI]")]
for c in MAIN:
    r = mean_row(c, "full_w5_e10", "H3_withinTF_auroc")
    cells.append(cell2(num(r.value), ci(r.ci_lo, r.ci_hi)))
line(cells)
cells = [lab(2, "scGPT, 12 layers (best layer)")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H3_withinTF_auroc")
    cells.append(f"{rng(p.scgpt.min(), p.scgpt.max())} ($L_{{{int(p.scgpt.idxmax())}}}$)")
line(cells)
cells = [lab(2, "Word2Vec minus scGPT at its best layer")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H3_withinTF_auroc")
    b = p.loc[p.scgpt.idxmax()]
    cells.append(cell2(num(b.delta, sign=True), ci(b.delta_ci_lo, b.delta_ci_hi)))
line(cells)
cells = [lab(2, "Layers with the difference's CI $<0$ / $>0$")]
for c in MAIN:
    p = paired(c, "full_w5_e10", "H3_withinTF_auroc")
    cells.append(f"{int((p.delta_ci_hi < 0).sum())} / {int((p.delta_ci_lo > 0).sum())}")
line(cells)

sub("H1a: assortativity of $\\SV{1}$--$\\SV{10}$ on the STRING network (descriptive; directions not index-matched)")
cells = [lab(2, "Word2Vec: directions with max-$T$ $p<0.05$ in all 5 seeds (largest mean $r$)")]
for c in MAIN:
    x = w[(w.context == c) & (w.corpus == "full_w5_e10") & (w.seed.astype(str) == "mean")
          & (w.task == "H1a_assortativity")]
    assert len(x) == 10
    cells.append(f"{int((x.frac_seeds_maxT_lt_005 == 1).sum())}/10 ({num(x.value.max())})")
line(cells)
cells = [lab(2, "scGPT: direction--layer tests with max-$T$ $p<0.05$ (largest $r$)")]
for c in MAIN:
    x = h1s[h1s.context == c]
    assert len(x) == 120
    cells.append(f"{int((x.p_maxT_120 < 0.05).sum())}/120 ({num(x.r.max())})")
line(cells)
cells = [lab(2, "Largest seed SD of Word2Vec (H1b, H2, H3)")]
for c in MAIN:
    sds = [mean_row(c, "full_w5_e10", t).sd_seeds for t in ["H1b_pair_cos_auroc", "H2_strat_tf_auroc",
                                                              "H3_withinTF_auroc"]]
    cells.append(num(max(sds)))
line(cells)

# ======================================================================================= B: saturation
panel("B. Saturation: full-corpus Word2Vec, seed 0, after 10 and 40 epochs")
line([lab(3, "Quantity"), "lung, 10 epochs", "lung, 40 epochs", "immune, 10 epochs", "immune, 40 epochs"])
rows.append("\\cmidrule(lr){1-7}")
SAT = [("lung", "full_w5_e10_s0"), ("lung", "full_w5_e40_s0"), ("immune", "full_w5_e10_s0"), ("immune", "full_w5_e40_s0")]
line([lab(3, "Pair updates")] + [sci(sat.loc[k].pair_updates_total) for k in SAT])
line([lab(3, "H1b [95\\% CI, delete-$d$]")] +
     [cell2(num(sat.loc[k].H1b_pair_cos_auroc), ci(sat.loc[k].H1b_pair_cos_auroc_ci_lo_deleted,
                                                   sat.loc[k].H1b_pair_cos_auroc_ci_hi_deleted)) for k in SAT])
line([lab(3, "H2 [95\\% CI]")] +
     [cell2(num(sat.loc[k].H2_strat_tf_auroc), ci(sat.loc[k].H2_strat_tf_auroc_ci_lo, sat.loc[k].H2_strat_tf_auroc_ci_hi))
      for k in SAT])
line([lab(3, "H2 permutation $p$ (null mean)")] +
     [f"{pval(sat.loc[k].H2_p)} ({num(sat.loc[k].H2_null_mean)})" for k in SAT])
line([lab(3, "H3 [95\\% CI]")] +
     [cell2(num(sat.loc[k].H3_withinTF_auroc), ci(sat.loc[k].H3_withinTF_auroc_ci_lo, sat.loc[k].H3_withinTF_auroc_ci_hi))
      for k in SAT])
line([lab(3, "H1a: largest $r$ (directions with max-$T$ $p<0.05$)")] +
     [f"{num(sat.loc[k].H1a_max_r)} ({int(sat.loc[k].H1a_n_sv_maxT_lt_005)}/10)" for k in SAT])
ck = {}
for c in ["lung", "immune"]:
    x = cka[(cka.context == c) & (cka.kind == "saturation_e40_vs_e10")]
    assert len(x) == 1
    ck[c] = x.iloc[0].cka
line([lab(3, "Linear CKA between the 10- and 40-epoch models"),
      f"\\multicolumn{{2}}{{>{{\\scriptsize}}c}}{{{num(ck['lung'])}}}",
      f"\\multicolumn{{2}}{{>{{\\scriptsize}}c@{{}}}}{{{num(ck['immune'])}}}"])
cells40 = {}
for c in ["lung", "immune"]:
    p = paired(c, "full_w5_e40", "H2_strat_tf_auroc")
    b = p.loc[p.scgpt.idxmax()]
    l11 = p.loc[11]
    ncov = int(((p.delta_nullcentred_ci_lo <= 0) & (p.delta_nullcentred_ci_hi >= 0)).sum())
    cells40[c] = (cell2(f"{num(b.delta_nullcentred, sign=True)} ($L_{{{int(p.scgpt.idxmax())}}}$)",
                        ci(b.delta_nullcentred_ci_lo, b.delta_nullcentred_ci_hi)),
                  cell2(num(l11.delta_nullcentred, sign=True), ci(l11.delta_nullcentred_ci_lo, l11.delta_nullcentred_ci_hi)),
                  f"{ncov}/12")
line([lab(3, "H2, Word2Vec minus scGPT at its best layer, null-centred"), "$^{a}$", cells40["lung"][0], "$^{a}$",
      cells40["immune"][0]])
line([lab(3, "H2, Word2Vec minus scGPT at $L_{11}$, null-centred"), "$^{a}$", cells40["lung"][1], "$^{a}$",
      cells40["immune"][1]])
line([lab(3, "H2, layers whose null-centred CI includes 0"), "$^{a}$", cells40["lung"][2], "$^{a}$", cells40["immune"][2]])

# ======================================================================================= C: window calibration
panel("C. Window calibration: Word2Vec trained on the extraction cells only (256 cells), seed 0")
line(["Quantity", "lung: window 5, 10 epochs", "lung: window 5, matched updates", "lung: window 512, 10 epochs",
      "immune: window 5, 10 epochs", "immune: window 5, matched updates", "immune: window 512, 10 epochs"])
rows.append("\\cmidrule(lr){1-7}")
CAL = [("lung", 5.0, 10.0, "extr_w5_e10_s0", "lung__extr__w5__e10__s0"),
       ("lung", 5.0, 571.0, "extr_w5_e571_s0", "lung__extr__w5__e571__s0"),
       ("lung", 512.0, 10.0, "extr_w512_e10_s0", "lung__extr__w512__e10__s0"),
       ("immune", 5.0, 10.0, "extr_w5_e10_s0", "immune__extr__w5__e10__s0"),
       ("immune", 5.0, 365.0, "extr_w5_e365_s0", "immune__extr__w5__e365__s0"),
       ("immune", 512.0, 10.0, "extr_w512_e10_s0", "immune__extr__w512__e10__s0")]
line(["Epochs; pair updates"] + [f"{int(e)}; {sci(models[m]['pair_updates_total'])}" for c, wd, e, k, m in CAL])
cells = ["H1b [95\\% CI, delete-$d$]"]
for c, wd, e, k, m in CAL:
    r = seed_row(c, "extr", wd, e, 0, "H1b_pair_cos_auroc")
    assert abs(r.value - cal.loc[(c, k)].H1b_pair_cos_auroc) < 1e-12
    cells.append(cell2(num(r.value), ci(r.ci_lo_deleted, r.ci_hi_deleted)))
line(cells)
cells = ["H2 (permutation $p$)"]
for c, wd, e, k, m in CAL:
    r = cal.loc[(c, k)]
    cells.append(cell2(num(r.H2_strat_tf_auroc), f"({pval(r.H2_p)})"))
line(cells)
cells = ["H3 [95\\% CI]"]
for c, wd, e, k, m in CAL:
    r = cal.loc[(c, k)]
    cells.append(cell2(num(r.H3_withinTF_auroc), ci(r.H3_withinTF_auroc_ci_lo, r.H3_withinTF_auroc_ci_hi)))
line(cells)
cells = ["H1a: largest $r$ (directions with max-$T$ $p<0.05$)"]
for c, wd, e, k, m in CAL:
    r = cal.loc[(c, k)]
    cells.append(f"{num(r.H1a_max_r)} ({int(r.H1a_n_sv_maxT_lt_005)}/10)")
line(cells)
cells = ["Linear CKA with the window-512 model"]
for c, wd, e, k, m in CAL:
    if wd == 512.0:
        cells.append("--")
        continue
    x = cka[(cka.context == c) & (cka.kind == "calibration") & (cka.a == f"extr_w512_e10_s0") & (cka.b == k)]
    assert len(x) == 1, (c, k)
    cells.append(num(x.iloc[0].cka))
line(cells)

# ----------------------------------------------------------------------------------------------- caption
caption = (
    "\\textbf{A co-occurrence-only instrument: Word2Vec trained on the token sets scGPT receives, compared with scGPT's "
    "residual stream on H1--H3.} "
    "Word2Vec: skip-gram, 256 dimensions, 10 negative samples, trained on every cell of each extraction's population "
    "(lung 15{,}000, lung (ext.) 9{,}000, kidney~B 8{,}000, immune 20{,}000, kidney~A 11{,}376 cells), each tokenised as "
    "scGPT tokenises it (non-zero genes, the 512 most expressed), with token and cell order reshuffled on every pass. "
    "Every representation is mean-centred over the observed genes $V$ and scored on identical genes, pairs, labels and "
    "resampling draws; scGPT is the whole-human checkpoint. "
    "H1b, AUROC of the cosine of STRING physical pairs (combined score $\\geq 700$) against one degree-preserving rewired "
    "non-pair set, with vertex-subsampling intervals ($m = n/2$, $B = 2{,}000$) under delete-$d$ scaling; "
    "H2, detection-stratified AUROC of the LOO-cosine TF score (Lambert TFs against all other observed genes), "
    "within-stratum permutation null ($N = 2{,}000$; floor $5.0\\times10^{-4}$) and gene-bootstrap intervals "
    "($B = 2{,}000$, scores held fixed); H3, within-TF stratified AUROC of the TF--target cosine with two-way TF $\\times$ "
    "target bootstrap intervals ($B = 2{,}000$); H1a, assortativity of $\\SV{1}$--$\\SV{10}$ of the observed-gene SVD on "
    "the STRING network with a Maslov--Sneppen null ($N = 1{,}000$; max-$T$ over 10 directions for Word2Vec and over the "
    "120 direction--layer tests for scGPT). These $p$-values are not those of the primary H1/H2 tests. "
    "Differences are paired (Word2Vec minus scGPT on the same draws); the best layer is the scGPT layer with the highest "
    "value. Under the within-stratum null the LOO score is not centred at 0.5 (null means 0.27--0.50 for Word2Vec), so H2 "
    "is judged against its own null, and the null-centred difference subtracts each instrument's null mean (a descriptive "
    "adjustment). Seed variation is reported separately and is not propagated into the intervals. "
    "\\textbf{B}: one seed; $^{a}$the 10-epoch paired differences are in Panel A (mean of 5 seeds). "
    "\\textbf{C}: the five-token window used for cost against a window spanning the whole cell (512), and the five-token "
    "window with its number of pair updates matched to the 512-token window (571 epochs for lung, 365 for immune). "
    "Pair updates, the number of (centre, context) updates in training. "
    "Per-seed values, every layer, the extraction-cells-only models for all extractions and all CKA values are in the CSV. "
    "For H1b, the vertex-subsampling intervals with the pre-registered scaling factor $\\sqrt{m/n}$, which omits the "
    "finite-population term and so gives half-widths smaller by a factor $\\sqrt{2}$, are also in the CSV (columns without "
    "the suffix \\texttt{\\_deleted}, beside the reported \\texttt{\\_deleted} columns). "
    "Full values: \\texttt{Table\\_S11\\_cooccurrence.csv}."
)

colspec = ("@{}>{\\scriptsize\\raggedright\\arraybackslash}p{3.1cm}"
           "*{6}{>{\\scriptsize\\raggedright\\arraybackslash}p{2.1cm}}@{}")
tex = ["\\setcounter{table}{10}",
       "\\begingroup", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}",
       "\\renewcommand{\\arraystretch}{1.05}", f"\\begin{{longtable}}{{{colspec}}}",
       f"\\caption{{{caption}}}\\label{{tab:S11}}\\\\", "\\toprule", "\\endfirsthead",
       f"\\multicolumn{{{NC}}}{{@{{}}l}}{{\\footnotesize Table~\\thetable\\ (continued)}}\\\\", "\\toprule", "\\endhead",
       "\\midrule", f"\\multicolumn{{{NC}}}{{r@{{}}}}{{\\footnotesize continued on next page}}\\\\", "\\endfoot",
       "\\bottomrule", "\\endlastfoot"]
tex += rows
tex += ["\\end{longtable}", "\\endgroup"]
OUT_TEX.write_text("\n".join(tex) + "\n")

# ----------------------------------------------------------------------------------------------- CSV
parts = []


def add(df, panel_id, source):
    df = df.copy()
    df.insert(0, "source_file", source)
    df.insert(0, "panel", panel_id)
    parts.append(df)


add(w, "A_B_C_word2vec_results", "E8_word2vec/w2v_results.csv (H1b: ci_lo_deleted/ci_hi_deleted = reported intervals; ci_lo/ci_hi = pre-registered sqrt(m/n) scaling)")
add(pv, "A_B_C_paired_vs_scgpt", "E8_word2vec/paired_vs_scgpt.csv (H1b: *_deleted = reported intervals)")
add(h1s, "A_scgpt_H1a_assortativity", "E8_word2vec/h1_assortativity_scgpt.csv")
add(sat.reset_index(), "B_saturation", "E8_word2vec/saturation.csv")
cal2 = cal.reset_index()
lo_d, hi_d = [], []
for _, r in cal2.iterrows():
    wd = 512.0 if "w512" in r.model else 5.0
    ep = float(r.model.split("_e")[1].split("_")[0])
    s = seed_row(r.context, "extr", wd, ep, 0, "H1b_pair_cos_auroc")
    lo_d.append(s.ci_lo_deleted)
    hi_d.append(s.ci_hi_deleted)
cal2.insert(cal2.columns.get_loc("H1b_pair_cos_auroc_ci_hi") + 1, "H1b_pair_cos_auroc_ci_lo_deleted", lo_d)
cal2.insert(cal2.columns.get_loc("H1b_pair_cos_auroc_ci_lo_deleted") + 1, "H1b_pair_cos_auroc_ci_hi_deleted", hi_d)
add(cal2, "C_calibration", "E8_word2vec/calibration.csv (+ delete-d H1b interval from w2v_results.csv, same model)")
add(cka, "B_C_cka", "E8_word2vec/cka.csv")
md = pd.DataFrame([dict(model=k, context=v.get("context"), corpus=v.get("corpus"), seed=v.get("seed"),
                        window=v["params"]["window"], epochs=v["params"]["epochs"], n_cells=v.get("n_cells"),
                        n_tokens=v.get("n_tokens"), pair_updates_total=v.get("pair_updates_total"),
                        seconds=v.get("seconds")) for k, v in models.items()])
add(md, "B_C_model_metadata", "E8_word2vec/models/*.json")
out = pd.concat(parts, ignore_index=True, sort=False)
_txt = " ".join(out.astype(str).values.ravel()).lower() + " " + " ".join(out.columns).lower()
for bad in ("legacy", "submitted", "revision", "verifier", "round-", "manuscript"):
    assert bad not in _txt, bad
out.to_csv(OUT_CSV, index=False)
print("wrote", OUT_TEX, OUT_CSV, out.shape)
