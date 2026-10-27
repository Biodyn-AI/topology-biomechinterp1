"""Build Supplementary Table S8 (rank reduction, H4) from revision3_experiments/E7_rank.

Writes
  sagmb_revision/supplement/S8.tex                     (LaTeX fragment, one longtable)
  sagmb_revision/supplementary_tables/Table_S8_rank.csv

All numbers are read from the E7 outputs. Interaction (PPI) intervals are the delete-d vertex-subsampling intervals
(columns ci_lo_ms / ci_hi_ms; CONVENTION3_ERRATA item 3). The only quantities computed here are counts of layers whose
interval excludes 0 and the range of the reported values.
Run:  /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python build_S8.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
E7 = ROOT / "revision3_experiments" / "E7_rank"
OUT_TEX = HERE / "S8.tex"
OUT_CSV = HERE.parent / "supplementary_tables" / "Table_S8_rank.csv"

EVAL = ["lung", "lung_external", "kidney_B", "immune"]
FIVE = EVAL + ["kidney_A"]
SIX = FIVE + ["immune_seed43"]
NAME = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune",
        "kidney_A": "kidney~A", "immune_seed43": "immune (seed 43)"}
TASK = {"PPI": "Interaction", "NODE": "TF identity", "EDGE": "Regulatory edges"}
STATIC = [("scgpt_static", "scGPT input embedding"), ("gf_v2_316m", "Geneformer V2-316M"),
          ("gf_v1_10m", "Geneformer V1-10M"), ("gf_v2_104m", "Geneformer V2-104M"),
          ("aido_100m", "AIDO.Cell-100M"), ("scfoundation", "scFoundation"), ("tgpt", "tGPT"),
          ("esm2_uce", "ESM2-15B (UCE per-protein)"), ("esm2_650m", "ESM2-650M"),
          ("seq_features", "Sequence descriptors"), ("hpa_coexp", "HPA cell-type expression"),
          ("scgpt_lung_contextual", "scGPT lung residual stream, layer mean$^{c}$")]


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
        return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def cell2(top, bottom):
    return f"{top}\\newline {bottom}"


def intc(n):
    return f"{int(n):,}".replace(",", "{,}")


def band_tex(b):
    b = str(b).replace("SV", "")
    i, j = b.split("-")
    return f"$\\SV{{{i}}}$--$\\SV{{{j}}}$"


# ----------------------------------------------------------------------------------------------- load
kstar = json.loads((E7 / "kstar.json").read_text())
delta = pd.read_csv(E7 / "delta.csv")
curves = {c: pd.read_csv(E7 / f"curves_{c}.csv") for c in SIX}
rtj = pd.read_csv(E7 / "remove_top_j.csv")
sens = pd.read_csv(E7 / "sens_R_N_ci.csv")
bvt = pd.read_csv(E7 / "bands_vs_topk.csv")
anyk = pd.read_csv(E7 / "any_k_test.csv")
post = pd.read_csv(E7 / "posthoc_best_k.csv")
ssum = pd.read_csv(E7 / "static_summary.csv")
scur = pd.read_csv(E7 / "static_curves.csv")

NC = 7
rows = []


def panel(title):
    if rows:
        rows.append("\\midrule")
    rows.append(f"\\multicolumn{{{NC}}}{{@{{}}p{{16.4cm}}@{{}}}}{{\\footnotesize\\textbf{{{title}}}}}\\\\")


def header(cells):
    rows.append(" & ".join(cells) + "\\\\")
    rows.append("\\cmidrule(lr){1-7}")


def line(cells):
    assert len(cells) == NC, cells
    rows.append(" & ".join(cells) + "\\\\")


def curve_mean(c, task, rep):
    x = curves[c]
    r = x[(x.task == task) & (x.rep == rep) & (x.layer.astype(str) == "mean")]
    assert len(r) == 1, (c, task, rep, len(r))
    return r.iloc[0]


# ======================================================================================= A: selection of k*
panel("A. Truncation $k^{*}$ chosen in the selection draw (immune, seed 43), before any evaluation extraction")
header(["Task", "$k^{*}$", "$m(k^{*})$", "$m(512)$", "Selection gap", "$m(3) - m(512)$", "$m(10) - m(512)$"])
for t in ["PPI", "NODE", "EDGE"]:
    k = kstar["tasks"][t]
    line([TASK[t], str(k["k_star"]), num(k["layer_mean_at_kstar"], 4), num(k["layer_mean_full"], 4),
          num(k["selection_gap"], 4, sign=True), num(curve_mean("immune_seed43", t, "T3").diff_vs_full, sign=True),
          num(curve_mean("immune_seed43", t, "T10").diff_vs_full, sign=True)])

# ======================================================================================= B: Delta per extraction
panel("B. Pre-registered estimand $\\Delta = m(k^{*}) - m(512)$, layer mean, in the four evaluation extractions")
header(["Task, extraction", "$m(k^{*})$", "$m(512)$", "$\\Delta$", "\\multicolumn{2}{>{\\scriptsize}l}{95\\% CI}",
        "Layers with CI $<0$ / $>0$"])
for t in ["PPI", "NODE", "EDGE"]:
    kk = kstar["tasks"][t]["k_star"]
    rows.append(f"\\multicolumn{{{NC}}}{{@{{}}l}}{{\\scriptsize\\textit{{{TASK[t]}, $k^{{*}} = {kk}$}}}}\\\\")
    for c in EVAL:
        d = delta[(delta.ctx == c) & (delta.task == t)]
        m = d[d.layer.astype(str) == "mean"].iloc[0]
        lay = d[d.layer.astype(str) != "mean"]
        assert len(lay) == 12
        below = int((lay.ci_hi_ms < 0).sum())
        above = int((lay.ci_lo_ms > 0).sum())
        assert int(m.k_star) == kk
        lab = f"\\quad {NAME[c]}"
        if t == "NODE":
            rows.append(" & ".join([lab, num(m.metric_kstar, 4), num(m.metric_full, 4), "0 (identical)",
                                    "\\multicolumn{2}{>{\\scriptsize}l}{--}", "--"]) + "\\\\")
        else:
            rows.append(" & ".join([lab, num(m.metric_kstar, 4), num(m.metric_full, 4), num(m.delta, 4, sign=True),
                                    f"\\multicolumn{{2}}{{>{{\\scriptsize}}l}}{{{ci(m.ci_lo_ms, m.ci_hi_ms, 4)}}}", f"{below} / {above}"])
                        + "\\\\")

# ======================================================================================= C: k = 3 and k = 10
panel("C. The leading $k$ directions only: layer-mean value at $k = 3$ and $k = 10$ (difference from $m(512)$ in parentheses)")
header(["Task, $k$"] + [NAME[c] for c in SIX])
for t in ["PPI", "NODE", "EDGE"]:
    line([f"{TASK[t]}, $k = 512$ (full)"] + [num(curve_mean(c, t, "T512").value) for c in SIX])
    for k in [3, 10]:
        cells = [f"{TASK[t]}, $k = {k}$"]
        for c in SIX:
            r = curve_mean(c, t, f"T{k}")
            cells.append(f"{num(r.value)} ({num(r.diff_vs_full, sign=True)})")
        line(cells)

# ======================================================================================= D: remove top j
panel("D. Removing the $j$ leading directions (keeping $\\SV{j+1}$--$\\SV{512}$): layer-mean difference from $m(512)$ [95\\% CI]")
header(["Task, directions removed"] + [NAME[c] for c in SIX])
for t in ["PPI", "NODE", "EDGE"]:
    for j in [1, 2, 3, 5, 10]:
        cells = [f"{TASK[t]}, $j = {j}$"]
        for c in FIVE:
            r = sens[(sens.ctx == c) & (sens.task == t) & (sens.rep == f"R{j}") & (sens.layer.astype(str) == "mean")]
            assert len(r) == 1
            r = r.iloc[0]
            cells.append(cell2(num(r["diff"], sign=True), ci(r.ci_lo_ms, r.ci_hi_ms)))
        r = rtj[(rtj.ctx == "immune_seed43") & (rtj.task == t) & (rtj.rep == f"R{j}")].iloc[0]
        cells.append(num(r.diff_vs_full, sign=True))
        line(cells)

# ======================================================================================= E: screen bands vs top-k
panel("E. The screen's bands against the leading directions of the same width: layer-mean difference [95\\% CI]")
header(["Task, band vs partner"] + [NAME[c] for c in SIX])
BANDS = [("Nnom_SV2-4", "T3", f"{band_tex('SV2-4')} vs top 3"), ("Nnom_SV5-7", "T3", f"{band_tex('SV5-7')} vs top 3"),
         ("Nnom_SV2-7", "T6", f"{band_tex('SV2-7')} vs top 6"),
         ("Nmap_SV5-7", "T3", f"{band_tex('SV4-6')}$^{{b}}$ vs top 3")]
for t in ["PPI", "NODE", "EDGE"]:
    for rep, partner, lab in BANDS:
        cells = [f"{TASK[t]}, {lab}"]
        for c in FIVE:
            r = sens[(sens.ctx == c) & (sens.task == t) & (sens.rep == rep) & (sens.layer.astype(str) == "mean")]
            assert len(r) == 1 and r.iloc[0].reference == partner
            r = r.iloc[0]
            cells.append(cell2(num(r["diff"], sign=True), ci(r.ci_lo_ms, r.ci_hi_ms)))
        r = bvt[(bvt.ctx == "immune_seed43") & (bvt.task == t) & (bvt.rep == rep)].iloc[0]
        assert r.partner_rep == partner
        cells.append(num(r.diff_vs_partner, sign=True))
        line(cells)

# ======================================================================================= F: static representations
panel("F. Eleven static gene representations on the common lung gene set (5{,}226 genes): full embedding and best truncation")
header(["Representation ($d$)", "Interaction: full", "Interaction: best $k$ (gain)", "TF identity: full",
        "TF identity: best $k$ (gain)", "Edges: full", "Edges: best $k$ (gain)"])
# shorter row labels in the PDF table (the caption defines them; the CSV keeps the long labels)
SHORT = {"esm2_uce": "ESM2-15B", "hpa_coexp": "HPA cell types",
         "scgpt_lung_contextual": "scGPT lung, layer mean$^{c}$"}
for key, lab in STATIC:
    if key == "scgpt_lung_contextual":
        rows.append("\\cmidrule(lr){1-7}")
    x = ssum[ssum.model == key].set_index("task")
    cells = [f"{SHORT.get(key, lab)} ({intc(x.loc['PPI'].dim_full)})"]
    for t in ["PPI", "NODE", "EDGE"]:
        r = x.loc[t]
        cells += [num(r.full_value), f"{int(r.best_k)} ({num(r.best_minus_full, sign=True)})"]
    line(cells)
s11 = ssum[ssum.model != "scgpt_lung_contextual"]
gain_max = s11.best_minus_full.max()
assert gain_max < 0.0225

# ----------------------------------------------------------------------------------------------- caption
caption = (
    "\\textbf{Rank reduction (H4): truncation to the leading spectral directions against the full representation.} "
    "$m$, layer-mean task metric (mean over the 12 layers) of the representation restricted to the first $k$ directions "
    "$U_k S_k$ of the observed-gene SVD at each layer ($m(512)$, the full representation). Tasks: interaction, AUROC of "
    "the cosine of STRING physical pairs (combined score $\\geq 700$) against one degree-preserving rewired set of non-pairs; "
    "TF identity, the H2 statistic (detection-stratified LOO-cosine AUROC, Lambert TFs against all other observed genes); "
    "regulatory edges, the H3 statistic (within-TF stratified AUROC of the TF--target cosine, TRRUST pairs against "
    "matched candidates). "
    "\\textbf{A}: $k^{*}$ maximises $m(k)$ over $k \\in \\{1, 2, 3, 5, 7, 10, 16, 32, 64, 128, 256, 512\\}$ in the selection draw "
    "(ties to the smallest $k$); selection gap $= m(k^{*}) - m(512)$ there. "
    "\\textbf{B}: $\\Delta$ in the evaluation extractions with paired 95\\% intervals ($B = 2{,}000$; the same resamples for "
    "$k^{*}$ and 512): interaction, vertex subsampling ($m = n/2$) with delete-$d$ scaling; TF identity, gene bootstrap; edges, "
    "two-way TF $\\times$ target bootstrap. For TF identity the selection draw chose $k^{*} = 512$, so $\\Delta = 0$ by "
    "construction. Last column: number of the 12 per-layer intervals entirely below / above 0. "
    "\\textbf{C}--\\textbf{E}: descriptive, not pre-registered; kidney~A and immune (seed 43, the selection draw) are shown "
    "for completeness; the last column has no interval. "
    "\\textbf{D}: removing the leading directions instead of keeping them; intervals as in B (interaction delete-$d$). "
    "\\textbf{E}: screen bands $\\SV{2}$--$\\SV{4}$, $\\SV{5}$--$\\SV{7}$, $\\SV{2}$--$\\SV{7}$ on the observed-gene basis "
    "against the first 3 or 6 directions (the screen used $\\SV{5}$--$\\SV{7}$ of the observed genes in its later kidney~A "
    "iterations; Table~S14); $^{b}$$\\SV{5}$--$\\SV{7}$ shifted down by one index, as when it is read in the reference "
    "decomposition of the full kidney~A matrix, zero rows included (there $\\SV{2}$--$\\SV{4}$ and $\\SV{2}$--$\\SV{7}$ "
    "become the first 3 and 6 directions, a difference of 0 by construction), shown in every extraction. "
    "\\textbf{F}: one vector per gene: input gene-token embeddings (scGPT, Geneformer, tGPT); gene outputs of one forward "
    "pass on a neutral input (AIDO.Cell-100M, scFoundation); ESM2-15B, the 5{,}120-dimensional per-protein embedding (UCE "
    "file), and residue-mean-pooled ESM2-650M; 32 sequence descriptors; and expression across 154 Human Protein Atlas cell "
    "types (128 principal components). Each is mean-centred on the 5{,}226 lung genes present in all eleven vocabularies and scored on identical "
    "genes, pairs and labels; $d$, dimension; best $k$ over the same grid, gain $=$ best minus full (no interval). "
    "$^{c}$Reference row: the scGPT lung residual stream on the same genes, layer mean. "
    f"No truncation exceeds the full embedding by more than {num(gain_max)} for any static representation and task. "
    "Per-layer values, all $k$, the descriptive any-$k$ permutation test and, for interaction, the intervals with the "
    "pre-registered scaling factor $\\sqrt{m/n}$, which omits the finite-population term and so gives half-widths smaller by "
    "a factor $\\sqrt{2}$ (columns \\texttt{ci\\_lo}/\\texttt{ci\\_hi}, beside the reported \\texttt{ci\\_lo\\_ms}/"
    "\\texttt{ci\\_hi\\_ms}), are in the CSV. "
    "Full values: \\texttt{Table\\_S8\\_rank.csv}."
)

colspec = ("@{}>{\\scriptsize\\raggedright\\arraybackslash}p{3.1cm}"
           "*{6}{>{\\scriptsize\\raggedright\\arraybackslash}p{2.1cm}}@{}")
tex = ["\\setcounter{table}{7}",
       "\\begingroup", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}",
       "\\renewcommand{\\arraystretch}{1.05}", f"\\begin{{longtable}}{{{colspec}}}",
       f"\\caption{{{caption}}}\\label{{tab:S8}}\\\\", "\\toprule", "\\endfirsthead",
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
    if "context" not in df.columns and "ctx" in df.columns:
        df.insert(0, "context", df["ctx"])
    elif "context" not in df.columns and "model" in df.columns:
        df.insert(0, "context", "lung (common 5,226-gene set)")
    df.insert(0, "source_file", source)
    df.insert(0, "panel", panel_id)
    parts.append(df)


ks = []
for t, k in kstar["tasks"].items():
    for kk, v in k["layer_mean_curve"].items():
        ks.append(dict(ctx=kstar["selection_context"], task=t, k=int(kk), value=v, k_star=k["k_star"],
                       layer_mean_at_kstar=k["layer_mean_at_kstar"], layer_mean_full=k["layer_mean_full"],
                       selection_gap=k["selection_gap"], note="selection record written before evaluation; its k = 1 entries use a 1e-12 norm offset that "
                            "orders exact +-1 cosines (exact k = 1 values in the C_D_E_curves rows); k = 1 is not the argmax "
                            "either way"))
add(pd.DataFrame(ks), "A_kstar_selection", "E7_rank/kstar.json")
add(delta, "B_delta", "E7_rank/delta.csv (ci_lo_ms/ci_hi_ms = reported intervals, PPI delete-d; PPI ci_lo/ci_hi = pre-registered sqrt(m/n) scaling)")
add(anyk, "B_any_k_descriptive", "E7_rank/any_k_test.csv")
for c in SIX:
    add(curves[c], "C_D_E_curves", f"E7_rank/curves_{c}.csv")
add(post, "C_posthoc_best_k", "E7_rank/posthoc_best_k.csv")
add(rtj, "D_remove_top_j", "E7_rank/remove_top_j.csv")
add(sens, "D_E_paired_CIs", "E7_rank/sens_R_N_ci.csv (ci_lo_ms/ci_hi_ms = reported intervals, PPI delete-d; PPI ci_lo/ci_hi = pre-registered sqrt(m/n) scaling)")
add(bvt, "E_bands_vs_topk", "E7_rank/bands_vs_topk.csv")
ssum2 = ssum.copy()
ssum2.insert(1, "label", ssum2.model.map(dict(STATIC)).str.replace("$^{c}$", "", regex=False))
add(ssum2, "F_static_summary", "E7_rank/static_summary.csv")
add(scur, "F_static_curves", "E7_rank/static_curves.csv")
out = pd.concat(parts, ignore_index=True, sort=False)
# neutral names for the public CSV (values unchanged): the screen's band mapping (Table S14)
for col_ in out.columns:
    if out[col_].dtype == object:
        out[col_] = out[col_].map(lambda x: x.replace("legacy_band_mapping.json", "screen_band_mapping.json")
                                  .replace("legacy_full_zero", "screen_full_zero") if isinstance(x, str) else x)
_txt = " ".join(out.astype(str).values.ravel()).lower() + " " + " ".join(out.columns).lower()
for bad in ("legacy", "submitted", "revision", "verifier", "round-", "manuscript"):
    assert bad not in _txt, bad
out.to_csv(OUT_CSV, index=False)
print("wrote", OUT_TEX, OUT_CSV, out.shape, "max static gain", gain_max)
