"""Build Table S14 (crosswalk between the screen's spectral directions and the observed-gene basis; band provenance)
from revision3_experiments/E9_screen.

Inputs (read only):
  E9_screen/crosswalk.csv            per layer x convention x direction: best-matching observed direction and |cos|
                                     (right singular vectors), plus band principal angles
  E9_screen/crosswalk_band_scan.csv  every contiguous observed band vs each screen band, mean cos^2 of principal angles
  E9_screen/screen_band_mapping.json the fixed band mapping (geometry only; copy of legacy_band_mapping.json)
  E9_screen/band_provenance.md       provenance timeline of the three bands, with verbatim quotes from the iteration files
  E9_screen/verify/v2_spectral.json  independent re-derivation of the pooled k -> k-1 statistics
  E9_screen/results.json             pooled crosswalk numbers (checked against)
Outputs:
  supplement/S14.tex
  supplementary_tables/Table_S14_band_crosswalk.csv
Only counts, means and ranges over the source rows are computed here (no SVD is recomputed).
"""
from __future__ import annotations
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
E9 = ROOT / "revision3_experiments" / "E9_screen"
OUT_TEX = ROOT / "sagmb_revision" / "supplement" / "S14.tex"
OUT_CSV = ROOT / "sagmb_revision" / "supplementary_tables" / "Table_S14_band_crosswalk.csv"
CITE = "Kendiukhov (2026b)"

# ----------------------------------------------------------------------------- text -> LaTeX
_UNI = {"–": "--", "—": "---", "→": "->", "≥": ">=", "≤": "<=", "−": "-",
        "′": "'", "²": "^2"}


def tex(s) -> str:
    """Escape plain text (or the light markdown of band_provenance.md) for LaTeX and typeset recurring notation."""
    s = str(s)
    store: list[str] = []

    def keep(latex: str) -> str:
        store.append(latex)
        return f"\x00{len(store) - 1}\x00"

    s = re.sub(r"\*\*(.+?)\*\*", lambda m: keep(r"\textbf{" + tex(m.group(1)) + "}"), s)
    s = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", lambda m: keep(r"\emph{" + tex(m.group(1)) + "}"), s)
    s = re.sub(r"`([^`]+)`", lambda m: keep(r"\texttt{" + tex(m.group(1)) + "}"), s)
    s = s.replace("\\|", "|").replace("kidney_A", keep("kidney~A"))
    s = s.replace("3×10⁻⁵", keep(r"$3\times10^{-5}$")).replace("×", keep(r"$\times$"))
    s = s.replace("…", keep(r"\ldots{}"))
    for k, v in _UNI.items():
        s = s.replace(k, v)
    s = re.sub(r'"([^"]*)"', lambda m: keep("``") + m.group(1) + keep("''"), s)
    s = s.replace("->", keep(r"$\to$")).replace(">=", keep(r"$\geq$")).replace("<=", keep(r"$\leq$"))
    s = re.sub(r"\bSV\(k, k\+2\)", lambda m: keep(r"$\SV{k}$--$\SV{k+2}$"), s)
    s = re.sub(r"\bSV(\d+)\s*-{1,2}\s*SV(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$--$\\SV{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bSV(\d+)-{1,2}(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$--$\\SV{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bSV(\d+)\+", lambda m: keep(f"$\\SV{{{m.group(1)}}}$+"), s)
    s = re.sub(r"\bSV(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$"), s)
    s = re.sub(r"\bL(\d{1,2})\s*-\s*L(\d{1,2})\b", lambda m: keep(f"$L_{{{m.group(1)}}}$--$L_{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bL(\d{1,2})\b", lambda m: keep(f"$L_{{{m.group(1)}}}$"), s)
    s = re.sub(r"(?<=\d)-(?=\d)", "--", s)
    s = re.sub(r"(?<![\w)\]\-])-(?=\d)", lambda m: keep("$-$"), s)
    s = re.sub(r"\brho\b", lambda m: keep(r"$\rho$"), s)
    s = re.sub(r"(?<![\w.$])([pNnzr])(?=\s*(?:=|<|>)\s*[-\d.(\x00])", lambda m: keep(f"${m.group(1)}$"), s)
    s = s.replace("~", keep(r"$\sim$")).replace("<", keep("$<$")).replace(">", keep("$>$"))
    s = re.sub(r"(?<![\w.,;:!?)\]])'", "`", s)
    s = s.replace("\\", r"\textbackslash{}")
    for ch in "&%$#_{}":
        s = s.replace(ch, "\\" + ch)
    s = s.replace("^", r"\^{}").replace("|", r"\textbar{}")
    s = s.replace("/", r"/\allowbreak{}").replace(r"\_", r"\_\allowbreak{}")
    while "\x00" in s:
        s = re.sub("\x00(\\d+)\x00", lambda m: store[int(m.group(1))], s)
    return s


def band_tex(b: str) -> str:
    a, c = b[2:].split("-")
    return f"$\\SV{{{a}}}$--$\\SV{{{c}}}$"


f3 = lambda x: f"{x:.3f}"

# ----------------------------------------------------------------------------- data
CW = pd.read_csv(E9 / "crosswalk.csv")
BS = pd.read_csv(E9 / "crosswalk_band_scan.csv")
MAP = json.loads((E9 / "screen_band_mapping.json").read_text())
V2 = json.loads((E9 / "verify" / "v2_spectral.json").read_text())
RES = json.loads((E9 / "results.json").read_text())
MD = (E9 / "band_provenance.md").read_text(encoding="utf-8")

v = CW[CW.kind == "vector"].copy()
v["k"] = v.legacy_sv.str[2:].astype(int)
v["b"] = v.best_observed_sv.str[2:].astype(int)

# ----------------------------------------------------------------------------- Panel A: per-layer vector correspondence
summ = []
for (ctx, conv, layer), g in v.groupby(["context", "legacy_convention", "layer"]):
    g2 = g[g.k >= 2]
    s1 = g[g.k == 1].iloc[0]
    summ.append(dict(context=ctx, legacy_convention=conv, layer=int(layer),
                     n_dirs_k2_12=len(g2), n_best_is_k_minus_1=int((g2.b == g2.k - 1).sum()),
                     n_dirs_k1_12=len(g), n_best_is_same_index=int((g.b == g.k).sum()),
                     mean_abs_cos_best_k2_12=float(g2.abs_cos_best.mean()), min_abs_cos_best_k2_12=float(g2.abs_cos_best.min()),
                     max_abs_cos_best_k2_12=float(g2.abs_cos_best.max()),
                     mean_abs_cos_best_k1_12=float(g.abs_cos_best.mean()),
                     mean_abs_cos_same_index_k1_12=float(g.abs_cos_same_index.mean()),
                     min_abs_cos_same_index_k1_12=float(g.abs_cos_same_index.min()),
                     legacy_SV1_best_observed_sv=s1.best_observed_sv, legacy_SV1_abs_cos_best=float(s1.abs_cos_best)))
SM = pd.DataFrame(summ)
kid = SM[(SM.context == "kidney_A") & (SM.legacy_convention == "legacy_full_zero")].set_index("layer")
k209 = SM[(SM.context == "kidney_A") & (SM.legacy_convention == "screen_209_true")].set_index("layer")
imm = SM[(SM.context == "immune") & (SM.legacy_convention == "screen_nonzero")].set_index("layer")
assert len(kid) == len(k209) == len(imm) == 12

# pooled checks against results.json and the independent verification
vk = v[(v.context == "kidney_A") & (v.legacy_convention == "legacy_full_zero")]
vk2 = vk[vk.k >= 2]
pooled_frac = float((vk2.b == vk2.k - 1).mean())
CWR = RES["crosswalk_kidney_A_legacy_full_zero"]
assert abs(pooled_frac - CWR["frac_layer_sv_pairs_with_best_match_k_minus_1"]) < 1e-12
assert abs(pooled_frac - V2["crosswalk_frac_legacy_k_best_is_observed_k-1"]) < 1e-12
assert abs(vk2.abs_cos_best.mean() - CWR["mean_abs_cos_best_SV2_12"]) < 1e-12
assert abs(vk[vk.k == 1].abs_cos_best.mean() - CWR["mean_abs_cos_best_SV1"]) < 1e-12
EXACT_KM1 = V2["crosswalk_mean_abs_cos_legacy_k_vs_obs_k-1"]
vi = v[(v.context == "immune") & (v.legacy_convention == "screen_nonzero")]
v209 = v[(v.context == "kidney_A") & (v.legacy_convention == "screen_209_true")]

rowsA = []
for l in range(12):
    a, b, c = kid.loc[l], k209.loc[l], imm.loc[l]
    rowsA.append(f"$L_{{{l}}}$ & {a.n_best_is_k_minus_1}/{a.n_dirs_k2_12} & {f3(a.mean_abs_cos_best_k2_12)} & "
                 f"{f3(a.min_abs_cos_best_k2_12)}--{f3(a.max_abs_cos_best_k2_12)} & "
                 f"{f3(a.legacy_SV1_abs_cos_best)} (${a.legacy_SV1_best_observed_sv.replace('SV', chr(92) + 'SV{')}}}$) & "
                 f"{f3(b.mean_abs_cos_best_k1_12)} & {c.n_best_is_same_index}/{c.n_dirs_k1_12} & {f3(c.mean_abs_cos_same_index_k1_12)} \\\\")
rowsA.append(r"\midrule")
rowsA.append(f"all layers & {int((vk2.b == vk2.k - 1).sum())}/{len(vk2)} ({100 * pooled_frac:.1f}\\%) & {f3(vk2.abs_cos_best.mean())} & "
             f"{f3(vk2.abs_cos_best.min())}--{f3(vk2.abs_cos_best.max())} & {f3(vk[vk.k == 1].abs_cos_best.mean())} (mean) & "
             f"{f3(v209.abs_cos_best.mean())} & {int((vi.b == vi.k).sum())}/{len(vi)} & {f3(vi.abs_cos_same_index.mean())} \\\\")

panelA = [r"\begin{tabular}{@{}lrrrrrrr@{}}", r"\toprule",
          r" & \multicolumn{4}{c}{kidney~A: full matrix, zero rows included} & \multicolumn{1}{c}{kidney~A: 209-gene} & "
          r"\multicolumn{2}{c}{immune: screen's basis} \\",
          r"\cmidrule(lr){2-5}\cmidrule(lr){6-6}\cmidrule(lr){7-8}",
          r"Layer & best $= \SV{k-1}$ & mean $|\cos|$ & range $|\cos|$ & $\SV{1}$: best $|\cos|$ & mean best $|\cos|$ & "
          r"best $= \SV{k}$ & mean $|\cos|$ \\",
          r"\midrule"] + rowsA + [r"\bottomrule", r"\end{tabular}"]

# ----------------------------------------------------------------------------- Panel B: band mapping
CONV = {
    ("kidney_A", "legacy_full_zero"): "full matrix (2,831 zero rows)",
    ("kidney_A", "screen_209_true"): "209 named genes (22 zero rows)",
    ("kidney_A", "screen_nonzero"): "observed (1,972), mean",
    ("kidney_A", "screen_nonzero_sum"): "observed (1,972), summed",
    ("immune", "screen_nonzero"): "observed (2,039)",
}
USED = {  # display of where_used (legacy_band_mapping.json), worded neutrally
    "iters 13-19: SV2, SV3, SV4 PPI co-pole axes on the 209-gene subset": "iterations 13--19 (co-pole axes)",
    "iter_0049 onwards: SV2-SV4 TRRUST proximity / TF-target tests": "iterations 49+ (TRRUST, TF/target)",
    "iter_0050 onwards": "iterations 50+",
    "iter_0056 onwards (joint 6-D classifier)": "iterations 56+ (joint classifier)",
    "iters 0059-0060 (layer_gene_embedding_sum.npy, observed rows)": "iterations 59--60",
    "phase-2 iterations that load cycle1_main as 'kidney' (mean embeddings, nonzero rows)": "later iterations on kidney~A",
    "submitted manuscript (D2 convention)": "--- (reference basis)",
}
USED_CSV = {"submitted manuscript (D2 convention)": "reference basis: band indices in the full-matrix decomposition (zero rows included; kidney_A)"}
ORDER = [("kidney_A", "legacy_full_zero"), ("kidney_A", "screen_209_true"), ("kidney_A", "screen_nonzero"),
         ("kidney_A", "screen_nonzero_sum"), ("immune", "screen_nonzero")]
bands = sorted(MAP["bands"], key=lambda b: (ORDER.index((b["context"], b["legacy_convention"])), ["SV2-4", "SV5-7", "SV2-7"].index(b["legacy_band"])))
assert len(bands) == 9
rowsB, mapcsv = [], []
for b in bands:
    s = BS[(BS.context == b["context"]) & (BS.legacy_convention == b["legacy_convention"]) & (BS.legacy_band == b["legacy_band"]) &
           (BS.observed_band == b["mapped_observed_band"])]
    assert len(s) == 12 and abs(s.mean_cos2.mean() - b["mean_cos2_layer_avg"]) < 1e-9
    won = b["n_layers_won"].get(b["mapped_observed_band"], 0)
    oth = {k: n for k, n in b["n_layers_won"].items() if k != b["mapped_observed_band"]}
    others = (", ".join(f"{band_tex(k)} {n}" for k, n in oth.items()) if len(oth) <= 2 else f"{len(oth)} other bands")
    ctx = "kidney~A" if b["context"] == "kidney_A" else "immune"
    rowsB.append(f"{band_tex(b['legacy_band'])} & {ctx} & {CONV[(b['context'], b['legacy_convention'])]} & {USED[b['where_used']]} & "
                 f"{band_tex(b['mapped_observed_band'])} & {f3(b['mean_cos2_layer_avg'])} & {f3(s.mean_cos2.min())}--{f3(s.mean_cos2.max())} & "
                 f"{f3(b['identity_band_mean_cos2'])} & {band_tex(b['runner_up'])} ({f3(b['runner_up_mean_cos2'])}) & "
                 f"{won}/12" + (f" ({others})" if others else "") + r" \\")
    mapcsv.append(dict(panel="B", kind="band_mapping", context=b["context"], legacy_convention=b["legacy_convention"],
                       legacy_band=b["legacy_band"], where_used=USED_CSV.get(b["where_used"], b["where_used"]),
                       mapped_observed_band=b["mapped_observed_band"], mean_cos2_layer_avg=b["mean_cos2_layer_avg"],
                       mean_cos2_min_over_layers=float(s.mean_cos2.min()), mean_cos2_max_over_layers=float(s.mean_cos2.max()),
                       identity_band_mean_cos2=b["identity_band_mean_cos2"], runner_up=b["runner_up"],
                       runner_up_mean_cos2=b["runner_up_mean_cos2"], n_layers_won=json.dumps(b["n_layers_won"]),
                       per_layer_winner=json.dumps(b["per_layer_winner"])))
# the three kidney_A full-matrix values quoted in the text
fz = {b["legacy_band"]: b for b in bands if b["legacy_convention"] == "legacy_full_zero"}
assert fz["SV2-4"]["mapped_observed_band"] == "SV1-3" and fz["SV5-7"]["mapped_observed_band"] == "SV4-6" and fz["SV2-7"]["mapped_observed_band"] == "SV1-6"
panelB = [r"\begin{tabular}{@{}l l >{\raggedright\arraybackslash}p{2.2cm} >{\raggedright\arraybackslash}p{2.25cm} l r r r l >{\raggedright\arraybackslash}p{1.85cm}@{}}",
          r"\toprule",
          r"Band & Matrix & Decomposition (rows) & Used in & Mapped & $\overline{\cos^2}$ & range & identity & runner-up & layers won \\",
          r"\midrule"] + rowsB + [r"\bottomrule", r"\end{tabular}"]

# ----------------------------------------------------------------------------- Panel C: provenance timeline (parsed from band_provenance.md)
lines = [ln for ln in MD.splitlines() if re.match(r"^\| \d+ \|", ln)]
assert len(lines) == 20
BAND = {1: "SV2-4", 2: "SV2-4", 3: "SV2-4", 4: "SV2-4", 5: "SV2-4", 6: "SV2-4", 7: "SV2-4", 8: "all", 9: "SV2-4", 10: "SV2-4",
        11: "SV5-7", 12: "SV5-7", 13: "SV2-4", 14: "SV5-7", 15: "SV5-7", 16: "SV5-7", 17: "SV2-7", 18: "SV2-7", 19: "SV2-7", 20: "all"}
DROP = [" (per-layer z in `iter_0015/h01_sv45_ppi_copole.json`)",
        " (recomputed exactly here: claims_recompute.json, `residualised_immune_iter51_52`)"]
steps, rowsC = [], []
for ln in lines:
    cells = [c.strip() for c in re.split(r"(?<!\\)\|", ln.strip())[1:-1]]
    assert len(cells) == 5, ln
    step, it, mat, what, res = cells
    step = int(step)
    res_disp = res
    for d in DROP:
        res_disp = res_disp.replace(d, "")
    it_disp = re.sub(r"iter_0*(\d+)", r"\1", it)
    band = BAND[step]
    band_disp = "all" if band == "all" else band_tex(band)
    rowsC.append(f"{step} & {tex(it_disp)} & {band_disp} & {tex(mat) if mat not in ('—', 'same') else ('same as above' if mat == 'same' else '---')} & "
                 f"{tex(what)} & {tex(res_disp) if res_disp != '—' else '---'} \\\\")
    steps.append(dict(panel="C", kind="provenance_step", step=step, iteration=it, band=band, matrix_convention=mat,
                      what_happened=what, result_or_criterion=res))
phases = []
for ln in MD.splitlines():
    m = re.match(r"^\| (A|A′|B|B′) \|", ln)
    if m:
        cells = [c.strip() for c in re.split(r"(?<!\\)\|", ln.strip())[1:-1]]
        phases.append(dict(panel="C", kind="phase", step=cells[0], iteration=cells[1], matrix_convention=cells[2],
                           what_happened=f"gene rows: {cells[3]}", result_or_criterion=f"SVD convention: {cells[4]}"))
assert len(phases) == 4

# ----------------------------------------------------------------------------- LaTeX
W = [0.4, 1.35, 1.3, 2.35, 5.1, 5.4]
spec = "@{}" + "".join(f">{{\\scriptsize\\raggedright\\arraybackslash}}p{{{w}cm}}" for w in W) + "@{}"
hdrC = r"\# & Iteration & Band & Matrix, gene rows, convention & What was done or decided & Result or criterion (quotations verbatim) \\"
caption = (
    r"\textbf{Correspondence of other spectral bases with the observed-gene basis, and provenance of the screen's three bands.} "
    r"\emph{Panel A}, correspondence of single directions, from the right singular vectors (512-dimensional feature space; $|\cos|$ is sign-free) "
    r"at each layer. kidney~A, full matrix: the decomposition of all 4,803 rows, mean-centred over all of them, including the 2,831 "
    r"never-observed zero rows, whose first direction is therefore the observation indicator; it is compared with the decomposition of the 1,972 "
    r"observed genes (this paper's convention). For its $\SV{k}$, $k = 2, \ldots, 12$: best $= \SV{k-1}$, number of directions whose "
    r"best-matching observed direction (largest $|\cos|$ among observed $\SV{1}$--$\SV{40}$) is $\SV{k-1}$; mean and range of that best $|\cos|$ over "
    r"the 11 directions; $\SV{1}$: best $|\cos|$ of its first direction with any observed direction (and that direction), showing that it "
    r"has no counterpart. Pooled over the 132 layer--direction cases, the mean $|\cos|$ between full-matrix $\SV{k}$ and observed $\SV{k-1}$ "
    r"exactly is " + f3(EXACT_KM1) + r" (independent re-derivation); at layers where all 11 best matches are $\SV{k-1}$ it equals the mean shown. "
    r"kidney~A, 209-gene: mean best $|\cos|$ over $\SV{1}$--$\SV{12}$ for the decomposition fitted to the 209 named genes at their true rows "
    r"(22 of them zero rows), on which the screen assembled $\SV{2}$--$\SV{4}$. immune, screen's basis: the screen's decomposition of the "
    r"2,039 observed rows (``centered nonzero embeddings'') against this paper's; best $= \SV{k}$, number of $\SV{1}$--$\SV{12}$ whose best match "
    r"is the same index, with the mean $|\cos|$ at the same index. "
    r"\emph{Panel B}, mapping of bands: band, the band's indices in that decomposition; matrix and decomposition (rows), the matrix and "
    r"the rows the decomposition was fitted to (mean, mean embeddings; summed, the summed-embedding matrix the screen loaded in iterations 59--60); "
    r"used in, the iterations of the screen whose records use the band in that decomposition (---, a reference basis: the full kidney~A "
    r"decomposition, never-observed rows included); mapped, the contiguous observed-gene band of the same width (ending at or "
    r"before $\SV{12}$) with the largest mean $\cos^2$ of the principal angles between the two subspaces, averaged over the 12 layers "
    r"($\overline{\cos^2}$; range, its range over the 12 layers); identity, the same quantity for the band with the same indices; runner-up, the "
    r"second-best band; layers won, number of layers at which the mapped band is the best (other winning bands in parentheses). The mapping uses the "
    r"geometry of the matrices only; no outcome data (interaction, regulatory or label sets) were read to fix it. "
    r"\emph{Panel C}, provenance of the bands from the iteration files (JSON records, iteration reports, the language-model brainstormer's "
    r"briefs and next-step notes): \#, step; iteration, iteration and hypothesis or file; band, the band the step bears on (all, every "
    r"direction); text in quotation marks is verbatim from those files. ``sig'' is the screen's count of layers with empirical "
    r"permutation $p < 0.05$; emp\_p and perm\_p, empirical permutation $p$; mean\_z, mean $z$ over layers; rbc, rank-biserial correlation; "
    r"TRRUST, STRING, reference networks of documented regulation and protein interaction; PPI, protein--protein interaction; TF, transcription "
    r"factor. None of the three bands was specified before the results that motivated it were seen. "
    r"No sampling intervals apply: the matrices are fixed and all values are exact geometric quantities. "
    r"Full values: \texttt{Table\_S14\_band\_crosswalk.csv}."
)
L = [r"\setcounter{table}{13}", r"\begingroup", r"\footnotesize", r"\setlength{\tabcolsep}{2.4pt}",
     r"\setlength{\LTpre}{6pt}\setlength{\LTpost}{6pt}",
     f"\\begin{{longtable}}{{{spec}}}",
     f"\\caption{{{caption}}}\\label{{tab:S14}}\\\\",
     r"\multicolumn{6}{@{}l}{\footnotesize\textbf{A.} $\SV{k}$ of the full-matrix, 209-gene and screen's immune decompositions against the observed-gene $\SV{k}$ and $\SV{k-1}$, per layer}\\[2pt]",
     r"\multicolumn{6}{@{}l@{}}{\scriptsize\setlength{\tabcolsep}{3.2pt}" + "\n".join(panelA) + r"}\\[6pt]",
     r"\multicolumn{6}{@{}l}{\footnotesize\textbf{B.} The three bands in each decomposition, mapped to the observed-gene basis}\\[2pt]",
     r"\multicolumn{6}{@{}l@{}}{\scriptsize\setlength{\tabcolsep}{2.4pt}" + "\n".join(panelB) + r"}\\[6pt]",
     r"\multicolumn{6}{@{}l}{\footnotesize\textbf{C.} How the bands were chosen}\\",
     r"\toprule", hdrC, r"\midrule", r"\endfirsthead",
     r"\multicolumn{6}{@{}l}{\footnotesize\textit{Table~\thetable{} (continued)}}\\",
     r"\toprule", hdrC, r"\midrule", r"\endhead",
     r"\midrule", r"\multicolumn{6}{@{}r@{}}{\scriptsize\textit{continued on next page}}\\", r"\endfoot",
     r"\bottomrule", r"\endlastfoot"]
L += rowsC
L += [r"\end{longtable}", r"\endgroup"]
OUT_TEX.write_text("\n".join(L) + "\n")

# ----------------------------------------------------------------------------- CSV (tidy union of the three panels)
a_vec = v.drop(columns=["k", "b", "legacy_band", "observed_band", "principal_angles_deg", "mean_cos2", "max_angle_deg"]).copy()
a_vec.insert(0, "panel", "A")
a_sum = SM.copy()
a_sum.insert(0, "kind", "layer_summary")
a_sum.insert(0, "panel", "A")
pooled = []
for (ctx, conv), g in v.groupby(["context", "legacy_convention"]):
    g2 = g[g.k >= 2]
    pooled.append(dict(panel="A", kind="pooled_summary", context=ctx, legacy_convention=conv, layer="all",
                       n_dirs_k2_12=len(g2), n_best_is_k_minus_1=int((g2.b == g2.k - 1).sum()),
                       n_dirs_k1_12=len(g), n_best_is_same_index=int((g.b == g.k).sum()),
                       mean_abs_cos_best_k2_12=float(g2.abs_cos_best.mean()), min_abs_cos_best_k2_12=float(g2.abs_cos_best.min()),
                       max_abs_cos_best_k2_12=float(g2.abs_cos_best.max()), mean_abs_cos_best_k1_12=float(g.abs_cos_best.mean()),
                       mean_abs_cos_same_index_k1_12=float(g.abs_cos_same_index.mean()),
                       min_abs_cos_same_index_k1_12=float(g.abs_cos_same_index.min()),
                       legacy_SV1_abs_cos_best=float(g[g.k == 1].abs_cos_best.mean())))
pooled.append(dict(panel="A", kind="pooled_summary_exact_k_minus_1", context="kidney_A", legacy_convention="legacy_full_zero", layer="all",
                   n_dirs_k2_12=len(vk2), mean_abs_cos_exact_k_minus_1=EXACT_KM1,
                   source="E9_screen/verify/v2_spectral.json (independent re-derivation; observed SV1-SV13 searched)"))
b_ang = CW[CW.kind == "band_principal_angles"].drop(columns=["legacy_sv", "best_observed_sv", "abs_cos_best", "abs_cos_same_index",
                                                              "second_best_abs_cos", "norm_in_observed_SV1_12"]).copy()
b_ang.insert(0, "panel", "B")
b_scan = BS.copy()
b_scan.insert(0, "panel", "B")
c_df = pd.DataFrame(steps + phases)
out = pd.concat([a_sum, pd.DataFrame(pooled), a_vec, pd.DataFrame(mapcsv), b_ang, b_scan, c_df], ignore_index=True)
for col_ in out.columns:
    if out[col_].dtype == object:
        out[col_] = out[col_].map(lambda x: x.replace("called \"immune\" in the preprint", f"called \"immune\" in {CITE}")
                                  .replace("the preprint", CITE) if isinstance(x, str) else x)
lead = ["panel", "kind", "context", "legacy_convention", "layer"]
out = out[lead + [c for c in out.columns if c not in lead]]
# neutral names for the public CSV (values unchanged): the compared bases, bands and directions are called "basis"; the
# full-matrix decomposition (zero rows included) is a reference basis, not one attributed to the screen
out = out.rename(columns={c: c.replace("legacy_", "basis_") for c in out.columns})
for col_ in out.columns:
    if out[col_].dtype == object:
        out[col_] = out[col_].map(lambda x: x.replace("legacy_full_zero", "full_matrix_zero_rows") if isinstance(x, str) else x)
txt = " ".join(out.astype(str).values.ravel()).lower() + " " + " ".join(out.columns).lower()
for bad in ("preprint", "submitted", "revision", "review", "verifier", "round-3", "legacy"):
    assert bad not in txt, bad
out.to_csv(OUT_CSV, index=False)
print("S14:", len(rowsA), "panel A rows;", len(rowsB), "bands;", len(rowsC), "steps;", len(out), "CSV rows")
