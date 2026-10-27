"""Build Supplementary Table S1 (spectral concentration) from revision3_experiments/E1_compression.

Reads only CSV/JSON outputs written by the E1 scripts; the only quantities computed here are trivial derived ones
(L0/L11 folds, counts of layer steps with an increase, counts of near-degenerate pairs, ranges over extractions).
Writes  sagmb_revision/supplement/S1.tex                     (LaTeX fragment, one longtable)
        sagmb_revision/supplementary_tables/Table_S1_compression.csv   (long-format values)
Run:    /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python sagmb_revision/supplement/build_S1.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REV = HERE.parent
ROOT = REV.parent
E1 = ROOT / "revision3_experiments" / "E1_compression"
OUT_TEX = HERE / "S1.tex"
OUT_CSV = REV / "supplementary_tables" / "Table_S1_compression.csv"

CORE = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
LABEL = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune",
         "kidney_A": "kidney~A"}
# 16 extractions: (context, group label, design label) in display order
ALL16 = [
    ("lung", "lung", "primary"),
    ("lung_external", "lung (ext.)", "replication"),
    ("kidney_B", "kidney~B", "second tissue"),
    ("immune", "immune", "continuity"),
    ("kidney_A", "kidney~A", "continuity"),
    ("immune_seed43", "immune", "cell draw, seed 43"),
    ("immune_seed44", "immune", "cell draw, seed 44"),
    ("immune_TSP14", "immune", "donor TSP14"),
    ("immune_TSP21", "immune", "donor TSP21"),
    ("immune_TSP25", "immune", "donor TSP25"),
    ("immune_Bcell", "immune", "B cells"),
    ("immune_CD4T", "immune", "CD4 T cells"),
    ("immune_CD8T", "immune", "CD8 T cells"),
    ("kidney_ckpt_whole_human", "kidney (192 cells)", "whole-human ckpt"),
    ("kidney_ckpt_kidney", "kidney (192 cells)", "kidney ckpt"),
    ("kidney_ckpt_brain", "kidney (192 cells)", "brain ckpt"),
]
DESIGN = {c: d for c, _, d in ALL16}
PLAIN = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney B", "immune": "immune",
         "kidney_A": "kidney A", "immune_seed43": "immune, cell draw seed 43", "immune_seed44": "immune, cell draw seed 44",
         "immune_TSP14": "immune, donor TSP14", "immune_TSP21": "immune, donor TSP21", "immune_TSP25": "immune, donor TSP25",
         "immune_Bcell": "immune, B cells", "immune_CD4T": "immune, CD4 T cells", "immune_CD8T": "immune, CD8 T cells",
         "kidney_ckpt_whole_human": "kidney (192 cells), whole-human checkpoint",
         "kidney_ckpt_kidney": "kidney (192 cells), kidney checkpoint",
         "kidney_ckpt_brain": "kidney (192 cells), brain checkpoint"}
SHOW_LAYERS = [0, 3, 6, 9, 11]

# ------------------------------------------------------------------------------------------------ inputs
rk = pd.read_csv(E1 / "rank_definitions.csv", float_precision="round_trip")
per = pd.read_csv(E1 / "per_context_layer.csv", float_precision="round_trip")
dt = pd.read_csv(E1 / "depth_trend.csv", float_precision="round_trip").set_index("context")
gaps = pd.read_csv(E1 / "gaps.csv", float_precision="round_trip")
gb1000 = pd.read_csv(E1 / "verify" / "gap_B1000.csv", float_precision="round_trip")
sh = pd.read_csv(E1 / "shuffle.csv", float_precision="round_trip")
mt = pd.read_csv(E1 / "matched.csv", float_precision="round_trip")
thr = pd.read_csv(E1 / "detection_threshold.csv", float_precision="round_trip")
nd = json.loads((E1 / "null_direction.json").read_text())
zr = pd.read_csv(E1 / "zero_row_artefact.csv", float_precision="round_trip")
dsum = pd.read_csv(E1 / "detection_summary.csv", float_precision="round_trip")
dcor = pd.read_csv(E1 / "detection_corr.csv", float_precision="round_trip")
assert gb1000.flag_B1000.astype(bool).all() and len(gb1000) == 9

spec = rk.merge(per[["context", "layer", "sv1_var", "sv1_var_boot_lo", "sv1_var_boot_hi", "sv2_var", "sv3_var",
                     "cum10_var", "n_sv_below_1e-3_rel"]], on=["context", "layer"], validate="1:1")
assert len(spec) == 192
# consistency of the two sources (same SVD)
chk = rk.merge(per[["context", "layer", "er", "pr"]], on=["context", "layer"], suffixes=("", "_per"))
assert (chk.er_sv - chk.er).abs().max() < 1e-9 and (chk.pr - chk.pr_per).abs().max() < 1e-9


def get(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).to_numpy()
    r = df[m]
    assert len(r) == 1, (kw, len(r))
    return r.iloc[0]


# ------------------------------------------------------------------------------------------------ derived
trend = {}
for c, _, _ in ALL16:
    s = spec[spec.context == c].sort_values("layer")
    d = {}
    for m in ("er_sv", "er_var", "pr"):
        v = s[m].to_numpy()
        d[f"{m}_L0"], d[f"{m}_L11"] = v[0], v[-1]
        d[f"{m}_fold_L0_over_L11"] = v[0] / v[-1]
        d[f"{m}_n_increases"] = int((np.diff(v) > 0).sum())
    d["sv1_var_L0"], d["sv1_var_L11"] = s.sv1_var.iloc[0], s.sv1_var.iloc[-1]
    d["sv1_var_fold_L11_over_L0"] = d["sv1_var_L11"] / d["sv1_var_L0"]
    d["n_V"] = int(s.n_V.iloc[0])
    # cross-check against depth_trend.csv (written by e1_compression.py)
    t = dt.loc[c]
    assert abs(d["er_sv_fold_L0_over_L11"] - t.er_fold_L0_over_L11) < 1e-9
    assert abs(d["pr_fold_L0_over_L11"] - t.pr_fold_L0_over_L11) < 1e-6
    assert d["er_sv_n_increases"] == int(t.er_n_increases) and d["pr_n_increases"] == int(t.pr_n_increases)
    d["spearman_er_layer"], d["spearman_pr_layer"] = t.spearman_er_layer, t.spearman_pr_layer
    d["spearman_sv1var_layer"] = t.spearman_sv1var_layer
    trend[c] = d

gaps["nd"] = gaps.near_degenerate.astype(bool)
assert ((gaps.basic_lo <= 1.0) == gaps.nd).all()
gk = gaps.pivot_table(index="k", columns="context", values="nd", aggfunc="sum")
gk_all = gaps.groupby("k").nd.agg(["sum", "count"])
g_tot = gaps.groupby("context").nd.agg(["sum", "count"])

# ------------------------------------------------------------------------------------------------ CSV (long)
rows = []


def add(section, context, statistic, value, layer=None, k=None, scheme=None, min_nobs=None, ci_lo=None, ci_hi=None,
        sd=None, n=None, interval="", source="", note=""):
    rows.append(dict(section=section, context=context, extraction=PLAIN.get(context, context),
                     design=DESIGN.get(context, ""), layer=layer, k=k, scheme=scheme, min_nobs=min_nobs,
                     statistic=statistic, value=value, ci_lo=ci_lo, ci_hi=ci_hi, sd=sd, n=n, interval=interval,
                     source=source, note=note))


SRC_SPEC = "E1_compression/rank_definitions.csv; E1_compression/per_context_layer.csv"
for r in spec.itertuples():
    add("spectrum_by_layer", r.context, "er_sv (effective rank, Roy-Vetterli, p=S/sum S)", r.er_sv, layer=r.layer,
        n=r.n_V, source=SRC_SPEC)
    add("spectrum_by_layer", r.context, "er_var (exp of the Shannon entropy of the variance shares, q=S^2/sum S^2)", r.er_var, layer=r.layer, n=r.n_V,
        source=SRC_SPEC)
    add("spectrum_by_layer", r.context, "pr (participation ratio)", r.pr, layer=r.layer, n=r.n_V, source=SRC_SPEC)
    add("spectrum_by_layer", r.context, "sv1_var (share of variance on SV1)", r.sv1_var, layer=r.layer, n=r.n_V,
        ci_lo=r.sv1_var_boot_lo, ci_hi=r.sv1_var_boot_hi,
        interval="95% gene-bootstrap percentile interval, B=200", source=SRC_SPEC)
    add("spectrum_by_layer", r.context, "sv2_var", r.sv2_var, layer=r.layer, n=r.n_V, source=SRC_SPEC)
    add("spectrum_by_layer", r.context, "sv3_var", r.sv3_var, layer=r.layer, n=r.n_V, source=SRC_SPEC)
    add("spectrum_by_layer", r.context, "cum10_var (cumulative share SV1-10)", r.cum10_var, layer=r.layer,
        n=r.n_V, source=SRC_SPEC)
    add("spectrum_by_layer", r.context, "n_sv_below_1e-3_rel (singular values < 1e-3 S1)",
        int(spec.loc[r.Index, "n_sv_below_1e-3_rel"]), layer=r.layer,
        n=r.n_V, source=SRC_SPEC)
for c, d in trend.items():
    for key, v in d.items():
        if key == "n_V":
            continue
        add("depth_trend_L0_to_L11", c, key, v, n=d["n_V"],
            source="derived from E1_compression/rank_definitions.csv; spearman_* from depth_trend.csv",
            note="n_increases = number of the 11 consecutive layer steps with an increase")
for r in gaps.itertuples():
    add("gap_ratio", r.context, "gap_ratio S_k/S_k+1", r.ratio, layer=r.layer, k=r.k, ci_lo=r.basic_lo,
        ci_hi=r.basic_hi, sd=r.boot_sd, n=None, interval="95% basic gene-bootstrap interval, B=200",
        source="E1_compression/gaps.csv", note=f"near_degenerate={bool(r.nd)} (rule: basic_lo <= 1)")
for c in [x for x, _, _ in ALL16]:
    for k in range(1, 12):
        add("near_degenerate_count_by_k", c, "n_near_degenerate (of 12 layers)", int(gk.loc[k, c]), k=k, n=12,
            source="derived from E1_compression/gaps.csv")
    add("near_degenerate_count_by_k", c, "n_near_degenerate (all k=1..11, of 132)", int(g_tot.loc[c, "sum"]),
        n=int(g_tot.loc[c, "count"]), source="derived from E1_compression/gaps.csv")
for k in range(1, 12):
    add("near_degenerate_count_by_k", "all_16", "n_near_degenerate (of 192 extraction-layers)",
        int(gk_all.loc[k, "sum"]), k=k, n=int(gk_all.loc[k, "count"]), source="derived from E1_compression/gaps.csv")
add("near_degenerate_count_by_k", "all_16", "n_near_degenerate (all pairs)", int(gaps.nd.sum()), n=len(gaps),
    source="derived from E1_compression/gaps.csv")
for r in gb1000.itertuples():
    add("gap_ratio_B1000_check", r.context, "gap_ratio S_k/S_k+1 (B=1000)", r.ratio, layer=r.layer, k=r.k,
        ci_lo=r.basic_lo_B1000, ci_hi=r.basic_hi_B1000, sd=r.boot_sd,
        interval="95% basic gene-bootstrap interval, B=1000", source="E1_compression/verify/gap_B1000.csv",
        note=f"flag_B1000={bool(r.flag_B1000)}")
for r in sh.itertuples():
    s = "E1_compression/shuffle.csv"
    add("feature_shuffle", r.context, "er_obs", r.er_obs, layer=r.layer, n=r.n_V, source=s)
    add("feature_shuffle", r.context, "er_shuf_mean", r.er_shuf_mean, layer=r.layer, sd=r.er_shuf_sd, n=r.n_draws,
        source=s, note=f"min {r.er_shuf_min:.4f}, max {r.er_shuf_max:.4f} over draws")
    add("feature_shuffle", r.context, "er_shuf_over_obs", r.er_shuf_over_obs, layer=r.layer, n=r.n_draws, source=s)
    add("feature_shuffle", r.context, "pr_shuf_mean", r.pr_shuf_mean, layer=r.layer, sd=r.pr_shuf_sd, n=r.n_draws,
        source=s)
    add("feature_shuffle", r.context, "sv1_var_shuf_mean", r.sv1_var_shuf_mean, layer=r.layer, sd=r.sv1_var_shuf_sd,
        n=r.n_draws, source=s)
for r in mt.itertuples():
    s = "E1_compression/matched.csv"
    note = (f"n_sub={r.n_sub}, {r.n_draws} draws without replacement; median n_obs of subsample "
            f"{r.median_nobs_sub}; scheme 'pooled' = common deciles of log n_obs pooled over the five extractions")
    for m in ("er", "pr", "sv1_var"):
        add("detection_matched_subsample", r.context, f"{m}_mean", getattr(r, f"{m}_mean"), layer=r.layer,
            scheme=r.scheme, ci_lo=getattr(r, f"{m}_q025"), ci_hi=getattr(r, f"{m}_q975"), sd=getattr(r, f"{m}_sd"),
            n=r.n_sub, interval="2.5-97.5% quantiles over 200 subsamples", source=s, note=note)
for r in thr.itertuples():
    s = "E1_compression/detection_threshold.csv"
    for m in ("n_genes", "er", "pr", "sv1_var", "cum10_var"):
        add("nobs_threshold", r.context, m, getattr(r, m), layer=r.layer, min_nobs=r.min_nobs, n=r.n_genes, source=s)
for c, dd in nd.items():
    for L, v in dd.items():
        for m, x in v.items():
            add("null_direction", c, m, x, layer=int(L), source="E1_compression/null_direction.json")
for r in zr.itertuples():
    s = "E1_compression/zero_row_artefact.csv"
    for m in ("er", "pr", "sv1_var", "cum10_var", "pearson_sv1_observed_indicator", "spearman_sv1_observed_indicator",
              "er_feature_shuffle_mean"):
        x = getattr(r, m)
        if pd.notna(x):
            add("zero_row_matrix", r.context, m, x, layer=r.layer, n=r.n_rows, source=s,
                note=f"matrix={r.matrix}; zero rows={r.n_zero_rows}")
for r in dsum.itertuples():
    s_ = "E1_compression/detection_summary.csv"
    add("detection_tracking", r.context, "sv1_spearman (SV1 coordinate vs log n_obs)", r.sv1_spearman, layer=r.layer,
        source=s_)
    add("detection_tracking", r.context, "top_abs_spearman (max |rho| over SV1-10)", r.top_abs_spearman,
        layer=r.layer, k=r.top_k, source=s_, note="k = index of the most detection-correlated of SV1-10")
# median |rho| of the other nine of SV1-10 (all but the most detection-correlated), per extraction and layer
for (c, Lr), g in dcor.groupby(["context", "layer"], sort=False):
    assert sorted(g.k) == list(range(1, 11))
    g = g.sort_values("abs_spearman", ascending=False)
    top = get(dsum, context=c, layer=Lr)
    assert int(g.k.iloc[0]) == int(top.top_k) and abs(g.abs_spearman.iloc[0] - top.top_abs_spearman) < 1e-9
    add("detection_tracking", c, "median_abs_spearman_other_nine (median |rho| over the other nine of SV1-10)",
        float(g.abs_spearman.iloc[1:].median()), layer=int(Lr), n=9, source="E1_compression/detection_corr.csv",
        note="the nine directions of SV1-10 other than the most detection-correlated one")
csv = pd.DataFrame(rows)
for col in ("layer", "k", "min_nobs", "n"):
    v = pd.to_numeric(csv[col], errors="coerce")
    assert ((v.dropna() % 1) == 0).all()
    csv[col] = v.astype("Int64")
csv.loc[csv.context == "all_16", "extraction"] = "all 16 extractions"
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
csv.to_csv(OUT_CSV, index=False)

# ------------------------------------------------------------------------------------------------ LaTeX
def f(x, d=1):
    return f"{x:,.{d}f}".replace(",", "{,}")


def rng_(vals, d=1):
    return f"{f(min(vals), d)}--{f(max(vals), d)}"


NC = 8
L = []
A = L.append
cap = (
    r"\caption[Spectral concentration with depth]{\textbf{Spectral concentration with depth, in every extraction.} "
    r"Singular values $S_1 \geq S_2 \geq \dots$ of the centred observed-gene matrix ($|V| \times 512$) at each layer. "
    r"ER, effective rank $\exp(-\sum_k p_k \log p_k)$ with $p_k = S_k/\sum_j S_j$ (Roy and Vetterli, 2007); "
    r"$\mathrm{ER}_{\mathrm{var}}$, exponential of the Shannon entropy of the variance shares $q_k = S_k^2/\sum_j S_j^2$; "
    r"PR, participation ratio "
    r"$(\sum_k S_k^2)^2/\sum_k S_k^4$; $\SV{1}$, share of variance $q_1$ on the first direction (\%); "
    r"$\SV{1\text{--}10}$, cumulative share of the first ten (\%). "
    r"\textbf{A}, the five main extractions at five layers ($|V|$ = 5{,}648, 5{,}616, 4{,}780, 2{,}039 and 1{,}972 genes "
    r"for lung, lung (ext.), kidney~B, immune and kidney~A). Fold, $L_0$ value over $L_{11}$ value for ER, $\mathrm{ER}_{\mathrm{var}}$ and "
    r"PR, $L_{11}$ over $L_0$ for the variance shares. Intervals for $\SV{1}$: 95\% gene-bootstrap percentile intervals "
    r"($B = 200$; genes resampled with replacement). No bootstrap interval is given for ER: resampling genes with replacement "
    r"duplicates rows and biases it downwards, so its sampling variability is taken from the without-replacement subsamples "
    r"of panel~D. "
    r"\textbf{B}, all sixteen extractions. Cell draw, an independent draw of 256 immune cells; donor, 256 cells of one donor; "
    r"B/CD4~T/CD8~T cells, 256 cells of one type; ckpt, one draw of 192 kidney cells (the 8{,}143-gene axis of kidney~B) "
    r"passed through the whole-human checkpoint used for every other extraction, the kidney-fine-tuned and the brain scGPT "
    r"checkpoint. Rises, number of the eleven consecutive layer steps at which ER / $\mathrm{ER}_{\mathrm{var}}$ / PR increased "
    r"(0, strictly decreasing). "
    r"\textbf{C}, identifiability and controls. Near-degenerate pairs: number of layers (of 12) at which the 95\% basic "
    r"bootstrap interval of the gap ratio $S_k/S_{k+1}$ ($B = 200$ gene resamples) reaches 1; last column, count and "
    r"percentage over all 192 extraction-layers of the sixteen extractions (all 2{,}112 pairs in the last of these rows). "
    r"The nine such pairs among $\SV{1}$--$\SV{4}$ ($k \leq 3$) of the five main extractions are flagged again at "
    r"$B = 1{,}000$. Feature shuffle: ER after each of the 512 features is permuted independently across genes (mean of "
    r"20 draws; SD $\leq 0.12$), divided by the observed ER, and the shuffled ER itself at $L_{11}$; last column, range over "
    r"the sixteen extractions. Null direction: singular values below $10^{-3} S_1$ per layer, and the absolute cosine between "
    r"the last right singular vector and the uniform vector $\mathbf{1}/\sqrt{512}$ at $L_0$ / $L_{11}$, with the "
    r"relative size of the last singular value. "
    r"\textbf{D}, detection controls. Matched: 1{,}300 genes subsampled without replacement within common deciles of "
    r"$\log n_{\mathrm{obs}}$ pooled over the five extractions, mean $\pm$ SD over 200 draws; $n_{\mathrm{obs}} \geq 3$: "
    r"the spectrum at $L_{11}$ after excluding genes contributed by fewer than three sampled cells (thresholds 5 and 10 "
    r"in the CSV). Last column, range over the five extractions. "
    r"Full values: \texttt{Table\_S1\_compression.csv}.}\label{tab:S1}\\"
)


def panel(title):
    A(r"\multicolumn{%d}{@{}l}{\scriptsize\textbf{%s}}\\" % (NC, title))


def mc(n, s, align="c"):
    return r"\multicolumn{%d}{%s}{\scriptsize %s}" % (n, align, s)


A(r"\setcounter{table}{0}")
A(r"\begingroup")
A(r"\footnotesize")
A(r"\setlength{\tabcolsep}{2.6pt}")
A(r"\renewcommand{\arraystretch}{0.97}")
A(r"\setlength{\LTleft}{0pt}\setlength{\LTright}{0pt}")
A(r"\begin{longtable}{@{\extracolsep{\fill}}>{\scriptsize}l>{\scriptsize}l*{6}{>{\scriptsize}r}@{}}")
A(cap)
A(r"\toprule")
A(r"\endfirsthead")
A(r"\multicolumn{%d}{@{}l}{\scriptsize\itshape Table~\thetable\ (continued)}\\" % NC)
A(r"\toprule")
A(r"\endhead")
A(r"\midrule")
A(r"\multicolumn{%d}{r@{}}{\scriptsize\itshape continued on next page}\\" % NC)
A(r"\endfoot")
A(r"\bottomrule")
A(r"\endlastfoot")
# ---- Panel A
panel(r"A\quad Spectral summaries at five layers, five main extractions")
A(r"\midrule")
A(r"Extraction & Measure & $L_0$ & $L_3$ & $L_6$ & $L_9$ & $L_{11}$ & Fold \\")
A(r"\midrule")
MEAS = [("ER", "er_sv", 1), (r"$\mathrm{ER}_{\mathrm{var}}$", "er_var", 1), ("PR", "pr", 1), (r"$\SV{1}$ (\%)", "sv1_var", 1),
        (r"$\SV{1\text{--}10}$ (\%)", "cum10_var", 1)]
for ci, c in enumerate(CORE):
    s = spec[spec.context == c].set_index("layer")
    for mi, (lab, col, d) in enumerate(MEAS):
        cells, ci_cells = [], []
        for Lr in SHOW_LAYERS:
            r = s.loc[Lr]
            if col == "sv1_var":
                cells.append(f(100 * r.sv1_var))
                ci_cells.append(f"{f(100 * r.sv1_var_boot_lo)}--{f(100 * r.sv1_var_boot_hi)}")
            elif col == "cum10_var":
                cells.append(f(100 * r.cum10_var))
            else:
                cells.append(f(r[col], d))
        if col in ("sv1_var", "cum10_var"):
            fold = s.loc[11, col] / s.loc[0, col]
        else:
            fold = s.loc[0, col] / s.loc[11, col]
        name = LABEL[c] if mi == 0 else ""
        end = r"\\*" if mi < len(MEAS) - 1 else r"\\"   # no page break inside one extraction's block
        A(f"{name} & {lab} & " + " & ".join(cells) + f" & {f(fold, 2)}$\\times$ " + end)
        if ci_cells:
            A(r" & \quad 95\% CI & " + " & ".join(ci_cells) + r" & \\*")
    if ci < len(CORE) - 1:
        A(r"\addlinespace[2pt]")
A(r"\midrule")
# ---- Panel B
panel(r"B\quad From the first to the last layer, all sixteen extractions")
A(r"\midrule")
A(r" & & & " + mc(3, r"Fold $L_0/L_{11}$") + r" & $\SV{1}$ (\%) & Rises \\")
A(r"\cmidrule(lr){4-6}")
A(r"Extraction & Design & $|V|$ & ER & $\mathrm{ER}_{\mathrm{var}}$ & PR & $L_0 \to L_{11}$ & ER/$\mathrm{ER}_{\mathrm{var}}$/PR \\")
A(r"\midrule")
prev = None
for i, (c, grp, des) in enumerate(ALL16):
    d = trend[c]
    if i in (5, 13):
        A(r"\addlinespace[2pt]")
    name = grp if grp != prev else ""
    prev = grp
    A(f"{name} & {des} & {f(d['n_V'], 0)} & {f(d['er_sv_fold_L0_over_L11'], 2)} & {f(d['er_var_fold_L0_over_L11'], 2)} & "
      f"{f(d['pr_fold_L0_over_L11'], 2)} & {f(100 * d['sv1_var_L0'])} $\\to$ {f(100 * d['sv1_var_L11'])} & "
      f"{d['er_sv_n_increases']}/{d['er_var_n_increases']}/{d['pr_n_increases']} \\\\")
A(r"\midrule")
# ---- Panel C
panel(r"C\quad Near-degenerate adjacent singular values, feature shuffle and the null direction")
A(r"\midrule")
hdr = " & ".join(LABEL[c] for c in CORE)
A(mc(2, "", "@{}l") + f" & {hdr} & All 16 \\\\")
A(r"\midrule")
A(mc(2, r"\itshape Near-degenerate pairs, $S_k/S_{k+1}$", "@{}l") + r" & & & & & & \\")
for k in range(1, 12):
    cells = [str(int(gk.loc[k, c])) for c in CORE]
    s_, n_ = int(gk_all.loc[k, "sum"]), int(gk_all.loc[k, "count"])
    A(mc(2, rf"\quad $k = {k}$ ($\SV{{{k}}}/\SV{{{k + 1}}}$)", "@{}l") + " & " + " & ".join(cells)
      + f" & {s_}/{n_} ({f(100 * s_ / n_)}\\%) \\\\")
cells = [f"{int(g_tot.loc[c, 'sum'])}/{int(g_tot.loc[c, 'count'])}" for c in CORE]
A(mc(2, r"\quad all pairs, $k = 1$--$11$", "@{}l") + " & " + " & ".join(cells)
  + f" & {f(int(gaps.nd.sum()), 0)}/{f(len(gaps), 0)} ({f(100 * gaps.nd.mean())}\\%) \\\\")
A(r"\addlinespace[2pt]")
A(mc(2, r"\itshape Feature shuffle", "@{}l") + r" & & & & & & \\")
all16 = [c for c, _, _ in ALL16]
for Lr in (0, 11):
    cells = [f(get(sh, context=c, layer=Lr).er_shuf_over_obs, 2) + r"$\times$" for c in CORE]
    allv = [get(sh, context=c, layer=Lr).er_shuf_over_obs for c in all16]
    A(mc(2, rf"\quad shuffled/observed ER, $L_{{{Lr}}}$", "@{}l") + " & " + " & ".join(cells)
      + f" & {rng_(allv, 2)} \\\\")
cells = [f(get(sh, context=c, layer=11).er_shuf_mean) for c in CORE]
allv = [get(sh, context=c, layer=11).er_shuf_mean for c in all16]
A(mc(2, r"\quad shuffled ER, $L_{11}$", "@{}l") + " & " + " & ".join(cells) + f" & {rng_(allv)} \\\\")
assert sh.er_shuf_sd.max() <= 0.12
A(r"\addlinespace[2pt]")
A(mc(2, r"\itshape Null direction", "@{}l") + r" & & & & & & \\")
nsv = spec.groupby("context")["n_sv_below_1e-3_rel"].agg(["min", "max"])
cells = [str(int(nsv.loc[c, "min"])) if nsv.loc[c, "min"] == nsv.loc[c, "max"] else
         f"{int(nsv.loc[c, 'min'])}--{int(nsv.loc[c, 'max'])}" for c in CORE]
assert nsv["min"].min() == nsv["max"].max() == 1
A(mc(2, r"\quad singular values $< 10^{-3} S_1$", "@{}l") + " & " + " & ".join(cells)
  + r" & 1 at all 192 \\")
cos0 = [nd[c]["0"]["abs_cos_last_dir_uniform"] for c in CORE]
cos11 = [nd[c]["11"]["abs_cos_last_dir_uniform"] for c in CORE]
assert len({f(x, 4) for x in cos0}) == 1 and len({f(x, 4) for x in cos11}) == 1   # identical to 4 decimals
relsv = [nd[c][L]["rel_last_sv"] for c in CORE for L in ("0", "11")]
A(mc(2, r"\quad $|\cos|$ with uniform, $L_0$/$L_{11}$", "@{}l") + " & "
  + mc(5, f"{f(cos0[0], 4)} / {f(cos11[0], 4)} in each of the five") + r" & \\")
mant = lambda x: f"{x / 1e-8:.1f}"
A(mc(2, r"\quad $S_{512}/S_1$ ($L_0$, $L_{11}$)", "@{}l") + " & "
  + mc(5, rf"${mant(min(relsv))}$--${mant(max(relsv))} \times 10^{{-8}}$") + r" & \\")
A(r"\midrule")
# ---- Panel D
panel(r"D\quad Detection controls, five main extractions")
A(r"\midrule")
A(mc(2, "", "@{}l") + f" & {hdr} & Range \\\\")
A(r"\midrule")
er11 = [get(spec, context=c, layer=11).er_sv for c in CORE]
A(mc(2, r"ER, all $V$, $L_{11}$ (panel A)", "@{}l") + " & " + " & ".join(f(x) for x in er11) + f" & {rng_(er11)} \\\\")
for Lr in (0, 11):
    vals = [get(mt, context=c, layer=Lr, scheme="pooled") for c in CORE]
    A(mc(2, rf"ER, matched 1{{,}}300 genes, $L_{{{Lr}}}$", "@{}l") + " & "
      + " & ".join(f"{f(v.er_mean)} $\\pm$ {f(v.er_sd)}" for v in vals) + f" & {rng_([v.er_mean for v in vals])} \\\\")
folds = [get(mt, context=c, layer=0, scheme="pooled").er_mean / get(mt, context=c, layer=11, scheme="pooled").er_mean
         for c in CORE]
A(mc(2, r"ER, matched, fold $L_0/L_{11}$", "@{}l") + " & " + " & ".join(f(x, 2) + r"$\times$" for x in folds)
  + f" & {rng_(folds, 2)} \\\\")
A(r"\addlinespace[2pt]")
for lab, col, scale, d in [(r"genes with $n_{\mathrm{obs}} \geq 3$", "n_genes", 1, 0),
                           (r"$n_{\mathrm{obs}} \geq 3$: ER, $L_{11}$", "er", 1, 1),
                           (r"$n_{\mathrm{obs}} \geq 3$: PR, $L_{11}$", "pr", 1, 1),
                           (r"$n_{\mathrm{obs}} \geq 3$: $\SV{1}$ (\%), $L_{11}$", "sv1_var", 100, 1)]:
    vals = [get(thr, context=c, layer=11, min_nobs=3)[col] * scale for c in CORE]
    A(mc(2, lab, "@{}l") + " & " + " & ".join(f(x, d) for x in vals) + f" & {rng_(vals, d)} \\\\")
A(r"\end{longtable}")
A(r"\endgroup")
OUT_TEX.write_text("\n".join(L) + "\n")
print(f"wrote {OUT_TEX} and {OUT_CSV} ({len(csv)} rows)")
print("null-direction relative last singular value range:", min(relsv), max(relsv))
