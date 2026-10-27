"""Build Supplementary Table S7 (regulatory edges, H3, and every edge-level sensitivity analysis).

Sources (each experiment's verifier-corrected SUMMARY.md is authoritative; numbers are read from its CSV/JSON outputs):
  revision3_experiments/E4a_edges     h3_summary.csv, h3_primary_<ctx>.csv, power.csv, sens_*.csv, reference_splits.csv,
                                      reference_split_reproduction.csv, verify/verify_results.json (reference-split percentiles)
  revision3_experiments/E4a_chipfix   sens_chip_fixed.csv, sens_chip_fixed_matchlevel0.csv (ChIP-seq standard, Panel B;
                                      supersedes E4a_edges/sens_chip.csv)
  revision3_experiments/E5_confounds  clogit_<ctx>.csv, residual_auroc.csv, stratum_weights.csv
  revision3_experiments/E4b_grn       grn_summary.csv, grn_deltas.csv, direction_test.csv
  revision3_experiments/E4c_probes    probe_results.csv, relational_signal.csv, verify/v1_rewired_did.csv,
                                      verify/v6_direction_cluster.csv, verify/v6_nb_sensitivity.csv, verify/v6_nb_counts.csv,
                                      verify/v6_auc_upper_ratio025.csv, verify/v3_paired_sources.csv, verify/v8_did_sources.csv,
                                      direction_test.csv

Writes
  sagmb_revision/supplement/S7.tex                       (LaTeX fragment, one longtable)
  sagmb_revision/supplementary_tables/Table_S7_edges.csv (complete values; one row per source row, tagged by panel)

No network statistic (vertex subsampling) enters this table: every interval is a two-way (factor x target) pigeonhole
bootstrap, a factor bootstrap (macro AUROC), a corrected resampled t (probes) or a gene-cluster bootstrap (probe direction),
so CONVENTION3_ERRATA item 3 does not apply. The only quantities computed here are ranges over layers, counts of layers
whose interval excludes the null value, and maxima/minima of reported values.
Run:  /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python build_S7.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
R3 = ROOT / "revision3_experiments"
E4A, E5, E4B, E4C = R3 / "E4a_edges", R3 / "E5_confounds", R3 / "E4b_grn", R3 / "E4c_probes"
OUT_TEX = HERE / "S7.tex"
OUT_CSV = HERE.parent / "supplementary_tables" / "Table_S7_edges.csv"

ALL16 = ["lung", "lung_external", "kidney_B", "kidney_ckpt_whole_human", "kidney_ckpt_kidney", "kidney_ckpt_brain",
         "immune", "immune_seed43", "immune_seed44", "immune_TSP14", "immune_TSP21", "immune_TSP25",
         "immune_Bcell", "immune_CD4T", "immune_CD8T", "kidney_A"]
LUNGS = ["lung", "lung_external"]
E5CTX = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
GRN5 = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
PROBE3 = ["lung", "lung_external", "kidney_B"]
CELLTYPE = ["immune_Bcell", "immune_CD4T", "immune_CD8T"]
NAME = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune",
        "kidney_A": "kidney~A", "immune_seed43": "immune (seed 43)", "immune_seed44": "immune (seed 44)",
        "immune_TSP14": "immune (donor TSP14)", "immune_TSP21": "immune (donor TSP21)",
        "immune_TSP25": "immune (donor TSP25)", "immune_Bcell": "immune (B cells)",
        "immune_CD4T": "immune (CD4 T cells)", "immune_CD8T": "immune (CD8 T cells)",
        "kidney_ckpt_whole_human": "kidney (whole-human ckpt)", "kidney_ckpt_kidney": "kidney (kidney ckpt)",
        "kidney_ckpt_brain": "kidney (brain ckpt)"}
COS = "signed cosine full-512"


# ----------------------------------------------------------------------------------------------- formatting
def isnan(x):
    return x is None or (isinstance(x, (float, np.floating)) and np.isnan(x))


def num(x, d=3, sign=False):
    """Number with a true minus sign; sign=True prints '+' for positive values."""
    if isnan(x):
        return "--"
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    if s.startswith("-"):
        if float(s) == 0.0:
            s = s[1:] if not sign else "+" + s[1:]
        else:
            return "$-$" + s[1:]
    if s.startswith("+"):
        return "$+$" + s[1:]
    return s


def ci(lo, hi, d=3):
    """Interval with an unbreakable thin space, so that a cell breaks between the estimate and its interval."""
    if lo < 0 or hi < 0:   # too wide for one line of a 1.57-cm column: allow a break after the comma
        return f"[{num(lo, d)}, {num(hi, d)}]"
    return f"[{num(lo, d)},\\,{num(hi, d)}]"


def rng(a, b, d=3, sign=False):
    if np.isclose(a, b) and f"{a:.{d}f}" == f"{b:.{d}f}":
        return num(a, d, sign)
    if a >= 0 and b >= 0 and not sign:
        return f"{num(a, d)}--{num(b, d)}"
    return f"{num(a, d, sign)} to {num(b, d, sign)}"


def pval(p):
    if isnan(p):
        return "--"
    if p >= 0.001:
        return f"{p:.4f}" if p < 0.1 else f"{p:.2f}"
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def intc(n):
    return f"{int(round(n)):,}".replace(",", "{,}")


def lay(layer):
    return f"$L_{{{int(layer)}}}$"


def mlab(label):
    """Analysis label spanning the two rows (lung, lung (ext.)) of a Panel-B block."""
    return f"\\multirow[t]{{4}}{{{W[0]}cm}}{{\\scriptsize\\raggedright {label}}}"


def frac12(n, tot=12):
    return f"{int(n)}/{tot}"


def nlayers_ranges(layers):
    """Compact list of layer indices, e.g. [1,2,3,4,5,6,7] -> L1--L7."""
    layers = sorted(int(x) for x in layers)
    if not layers:
        return ""
    out, start, prev = [], layers[0], layers[0]
    for x in layers[1:] + [None]:
        if x is not None and x == prev + 1:
            prev = x
            continue
        out.append(lay(start) if start == prev else f"{lay(start)}--{lay(prev)}")
        if x is not None:
            start = prev = x
    return ", ".join(out)


# ----------------------------------------------------------------------------------------------- grid
W = [1.9, 1.5] + [1.572] * 8             # column widths (cm)
SEP_PT = 1.5                               # \tabcolsep
SEPCM = 2 * SEP_PT * 2.54 / 72.27
NC = len(W)


def spanw(i, j):
    """Width of a cell spanning columns i..j, rounded down so that no spanning cell widens the columns it covers."""
    return np.floor((sum(W[i:j + 1]) + (j - i) * SEPCM) * 1000 - 1) / 1000


TOTAL = spanw(0, NC - 1)
COLSPEC = "@{}" + "".join(f">{{\\scriptsize\\raggedright\\arraybackslash}}p{{{w}cm}}" for w in W) + "@{}"
rows = []


def panel(title):
    if rows:
        rows.append("\\midrule")
    rows.append(f"\\multicolumn{{{NC}}}{{@{{}}p{{{TOTAL:.3f}cm}}@{{}}}}{{\\footnotesize\\textbf{{{title}}}}}\\\\*")


def sub(title):
    rows.append(f"\\multicolumn{{{NC}}}{{@{{}}p{{{TOTAL:.3f}cm}}@{{}}}}{{\\scriptsize\\textit{{{title}}}}}\\\\*")


def line(cells, rule=False, nobreak=False):
    out, i = [], 0
    for c in cells:
        text, span = (c, 1) if isinstance(c, str) else c
        if span == 1:
            out.append(text)
        else:
            j = i + span - 1
            spec = ("@{}" if i == 0 else "") + \
                f">{{\\scriptsize\\raggedright\\arraybackslash}}p{{{spanw(i, j):.3f}cm}}" + ("@{}" if j == NC - 1 else "")
            out.append(f"\\multicolumn{{{span}}}{{{spec}}}{{{text}}}")
        i += span
    assert i == NC, (i, cells)
    rows.append(" & ".join(out) + ("\\\\*" if nobreak else "\\\\"))
    if rule:
        # header rule that admits no page break between a header and its first data row (a booktabs \cmidrule
        # leaves a legal breakpoint at the glue below the rule)
        rows.append("\\noalign{\\nobreak\\vskip 1.5pt\\hrule height \\cmidrulewidth\\nobreak\\vskip 2pt\\nobreak}")


def header(cells):
    line(cells, rule=True, nobreak=True)


csv_parts = []


def keep(df, panel_id, source, **extra):
    d = df.copy()
    d.insert(0, "source_file", source)
    d.insert(0, "panel", panel_id)
    for k, v in extra.items():
        d[k] = v
    csv_parts.append(d)


# ======================================================================================= load E4a
summ = pd.read_csv(E4A / "h3_summary.csv").set_index("context")
prim = {c: pd.read_csv(E4A / f"h3_primary_{c}.csv") for c in ALL16}
power = pd.read_csv(E4A / "power.csv")
match = pd.read_csv(E4A / "sens_matching.csv")
short = pd.read_csv(E4A / "sens_shortcut.csv")
strata = pd.read_csv(E4A / "sens_strata.csv")
chip = pd.read_csv(R3 / "E4a_chipfix" / "sens_chip_fixed.csv")
chipm = pd.read_csv(R3 / "E4a_chipfix" / "sens_chip_fixed_matchlevel0.csv")  # fully matched candidates only
noise = pd.read_csv(E4A / "sens_label_noise.csv")
macro = pd.read_csv(E4A / "sens_macro.csv")
ctxl = pd.read_csv(E4A / "sens_contextual.csv")
pax = pd.read_csv(E4A / "sens_example_pax5_cd79a.csv")
legacy = pd.read_csv(E4A / "reference_splits.csv")   # copy of legacy_splits.csv
legrep = pd.read_csv(E4A / "reference_split_reproduction.csv")   # copy of legacy_split_reproduction.csv
legver = json.loads((E4A / "verify" / "verify_results.json").read_text())["legacy"]
repl = pd.read_csv(E4A / "h3_replication_lung_external.csv")
mc3 = json.loads((E4A / "verify" / "verify_part3.json").read_text())   # Monte-Carlo check of the borderline C label

for c in ALL16:
    assert len(prim[c]) == 12 and sorted(prim[c].layer) == list(range(12)), c
    assert np.isclose(prim[c].auroc.min(), summ.loc[c, "auroc_min"]) and np.isclose(prim[c].auroc.max(), summ.loc[c, "auroc_max"])

# ======================================================================================= Panel A
panel("A. Primary benchmark: within-factor AUROC of the full-dimensional signed cosine in all 16 extractions")
header(["Extraction", "Factors; positives; candidates", "AUROC, range $L_0$--$L_{11}$", "Best layer: AUROC [95\\% CI]",
        "Max.\\ UB", "Min.\\ max-$T$ $p$", "MDE: AUROC; $\\Delta$", "$|r|$ AUROC: extr.; 3k",
        "Max.\\ $\\Delta$ UB: extr.; 3k", "Verdict"])
VMAP = {"detected": "D", "AUROC>=0.60 excluded": "U", "not above co-expression": "C",
        "not above co-expression (3k tier)": "C"}


def verdict_codes(v):
    return [VMAP[p.strip()] for p in v.split("+")]


for c in ALL16:
    s = summ.loc[c]
    codes = []
    if s.detected:
        codes.append("D")
    if s.auroc_ge_060_excluded:
        codes.append("U")
    if s.not_above_coexpression:
        codes.append("C")
    assert codes == verdict_codes(s.verdict), (c, codes, s.verdict)
    codes3k = verdict_codes(s.verdict_3k_comparator)
    v = ", ".join(codes)
    extra3k = [x for x in codes3k if x not in codes]
    if extra3k:
        v += "; " + ", ".join(f"{x}$_{{\\mathrm{{3k}}}}$" for x in extra3k)
    p = prim[c]
    b = p.loc[p.auroc.idxmax()]
    assert int(b.layer) == int(s.best_layer) and np.isclose(b.ci_lo, s.best_ci_lo)
    line([NAME[c], f"{intc(s.n_tf)}; {intc(s.n_pos)}; {intc(s.n_neg)}", rng(s.auroc_min, s.auroc_max),
          f"{lay(s.best_layer)} {num(s.best_auroc)} {ci(s.best_ci_lo, s.best_ci_hi)}", num(s.max_ub95), pval(s.min_p_maxT),
          f"{num(s.mde_auroc_sim_two_way_bootstrap)}; {num(s.mde_delta_sim_two_way_bootstrap)}",
          f"{num(s.auroc_cx_ext)}; {num(s.auroc_cx_3k)}",
          f"{num(s.max_delta_cxext_ub95)}{'$^{*}$' if c == mc3['context'] else ''}; {num(s.max_delta_cx3k_ub95)}", v])
keep(summ.reset_index(), "A", "E4a_edges/h3_summary.csv")
for c in ALL16:
    keep(prim[c], "A", f"E4a_edges/h3_primary_{c}.csv")
keep(power, "A;D", "E4a_edges/power.csv")
keep(repl, "A", "E4a_edges/h3_replication_lung_external.csv")
assert np.isclose(summ.loc[mc3["context"], "max_delta_cxext_ub95"], mc3["replay_seed"]["delta_ub95"])
keep(pd.DataFrame([{"context": mc3["context"], "layer": mc3["layer"], "delta_ub95_replay": mc3["replay_seed"]["delta_ub95"],
                    "delta_ub95_indep_seed_1": mc3["indep_seed_1"]["delta_ub95"],
                    "delta_ub95_indep_seed_2": mc3["indep_seed_2"]["delta_ub95"]}]), "A", "E4a_edges/verify/verify_part3.json")

# ======================================================================================= Panel B
panel("B. Sensitivity analyses in the two lung extractions (cosine at $L_0$, and over the twelve layers)")
header(["Analysis", "Extraction", "Factors; positives", "$L_0$ [95\\% CI]", "Range $L_0$--$L_{11}$", "Max.\\ UB",
        "Layers with CI $>0.5$", "$|r|$ AUROC, same pairs", "$\\Delta$ vs $|r|$, range", "Max.\\ $\\Delta$ UB"])


def layer_block(df):
    d = df[df.layer >= 0].sort_values("layer")
    assert sorted(d.layer.astype(int)) == list(range(12)), (len(d),)
    return d


def cos_row(label, c, d, n_txt, cx=None, delta=True, flag=""):
    nb = c == "lung"
    d = layer_block(d)
    l0 = d[d.layer == 0].iloc[0]
    cells = [label, NAME[c], n_txt, f"{num(l0.auroc)} {ci(l0.ci_lo, l0.ci_hi)}", rng(d.auroc.min(), d.auroc.max()),
             num(d.ub95.max()), frac12((d.ci_lo > 0.5).sum()), "--" if cx is None else num(cx) + flag]
    if delta:
        cells += [rng(d.delta_cxext.min(), d.delta_cxext.max(), sign=True) + flag, num(d.delta_cxext_ub95.max()) + flag]
    else:
        cells += ["--", "--"]
    line(cells, nobreak=nb)


# B1 primary (reference)
for c in LUNGS:
    s = summ.loc[c]
    cos_row(mlab("All pairs (primary)") if c == "lung" else "", c, prim[c], f"{intc(s.n_tf)}; {intc(s.n_pos)}", s.auroc_cx_ext)
# B2 shortcut scores (layer-free)
for scorer, lab in [("shortcut: target TRRUST in-degree (total, incl. own edge)", "Shortcut: target in-degree, incl.\\ own edge"),
                    ("shortcut: target other-TF in-degree (excl. own edge)", "Shortcut: target's other-factor in-degree")]:
    for c in LUNGS:
        r = match[(match.context == c) & (match.block == "shortcut") & (match.scorer == scorer)]
        assert len(r) == 1
        r = r.iloc[0]
        line([mlab(lab) if c == "lung" else "", NAME[c], f"{intc(summ.loc[c, 'n_tf'])}; {intc(summ.loc[c, 'n_pos'])}",
              f"{num(r.auroc)} {ci(r.ci_lo, r.ci_hi)}", "(no layer)", num(r.ub95), "--", "--", "--", "--"], nobreak=c == "lung")
# B3 strata TF x other-TF in-degree
for c in LUNGS:
    blk = match[(match.context == c) & (match.block == "strata TF x other-TF in-degree")]
    d = blk[blk.scorer == COS]
    cxr = blk[blk.scorer == "|co-expression| extraction cells"].iloc[0]
    npos = d.n_pos_in_informative_strata.dropna().iloc[0]
    assert d.delta_cx_basis.iloc[0].startswith("original") and cxr.delta_cx_basis.startswith("original")
    cos_row(mlab("Strata: factor $\\times$ other-factor in-degree") if c == "lung" else "", c, d,
            f"--; {intc(npos)}$^{{a}}$", cxr.auroc, flag="$^{b}$")
# B4 ChIP-seq positives
for c in LUNGS:
    blk = chip[chip.context == c]
    d = blk[blk.scorer == COS]
    cx = blk[blk.scorer == "|co-expression| extraction cells"].iloc[0].auroc
    cos_row(mlab("ChIP-seq positives (not in TRRUST)") if c == "lung" else "", c, d,
            f"{intc(d.n_tf.iloc[0])}; {intc(d.n_pos.iloc[0])}", cx)
for c in LUNGS:
    blk = chipm[chipm.context == c]
    d = blk[blk.scorer == COS]
    cx = blk[blk.scorer == "|co-expression| extraction cells"].iloc[0].auroc
    cos_row(mlab("ChIP-seq positives, fully matched candidates") if c == "lung" else "", c, d,
            f"{intc(d.n_tf.iloc[0])}; {intc(d.n_pos.iloc[0])}", cx)
# B5 strata
STRATA = [("PMID count >= 2", "$\\geq 2$ publications"),
          ("PMID count = 1", "1 publication"),
          ("TF n_obs (scGPT input) >= 10% of extraction cells [CONVENTION3 detection]", "Factor is input in $\\geq 10\\%$ of cells"),
          ("TF n_obs (scGPT input) < 10% of extraction cells [CONVENTION3 detection]", "Factor is input in $<10\\%$ of cells"),
          ("TF detected in >= 10% of extraction cells", "Factor expressed in $\\geq 10\\%$ of cells"),
          ("TF detected in < 10% of extraction cells", "Factor expressed in $<10\\%$ of cells")]
for st, lab in STRATA:
    for c in LUNGS:
        blk = strata[(strata.context == c) & (strata.stratum == st)]
        d = blk[blk.scorer == COS]
        cxr = blk[blk.scorer == "|co-expression| extraction cells"].iloc[0]
        flag = "$^{b}$" if d.delta_cx_basis.iloc[0].startswith("original") else ""
        assert cxr.delta_cx_basis.startswith("original") == bool(flag), (c, st)
        cos_row(mlab(lab) if c == "lung" else "", c, d, f"{intc(d.n_tf.iloc[0])}; {intc(d.n_pos.iloc[0])}", cxr.auroc, flag=flag)
# B6 label contamination
for pi in [0.05, 0.10, 0.20]:
    for c in LUNGS:
        d = noise[(noise.context == c) & np.isclose(noise.pi, pi)].sort_values("layer")
        assert len(d) == 12
        l0 = d[d.layer == 0].iloc[0]
        line([mlab(f"Contamination $\\pi = {pi:.2f}$") if c == "lung" else "", NAME[c], "--", num(l0.auroc_corrected),
              rng(d.auroc_corrected.min(), d.auroc_corrected.max()), num(d.ub95_corrected.max()),
              f"UB $<0.60$: {frac12(d.ub95_corrected_below_060.sum())}", "--", "--", "--"], nobreak=c == "lung")
# B7 macro per-factor AUROC
for c in LUNGS:
    blk = macro[macro.context == c]
    d = blk[blk.scorer == COS].sort_values("layer")
    assert len(d) == 12
    l0 = d[d.layer == 0].iloc[0]
    cx = blk[blk.scorer == "cx_ext"].iloc[0]
    line([mlab("Macro AUROC$^{c}$") if c == "lung" else "", NAME[c], f"{intc(l0.n_tf)}; --",
          f"{num(l0.macro_auroc)} {ci(l0.ci_lo, l0.ci_hi)}", rng(d.macro_auroc.min(), d.macro_auroc.max()), num(d.ub95.max()),
          frac12((d.ci_lo > 0.5).sum()), num(cx.macro_auroc),
          (f"factors $>0.5$: {100 * d.frac_tf_above_half.min():.0f}--{100 * d.frac_tf_above_half.max():.0f}\\%", 2)], nobreak=c == "lung")
# B8 residualised AUROC (E5)
res = pd.read_csv(E5 / "residual_auroc.csv")
for c in LUNGS:
    d = res[res.context == c].sort_values("layer")
    assert len(d) == 12
    l0 = d[d.layer == 0].iloc[0]
    dd = d.delta_resid_minus_raw
    line([mlab("Residualised$^{d}$") if c == "lung" else "", NAME[c], f"{intc(l0.n_tf)}; {intc(l0.n_pos)}",
          f"{num(l0.auroc_resid)} {ci(l0.auroc_resid_lo, l0.auroc_resid_hi)}", rng(d.auroc_resid.min(), d.auroc_resid.max()),
          num(d.auroc_resid_upper95_onesided.max()), frac12((d.auroc_resid_lo > 0.5).sum()), "--",
          (f"resid.\\ $-$ raw: {rng(dd.min(), dd.max(), sign=True)}", 2)], nobreak=c == "lung")
# B9 conditional logit (E5)
clog = {c: pd.read_csv(E5 / f"clogit_{c}.csv") for c in E5CTX}
for c in LUNGS:
    x = clog[c]
    d = x[x.layer != "joint_mean_over_layers"].copy()
    d["layer"] = d.layer.astype(int)
    d = d.sort_values("layer")
    assert len(d) == 12
    l0 = d[d.layer == 0].iloc[0]
    j = x[x.layer == "joint_mean_over_layers"].iloc[0]
    line([mlab("Conditional logit OR$^{d}$") if c == "lung" else "", NAME[c], f"{intc(l0.n_tf)}; {intc(l0.n_pos)}",
          f"{num(l0.or_adj)} {ci(l0.or_adj_boot_lo, l0.or_adj_boot_hi)}", rng(d.or_adj.min(), d.or_adj.max()), "--",
          f"CI $>1$: {frac12((d.or_adj_boot_lo > 1).sum())}", f"unadj.\\ {rng(d.or_unadj.min(), d.or_unadj.max())}",
          (f"layer mean {num(j.or_adj)} {ci(j.or_adj_boot_lo, j.or_adj_boot_hi)}; unadj.\\ {num(j.or_unadj)}", 2)], nobreak=c == "lung")
sw = pd.read_csv(E5 / "stratum_weights.csv").set_index("context")
for c in LUNGS:
    s = sw.loc[c]
    line([mlab("Effective no.\\ of factors$^{e}$") if c == "lung" else "", NAME[c], f"{intc(s.n_tf)}; --",
          (f"pooled AUROC {s.eff_n_tf_auroc:.1f} (top 5: {100 * s.auroc_weight_top5:.0f}\\%)", 2),
          (f"cond.\\ logit {s.eff_n_tf_clogit:.1f} (top 5: {100 * s.clogit_info_top5:.0f}\\%)", 2), ("", 3)], nobreak=c == "lung")

keep(match, "B", "E4a_edges/sens_matching.csv")
keep(short, "B", "E4a_edges/sens_shortcut.csv")
keep(strata, "B", "E4a_edges/sens_strata.csv")
keep(chip, "B", "E4a_chipfix/sens_chip_fixed.csv", chip_candidates="all")
keep(chipm, "B", "E4a_chipfix/sens_chip_fixed_matchlevel0.csv", chip_candidates="fully matched only")
keep(noise, "B", "E4a_edges/sens_label_noise.csv")
keep(macro, "B", "E4a_edges/sens_macro.csv")
keep(res, "B", "E5_confounds/residual_auroc.csv")
for c in E5CTX:
    x = clog[c].copy()
    x["layer"] = x["layer"].astype(str)
    keep(x, "B", f"E5_confounds/clogit_{c}.csv")
keep(sw.reset_index(), "B", "E5_confounds/stratum_weights.csv")

# ======================================================================================= Panel C
panel("C. Cell-type extractions restricted to factors detected in $\\geq 10\\%$ of their cells, and the PAX5$\\to$CD79A example")
header(["Subset", "Factors; positives; candidates", "$L_0$ [95\\% CI]", "Best layer [95\\% CI]", "Range $L_0$--$L_{11}$",
        "Max.\\ UB", "$|r|$ extr.\\ [95\\% CI]", "$\\Delta$ at $L_0$ [95\\% CI]", "Layers with $\\Delta$ CI $>0$",
        "$\\Delta_{\\mathrm{3k}}$ at $L_0$ [95\\% CI]"])
for c in CELLTYPE:
    blk = ctxl[ctxl.context == c]
    cx = blk[blk.layer == -1].iloc[0]
    d = layer_block(blk[blk.scorer == COS])
    l0 = d[d.layer == 0].iloc[0]
    b = d.loc[d.auroc.idxmax()]
    pos_layers = d[d.delta_cxext_ci_lo > 0].layer.tolist()
    nl = frac12(len(pos_layers)) + (f" ({nlayers_ranges(pos_layers)})" if pos_layers else "")
    line([NAME[c], f"{intc(l0.n_tf)}; {intc(l0.n_pos)}; {intc(l0.n_neg)}", f"{num(l0.auroc)} {ci(l0.ci_lo, l0.ci_hi)}",
          f"{lay(b.layer)} {num(b.auroc)} {ci(b.ci_lo, b.ci_hi)}", rng(d.auroc.min(), d.auroc.max()), num(d.ub95.max()),
          f"{num(cx.auroc)} {ci(cx.ci_lo, cx.ci_hi)}",
          f"{num(l0.delta_cxext, sign=True)} {ci(l0.delta_cxext_ci_lo, l0.delta_cxext_ci_hi)}", nl,
          f"{num(l0.delta_cx3k, sign=True)} {ci(l0.delta_cx3k_ci_lo, l0.delta_cx3k_ci_hi)}"])
keep(ctxl, "C", "E4a_edges/sens_contextual.csv")

sub("PAX5$\\to$CD79A: rank of the documented pair among all of PAX5's pairs (documented and candidate), over $L_0$--$L_{11}$")
header(["Extraction", "Detected: PAX5; CD79A", "Rank ($n$ pairs)", "Candidates below", "$|r|$ rank: extr.; 3k",
        "Extraction", "Detected: PAX5; CD79A", "Rank ($n$ pairs)", "Candidates below", "$|r|$ rank: extr.; 3k"])
L = pax[~pax.layer.isin(["cx_ext", "cx_3k"])].copy()
L["layer"] = L.layer.astype(int)
PAXORDER = ["immune_Bcell", "immune_CD4T", "immune_CD8T", "immune", "immune_seed43", "immune_seed44", "immune_TSP14",
            "immune_TSP21", "immune_TSP25", "lung", "lung_external", "kidney_B", "kidney_A", "kidney_ckpt_whole_human",
            "kidney_ckpt_kidney", "kidney_ckpt_brain"]
assert sorted(PAXORDER) == sorted(L.context.unique())


def pax_cells(c):
    d = L[L.context == c]
    assert len(d) == 12
    cxe = pax[(pax.context == c) & (pax.layer == "cx_ext")].iloc[0].rank_among_all_candidates
    cx3 = pax[(pax.context == c) & (pax.layer == "cx_3k")].iloc[0].rank_among_all_candidates
    r = d.rank_among_all_candidates
    f = d.frac_negatives_below
    rk = f"{int(r.min())}" if r.min() == r.max() else f"{int(r.min())}--{int(r.max())}"
    fb = (f"{100 * f.min():.0f}\\%" if round(100 * f.min()) == round(100 * f.max())
          else f"{100 * f.min():.0f}--{100 * f.max():.0f}\\%")
    return [NAME[c], f"{100 * d.PAX5_det_frac.iloc[0]:.1f}\\%; {100 * d.CD79A_det_frac.iloc[0]:.1f}\\%",
            f"{rk} ({int(d.n_candidates.iloc[0])})", fb, f"{int(cxe)}; {int(cx3)}"]


for i in range(0, len(PAXORDER), 2):
    line(pax_cells(PAXORDER[i]) + pax_cells(PAXORDER[i + 1]))
imm9 = [c for c in PAXORDER if c.startswith("immune")]
ri = L[L.context.isin(imm9)].rank_among_all_candidates
line([(f"All nine immune extractions: rank {int(ri.min())}--{int(ri.max())} at every layer.", NC)])
keep(pax, "C", "E4a_edges/sens_example_pax5_cd79a.csv")

# ======================================================================================= Panel D
panel("D. An 80/20 held-out split of the kidney~A pairs (49 positives, 147 candidates, 27 factors, 117 targets): 500 random splits")
header([("Scorer ($L_0$)", 2), "Quantity", "Mean", "SD", "2.5--97.5\\%", "Min--max", "Fraction $>0.5$ ($>0$ for $\\Delta$)",
        ("Reference split (percentile)", 2)])
LEG = [("full-512 L0 cosine [submitted: uncentred]", "Full 512-d, uncentred (screen convention)", "full-512 L0 cosine [submitted: uncentred]"),
       ("full-512 L0 cosine [paper: V-centred]", "Full 512-d, $V$-centred", None),
       ("SV2-SV7 L0 cosine [submitted: zero-row SVD]", "$\\SV{2}$--$\\SV{7}$, zero-row SVD (reference basis)", "SV2-SV7 L0 cosine [submitted: zero-row SVD]"),
       ("SV2-SV7 L0 cosine [paper: observed-gene SVD]", "$\\SV{2}$--$\\SV{7}$, observed-gene SVD", None),
       ("SV5-SV7 L0 cosine [submitted: zero-row SVD]", "$\\SV{5}$--$\\SV{7}$, zero-row SVD (reference basis)", "SV5-SV7 L0 cosine [submitted: zero-row SVD]"),
       ("SV5-SV7 L0 cosine [paper: observed-gene SVD]", "$\\SV{5}$--$\\SV{7}$, observed-gene SVD", None),
       ("|co-expression| (3,000 cells, run_trees definition)", "$|r|$, 3{,}000 cells", "|co-expression| (3,000 cells, run_trees definition)"),
       ("random", "Random scores", "random"),
       ("stratified_full512_submitted", "Full 512-d, uncentred (screen convention), within-factor AUROC", None)]
obs_delta = power[power.design.str.startswith("kidney_A legacy")].observed_delta_cxext_min.iloc[0]
for key, lab, repkey in LEG:
    blk = legacy[legacy.scorer == key]
    assert len(blk) >= 1, key
    for qi, (_, r) in enumerate(blk.iterrows()):
        isd = r.quantity.startswith("delta")
        sub_txt = "--"
        if repkey is not None and not isd:
            rr = legrep[legrep.scorer == repkey].iloc[0]
            sub_txt = num(rr.auroc)
            if key.startswith("full-512") and "uncentred" in key:
                assert np.isclose(rr.auroc, 0.5914202415660141, atol=1e-6)
                sub_txt += f" {ci(rr.ci_lo, rr.ci_hi)} ({legver['pct_0591']:.1f}th)"
            if key.startswith("|co-expression|"):
                sub_txt += f" ({legver['pct_cx_0556']:.1f}th)"
        if isd and key.startswith("full-512") and "uncentred" in key:
            sub_txt = f"{num(obs_delta, sign=True)} ({legver['pct_delta_submitted']:.1f}th)"
            assert np.isclose(r.frac_above_0_5, legver["frac_delta_pos"])
        line([(lab if qi == 0 else "", 2), "$\\Delta$ vs $|r|$" if isd else "AUROC", num(r["mean"], sign=isd), num(r.sd),
              f"{num(r.p2_5, sign=isd)} to {num(r.p97_5, sign=isd)}" if isd else f"{num(r.p2_5)}--{num(r.p97_5)}",
              f"{num(r['min'], sign=isd)} to {num(r['max'], sign=isd)}" if isd else f"{num(r['min'])}--{num(r['max'])}",
              f"{r.frac_above_0_5:.3f}", (sub_txt, 2)])
keep(legacy, "D", "E4a_edges/reference_splits.csv")
keep(legrep, "D", "E4a_edges/reference_split_reproduction.csv")
keep(pd.DataFrame([legver]), "D", "E4a_edges/verify/verify_results.json (reference-split percentiles)")

pw = power[power.design.str.startswith("kidney_A legacy")]
r0 = pw.iloc[0]
sub(f"Minimum detectable effect of that design (80\\% power, two-sided $\\alpha = 0.05$; target ICC {r0.rho_target_icc:.2f}, "
    f"factor ICC {r0.rho_tf_icc:.2f}; {int(r0.nsim_per_grid)} simulations per grid point)")
header([("Statistic", 2), "Test", "MDE AUROC", "Size at null", "MDE $\\Delta$", "Size at null ($\\Delta$)", ("", 3)])
TESTNAME = {"two_way_bootstrap": "two-way bootstrap", "tf_only_bootstrap": "factor-only bootstrap",
            "iid_delong": "i.i.d.\\ DeLong", "permutation": "Wald, perm.\\ SD"}
for stat, slab in [("pooled AUROC (as submitted)", "Pooled AUROC (positives against all held-out non-edges)"),
                   ("within-TF stratified AUROC", "Within-factor AUROC")]:
    blk = pw[pw.statistic == stat]
    first = True
    for t in ["two_way_bootstrap", "tf_only_bootstrap", "iid_delong", "permutation"]:
        r = blk[blk.test == t]
        if not len(r):
            continue
        r = r.iloc[0]
        line([(slab if first else "", 2), TESTNAME[t], num(r.mde_auroc_sim), f"{r.size_auroc_sim:.3f}", num(r.mde_delta_sim),
              "--" if isnan(r.size_delta_sim) else f"{r.size_delta_sim:.3f}", ("", 3)])
        first = False
    r = blk.iloc[0]
    line([("", 2), "oracle SD", num(r.mde_auroc_oracle_sd), "--", num(r.mde_delta_oracle_sd), "--", ("", 3)])
r = pw.iloc[0]
line([("Either statistic", 2), "Hanley--McNeil", num(r.mde_auroc_hanley_mcneil), "--", num(r.mde_delta_hanley_mcneil), "--",
      ("", 3)])

# ======================================================================================= Panel E
gs = pd.read_csv(E4B / "grn_summary.csv")
gd = pd.read_csv(E4B / "grn_deltas.csv")
panel("E. GENIE3 (extra-trees) and GRNBoost2-style boosting on the same benchmark (T1, the extraction's own cells; T2, 3{,}000 cells)")
header(["Extraction (cosine $L_0$)", "Method, tier (cells)", "AUROC, seed ensemble [95\\% CI]", "UB", "Seed mean (SD; seeds)",
        "$|r|$, same cells", "$\\Delta$ vs $|r|$ [95\\% CI]", "$\\Delta$ vs cosine $L_0$ [95\\% CI]",
        "$\\Delta$ vs cosine, range $L_0$--$L_{11}$", "T2 $-$ T1 [95\\% CI]"])
MNAME = {"genie3": "GENIE3", "grnboost2": "GRNBoost2"}
n_cos_ex0 = 0
n_cos_tot = 0
n_cx_ex0 = 0
n_cos_ub_lt0 = []
cx_ub_large = []
for c in GRN5:
    cos0 = gs[(gs.context == c) & (gs.method == "scgpt_cos") & (gs.layer == 0)].iloc[0].auroc
    assert np.isclose(cos0, prim[c][prim[c].layer == 0].iloc[0].auroc, atol=1e-12)
    first = True
    for m in ["genie3", "grnboost2"]:
        for t in ["T1", "T2"]:
            e = gs[(gs.context == c) & (gs.method == m) & (gs.tier == t) & (gs.seed == "ens")]
            ms = gs[(gs.context == c) & (gs.method == m) & (gs.tier == t) & (gs.seed == "seed_mean_sd")]
            cx = gs[(gs.context == c) & (gs.method == "abs_coexpr") & (gs.tier == t)]
            assert len(e) == 1 and len(ms) == 1 and len(cx) == 1
            e, ms, cx = e.iloc[0], ms.iloc[0], cx.iloc[0]
            dcx = gd[(gd.context == c) & (gd.method == m) & (gd.tier == t) & (gd.comparison == "grn_minus_abs_coexpr_same_cells")]
            dcos = gd[(gd.context == c) & (gd.method == m) & (gd.tier == t) & (gd.comparison == "grn_minus_scgpt_cos")]
            assert len(dcx) == 1 and len(dcos) == 12
            dcx = dcx.iloc[0]
            dcos = dcos.sort_values("comparator_layer")
            d0 = dcos[dcos.comparator_layer == 0].iloc[0]
            n_cos_ex0 += int(((dcos.ci_lo > 0) | (dcos.ci_hi < 0)).sum())
            n_cos_tot += len(dcos)
            n_cx_ex0 += int((dcx.ci_lo > 0) or (dcx.ci_hi < 0))
            n_cos_ub_lt0 += [(c, m, t, int(r.comparator_layer), r.upper95) for _, r in dcos.iterrows() if r.upper95 < 0]
            if c in ["lung", "lung_external", "kidney_B"]:
                cx_ub_large.append(dcx.upper95)
            t21 = "--"
            if t == "T2":
                dt = gd[(gd.context == c) & (gd.method == m) & (gd.comparison == "T2_minus_T1_same_method")]
                assert len(dt) == 1
                dt = dt.iloc[0]
                t21 = f"{num(dt.delta, sign=True)} {ci(dt.ci_lo, dt.ci_hi)}"
            lab = f"{NAME[c]} ({num(cos0)})" if first else ""
            first = False
            line([lab, f"{MNAME[m]} {t} ({intc(e.n_cells)})", f"{num(e.auroc)} {ci(e.ci_lo, e.ci_hi)}", num(e.upper95),
                  f"{num(ms.auroc)} ({num(ms.seed_sd, 4)}; {int(ms.n_seeds)})", num(cx.auroc),
                  f"{num(dcx.delta, sign=True)} {ci(dcx.ci_lo, dcx.ci_hi)}", f"{num(d0.delta, sign=True)} {ci(d0.ci_lo, d0.ci_hi)}",
                  rng(dcos.delta.min(), dcos.delta.max(), sign=True), t21])
assert len(n_cos_ub_lt0) == 1, n_cos_ub_lt0
oc = n_cos_ub_lt0[0]
line([(f"Two-sided intervals excluding 0: $\\Delta$ vs $|r|$, {n_cx_ex0} of 20; $\\Delta$ vs cosine, {n_cos_ex0} of {n_cos_tot} "
       f"(extraction $\\times$ method $\\times$ tier $\\times$ layer); one one-sided bound lies below 0 ({NAME[oc[0]]} "
       f"{MNAME[oc[1]]} {oc[2]} vs {lay(oc[3])}, UB {num(oc[4], 4)}; unadjusted). One-sided UB of $\\Delta$ vs $|r|$ in lung, "
       f"lung (ext.) and kidney~B: at most {num(max(cx_ub_large))}.", NC)])
keep(gs, "E", "E4b_grn/grn_summary.csv")
keep(gd, "E", "E4b_grn/grn_deltas.csv")
keep(pd.read_csv(E4B / "verify" / "v5_h3_margins.csv"), "E", "E4b_grn/verify/v5_h3_margins.csv")

# ======================================================================================= Panel F
pr = pd.read_csv(E4C / "probe_results.csv")
rs = pd.read_csv(E4C / "relational_signal.csv")
did = pd.read_csv(E4C / "verify" / "v1_rewired_did.csv")
nbs = pd.read_csv(E4C / "verify" / "v6_nb_sensitivity.csv")
nbc = pd.read_csv(E4C / "verify" / "v6_nb_counts.csv")
ubg = pd.read_csv(E4C / "verify" / "v6_auc_upper_ratio025.csv")
paired = pd.read_csv(E4C / "verify" / "v3_paired_sources.csv")
dids = pd.read_csv(E4C / "verify" / "v8_did_sources.csv")
FEAT = {"bilinear": "bilinear", "tg_bilinear": "$[x_j$, bilinear$]$", "concat": "concatenation", "hadamard": "Hadamard",
        "absdiff": "$|x_i - x_j|$", "tg_only": "target only", "tg_hadamard": "$[x_j$, Hadamard$]$"}


def m1(c, src="scgpt", feat=None, regime="C3"):
    d = pr[(pr.context == c) & (pr.model == "M1") & (pr.regime == regime) & (pr.source == src)]
    return d if feat is None else d[d.feature == feat]


def rsl(c, src, feat, regime="C3", model="M1"):
    return rs[(rs.context == c) & (rs.model == model) & (rs.regime == regime) & (rs.source == src) & (rs.feature == feat)]


def nbrow(c, src, feat, comp, kind, layer):
    d = nbs[(nbs.context == c) & (nbs.source == src) & (nbs.feature == feat) & (nbs.comparison == comp) &
            (nbs.ratio_kind == kind) & (nbs.layer == layer)]
    assert len(d) == 1, (c, src, feat, comp, kind, layer, len(d))
    return d.iloc[0]


def nbcount(c, src, feat, comp, kind):
    d = nbc[(nbc.context == c) & (nbc.source == src) & (nbc.feature == feat) & (nbc.comparison == comp) & (nbc.ratio_kind == kind)]
    assert len(d) == 1, (c, src, feat, comp, kind)
    return d.iloc[0]


panel("F. Supervised probes on pair-disjoint folds (neither gene of a test pair seen in training; 5 $\\times$ 5 folds)")
header([("Quantity", 4), ("lung", 2), ("lung (ext.)", 2), ("kidney~B", 2)])


def frow(label, vals):
    line([(label, 4)] + [(v, 2) for v in vals])


vals = []
for c in PROBE3:
    s = summ.loc[c]
    vals.append(f"{intc(s.n_pos + s.n_neg)} pairs; {intc(s.n_pos)} positives")
frow("Benchmark (each test fold holds about 4\\% of the pairs)", vals)

v1, v2, v3 = [], [], []
for c in PROBE3:
    d = m1(c)
    d = d[~d.feature.isin(["tf_only"])]
    b = d.loc[d.auc_mean.idxmax()]
    u = ubg[(ubg.context == c) & (ubg.source == "scgpt") & (ubg.feature == b.feature) & (ubg.layer == b.layer)]
    assert len(u) == 1 and np.isclose(u.iloc[0].auc, b.auc_mean)
    v1.append(f"{FEAT[b.feature]}, {lay(b.layer)}")
    v2.append(f"{num(b.auc_mean)} {ci(b.auc_ci_lo, b.auc_ci_hi)}; {num(u.iloc[0].hi_025)}")
    t = m1(c, feat="tg_only")
    assert len(t) == 12
    v3.append(rng(t.auc_mean.min(), t.auc_mean.max()))
frow("Best scGPT probe: feature, layer", v1)
frow("\\quad AUROC [95\\% CI, pair ratio]; upper end of the 95\\% CI, gene ratio", v2)
frow("Target-only probe: AUROC, range over layers (factor-only: 0.500 by construction)", v3)

v = []
for c in PROBE3:
    d = rsl(c, "scgpt", "concat")
    v.append(rng(d.delta_node.min(), d.delta_node.max(), sign=True))
frow("Concatenation (additive) $-$ target-only: range", v)
for feat in ["hadamard", "absdiff"]:
    v = []
    for c in PROBE3:
        d = rsl(c, "scgpt", feat)
        cp, cg = nbcount(c, "scgpt", feat, "delta_node", "pairs"), nbcount(c, "scgpt", feat, "delta_node", "genes_0.25")
        assert int(cp.n_lo_gt0) == 0
        v.append(f"{rng(d.delta_node.min(), d.delta_node.max(), sign=True)}; {int(cp.n_hi_lt0)}; {int(cg.n_hi_lt0)}")
    frow(f"{FEAT[feat]} $-$ target-only: range; layers with CI $<0$ (pair; gene ratio)", v)

va, vb, vc, vd = [], [], [], []
for c in PROBE3:
    d = rsl(c, "scgpt", "bilinear").sort_values("layer")
    assert len(d) == 12
    va.append(rng(d.delta_node.min(), d.delta_node.max(), sign=True))
    b = d.loc[d.delta_node.idxmax()]
    g = nbrow(c, "scgpt", "bilinear", "delta_node", "genes_0.25", int(b.layer))
    p = nbrow(c, "scgpt", "bilinear", "delta_node", "pairs", int(b.layer))
    assert np.isclose(p.lo, b.delta_node_lo) and np.isclose(p.hi, b.delta_node_hi)
    vb.append(f"{lay(b.layer)}: {num(b.delta_node, sign=True)} {ci(b.delta_node_lo, b.delta_node_hi)} {ci(g.lo, g.hi)}")
    cp, cg = nbcount(c, "scgpt", "bilinear", "delta_node", "pairs"), nbcount(c, "scgpt", "bilinear", "delta_node", "genes_0.25")
    assert int(cp.n_lo_gt0) == int((d.delta_node_lo > 0).sum())
    vc.append(f"{int(cp.n_lo_gt0)}; {int(cg.n_lo_gt0)}")
    t = rsl(c, "scgpt", "tg_bilinear")
    tp, tg = nbcount(c, "scgpt", "tg_bilinear", "delta_node", "pairs"), nbcount(c, "scgpt", "tg_bilinear", "delta_node", "genes_0.25")
    vd.append(f"{rng(t.delta_node.min(), t.delta_node.max(), sign=True)}; {int(tp.n_lo_gt0)}; {int(tg.n_lo_gt0)}")
frow("Bilinear $-$ target-only: range over layers", va)
frow("\\quad largest gain [95\\% CI, pair ratio] [gene ratio]", vb)
frow("\\quad layers with CI $>0$: pair; gene ratio", vc)
frow("$[x_j$, bilinear$]$ $-$ target-only: range; layers with CI $>0$ (pair; gene ratio)", vd)

va, vb, vc, vd = [], [], [], []
for c in PROBE3:
    x = did[(did.context == c) & (did.test == "bilinear_vs_rewired")].sort_values("layer")
    rr = rsl(c, "scgpt", "bilinear").sort_values("layer")
    assert len(x) == 12 and np.allclose(x.p_emp.values, rr.p_emp_rewired.values) and np.allclose(x.p_maxT.values, rr.p_maxT_layers.values)
    va.append(f"{frac12((x.p_emp < 0.05).sum())}; {frac12((x.p_maxT < 0.05).sum())}")
    y = did[(did.context == c) & (did.test == "concat_vs_rewired")].sort_values("layer")
    vb.append(f"{frac12((y.p_emp < 0.05).sum())}; {frac12((y.p_maxT < 0.05).sum())}; true {rng(y.obs.min(), y.obs.max())} vs "
              f"{rng(y.null_mean.min(), y.null_mean.max())}")
    z = did[(did.context == c) & (did.test == "DiD_bilinear_minus_concat")].sort_values("layer")
    sig = z[z.p_maxT < 0.05].layer.tolist()
    vc.append(f"{rng(z.obs.min(), z.obs.max(), sign=True)} (rewired {rng(z.null_mean.min(), z.null_mean.max(), sign=True)}); "
              f"{frac12(len(sig))}" + (f" ({nlayers_ranges(sig)})" if sig else ""))
    cp = nbcount(c, "scgpt", "bilinear", "relational_signal", "pairs")
    cg = nbcount(c, "scgpt", "bilinear", "relational_signal", "genes_0.25")
    vd.append(f"{int(cp.n_lo_gt0)}; {int(cg.n_lo_gt0)}")
frow("Bilinear vs the same probe on rewired labels ($N = 100$): layers $p<0.05$; max-$T<0.05$", va)
frow("Concatenation vs rewired labels: layers $p<0.05$; max-$T<0.05$; AUROC true vs rewired (mean)", vb)
frow("Difference in differences, bilinear $-$ concatenation: true labels, range (rewired-label mean); layers max-$T<0.05$", vc)
frow("Relational signal (probe $-$ max[rewired, target-only]): layers with CI $>0$, pair; gene ratio", vd)

v = []
for c in PROBE3:
    d = rsl(c, "scgpt", "bilinear")
    cp, cg = nbcount(c, "scgpt", "bilinear", "delta_cosine", "pairs"), nbcount(c, "scgpt", "bilinear", "delta_cosine", "genes_0.25")
    v.append(f"{rng(d.delta_cosine.min(), d.delta_cosine.max(), sign=True)}; {int(cp.n_lo_gt0)}; {int(cg.n_lo_gt0)}")
frow("Bilinear $-$ signed cosine on the same folds: range; layers with CI $>0$ (pair; gene ratio)", v)

for src, slab in [("static", "scGPT input gene-token embedding (no cell context)"), ("esm2", "ESM2-15B per-protein embedding (UCE)")]:
    va, vb = [], []
    for c in PROBE3:
        d = rsl(c, src, "bilinear")
        assert len(d) == 1
        d = d.iloc[0]
        g = nbrow(c, src, "bilinear", "delta_node", "genes_0.25", -1)
        va.append(f"{num(d.probe_auc)}; {num(d.delta_node, sign=True)} {ci(d.delta_node_lo, d.delta_node_hi)} {ci(g.lo, g.hi)}")
        dd = dids[(dids.context == c) & (dids.source == src)].iloc[0]
        vb.append(f"{pval(d.p_emp_rewired)}; {pval(dd.p_emp)}")
    frow(f"{slab}, bilinear: AUROC; $-$ target-only [pair ratio] [gene ratio]", va)
    frow("\\quad $p$ vs rewired labels; difference-in-differences $p$", vb)
v, v2 = [], []
for c in PROBE3:
    x = paired[(paired.context == c) & paired.comparison.str.contains("static_bilinear")]
    assert len(x) == 12
    e = paired[(paired.context == c) & paired.comparison.str.contains("esm2_bilinear")]
    v.append(f"{(x.lo > 0).sum()}; {(x.hi < 0).sum()}")
    v2.append(f"{(x.lo_excess > 0).sum()}; {(e.lo_excess > 0).sum()}")
frow("scGPT bilinear $-$ input-embedding bilinear (paired): layers with CI $>0$; $<0$", v)
frow("scGPT gain over its target-only control $-$ the same gain of the input embedding; of ESM2-15B: layers with CI $>0$", v2)

va, vb = [], []
for c in PROBE3:
    d = rs[(rs.context == c) & rs.model.isin(["M2", "M3", "M2_tuned", "M3_tuned"]) & (rs.source == "scgpt")]
    d3, d2 = d[d.regime == "C3"], d[d.regime == "C2"]
    va.append(f"{rng(d3.probe_auc.min(), d3.probe_auc.max())}; {(d3.delta_node_lo > 0).sum()}; {(d3.delta_node_hi < 0).sum()} of {len(d3)}")
    vb.append(f"{(d2.delta_node_lo > 0).sum()} of {len(d2)}")
frow("Gradient boosting and MLP (fixed and tuned; $L_0$, $L_5$, $L_{11}$; concatenation, Hadamard): AUROC range; rows above; "
     "below own target-only", va)
frow("\\quad held-out factors (targets seen): rows above own target-only", vb)

keep(pr, "F", "E4c_probes/probe_results.csv")
keep(rs, "F", "E4c_probes/relational_signal.csv")
keep(did, "F", "E4c_probes/verify/v1_rewired_did.csv")
keep(nbs, "F", "E4c_probes/verify/v6_nb_sensitivity.csv")
keep(nbc, "F", "E4c_probes/verify/v6_nb_counts.csv")
keep(ubg, "F", "E4c_probes/verify/v6_auc_upper_ratio025.csv")
keep(paired, "F", "E4c_probes/verify/v3_paired_sources.csv")
keep(dids, "F", "E4c_probes/verify/v8_did_sources.csv")

# ======================================================================================= Panel G
gdir = pd.read_csv(E4B / "direction_test.csv")
panel("G. Direction: documented pairs $i \\to j$ of two factors whose reverse $j \\to i$ is in no regulatory resource")
sub("GENIE3 and GRNBoost2 importances, seed ensemble: fraction of untied pairs with importance($i\\to j$) $>$ importance($j\\to i$)")
header(["Extraction", "Method", "Pairs; untied T1; T2", "T1 fraction", "T1 sign-test $p$", "T1 regulator-cluster CI",
        "T2 fraction", "T2 sign-test $p$", "T2 regulator-cluster CI", "T1 per-seed range"])
n_p = 0
n_p_lt = 0
for c in GRN5:
    first = True
    for m in ["genie3", "grnboost2"]:
        e1 = gdir[(gdir.context == c) & (gdir.method == m) & (gdir.tier == "T1") & (gdir.seed == "ens")].iloc[0]
        e2 = gdir[(gdir.context == c) & (gdir.method == m) & (gdir.tier == "T2") & (gdir.seed == "ens")].iloc[0]
        ps = gdir[(gdir.context == c) & (gdir.method == m) & (gdir.tier == "T1") & ~gdir.seed.isin(["ens"])]
        assert len(ps) == 5
        for e in (e1, e2):
            n_p += 1
            n_p_lt += int(e.sign_test_p_two_sided < 0.05)
        line([NAME[c] if first else "", MNAME[m], f"{intc(e1.n_pairs)}; {intc(e1.n_untied)}; {intc(e2.n_untied)}",
              num(e1.frac_forward_gt), pval(e1.sign_test_p_two_sided), ci(e1.src_cluster_lo, e1.src_cluster_hi),
              num(e2.frac_forward_gt), pval(e2.sign_test_p_two_sided), ci(e2.src_cluster_lo, e2.src_cluster_hi),
              rng(ps.frac_forward_gt.min(), ps.frac_forward_gt.max())])
        first = False
line([(f"Sign tests with $p<0.05$: {n_p_lt} of {n_p} (none below the Bonferroni threshold 0.05/{n_p} = {0.05 / n_p:.4f}). "
       "$|r|$ and the cosine are symmetric and cannot score direction.", NC)])
keep(gdir, "G", "E4b_grn/direction_test.csv")

pdir = pd.read_csv(E4C / "verify" / "v6_direction_cluster.csv")
pdt = pd.read_csv(E4C / "direction_test.csv")
sub("Logistic probe on the ordered concatenation, pair-disjoint folds: fraction of pairs scored forward $>$ reverse")
header(["Extraction", "Pairs; genes", "Fraction, range $L_0$--$L_{11}$", "Layers, binomial CI excl.\\ 0.5",
        "Layers, gene-cluster CI excl.\\ 0.5", ("Largest fraction [gene-cluster CI]", 2), "Strict subset: pairs; range",
        "Out-degree $i > j$", "In-degree $j > i$"])
n_cl = 0
for c in PROBE3:
    d = pdir[pdir.context == c].sort_values("layer")
    t = pdt[pdt.context == c].sort_values("layer")
    assert len(d) == 12 and len(t) == 12 and np.allclose(d.frac_forward.values, t.frac_forward_higher.values)
    nb = int(((t.ci_lo > 0.5) | (t.ci_hi < 0.5)).sum())
    ex = d[(d.ci_lo_cluster > 0.5) | (d.ci_hi_cluster < 0.5)]
    n_cl += len(ex)
    b = d.loc[d.frac_forward.idxmax()]
    line([NAME[c], f"{intc(d.n_pairs.iloc[0])}; {intc(d.n_genes.iloc[0])}", rng(d.frac_forward.min(), d.frac_forward.max()),
          frac12(nb), frac12(len(ex)) + (f" ({nlayers_ranges(ex.layer)})" if len(ex) else ""),
          (f"{lay(b.layer)} {num(b.frac_forward)} {ci(b.ci_lo_cluster, b.ci_hi_cluster)}", 2),
          f"{intc(t.n_strict.iloc[0])}; {rng(t.frac_forward_strict.min(), t.frac_forward_strict.max())}",
          num(t.frac_outdeg_i_gt_j.iloc[0]), num(t.frac_indeg_j_gt_i.iloc[0])])
line([(f"Gene-cluster intervals excluding 0.5: {n_cl} of 36 extraction $\\times$ layer combinations.", NC)])
keep(pdir, "G", "E4c_probes/verify/v6_direction_cluster.csv")
keep(pdt, "G", "E4c_probes/direction_test.csv")

# ======================================================================================= caption + assembly
s_l, s_x = summ.loc["lung"], summ.loc["lung_external"]
caption = (
    "\\textbf{Regulatory edges (H3): the benchmark in every extraction and every edge-level sensitivity analysis.} "
    "Positives: TRRUST pairs with both genes observed; candidates: up to ten per positive, pairs of the same factor whose "
    "target is a documented target of another factor and that are absent, in either direction, from TRRUST and "
    "from DoRothEA at any confidence level, and matched where possible to the positive on the target's decile of TRRUST in-degree and of "
    "$n_{\\mathrm{obs}}$ (on both for MATCH_RANGE of candidates, otherwise on in-degree or not at all; CSV column "
    "\\texttt{frac\\_cand\\_fully\\_matched}); factors with $\\geq 3$ of each. Statistic: within-factor stratified AUROC of the signed "
    "cosine of the $V$-centred 512-dimensional vectors at layer $L_\\ell$. CI, 95\\% interval; UB, one-sided 95\\% upper "
    "bound; both from the two-way (factor $\\times$ target) pigeonhole bootstrap, $B = 2{,}000$ unless stated (1{,}000 for "
    "some strata, CSV column \\texttt{B}); Max.\\ UB, the largest UB over the twelve layers. $|r|$, absolute Pearson "
    "co-expression on the same pairs, on the extraction's own cells (extr.; 256, or 192 for kidney~B and the checkpoints) or "
    "on 3{,}000 cells of the same population (3k); $\\Delta$, cosine minus $|r|$ AUROC on the same bootstrap draws. "
    "\\textbf{A}: max-$T$ $p$, within-factor label permutation applied jointly to the twelve layers ($N = 10{,}000$; floor "
    "$1.0\\times10^{-4}$); MDE, the AUROC (and $\\Delta$) detected with 80\\% power at two-sided $\\alpha = 0.05$ in a "
    "simulation that resamples factors as clusters with the observed between-factor heterogeneity and tests with the two-way "
    "bootstrap. Verdict: D, detected (max-$T$ $p < 0.05$); U, every layer's UB $< 0.60$ (AUROC $\\geq 0.60$ excluded); C, "
    "every layer's $\\Delta$ UB $< 0.05$ (not above co-expression), against the extraction cells (C$_{\\mathrm{3k}}$: "
    "against the 3k tier only). $^{*}$MC_TXT Pre-declared replication (same sign, and a lung (ext.) CI excluding 0.5): the sign agrees "
    "at REPL_SIGN/12 layers and the CI excludes 0.5 at REPL_EXCL/12, so it is not met. "
    "\\textbf{B}: Shortcut rows score a pair by its target's TRRUST in-degree (no "
    "layer). ChIP-seq rows: positives, DoRothEA ChIP-seq pairs (level D) absent from TRRUST and from DoRothEA levels A--C; candidates constructed as for TRRUST, "
    "excluding in either direction every pair documented in any regulatory resource at any level, the ChIP-seq positives "
    "included; about half of these candidates (49--55\\%) could not be matched on in-degree, so the fully matched rows keep only candidates matched on "
    "both deciles (and the positives of factors that have at least one; CSV column \\texttt{chip\\_candidates}). $^{a}$Positives in strata that hold both classes. $^{b}$$|r|$ AUROC, and $\\Delta$, computed with exactly tied "
    "$|r|$ values ordered by floating-point noise rather than counted as one half (effect on the $|r|$ AUROC, and hence on "
    "$\\Delta$, $\\leq 0.00012$ where checked). Contamination: AUROC and UB corrected for a fraction $\\pi$ "
    "of undocumented edges among the candidates, $0.5 + (A - 0.5)/(1 - \\pi)$; seventh column, layers whose corrected UB is "
    "below 0.60. $^{c}$Mean of the per-factor AUROCs with a factor-bootstrap CI; last column, fraction of factors above "
    "0.5 (range over layers). $^{d}$Covariates: the pair's log $|r|$ and the target's log mean expression, log "
    "$n_{\\mathrm{obs}}$, STRING physical degree and GO annotation count (factor-level covariates cancel within the factor "
    "stratum); two-way bootstrap with refitting, $B = 500$. Residualised: AUROC of the cosine after OLS residualisation on "
    "the covariates. Conditional logit: odds ratio (OR) per pooled SD of the cosine, adjusted (unadjusted in the eighth "
    "column); layer mean, mean of the twelve coefficients. $^{e}$Effective number of factors under each statistic's "
    "weights ($n_{\\mathrm{pos}} n_{\\mathrm{cand}}$ per factor; Fisher information) and the share of the five largest. "
    "\\textbf{C}: post hoc subsets, no multiplicity adjustment; $\\Delta_{\\mathrm{3k}}$, difference from $|r|$ on 3{,}000 "
    "cells. PAX5$\\to$CD79A (illustration): detected, fraction of the extraction's cells expressing the gene; rank 1 is the "
    "highest cosine among PAX5's documented and candidate pairs ($n$); candidates below, share of PAX5's candidates scored "
    "lower. "
    "\\textbf{D}: each split holds out a random fifth of the 245 kidney~A TRRUST pairs (49 positives, three candidate "
    "non-edges each); nothing is trained, so the spread over splits is the sampling distribution of the test statistic, not a CI. "
    "Uncentred, cosine of raw vectors, as the screen scored pairs; zero-row SVD, the reference decomposition of the full "
    "kidney~A matrix including the unobserved all-zero gene rows (Table~S14). Reference split: one fixed 80/20 split of these pairs; its value (first row with its i.i.d.\\ CI) and its "
    "percentile among the 500 splits. MDE as in A, for the pooled AUROC (positives against all held-out non-edges) and for the within-factor AUROC; size, rejection rate at a population AUROC "
    "of 0.5 (nominal 0.05); perm.\\ SD, Wald test with the permutation SD; oracle SD, the simulated sampling SD. "
    "\\textbf{E}: regulator importance, normalised per target, of each factor for each target: GENIE3, the extra-trees "
    "variant (300 trees, scikit-learn); GRNBoost2, GRNBoost2-style gradient boosting (scikit-learn, early stopping on a "
    "tenth of the cells); T1, five seeds on the extraction cells; T2, two seeds on 3{,}000 cells; AUROC "
    "of the seed-averaged score, and mean (SD) of the per-seed AUROCs; $\\Delta$, the method's AUROC minus that of $|r|$ on "
    "the same cells or of the cosine at the stated layer, on the same bootstrap draws. Cosine $L_0$ (first column): scGPT on "
    "the same pairs. "
    "\\textbf{F}: L2-penalised logistic probes on features of the factor's and the target's vectors $x_i$, $x_j$ "
    "(bilinear: outer product of 16 principal components fitted on training genes; $[x_j, \\cdot]$: target block added); "
    "AUROC, mean over 25 test folds; CI, corrected resampled $t$ (Nadeau--Bengio) with the test-to-training ratio counted in "
    "pairs (about 0.0625) or in genes (0.25, the cautious reading). $-$ target-only: minus, fold by fold, the better of the "
    "target-only and the factor-only probe (the latter scores 0.5 within factor). Rewired labels: degree-preserving swaps of targets "
    "between factors' positives, 100 draws, each retrained on one repetition's folds; $p = (b + 1)/(N + 1)$, floor 0.0099; "
    "max-$T$ over layers. Difference in differences: bilinear minus concatenation AUROC on the true labels against the same "
    "difference on each rewired draw. Input gene-token embedding: scGPT whole-human checkpoint; ESM2-15B: 5{,}120-dimensional "
    "per-protein embedding (UCE); both without layers. Non-linear: histogram gradient boosting and a two-layer perceptron, "
    "one repetition; rows, probe $\\times$ layer $\\times$ feature rows whose interval (pair ratio) lies above or below the "
    "model's own target-only control. "
    "\\textbf{G}: fractions among untied pairs; sign test, exact two-sided; regulator-cluster CI, bootstrap over source "
    "factors; probe CIs, binomial (pairs independent) or gene-cluster bootstrap ($B = 2{,}000$); strict subset, reverse "
    "absent from every regulatory set; out-/in-degree, the same fraction for TRRUST degrees (no model). "
    "Full values: \\texttt{Table\\_S7\\_edges.csv}."
).replace("MC_TXT", f"{mc3['replay_seed']['delta_ub95']:.4f} at $L_{{{mc3['layer']}}}$; two independent bootstrap seeds give "
          f"{mc3['indep_seed_1']['delta_ub95']:.4f} and {mc3['indep_seed_2']['delta_ub95']:.4f}, so C fails by less than the "
          f"Monte-Carlo error.").replace("REPL_SIGN", str(int(repl.same_sign.sum()))).replace("REPL_EXCL", str(int(repl.lung_ext_ci_excludes_05.sum())))
fm = summ.frac_neg_fully_matched
caption = caption.replace("MATCH_RANGE", f"{100 * fm.min():.0f}--{100 * fm.max():.0f}\\%")

ncols = NC
tex = [
    f"\\setcounter{{table}}{{6}}",
    "\\begingroup",
    "% The caption is longer than one line of 16,383pt, the limit of longtable's centring test; within this group the",
    "% caption is always set as a paragraph (the branch longtable itself takes for any caption wider than the line).",
    "\\makeatletter",
    "\\def\\LT@makecaption#1#2#3{\\LT@mcol\\LT@cols c{\\hbox to\\z@{\\hss\\parbox[t]\\LTcapwidth{\\reset@font",
    "  #1{#2: }#3\\endgraf\\vskip\\baselineskip}\\hss}}}",
    "\\makeatother",
    "\\footnotesize",
    f"\\setlength{{\\tabcolsep}}{{{SEP_PT}pt}}",
    "\\renewcommand{\\arraystretch}{1.05}",
    f"\\begin{{longtable}}{{{COLSPEC}}}",
    f"\\caption{{{caption}}}\\label{{tab:S7}}\\\\",
    "\\toprule",
    "\\endfirsthead",
    # zero-width boxes: the head and foot are set before any column width is known, and a natural-width span would
    # widen the last column on a first run (longtable would then need a third run to settle)
    f"\\multicolumn{{{ncols}}}{{@{{}}l}}{{\\hbox to 0pt{{\\footnotesize Table~\\thetable\\ (continued)\\hss}}}}\\\\",
    "\\toprule",
    "\\endhead",
    "\\midrule",
    f"\\multicolumn{{{ncols}}}{{r@{{}}}}{{\\hbox to 0pt{{\\hss\\footnotesize continued on next page}}}}\\\\",
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot",
] + rows + ["\\end{longtable}", "\\endgroup", ""]
_ttxt = "\n".join(tex[1:]).lower()
for bad in ("submitted", "revision", "review", "verifier", "implementer", "run_trees", "convention3", "legacy", "esm2-3b"):
    assert bad not in _ttxt, bad
OUT_TEX.write_text("\n".join(tex))

# ----------------------------------------------------------------------------------------------- CSV
out = pd.concat(csv_parts, ignore_index=True, sort=False)
if "layer" in out.columns:
    out["layer"] = out["layer"].astype(str).replace({"nan": ""})
# neutral labels for the public CSV (values unchanged); the internal names of the analysis-convention document, of the code
# audit and of the earlier analysis code are dropped
# fragment, and the internal names of the analysis-convention document and of the code audit are dropped
NEUTRAL = [("[submitted: ", "["), ("stratified_full512_submitted", "stratified_full512_uncentred"),
           ("(3,000 cells, run_trees definition)", "(3,000 cells, reference-split computation)"),
           ("kidney_A legacy (run_trees.py split, 49/147)", "kidney_A 80/20 reference split (49/147)"),
           ("pooled AUROC (as submitted)", "pooled AUROC (positives against all held-out non-edges)"),
           ("pooled AUROC, all-V non-edges (legacy)", "pooled AUROC, all-V non-edges"),
           ("h3_legacy_pairs (all-V non-edges)", "h3_allV_pairs (all-V non-edges)"),
           ("original float-tie-broken |co-expression| (not recomputed)", "ties ordered by floating-point noise"),
           (" [CONVENTION3 detection]", " [input-based detection]"),
           ("verifier (CONVENTION3 detection definition)", "sensitivity analysis (detection defined by scGPT input n_obs)")]
DELTA_BASIS = {"tie-corrected": "ties counted as 1/2"}
for col_ in out.columns:
    if out[col_].dtype == object and col_ != "source_file":
        def _neutral(x):
            if not isinstance(x, str):
                return x
            for a, b in NEUTRAL:
                x = x.replace(a, b)
            return x
        out[col_] = out[col_].map(_neutral)
if "delta_cx_basis" in out.columns:
    out["delta_cx_basis"] = out["delta_cx_basis"].replace(DELTA_BASIS)
out = out.rename(columns={"pct_delta_submitted": "pct_delta_reference_split"})
# undocumented pairs are candidate non-edges, not negatives: column names follow the caption
CAND_COLS = {"n_neg": "n_cand", "frac_neg_fully_matched": "frac_cand_fully_matched", "frac_neg": "frac_cand",
             "frac_negatives_below": "frac_candidates_below", "frac_tied_pos_neg": "frac_tied_pos_cand",
             "mean_n_neg": "mean_n_cand", "n_candidates": "n_pairs_of_factor",
             "rank_among_all_candidates": "rank_among_factor_pairs"}
assert not set(CAND_COLS.values()) & set(out.columns)
out = out.rename(columns=CAND_COLS)
assert not [c_ for c_ in out.columns if "neg" in c_.lower()], [c_ for c_ in out.columns if "neg" in c_.lower()]
_txt = " ".join(out.drop(columns=["source_file"]).astype(str).values.ravel()).lower() + " " + " ".join(out.columns).lower()
for bad in ("submitted", "revision", "review", "verifier", "implementer", "run_trees", "convention3", "round-", "preprint",
            "legacy", "not recomputed", "tie-corrected", "screen split", "screen definition", "screen's statistic"):
    assert bad not in _txt, bad
assert not out.source_file.astype(str).str.contains("legacy", case=False).any()
out.to_csv(OUT_CSV, index=False)
print(f"wrote {OUT_TEX} ({len(rows)} rows) and {OUT_CSV} ({len(out)} rows, {out.shape[1]} columns)")
print(f"E: delta-vs-cosine intervals excluding 0: {n_cos_ex0}/{n_cos_tot}; delta-vs-|r| excluding 0: {n_cx_ex0}/20")
print(f"G: GRN sign tests p<0.05: {n_p_lt}/{n_p}; probe gene-cluster CIs excluding 0.5: {n_cl}/36")
