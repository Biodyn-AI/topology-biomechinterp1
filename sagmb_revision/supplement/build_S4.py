"""Build Supplementary Table S4 (kidney A bridge from the screen's co-pole statistic to H1) and its CSV.

Sources (read only): revision3_experiments/E2_ppi/bridge_kidneyA.csv, bridge_kidneyA_cells.csv
Writes sagmb_revision/supplement/S4.tex and sagmb_revision/supplementary_tables/Table_S4_bridge.csv.
Derived quantities (own computation from bridge_kidneyA_cells.csv): number of the 36 cells with max-T adjusted
p < 0.05, number at the permutation floor, and the observed / null-mean ratio of the tracked direction.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
E2 = ROOT / "revision3_experiments" / "E2_ppi"
OUT_TEX = ROOT / "sagmb_revision" / "supplement" / "S4.tex"
OUT_CSV = ROOT / "sagmb_revision" / "supplementary_tables" / "Table_S4_bridge.csv"


def f(x, d=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    s = f"{float(x):.{d}f}"
    if s.startswith("-"):
        s = s[1:] if float(s) == 0 else "$-$" + s[1:]
    return s


def cell(L, k):
    return f"$L_{{{int(L)}}}\\,\\SV{{{int(k)}}}$"


st = pd.read_csv(E2 / "bridge_kidneyA.csv")
cl = pd.read_csv(E2 / "bridge_kidneyA_cells.csv")
ORDER = ["S0_submitted", "S1_observed_SVD", "S1p_offpath_physical_on_187", "S2_V_universe", "S3_physical_STRING",
         "S4_assortativity", "S5_Maslov_Sneppen"]
assert list(st.step) == ORDER

# tracked direction: L5 'SV3' of the zero-row full SVD (a reference decomposition, not the screen's own 209-gene one)
# = L5 SV2 of the observed-gene SVD from S1 on
TRACK = {s: (5, 3) if s == "S0_submitted" else (5, 2) for s in ORDER}
LABEL = {"S0_submitted": "S0", "S1_observed_SVD": "S1", "S1p_offpath_physical_on_187": "S1$'$",
         "S2_V_universe": "S2", "S3_physical_STRING": "S3", "S4_assortativity": "S4", "S5_Maslov_Sneppen": "S5"}
CHANGE = {
    "S0_submitted": "The screen's ingredients (the 187 observed genes of its curated list; functional STRING pairs, score "
                    "$\\geq 0.7$, among them), with coordinates from the zero-row SVD",
    "S1_observed_SVD": "SVD of the observed-gene matrix ($V$, 1,972 genes); coordinates still restricted to the 187 genes",
    "S1p_offpath_physical_on_187": "Off-path check: S1 with physical STRING $\\geq 700$ pairs among the 187 genes",
    "S2_V_universe": "Universe $V$: poles and shuffle over all 1,972 observed genes; same 809 functional pairs",
    "S3_physical_STRING": "Physical STRING $\\geq 700$ edges within $V$ instead of functional pairs",
    "S4_assortativity": "Scalar assortativity $r$ instead of the co-pole rate",
    "S5_Maslov_Sneppen": "Maslov--Sneppen rewiring instead of the label shuffle (the primary H1 run of kidney~A)",
}
STAT = {s: "co-pole rate; label shuffle, one-sided" for s in ORDER[:5]}
STAT["S4_assortativity"] = "$r$; label shuffle, two-sided"
STAT["S5_Maslov_Sneppen"] = "$r$; MS rewiring, two-sided"

rows, csv_steps = [], []
for s in ORDER:
    t = st[st.step == s].iloc[0]
    c = cl[cl.step == s]
    assert len(c) == 36
    L, k = TRACK[s]
    tr = c[(c.layer == L) & (c.sv == k)].iloc[0]
    if s == "S0_submitted":
        assert abs(tr.stat - t.SV3_L5_stat) < 1e-12
    n_maxt = int((c.p_adj36 < 0.05).sum())
    n_floor = int(np.isclose(c.p_adj36, t.p_floor, rtol=0, atol=1e-12).sum())
    best = c.iloc[int(np.argmax(c.z.values if t["tail"] == "greater" else np.abs(c.z.values)))]
    assert int(best.layer) == int(t.best_layer) and int(best.sv) == int(t.best_sv)
    co = t["tail"] == "greater"
    bonf = f"{int(t.n_cells_p_lt_bonf183)}" if t.bonf183_resolvable else "n.r."
    ratio = tr.stat / tr.null_mean if co else np.nan
    rows.append(" & ".join([
        LABEL[s], CHANGE[s], f"{int(t.n_genes):,}/{int(t.n_pairs):,}", STAT[s], f"{int(t.N):,}",
        f"{int(t.n_cells_p_lt_05)}/{bonf}/{n_maxt}",
        cell(L, k), f(tr.stat), f(tr.null_mean), f(ratio, 2) if co else "--", f(tr.z, 1),
        cell(t.best_layer, t.best_sv), f(t.best_stat), f(t.best_z, 1)]) + " \\\\")
    d = t.to_dict()
    d.update(n_cells_p_adj36_lt_05=n_maxt, n_cells_p_adj36_at_floor=n_floor, tracked_layer=L, tracked_sv=k,
             tracked_stat=tr.stat, tracked_null_mean=tr.null_mean, tracked_null_sd=tr.null_sd, tracked_z=tr.z,
             tracked_p_raw=tr.p_raw, tracked_p_adj36=tr.p_adj36, tracked_ratio_obs_null=ratio)
    csv_steps.append(d)
    if s == "S1p_offpath_physical_on_187":
        rows[-1] = "\\itshape " + rows[-1].replace(" & ", " & \\itshape ")

s5 = [d for d in csv_steps if d["step"] == "S5_Maslov_Sneppen"][0]
s1 = st[st.step == "S1_observed_SVD"].iloc[0]
s1p = st[st.step == "S1p_offpath_physical_on_187"].iloc[0]

cap = (
    "\\textbf{Step-by-step bridge on kidney~A from the screen's interaction co-pole statistic to the H1 statistic.} "
    "Each step changes one factor of the previous step (cumulatively, except the off-path check S1$'$, shown in italics, which "
    "branches from S1). Every step evaluates 36 cells, $\\SV{2}$--$\\SV{4}$ at $L_0$--$L_{11}$ of that step's decomposition (at S0 "
    "the zero-row SVD, whose $\\SV{2}$--$\\SV{4}$ correspond to observed $\\SV{1}$--$\\SV{3}$; from S1 the observed-gene "
    "SVD), with one null draw applied jointly to all 36. "
    "Co-pole rate, fraction of pairs whose two genes fall in the same pole (the top or bottom $K = \\lfloor 0.27\\,n \\rfloor$ genes "
    "of the coordinate, $n$ the number of genes in the universe); label shuffle, permutation of coordinates (or pole labels) across "
    "the universe; MS, Maslov--Sneppen degree-preserving rewiring; $r$, scalar assortativity of the coordinate over the pairs. "
    "Zero-row SVD, SVD of the full 4,803-gene matrix including the all-zero rows of unobserved genes. "
    "Genes/pairs, size of the gene universe and number of gene pairs; $N$, null draws. "
    "Cells, number of the 36 with raw $p < 0.05$ / raw $p < 0.05/183$ (the screen's Bonferroni threshold over its 183 "
    "hypotheses; n.r., not resolvable because $1/(N+1) > 0.05/183$) / max-$T$ adjusted $p < 0.05$ over the 36 cells "
    "(descriptive). The minimum adjusted $p$ is at the permutation floor $1/(N+1)$ at every step except S1$'$ "
    f"({f(s1p.min_p_adj36)}); at S5, {s5['n_cells_p_adj36_at_floor']} of the 36 cells are at the floor. "
    "Tracked direction: $L_5\\,\\SV{3}$ of the zero-row SVD at S0, "
    "which is $L_5\\,\\SV{2}$ of the observed-gene SVD (the zero-row SVD's leading direction separates unobserved from observed "
    "genes and shifts the low indices by one; observed-gene $L_5\\,\\SV{3}$ at S1 has a co-pole rate of "
    f"{f(s1.SV3_L5_stat)}, $z = {s1.SV3_L5_z:.2f}$). "
    "Value, statistic of that direction; null, null mean; ratio, value / null mean (co-pole steps); "
    "$z = (\\text{value} - \\text{null mean})/\\text{null SD}$; best cell, the cell of the 36 with the largest $z$ "
    "(largest $|z|$ for two-sided steps). "
    "S0 does not use the screen's own decomposition, of a 209-gene subset that includes 22 zero rows and maps only weakly onto "
    "the observed-gene basis (Table~S14); the screen's values on that subset are in Table~S13 (claim C21). "
    "At S3 the co-pole null mean of the tracked direction under "
    "MS rewiring is about twice its label-shuffle value, so the label-shuffle ratio overstates the enrichment about twofold "
    "(Table~S3, panel E); S4 and S5 use assortativity, whose null mean is close to zero under both nulls. "
    "Full values: \\texttt{Table\\_S4\\_bridge.csv}."
)

lines = [
    "\\setcounter{table}{3}",
    "\\begin{table}[!htbp]",
    "\\centering\\footnotesize",
    "\\caption{" + cap + "}\\label{tab:S4}",
    "\\scriptsize\\setlength{\\tabcolsep}{2pt}",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}l>{\\raggedright\\arraybackslash}p{3.65cm}r"
    ">{\\raggedright\\arraybackslash}p{1.6cm}rclrrrrlrr@{}}",
    "\\toprule",
    "& & & & & & \\multicolumn{5}{c}{Tracked direction} & \\multicolumn{3}{c}{Best cell of 36} \\\\",
    "\\cmidrule(lr){7-11}\\cmidrule(l){12-14}",
    "Step & Change & Genes/pairs & Statistic; null & $N$ & Cells & Cell & Value & Null & Ratio & $z$ & Cell & Value & $z$ \\\\",
    "\\midrule",
] + rows + ["\\bottomrule", "\\end{tabular*}", "\\end{table}", ""]
OUT_TEX.write_text("\n".join(lines))

a = pd.DataFrame(csv_steps)
a.insert(0, "source_file", "E2_ppi/bridge_kidneyA.csv (+ derived from bridge_kidneyA_cells.csv)")
a.insert(0, "panel", "steps")
b = cl.copy()
b.insert(0, "source_file", "E2_ppi/bridge_kidneyA_cells.csv")
b.insert(0, "panel", "cells")
b.insert(2, "context", "kidney_A")
a.insert(2, "context", "kidney_A")
out = pd.concat([a, b], ignore_index=True, sort=False)
# neutral public labels for the screen step (S0) and its gene set
out["step"] = out["step"].replace({"S0_submitted": "S0_screen_ingredients"})
out["description"] = out["description"].str.replace("187 legacy genes", "the screen's 187 genes", regex=False)
assert not out.astype(str).apply(lambda s: s.str.contains("submitted|legacy", case=False)).any().any()
out.to_csv(OUT_CSV, index=False)
print("wrote", OUT_TEX, OUT_CSV, out.shape)
