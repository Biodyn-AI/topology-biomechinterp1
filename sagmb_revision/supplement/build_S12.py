"""Build Table S12 (the 183 screened hypothesis records) from revision3_experiments/E9_screen.

Inputs (read only):
  E9_screen/S2_hypotheses.csv        one row per JSON hypothesis record (409; 183 flagged in_183)
  E9_screen/S1_family_table.csv      per-JSON-family counts (directions involved)
  E9_screen/audit_409_records.csv    p-value extraction details per record
  E9_screen/audit_409.json           headline-record keys; all-record audit
  E9_screen/results.json             re-derived multiplicity counts (checked against)
Outputs:
  supplement/S12.tex
  supplementary_tables/Table_S12_screen_records.csv
Only counts and display rounding are computed here; every p-value and pass/fail flag is copied.
"""
from __future__ import annotations
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
E9 = ROOT / "revision3_experiments" / "E9_screen"
OUT_TEX = ROOT / "sagmb_revision" / "supplement" / "S12.tex"
OUT_CSV = ROOT / "sagmb_revision" / "supplementary_tables" / "Table_S12_screen_records.csv"
CITE = "Kendiukhov (2026b)"
THR = 0.05 / 183

FAM = {  # JSON family -> display code
    "module_structure": "MOD", "manifold_distance": "MAN", "intrinsic_dimensionality": "IDM",
    "graph_topology": "GRT", "null_sensitivity": "NUL", "topology_stability": "TST",
    "cross_model_alignment": "CMA", "persistent_homology": "PHO",
}

# ----------------------------------------------------------------------------- text -> LaTeX
_UNI = {"–": "--", "—": "---", "→": "->", "×": "\x01", "≥": ">=", "≤": "<=",
        "−": "-", "…": "...", "′": "'", "²": "^2"}


def tex(s) -> str:
    """Escape plain text for LaTeX and typeset recurring notation (SVk, L_l, ranges, arrows, sci. numbers)."""
    s = str(s)
    for k, v in _UNI.items():
        s = s.replace(k, v)
    store: list[str] = []

    def keep(latex: str) -> str:
        store.append(latex)
        return f"\x00{len(store) - 1}\x00"

    def sci(m):
        mant, sign, ex = m.group(1), m.group(2), m.group(3).lstrip("0") or "0"
        return keep(f"${mant}{{\\times}}10^{{{'-' if sign == '-' else ''}{ex}}}$")

    s = s.replace("\x01", keep(r"$\times$"))
    s = re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)e([-+]?)(\d+)\b", sci, s)
    s = s.replace("->", keep(r"$\to$")).replace(">=", keep(r"$\geq$")).replace("<=", keep(r"$\leq$"))
    s = s.replace("+-", keep(r"$\pm$"))
    s = re.sub(r"\bSV\(K,K\+2\)", lambda m: keep(r"$\SV{k}$--$\SV{k+2}$"), s)
    s = re.sub(r"\bSV(\d+)\s*-{1,2}\s*SV(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$--$\\SV{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bSV(\d+)-{1,2}(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$--$\\SV{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bSV(57|24)\b", lambda m: keep(f"SV{m.group(1)}"), s)          # record labels for bands, kept verbatim
    s = re.sub(r"\bSV(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$"), s)
    s = re.sub(r"\bL(\d{1,2})\s*-\s*L?(\d{1,2})\b", lambda m: keep(f"$L_{{{m.group(1)}}}$--$L_{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bL(\d{1,2})\b", lambda m: keep(f"$L_{{{m.group(1)}}}$"), s)
    s = re.sub(r"(?<=\d)-(?=\d)", "--", s)
    s = re.sub(r"(?<![\w)\]\-])-(?=\d)", lambda m: keep("$-$"), s)
    s = re.sub(r"(?<=\d)x\b", lambda m: keep(r"$\times$"), s)
    s = re.sub(r"(?<=\d) x (?=\d)", lambda m: keep(r"$\times$"), s)
    s = re.sub(r"\brho\b", lambda m: keep(r"$\rho$"), s)
    s = re.sub(r"(\d)\s*deg\b", lambda m: m.group(1) + keep(r"$^\circ$"), s)
    s = s.replace("~", keep(r"$\sim$")).replace("<", keep("$<$")).replace(">", keep("$>$"))
    s = re.sub(r"(?<![\w.,;:!?)\]])'", "`", s)
    s = s.replace("\\", r"\textbackslash{}")
    for ch in "&%$#_{}":
        s = s.replace(ch, "\\" + ch)
    s = s.replace("^", r"\^{}").replace("|", r"\textbar{}")
    s = re.sub("\x00(\\d+)\x00", lambda m: store[int(m.group(1))], s)
    return s


def fmt_p(p) -> str:
    """Display rounding of a p-value (3 significant digits)."""
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "---"
    p = float(p)
    if p == 0:
        return "0"
    if p >= 1e-3:
        return f"{p:.3g}"
    m, e = f"{p:.2e}".split("e")
    m = m.rstrip("0").rstrip(".")
    return f"${m}{{\\times}}10^{{{int(e)}}}$"


def matrix_display(m: str) -> str:
    parts = m.split(";")
    if m == "kidney_A;kidney_A_seed43;kidney_A_seed44":
        return "kidney~A (3 draws)"
    if m == "immune;immune_seed43;immune_seed44":
        return "immune (3 draws)"
    names = {"kidney_A": "kidney~A", "immune": "immune", "lung_attention": "lung attention",
             "random_gaussian": "Gaussian null", "geneformer_static": "Geneformer"}
    return " + ".join(names.get(x, tex(x)) for x in parts)


def short_title(name: str, n: int = 55) -> str:
    name = " ".join(str(name).split())
    if len(name) <= n:
        return tex(name)
    cut = name[:n].rsplit(" ", 1)[0].rstrip(" ,;:-—–(")
    return tex(cut) + r"\ldots{}"


def norm_dir(t: str) -> str:
    """'SV2-SV4' and 'SV2-4' are two spellings of one band in the records; count and show them once."""
    return re.sub(r"^SV(\d+)-SV(\d+)$", r"SV\1-\2", t)


def directions_display(v: str) -> str:
    if v == "none":
        return "---"
    toks = list(dict.fromkeys(norm_dir(t) for t in v.split(";") if t))
    return ", ".join(tex(t) for t in toks)


def as_bool(x) -> bool:
    return str(x) == "True"


# ----------------------------------------------------------------------------- data
S2 = pd.read_csv(E9 / "S2_hypotheses.csv", keep_default_na=False, na_values=[""])   # keep the literal 'n/a' entries
S1 = pd.read_csv(E9 / "S1_family_table.csv")
AU = pd.read_csv(E9 / "audit_409_records.csv")
A = json.loads((E9 / "audit_409.json").read_text())
RES = json.loads((E9 / "results.json").read_text())

assert len(S2) == 409 and len(AU) == 409
assert (S2.iteration.values == AU.iteration.values).all() and (S2.record_index.values == AU.record_index.values).all()
headline = {(h["iteration"], h["hypothesis_id"]): h["headline"] for h in A["headline_findings_all_records"]}

d = S2.copy()
d["p_recorded"] = pd.to_numeric(d.p_round1, errors="coerce")
d["p_floored"] = pd.to_numeric(d.primary_p_round2, errors="coerce")
d183 = d[d.in_183].copy()
d183["bonf_recorded"] = d183.p_recorded <= THR          # recorded p, zeros kept, '<=' (as in the audit)
d183["bonf_floored"] = d183.bonferroni_pass_round2.map(as_bool)
d183["bh"] = d183.bh_pass_round2.map(as_bool)

# checks against the audit's own counts
R = RES["audit_183_rederived"]
assert int(d183.p_recorded.notna().sum()) == R["n_extractable"] == 112
assert int((d183.p_recorded == 0).sum()) == R["n_p_equal_zero"] == 18
assert int(d183.bonf_recorded.sum()) == R["round1_bonferroni_pass_0p05_over_183"] == 43
assert int(d183.bonf_floored.sum()) == R["round2_bonferroni_pass"] == 25
assert int(d183.bh.sum()) == R["round2_bh_pass"] == 70
assert int(d183.spectral.sum()) == 108
assert int(d183.matrix_loaded.str.contains("kidney_A").sum()) == 115
assert int(d183.matrix_loaded.str.contains("immune").sum()) == 67

# starred (headline) records: their R/F marks, and how the headline-level count of the audit differs (re-measured p for some claims)
star = d183[[(a, b) in headline for a, b in zip(d183.iteration, d183.hypothesis_id)]]
assert len(star) == 11
N_STAR_R, N_STAR_F = int(star.bonf_recorded.sum()), int(star.bonf_floored.sum())
SENS = R["headline_round2_sensitivity"]["per_headline"]
HL_R2 = R["headline_bonferroni_round2"]
assert N_STAR_R == R["headline_bonferroni_round1"] == 8 and HL_R2 == sum(v["survives_round2"] for v in SENS.values()) == 5
# the audit scores a headline by its record's p (a zero floored at 1/1001) unless it substitutes a re-measured p
rec_p = {headline[(a, b)]: (1 / 1001 if p == 0 else p) for a, b, p in zip(star.iteration, star.hypothesis_id, star.p_recorded)}
REMEAS = [k for k, v in SENS.items() if not math.isclose(v["p_round2"], rec_p[k], rel_tol=1e-9)]
N_REMEAS, N_REMEAS_PASS = len(REMEAS), int(sum(SENS[k]["survives_round2"] for k in REMEAS))
assert (N_STAR_F, N_REMEAS, N_REMEAS_PASS) == (3, 5, 2) and N_STAR_F + N_REMEAS_PASS == HL_R2

# ----------------------------------------------------------------------------- Panel A
s1 = S1.set_index("json_family")
famrows = []
for fam, g in d183.groupby("json_family"):
    kid = g.matrix_loaded.str.contains("kidney_A")
    imm = g.matrix_loaded.str.contains("immune")
    row = dict(code=FAM[fam], fam=fam, n=len(g), spectral=int(g.spectral.sum()), ext=int(g.p_recorded.notna().sum()),
               zero=int((g.p_recorded == 0).sum()), b1=int(g.bonf_recorded.sum()), b2=int(g.bonf_floored.sum()),
               bh=int(g.bh.sum()), kid=int(kid.sum()), imm=int(imm.sum()), oth=int((~kid & ~imm).sum()),
               mis=int(g.gene_row_mapping.str.startswith("MISALIGNED").sum()))
    # copied counts must agree with S1_family_table.csv
    r1 = s1.loc[fam]
    assert row["n"] == r1.n_183 == r1.n_all and row["ext"] == r1.n_extractable_p_183
    assert row["b2"] == r1.n_bonferroni_pass_round2_183 and row["bh"] == r1.n_bh_pass_round2_183
    assert row["kid"] == r1.n_183_on_kidney_A and row["imm"] == r1.n_183_on_immune and row["mis"] == r1.n_183_misaligned_labels
    # directions and bands: records naming each label, with the two spellings of a band merged at record level
    # (S1_family_table.csv lists 'SV2-4' and 'SV2-SV4' separately); unmerged labels must agree with S1
    sv = Counter()
    for vv in g[g.spectral].sv_indices_recorded:
        for t in dict.fromkeys(norm_dir(t) for t in vv.split(";") if t and t != "none"):
            sv[t] += 1
    s1toks = [re.match(r"^(.*)\((\d+)\)$", t).groups() for t in str(r1.sv_directions_bands).split(";")] \
        if isinstance(r1.sv_directions_bands, str) else []
    merged = {norm_dir(lab) for lab, _ in s1toks if norm_dir(lab) != lab}
    for lab, n in s1toks:
        if lab not in merged and norm_dir(lab) == lab:
            assert sv[lab] == int(n), (fam, lab, sv[lab], n)
    top = sv.most_common()
    row["dirs"] = ", ".join(f"{tex(k)} ({n})" for k, n in top[:5]) + (", \\ldots{}" if len(top) > 5 else "")
    famrows.append(row)
famrows.sort(key=lambda r: -r["n"])
tot = {k: sum(r[k] for r in famrows) for k in ("n", "spectral", "ext", "zero", "b1", "b2", "bh", "kid", "imm", "oth", "mis")}
AL, AX = RES["audit_all"], RES["audit_all_excl_aborted"]

panelA = []
panelA.append(r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{3.3cm}rrrrrrrrrr>{\raggedright\arraybackslash}p{4.3cm}@{}}")
panelA.append(r"\toprule")
panelA.append(r" & & & \multicolumn{2}{c}{$p$-values} & \multicolumn{3}{c}{passes} & \multicolumn{3}{c}{matrix used} & \\")
panelA.append(r"\cmidrule(lr){4-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}")
panelA.append(r"JSON family (code) & $n$ & spec. & extr. & $=0$ & Bonf.\ rec. & Bonf.\ fl. & BH & kid.~A (mis.) & imm. & other & directions and bands involved (records) \\")
panelA.append(r"\midrule")
for r in famrows:
    panelA.append(f"\\texttt{{{tex(r['fam'])}}} ({r['code']}) & {r['n']} & {r['spectral']} & {r['ext']} & {r['zero']} & {r['b1']} & "
                  f"{r['b2']} & {r['bh']} & {r['kid']} ({r['mis']}) & {r['imm']} & {r['oth']} & {r['dirs']} \\\\")
panelA.append(r"\midrule")
panelA.append(f"All 183 records (iterations 1--63) & {tot['n']} & {tot['spectral']} & {tot['ext']} & {tot['zero']} & {tot['b1']} & "
              f"{tot['b2']} & {tot['bh']} & {tot['kid']} ({tot['mis']}) & {tot['imm']} & {tot['oth']} & Bonferroni threshold $0.05/183 = "
              f"{fmt_p(THR).strip('$')}$ \\\\")
panelA.append(f"All {AL['n_total']} records (iterations 1--138) & {AL['n_total']} & --- & {AL['n_extractable']} & {AL['n_p_equal_zero']} & --- & "
              f"{AL['n_bonferroni_pass']} & {AL['n_bh_pass']} & --- & --- & --- & Bonferroni at $0.05/{AL['n_extractable']}$ \\\\")
panelA.append(f"\\quad excl.\\ aborted duplicate run & {AX['n_total']} & --- & {AX['n_extractable']} & {AX['n_p_equal_zero']} & --- & "
              f"{AX['n_bonferroni_pass']} & {AX['n_bh_pass']} & --- & --- & --- & Bonferroni at $0.05/{AX['n_extractable']}$ \\\\")
panelA.append(r"\bottomrule")
panelA.append(r"\end{tabular}")

# ----------------------------------------------------------------------------- Panel B rows
rowsB = []
for r in d183.itertuples():
    it = int(r.iteration[5:9])
    hid = tex(r.hypothesis_id).replace(r"\_", r"\_\allowbreak{}")
    if (r.iteration, r.hypothesis_id) in headline:
        hid += r"$^{\star}$"
    mat = matrix_display(r.matrix_loaded)
    if str(r.gene_row_mapping).startswith("MISALIGNED"):
        mat += r"$^{\mathrm{m}}$"
    p1 = fmt_p(r.p_recorded)
    p2 = fmt_p(r.p_floored)
    marks = " ".join(m for m, ok in (("R", r.bonf_recorded), ("F", r.bonf_floored), ("H", r.bh)) if ok)
    rowsB.append(f"{it} & {hid} & {short_title(r.name)} & {FAM[r.json_family]} & {mat} & {directions_display(r.sv_indices_recorded)} & "
                 f"{p1} & {p2} & {marks} & {tex(r.decision)} \\\\")

# ----------------------------------------------------------------------------- LaTeX
W = dict(it=0.42, id=1.3, ti=4.25, fa=0.7, ma=1.55, di=1.95, p1=1.45, p2=1.45, pa=0.8, de=1.5)
col = lambda w, al="raggedright": f">{{\\scriptsize\\{al}\\arraybackslash}}p{{{w}cm}}"
spec = "@{}" + "".join([col(W["it"], "raggedleft"), col(W["id"]), col(W["ti"]), col(W["fa"]), col(W["ma"]), col(W["di"]),
                        col(W["p1"], "raggedleft"), col(W["p2"], "raggedleft"), col(W["pa"]), col(W["de"])]) + "@{}"
hdrB = (r"It. & ID & Hypothesis (short title) & Fam. & Matrix & Directions recorded & $p$ recorded & $p$ floored & Passes & Decision \\")

caption = (
    r"\textbf{The 183 hypothesis records of the automated screen (iterations 1--63).} "
    r"Each record is one hypothesis of the two-agent screen as archived in the iteration's JSON file, parsed together with the "
    r"iteration's code. The first 63 iterations produced these 183 records, whose headline findings " + CITE + r" (arXiv:2602.22247) reported. "
    r"\emph{Panel A}, per JSON family (the 8 families the records use; " + CITE + r" regroups them into 13): $n$, records; spec., records "
    r"involving a spectral decomposition; extr., records whose result string yields a primary $p$-value; $=0$, of these, $p$ stored as exactly 0 "
    r"(permutation tests with $N = 50$--1,000 draws); Bonf.\ rec., records passing Bonferroni at $0.05/183$ with the recorded $p$ "
    r"($p \leq 0.05/183$, zeros kept); Bonf.\ fl., passing with the floored $p$ ($p < 0.05/183$), where a zero is replaced by its "
    r"resolution floor $1/(N+1)$, with $N$ from the record's method text where stated and 1,000 otherwise; BH, passing Benjamini--Hochberg at "
    r"$q = 0.05$ over the 112 extractable floored $p$-values; kid.~A, records using the kidney~A matrix (a Tabula Sapiens kidney extraction, "
    r"4,803 genes, described as immune in " + CITE + r"), with the number whose gene labels were misaligned in parentheses (mis.; "
    r"iterations 27 and 29--39, where blank lines of the gene list were dropped so that all 209 labels sat on rows 0--208); imm., records "
    r"using the immune extraction (one record uses both); other, neither (one Gaussian null matrix, one Geneformer token-embedding matrix); "
    r"directions and bands, the five most frequent spectral directions or bands named in the family's spectral records (SV24 and SV57 are "
    r"the records' own labels for $\SV{2}$--$\SV{4}$ and $\SV{5}$--$\SV{7}$; PC, principal component). The last two rows repeat the audit over "
    r"all 409 records of the loop (137 JSON files; 3 records belong to an aborted duplicate run), with the Bonferroni denominator equal to the "
    r"number of extractable $p$-values; this is a sensitivity analysis. "
    r"\emph{Panel B}, every record in iteration order: It., iteration; ID, hypothesis identifier within the iteration ($^{\star}$, record keyed to "
    r"one of the 11 headline findings in the multiplicity audit); Fam., JSON family code (Panel A); Matrix, embedding matrix loaded (3 draws, the "
    r"main extraction and two further 256-cell draws; lung attention, scGPT attention scores from a lung extraction; $^{\mathrm{m}}$, misaligned "
    r"gene labels); Directions recorded, spectral directions or bands "
    r"named in the record (---, none); $p$ recorded, primary $p$ as extracted from the result string (---, none extractable); $p$ floored, after "
    r"replacing zeros by $1/(N+1)$; Passes: R, Bonferroni at $0.05/183$ with the recorded $p$; F, Bonferroni at $0.05/183$ with the "
    r"floored $p$; H, Benjamini--Hochberg at $q = 0.05$ with the floored $p$; Decision, the screen's own decision for the record. "
    r"Of the 11 starred records, " + str(N_STAR_R) + r" carry R and " + str(N_STAR_F) + r" carry F; the headline-level count after flooring ("
    + str(HL_R2) + r" of the 11) scores " + str(N_REMEAS) + r" of the headline findings by a $p$-value recomputed for the claim itself instead of the "
    r"record's primary $p$ shown here (for the transcription-factor classifier, a Stouffer combination of its 12 per-layer $p$-values, "
    r"which is invalid because the layers share genes), and " + str(N_REMEAS_PASS) +
    r" of these pass; replacing the $t$-approximation $p = 1.4\times10^{-24}$ of the STRING confidence gradient by its exact value, "
    r"0.0167, and setting aside the Stouffer combination leaves 3 (Table~S13, C23 and C28). "
    r"Titles are truncated to about 55 characters. No intervals are reported: $p$-values are as computed by the screen and are displayed to "
    r"three significant digits. Because each iteration's hypotheses were chosen after seeing earlier results, no correction over the recorded "
    r"family bounds the error rate of the screen as a whole. "
    r"The file also gives, for every record, its gene set and gene-row mapping, decomposition convention, reference sets, "
    r"script provenance and the claims of Table~S13 it supports, and the same fields for the 226 later records. "
    r"Full values: \texttt{Table\_S12\_screen\_records.csv}."
)

L = []
L.append(r"\setcounter{table}{11}")
L.append(r"\begingroup")
L.append(r"\footnotesize")
L.append(r"\setlength{\tabcolsep}{2.4pt}")
L.append(r"\setlength{\LTpre}{6pt}\setlength{\LTpost}{6pt}")
L.append(f"\\begin{{longtable}}{{{spec}}}")
L.append(f"\\caption{{{caption}}}\\label{{tab:S12}}\\\\")
L.append(r"\multicolumn{10}{@{}l}{\footnotesize\textbf{A.} Records per JSON family and multiple-testing passes}\\[2pt]")
L.append(r"\multicolumn{10}{@{}l@{}}{\scriptsize\setlength{\tabcolsep}{2.4pt}" + "\n".join(panelA) + r"}\\[4pt]")
L.append(r"\multicolumn{10}{@{}l}{\footnotesize\textbf{B.} The 183 records}\\")
L.append(r"\toprule")
L.append(hdrB)
L.append(r"\midrule")
L.append(r"\endfirsthead")
L.append(r"\multicolumn{10}{@{}l}{\footnotesize\textit{Table~\thetable{} (continued)}}\\")
L.append(r"\toprule")
L.append(hdrB)
L.append(r"\midrule")
L.append(r"\endhead")
L.append(r"\midrule")
L.append(r"\multicolumn{10}{@{}r@{}}{\scriptsize\textit{continued on next page}}\\")
L.append(r"\endfoot")
L.append(r"\bottomrule")
L.append(r"\endlastfoot")
L.extend(rowsB)
L.append(r"\end{longtable}")
L.append(r"\endgroup")
OUT_TEX.write_text("\n".join(L) + "\n")

# ----------------------------------------------------------------------------- CSV (all 409 records; 183 flagged)
c = S2.copy()
c.insert(0, "record_set", np.where(c.in_183, "183 screened records (iterations 1-63)",
                                   np.where(c.aborted_duplicate_run, "aborted duplicate run (iteration 64)",
                                            "later loop records (iterations 64-138)")))
c.insert(1, "in_table_S12_panel_B", c.in_183)
c["family_code"] = c.json_family.map(FAM).fillna("")
c["headline_finding"] = [headline.get((a, b), "") if i else "" for a, b, i in zip(c.iteration, c.hypothesis_id, c.in_183)]
c["matrix_display"] = [matrix_display(m).replace("~", " ") for m in c.matrix_loaded]
p1 = pd.to_numeric(c.p_round1, errors="coerce")
c["bonferroni_pass_183_recorded_p"] = np.where(c.in_183, (p1 <= THR), "n/a")
for col_ in ("p_source", "n_p_in_string", "n_perm_method", "n_perm_used", "resolution_limited", "bh_q", "bonferroni_alpha"):
    c[f"{col_}_all_records_audit"] = AU[col_].values
c = c.rename(columns={
    "p_round1": "p_recorded", "primary_p_round2": "p_floored_183", "n_perm_round2": "n_perm_for_floor_183",
    "resolution_limited_round2": "p_stored_as_zero_183", "bh_q_round2": "bh_q_183", "bh_pass_round2": "bh_pass_183",
    "bonferroni_pass_round2": "bonferroni_pass_183_floored_p", "n_perm_round2_regex": "n_perm_from_method_text",
    "p_this_audit": "p_floored_all_records", "preprint_claim": "ledger_claim_id",
    "preprint_family_heuristic": "family_as_regrouped_in_Kendiukhov2026b",
})
for col_ in c.columns:
    if c[col_].dtype == object:
        c[col_] = c[col_].map(lambda v: v.replace("the preprint's", f"{CITE}'s").replace("The preprint", CITE)
                              .replace("the preprint", CITE).replace("Preprint", CITE).replace("preprint", CITE)
                              if isinstance(v, str) else v)
        c[col_] = c[col_].map(lambda v: v.replace("n/a (not in 183)", "n/a (not among the 183)") if isinstance(v, str) else v)
c.to_csv(OUT_CSV, index=False)
print("S12:", len(rowsB), "records in Panel B;", len(c), "rows in CSV;", OUT_TEX)
