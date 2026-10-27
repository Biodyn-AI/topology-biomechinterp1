"""Build Table S13 (claims ledger: the 61 numerical claims of the screen's report) from revision3_experiments/E9_screen.

Inputs (read only):
  E9_screen/E9_claims_ledger.csv     61 claims: reported value, provenance, recomputation, status, observed-gene contrast
  E9_screen/results.json             status counts (checked against)
  sagmb_revision/manuscript.tex      section / table numbers of the main text (for the pointers in the last column), counted
                                     from the source; compared with manuscript.aux and with the label order of drafts/*.tex
Outputs:
  supplement/S13.tex
  supplementary_tables/Table_S13_claims_ledger.csv
Every value is copied from the ledger; this script only classifies the ledger's own status text into five
categories, attaches a pointer to the main-text section that examines the same structure, and formats text.
"""
from __future__ import annotations
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
E9 = ROOT / "revision3_experiments" / "E9_screen"
AUX = ROOT / "sagmb_revision" / "manuscript.aux"
MAIN = ROOT / "sagmb_revision" / "manuscript.tex"
OUT_TEX = ROOT / "sagmb_revision" / "supplement" / "S13.tex"
OUT_CSV = ROOT / "sagmb_revision" / "supplementary_tables" / "Table_S13_claims_ledger.csv"
CITE = "Kendiukhov (2026b)"

# ----------------------------------------------------------------------------- text -> LaTeX
_UNI = {"\u2013": "--", "\u2014": "---", "\u2192": "->", "\u2265": ">=", "\u2264": "<=", "\u2212": "-",
        "\u2026": "...", "\u2032": "'", "\u00b2": "^2"}


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

    s = s.replace("\u00d7", keep(r"$\times$"))
    s = s.replace("589-270-141", keep("$589-270-141$")).replace("kidney_A", keep("kidney~A"))
    s = s.replace("66 - 1", keep("$66 - 1$")).replace("ER lumen", "endoplasmic-reticulum lumen")
    s = re.sub(r"\bDelta\b", lambda m: keep(r"$\Delta$"), s)
    s = re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)e([-+]?)(\d+)\b", sci, s)
    s = s.replace("->", keep(r"$\to$")).replace(">=", keep(r"$\geq$")).replace("<=", keep(r"$\leq$"))
    s = s.replace("+-", keep(r"$\pm$"))
    s = re.sub(r"\bSV(\d+)\s*-{1,2}\s*SV(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$--$\\SV{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bSV(\d+)-{1,2}(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$--$\\SV{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bSV(\d+)\b", lambda m: keep(f"$\\SV{{{m.group(1)}}}$"), s)
    s = re.sub(r"\bL(\d{1,2})\s*-\s*L?(\d{1,2})\b", lambda m: keep(f"$L_{{{m.group(1)}}}$--$L_{{{m.group(2)}}}$"), s)
    s = re.sub(r"\bL(\d{1,2})\b", lambda m: keep(f"$L_{{{m.group(1)}}}$"), s)
    s = re.sub(r"\br_rb\b", lambda m: keep(r"$r_{\mathrm{rb}}$"), s)
    s = re.sub(r"(?<=\d)-(?=\d)", "--", s)
    s = re.sub(r"(?<![\w)\]\-])-(?=\d)", lambda m: keep("$-$"), s)
    s = re.sub(r"(?<=\d)x\b", lambda m: keep(r"$\times$"), s)
    s = re.sub(r"(?<=\d) x (?=\d)", lambda m: keep(r"$\times$"), s)
    s = re.sub(r"\brho\b", lambda m: keep(r"$\rho$"), s)
    s = re.sub(r"(?<![\w.$])([pNnzr])(?=\s*(?:=|<|>)\s*[-\d.(\x00])", lambda m: keep(f"${m.group(1)}$"), s)
    s = re.sub(r"(\d)\s*deg\b", lambda m: m.group(1) + keep(r"$^\circ$"), s)
    s = s.replace("~", keep(r"$\sim$")).replace("<", keep("$<$")).replace(">", keep("$>$"))
    s = re.sub(r"(?<![\w.,;:!?)\]])'", "`", s)
    s = s.replace("\\", r"\textbackslash{}")
    for ch in "&%$#_{}":
        s = s.replace(ch, "\\" + ch)
    s = s.replace("^", r"\^{}").replace("|", r"\textbar{}")
    s = s.replace("/", r"/\allowbreak{}").replace(r"\_", r"\_\allowbreak{}")
    s = re.sub("\x00(\\d+)\x00", lambda m: store[int(m.group(1))], s)
    return s


def cite_neutral(v):
    """Refer to the screen's report by its citation."""
    if not isinstance(v, str):
        return v
    v = v.replace("the preprint's own", f"{CITE}'s own").replace("the preprint's", f"{CITE}'s")
    v = v.replace("The preprint", CITE).replace("the preprint", CITE).replace("Preprint", CITE).replace("preprint", CITE)
    return v


SHOUT = {"NOT": "not", "KIDNEY": "kidney", "ALL": "all", "SIGN REVERSED": "sign reversed", "DEcreases": "decreases",
         "INCLUDING": "including", "EMBEDDINGS": "embeddings", "FULL": "full", "ATTENTION": "attention", "LUNG": "lung",
         "PHYSICAL": "physical"}


def calm(v: str) -> str:
    for k, w in SHOUT.items():
        v = re.sub(rf"\b{k}\b", w, v)
    return v


def disp(v) -> str:
    if not isinstance(v, str) or v.strip() in ("", "nan", "n/a"):
        return "---"
    return tex(calm(cite_neutral(v)))


# ----------------------------------------------------------------------------- data
LED = pd.read_csv(E9 / "E9_claims_ledger.csv", keep_default_na=False, na_values=[""])   # keep the literal 'n/a' entries
RES = json.loads((E9 / "results.json").read_text())
assert len(LED) == 61 and LED.claim_id.tolist() == [f"C{i:02d}" for i in range(1, 62)]
assert LED.reproduces.value_counts().to_dict() == RES["ledger"]["reproduces_counts"]
# this paper's columns call undocumented pairs candidates, not negatives (wording only; the screen's reported value is
# quoted as reported)
for _col, _a, _b in [("gene_label_set", "1,523 neg pairs", "1,523 TRRUST-negative candidate pairs"),
                     ("recomputed_own_convention", "1,523 neg)", "1,523 candidates)"),
                     ("reproduces", "with 1,523 negatives", "with 1,523 candidate pairs")]:
    assert LED[_col].astype(str).str.count(re.escape(_a)).sum() == 1, (_col, _a)
    LED[_col] = LED[_col].map(lambda v: v.replace(_a, _b) if isinstance(v, str) else v)


def strip_outer_paren(q: str) -> str:
    q = q.strip()
    if q.startswith("("):
        depth = 0
        for i, ch in enumerate(q):
            depth += ch == "("
            depth -= ch == ")"
            if depth == 0:
                return (q[1:i] + q[i + 1:]).strip()
    return q


def status(s: str) -> tuple[str, str]:
    s = s.strip()
    if s.startswith("not recomputable"):
        return "not recomputable", "iteration folders overwritten by an earlier loop"
    if s.startswith("not recomputed"):
        return "provenance only", strip_outer_paren(s[len("not recomputed"):])
    if s.startswith("partly"):
        return "partly", strip_outer_paren(s[len("partly"):].lstrip(":").strip())
    if s.startswith("no"):
        return "does not reproduce", strip_outer_paren(s[2:].strip().lstrip("-").strip())
    if s.startswith("yes"):
        return "reproduces", strip_outer_paren(s[3:].strip())
    raise ValueError(s)


OWN_TABLES = [("STRING PPI co-pole by axis (Table 1, Table S4)", "STRING PPI co-pole by axis (Tables 1 and S4 of Kendiukhov 2026b)"),
              ("SV2 per-layer = Table S4 exactly", "SV2 per-layer = Table S4 of Kendiukhov (2026b) exactly")]
for col_ in ("claim", "recomputed_own_convention"):
    for a, b in OWN_TABLES:
        LED[col_] = LED[col_].str.replace(a, b, regex=False)
assert LED.claim.str.contains("Kendiukhov 2026b").sum() == 1 and LED.recomputed_own_convention.str.contains("Kendiukhov").sum() == 1
LED["status_category"], LED["status_detail"] = zip(*LED.reproduces.map(status))
cnt = LED.status_category.value_counts().to_dict()
assert cnt == {"reproduces": 25, "provenance only": 27, "does not reproduce": 5, "partly": 2, "not recomputable": 2}, cnt

# ----------------------------------------------------------------------------- main-text pointers
# Section and table numbers are counted from the current main-text source (manuscript.tex); manuscript.aux can lag behind an
# edit (e.g. a subsection added to Methods), so it is only compared with and reported, never trusted on its own.
def count_numbers(src: str) -> dict:
    src = re.sub(r"(?<!\\)%.*", "", src)
    nums, sec_n, sub_n, tab_n, cur, cur_tab = {}, 0, 0, 0, None, None
    for m in re.finditer(r"\\(section|subsection)(\*?)\{|\\begin\{table\}|\\label\{((?:sec|tab):[^}]+)\}", src):
        if m.group(1) == "section":
            if m.group(2):
                cur = None
            else:
                sec_n, sub_n = sec_n + 1, 0
                cur = str(sec_n)
        elif m.group(1) == "subsection":
            if m.group(2):
                cur = None
            else:
                sub_n += 1
                cur = f"{sec_n}.{sub_n}"
        elif m.group(0).startswith(r"\begin{table}"):
            tab_n += 1
            cur_tab = str(tab_n)
        else:
            lab = m.group(3)
            nums[lab] = cur_tab if lab.startswith("tab:") else cur
    return nums


NUM = count_numbers(MAIN.read_text())
aux = AUX.read_text() if AUX.exists() else ""
NUM_AUX = {m.group(1): m.group(2) for m in re.finditer(r"\\newlabel\{((?:sec|tab):[^}]+)\}\{\{([^}]*)\}", aux)}
for lab in sorted(set(NUM) | set(NUM_AUX)):
    if NUM.get(lab) != NUM_AUX.get(lab):
        print(f"note: {lab} is {NUM.get(lab)} in manuscript.tex but {NUM_AUX.get(lab)} in manuscript.aux (aux not yet refreshed)")
# the drafts carry the same labels in the same order
DRAFTS = ["intro", "results_a", "results_h1", "results_h2", "results_h3", "results_h4", "results_general", "results_screen",
          "discussion", "methods"]
draft_src = "\n".join((ROOT / "sagmb_revision" / "drafts" / f"{d}.tex").read_text() for d in DRAFTS)
draft_labels = re.findall(r"\\label\{((?:sec|tab):[^}]+)\}", re.sub(r"(?<!\\)%.*", "", draft_src))
main_labels = re.findall(r"\\label\{((?:sec|tab):[^}]+)\}", re.sub(r"(?<!\\)%.*", "", MAIN.read_text()))
assert draft_labels == main_labels, "main-text drafts and manuscript.tex differ in structure"


def sec(*labels):
    return "; ".join(("Sec.~" if l.startswith("sec:") else "Table~") + NUM[l] for l in labels)


# claim -> (main-text labels, extra pointer text); None = not analysed here
P = {
    "C01": (("sec:screen", "sec:m-screen"), "Table~S12"),
    "C02": (("sec:collapse", "tab:compression", "tab:ledger"), "Table~S1"),
    "C03": (("sec:collapse", "tab:compression", "tab:ledger"), ""),
    "C04": (("sec:collapse",), "all observed genes"),
    "C05": None,
    "C06": (("sec:collapse",), "participation ratio, all observed genes"),
    "C07": (("sec:collapse",), "feature-permutation control"),
    "C08": (("sec:compartment", "tab:ledger"), "Table~S2"),
    "C09": (("sec:compartment",), "Table~S2"),
    "C10": (("sec:compartment",), "max-$T$ over 288 tests"),
    "C11": (("sec:compartment",), "compartment gradient"),
    "C12": (("sec:compartment",), "first twelve directions tested"),
    "C13": (("sec:grn",), "within-factor edge benchmark"),
    "C14": None, "C15": None, "C16": None,
    "C17": (("sec:tf",), "single directions, Table~S6"),
    "C18": None,
    "C19": (("tab:ledger",), ""),
    "C20": None,
    "C21": (("sec:ppi", "tab:ledger"), "Tables~S3, S4"),
    "C22": (("sec:screen",), "Table~S14"),
    "C23": (("sec:ppi", "tab:ledger"), "STRING confidence bands"),
    "C24": None,
    "C25": (("sec:ppi",), "hub removal, degree-stratified null; Table~S3"),
    "C26": None,
    "C27": (("sec:grn",), "within-factor edge benchmark"),
    "C28": (("sec:tf", "tab:ledger"), "Table~S6"),
    "C29": (("sec:tf",), "all 78 contiguous bands"),
    "C30": (("sec:tf",), "all 78 contiguous bands"),
    "C31": (("sec:grn",), "co-expression comparator"),
    "C32": (("sec:grn", "sec:rank"), ""),
    "C33": (("sec:grn", "tab:ledger"), "Table~S7"),
    "C34": (("sec:grn",), "every layer, Table~S7"),
    "C35": (("sec:grn", "sec:rank"), ""),
    "C36": (("sec:contexts", "tab:contexts"), ""),
    "C37": None, "C38": None,
    "C39": (("tab:ledger",), ""),
    "C40": None, "C41": None, "C42": None, "C43": None,
    "C44": (("tab:ledger",), ""),
    "C45": None, "C46": None, "C47": None, "C48": None, "C49": None,
    "C50": (("sec:rank",), "static embeddings as comparators"),
    "C51": None, "C52": None, "C53": None, "C54": None,
    "C55": (("sec:crosscontext",), "variance across cell draws"),
    "C56": None, "C57": None, "C58": None,
    "C59": (("sec:contexts", "sec:screen"), "observed genes"),
    "C60": (("sec:grn", "tab:contexts"), "edge benchmark counts"),
    "C61": (("sec:collapse",), ""),
}
assert set(P) == set(LED.claim_id)


def pointer(cid: str) -> str:
    p = P[cid]
    if p is None:
        return ""
    labels, extra = p
    out = sec(*labels)
    return out + (f" ({extra})" if extra else "")


LED["this_paper_pointer"] = [pointer(c) for c in LED.claim_id]


def this_paper_cell(r) -> str:
    parts = []
    oc = r.observed_convention_contrast
    if isinstance(oc, str) and oc.strip():
        parts.append(r"\emph{Observed genes / true rows:} " + disp(oc))
    if r.this_paper_pointer:
        parts.append(r.this_paper_pointer)
    return "; ".join(parts) if parts else "---"


def source_cell(r) -> str:
    its = []
    for s in str(r.source_iterations).split(";"):
        if s.strip().startswith("iterations/iter_0001-0063"):
            its.append("all records, iterations 1--63")
            continue
        m = re.match(r"iter_(\d{4})/(.+)$", s.strip())
        its.append(f"{int(m.group(1))}/{tex(m.group(2))}" if m else tex(calm(s.strip())))
    mat = str(r.matrix_and_true_tissue)
    if mat in ("nan", "n/a"):
        data = ""
    elif mat.startswith("kidney_A;") or mat.startswith("kidney_A; "):
        data = "kidney~A (3 draws)"
    elif mat.startswith("kidney_A"):
        data = "kidney~A"
    elif mat.startswith("immune;") or mat.startswith("immune; "):
        data = "immune (3 draws)"
    elif mat.startswith("immune"):
        data = "immune"
    elif mat.startswith("historical"):
        data = "earlier loop"
    elif mat.startswith("Geneformer token"):
        data = "Geneformer + kidney~A"
    elif mat.startswith("Geneformer"):
        data = "Geneformer"
    elif mat.startswith("scGPT ATTENTION"):
        data = "lung attention, kidney~A labels"
    else:
        raise ValueError(mat)
    conv = str(r.svd_convention)
    lab = str(r.gene_label_set)
    if conv.startswith("subset-fitted SVD of the 209"):
        cv = "209-gene subset (22 zero rows)"
    elif conv.startswith("subset-fitted SVD of rows 0-208"):
        cv = "misaligned rows 0--208"
    elif conv.startswith("observed non-zero rows"):
        cv = "observed rows"
    elif conv.startswith("full matrix incl"):
        cv = "full matrix (zero rows)"
    elif conv.startswith("random 2,000-row") or conv.startswith("2,000-row"):
        cv = "2,000-row subsample (zero rows)"
    elif conv.startswith("zero rows included"):
        cv = "zero rows included"
    elif conv.startswith("not spectral"):
        cv = "not spectral" + ("; misaligned rows" if "misaligned" in lab else "")
    else:
        cv = ""
    tail = "; ".join(x for x in (data, cv) if x)
    return ", ".join(its) + (r"\newline " + tail if tail else "")


def status_cell(r) -> str:
    d = r.status_detail
    return r"\textit{" + r.status_category + "}" + (": " + disp(d) if d else "")


rows = []
for r in LED.itertuples():
    rows.append(f"{r.claim_id} & {disp(r.claim)} & {source_cell(r)} & {disp(r.preprint_value)} & "
                f"{disp(r.recomputed_own_convention) if r.status_category not in ('provenance only', 'not recomputable') else '---'} & "
                f"{status_cell(r)} & {this_paper_cell(r)} \\\\")

# ----------------------------------------------------------------------------- LaTeX
W = [0.55, 1.85, 1.9, 2.75, 3.5, 2.2, 2.95]
col = lambda w: f">{{\\scriptsize\\raggedright\\arraybackslash}}p{{{w}cm}}"
spec = "@{}" + "".join(col(w) for w in W) + "@{}"
hdr = (r"ID & Claim & Source records; data & As reported & Screen's procedure, recomputed & Status & "
       r"Observed-gene value; this paper \\")
caption = (
    r"\textbf{Claims ledger: the 61 numerical claims reported from the screen, traced and recomputed.} "
    r"One row per numerical claim of " + CITE + r" (abstract, results, methods and tables), numbered C01--C61. "
    r"ID, claim identifier (the file of Table~S12 lists, for each record, the claims it supports); Claim, what was claimed; Source records, iteration/hypothesis of "
    r"the record(s) behind the claim, with the matrix the record used (kidney~A, the Tabula Sapiens kidney extraction with 4,803 genes, "
    r"described as immune in " + CITE + r"; immune, the immune extraction; 3 draws, the main extraction and two further 256-cell draws; "
    r"earlier loop, iterations of an earlier run of the loop whose folders were overwritten) and its decomposition: 209-gene subset (22 "
    r"zero rows), decomposition fitted to the 209 named genes at their true rows, 22 of which are never-observed all-zero rows; misaligned rows "
    r"0--208, the same with every label on the wrong row (blank lines of the gene list dropped before numbering); observed rows, only genes "
    r"observed in at least one sampled cell, mean-centred over them (the convention of this paper); full matrix (zero rows), all "
    r"rows including the never-observed ones, whose first direction is then the observation indicator. As reported, the value in " + CITE + r". "
    r"Screen's procedure, recomputed, the value recomputed here with the screen's own code conventions on the matrix and gene rows the record "
    r"used (brackets in C02: minimum--maximum over 50 random subsamples). Status: \textit{reproduces}, exactly or within resampling error "
    r"(25 claims); \textit{partly} (2); \textit{does not reproduce} (5); \textit{provenance only}, traced to its record and code but not "
    r"recomputed (27); \textit{not recomputable}, records not preserved (2). Observed-gene value; this paper: the same quantity "
    r"recomputed on the observed genes or on the labels' true rows, where computed, and the section and table of the main text that "
    r"examine the corresponding structure under the convention of this paper (---, not examined here). "
    r"Abbreviations: OR, odds ratio; ER, effective rank; TwoNN, two-nearest-neighbour intrinsic dimension; "
    r"emp $p$, empirical permutation $p$; $r_{\mathrm{rb}}$ and rbc, rank-biserial correlation; OOV, out of vocabulary; GC, germinal centre; "
    r"FFL, feed-forward loop; CKA, centred kernel alignment; U, left singular vectors; S, singular values; TF, transcription factor; "
    r"PPI, protein--protein interaction; GO BP, Gene Ontology biological process; EV, extracellular vesicle; mito, mitochondrion; "
    r"DC, dendritic cell; CV, cross-validation; OOF, out-of-fold. "
    r"No new intervals are computed for this table. "
    r"Full values: \texttt{Table\_S13\_claims\_ledger.csv}."
)

L = [r"\setcounter{table}{12}", r"\begingroup", r"\footnotesize", r"\setlength{\tabcolsep}{2.4pt}",
     r"\setlength{\LTpre}{6pt}\setlength{\LTpost}{6pt}",
     f"\\begin{{longtable}}{{{spec}}}",
     f"\\caption{{{caption}}}\\label{{tab:S13}}\\\\",
     r"\toprule", hdr, r"\midrule", r"\endfirsthead",
     r"\multicolumn{7}{@{}l}{\footnotesize\textit{Table~\thetable{} (continued)}}\\",
     r"\toprule", hdr, r"\midrule", r"\endhead",
     r"\midrule", r"\multicolumn{7}{@{}r@{}}{\scriptsize\textit{continued on next page}}\\", r"\endfoot",
     r"\bottomrule", r"\endlastfoot"]
L += rows
L += [r"\end{longtable}", r"\endgroup"]
OUT_TEX.write_text("\n".join(L) + "\n")

# ----------------------------------------------------------------------------- CSV
NEUTRAL = [
    ("Round-2 audit: ", "Multiplicity audit (floored p): "),
    ("Round-2 re-measurement (revision2_experiments/headline_rebuilt.csv, N = 100,000 label shuffles): ",
     "Re-measurement with N = 100,000 label shuffles: "),
    ("Round 2: p=0.004 fails Bonferroni.", "Multiplicity audit (floored p): p=0.004 fails Bonferroni."),
    ("Round-3 local STRING file (contexts.STRING_LINKS = the v12 PHYSICAL subnetwork that CONVENTION3 section 3 specifies, 1,477,610 links)",
     "The local STRING file used in this paper (the v12 physical subnetwork, 1,477,610 links)"),
    ("Round-2 re-measurement (N = 100,000)", "Re-measurement (N = 100,000 label shuffles)"),
    ("Round-2 audit used the record's p", "The headline-level count of Table S12 uses the record's p"),
    ("round 2 replaced it by a Stouffer combination over layers (invalid: layers are dependent).",
     "the headline-level count of Table S12 scores this finding by a Stouffer combination of its 12 per-layer p-values, which "
     "is invalid because the layers share genes, so it is set aside in the final count of three."),
    ("dropping also the Stouffer-combined TF/target p leaves 3",
     "setting aside the invalid Stouffer combination for the TF/target classifier leaves 3"),
    ("(band_provenance.md)", "(Table S14)"),
    ("p_r2=", "p_floored="),
]
c = LED.drop(columns=["status_detail"]).copy()
for col_ in ("notes", "multiplicity_round2"):
    for a, b in NEUTRAL:
        c[col_] = c[col_].map(lambda v: v.replace(a, b) if isinstance(v, str) else v)
for col_ in c.columns:
    if c[col_].dtype == object:
        c[col_] = c[col_].map(cite_neutral)
c = c.rename(columns={"preprint_value": "reported_value", "preprint_location": "reported_location_in_Kendiukhov2026b",
                      "multiplicity_round2": "multiplicity_of_source_records_183_floored_p",
                      "this_paper_pointer": "this_paper_main_text_pointer"})
c["this_paper_main_text_pointer"] = c["this_paper_main_text_pointer"].str.replace("~", " ", regex=False).str.replace("$", "", regex=False)
txt = " ".join(c.astype(str).values.ravel())
for bad in ("round-2", "round 2", "round-3", "revision", "review", "submitted", "verifier", "preprint"):
    assert bad not in txt.lower(), bad
c.to_csv(OUT_CSV, index=False)
print("S13:", len(rows), "claims;", cnt)
