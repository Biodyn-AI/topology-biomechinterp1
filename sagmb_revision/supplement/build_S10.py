"""Build Supplementary Table S10 (cross-model comparison) from revision3_experiments/E6_crossmodel.

Writes
  sagmb_revision/supplement/S10.tex                               (LaTeX fragment, one longtable)
  sagmb_revision/supplementary_tables/Table_S10_cross_model.csv

All numbers are read from the E6 outputs: attrition.csv, tierA_results.csv, tierB_results.csv,
node_detstrat_results.csv, pairwise_deltas_vs_scgpt.csv, selection_bias.csv, projection_control.csv,
work/build_meta.json, work/selection_meta.json and the verifier files verify/centring_sensitivity.csv and
verify/centring_delta.csv. Network-statistic (PPI) intervals are always the delete-d columns (*_deleted).
Quantities computed here: ranges and counts over the 12 layers, the best layer (largest point estimate), the
labelled-minus-Tier-A centring differences, counts of target-only genes, and the Kendall tau between the 178-gene
values and the labelled-centred Tier A values (a point estimate the verifier printed but did not store; asserted
against the SUMMARY value 0.964).
Run:  /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python build_S10.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
E6 = ROOT / "revision3_experiments" / "E6_crossmodel"
OUT_TEX = HERE / "S10.tex"
OUT_CSV = HERE.parent / "supplementary_tables" / "Table_S10_cross_model.csv"

CTX = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
CTXNAME = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney~B", "immune": "immune",
           "kidney_A": "kidney~A"}
CELLS = {"lung": 256, "lung_external": 256, "kidney_B": 192, "immune": 256, "kidney_A": 256}
# display order: token tables, forward-pass outputs, sequence, expression atlas
GROUPS = [
    ("Single-cell models: input gene-token tables",
     ["scgpt_static", "gf_v1_10m", "gf_v2_104m", "gf_v2_316m", "tgpt"]),
    ("Single-cell models without a token table: one forward pass on a neutral input",
     ["aido_100m", "scfoundation"]),
    ("Protein sequence", ["esm2_uce", "esm2_650m", "seq_features"]),
    ("Expression atlas", ["hpa_coexp"]),
]
MODELS = [m for _, ms in GROUPS for m in ms]
NAME = {"scgpt_static": "scGPT (whole-human)", "gf_v1_10m": "Geneformer V1-10M", "gf_v2_104m": "Geneformer V2-104M",
        "gf_v2_316m": "Geneformer V2-316M", "tgpt": "tGPT", "aido_100m": "AIDO.Cell-100M",
        "scfoundation": "scFoundation", "esm2_uce": "ESM2-15B", "esm2_650m": "ESM2-650M",
        "seq_features": "32 sequence descriptors", "hpa_coexp": "HPA co-expression"}
TYPE = {"scgpt_static": "token table", "gf_v1_10m": "token table$^{a}$", "gf_v2_104m": "token table$^{a}$",
        "gf_v2_316m": "token table$^{a}$", "tgpt": "token table", "aido_100m": "forward pass$^{b}$",
        "scfoundation": "forward pass$^{b}$", "esm2_uce": "protein LM (UCE file)",
        "esm2_650m": "protein LM, mean-pooled", "seq_features": "hand-computed",
        "hpa_coexp": "154 cell types, PCA"}
TASKS = ["node_lambert", "node_trrust", "ppi", "edge"]
TASKNAME = {"node_lambert": "Lambert", "node_trrust": "TRRUST", "ppi": "PPI", "edge": "edge"}


# ---------------------------------------------------------------- formatting
def isnan(x):
    return x is None or (isinstance(x, float) and np.isnan(x))


def num(x, d=3, sign=False):
    if isnan(x):
        return "--"
    if sign and x != 0 and abs(x) < 0.5 * 10 ** (-d):
        d += 1  # show a tiny signed difference with one more decimal rather than as -0.000
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    if s.startswith("-"):
        return "$-$" + s[1:]
    if s.startswith("+"):
        return "$+$" + s[1:]
    return s


def ci(lo, hi, d=3):
    return f"[{num(lo, d)}, {num(hi, d)}]"


def rng(a, b, d=3):
    if round(a, d) >= 0 and round(b, d) >= 0:
        return f"{num(a, d)}--{num(b, d)}"
    return f"{num(a, d, True)} to {num(b, d, True)}"


def pval(p):
    if p >= 0.01:
        return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def intc(n):
    return f"{int(round(n)):,}".replace(",", "{,}")


def nl(*parts):
    return "\\newline ".join(parts)


# ---------------------------------------------------------------- load
att = pd.read_csv(E6 / "attrition.csv")
tA = pd.read_csv(E6 / "tierA_results.csv")
tB = pd.read_csv(E6 / "tierB_results.csv")
det = pd.read_csv(E6 / "node_detstrat_results.csv")
pdl = pd.read_csv(E6 / "pairwise_deltas_vs_scgpt.csv")
sel = pd.read_csv(E6 / "selection_bias.csv")
proj = pd.read_csv(E6 / "projection_control.csv")
csens = pd.read_csv(E6 / "verify" / "centring_sensitivity.csv")
cdel = pd.read_csv(E6 / "verify" / "centring_delta.csv")
bmeta = json.loads((E6 / "work" / "build_meta.json").read_text())
smeta = json.loads((E6 / "work" / "selection_meta.json").read_text())
res = json.loads((E6 / "results.json").read_text())

assert res["tierA_size"] == 17693
assert set(tA.embedding) == set(MODELS)


def A_row(model, task, est="full"):
    r = tA[(tA.embedding == model) & (tA.task == task) & (tA.estimator == est)]
    assert len(r) == 1, (model, task, est)
    return r.iloc[0]


def B_row(s, emb, task, est="full"):
    r = tB[(tB.set == s) & (tB.embedding == emb) & (tB.task == task) & (tB.estimator == est)]
    assert len(r) == 1, (s, emb, task, est)
    return r.iloc[0]


def det_row(s, emb):
    r = det[(det.set == s) & (det.embedding == emb)]
    assert len(r) == 1, (s, emb)
    return r.iloc[0]


def ival(row, task):
    """(value, lo, hi) with the delete-d interval for PPI."""
    if task == "ppi":
        assert not isnan(row["ci_lo_deleted"])
        return row["value"], row["ci_lo_deleted"], row["ci_hi_deleted"]
    return row["value"], row["ci_lo"], row["ci_hi"]


def att_n(stage, setname):
    r = att[(att.stage == stage) & (att.set == setname)]
    assert len(r) == 1, (stage, setname)
    return int(r.n_genes.iloc[0])


# ---------------------------------------------------------------- CSV
rows = []


def add(**kw):
    rows.append(kw)


# Panel A: attrition (every row of attrition.csv); internal set names reworded neutrally for the public CSV
NEUTRAL = {
    "8 GPL-rollout models (scGPT..ESM2-UCE)":
        "8 vocabularies: the seven single-cell models and ESM2-15B (UCE)",
    "Geneformer V1-10M & V2-104M & ESM2-UCE (submitted models)":
        "Geneformer V1-10M & V2-104M & ESM2-15B (UCE) (the vocabularies that define the 186-gene set)",
    "... & Geneformer V1-10M & V2-104M & ESM2-UCE vocabularies (submitted Table 1)":
        "... & Geneformer V1-10M & V2-104M & ESM2-15B (UCE) vocabularies (186-gene set)",
    "1152 rows (label 1/0) of the round-1 TRRUST-filtered kidney_A edge list":
        "1152 rows (label 1/0) of the automated screen's TRRUST-filtered kidney_A edge list",
    "kidney_A V & Tier A, TRRUST-participating (what the vocabularies would have allowed)":
        "kidney_A V & Tier A, TRRUST-participating (the most a kidney_A-only TRRUST analysis could use)",
}
for _, r in att.iterrows():
    note = r["note"] if isinstance(r["note"], str) else ""
    add(panel="A", source_file="attrition.csv", block=r["stage"], set=NEUTRAL.get(r["set"], r["set"]),
        statistic="n_genes", value=r["n_genes"], note=NEUTRAL.get(note, note))
for k in ["tierA", "legacy186", "legacy178"] + [f"tierB_{c}" for c in CTX]:
    for q, v in bmeta[k].items():
        add(panel="A", source_file="work/build_meta.json", block=k, set=k, statistic=q, value=v,
            note="n_imputed: Tier A genes not in lung V, assigned the median lung detection decile" if q == "n_imputed" else "")

COLS_RES = ["set", "embedding", "model", "label", "kind", "ctx", "layer", "d", "n_genes", "task", "estimator",
            "metric", "n_pos", "n_neg", "n_tf", "value", "p", "p_N", "p_method", "best_band", "value_vocab_centred",
            "value_detection_stratified"]


def add_result_rows(df, panel, src):
    for _, r in df.iterrows():
        d = {c: r[c] for c in COLS_RES}
        if r["task"] == "ppi" and r["estimator"] != "bestband" and r["estimator"] != "bestband_alllayers":
            d.update(ci_lo=r["ci_lo_deleted"], ci_hi=r["ci_hi_deleted"],
                     ci_method="vertex subsampling m=n/2 B=200, delete-d sqrt(m/(n-m))",
                     ci_lo_preregistered=r["ci_lo"], ci_hi_preregistered=r["ci_hi"])
        else:
            d.update(ci_lo=r["ci_lo"], ci_hi=r["ci_hi"], ci_method=r["ci_method"])
        d.update(panel=panel, source_file=src, statistic="AUROC")
        add(**d)


add_result_rows(tA, "B", "tierA_results.csv")
for _, r in det.iterrows():
    add(panel="B" if r["set"] == "tierA" else ("C" if r["set"].startswith("tierB") else "D"),
        source_file="node_detstrat_results.csv", set=r["set"], embedding=r["embedding"], label=r["label"],
        n_genes=r["n_genes"], n_pos=r["n_pos"], n_neg=r["n_neg"], task="node_lambert_detstrat",
        estimator=r["estimator"], metric=r["metric"], statistic="AUROC", value=r["value"], ci_lo=r["ci_lo"],
        ci_hi=r["ci_hi"], ci_method=r["ci_method"], p=r["p"], p_N=r["p_N"], p_method=r["p_method"],
        note=f"unstratified {r['unstratified_value']:.6f} [{r['unstratified_ci_lo']:.6f}, {r['unstratified_ci_hi']:.6f}]")
for _, r in csens.iterrows():
    for conv, col in (("centred over all Tier A genes", "trrust_centre_tierA"),
                      ("centred over the 2,753 TRRUST-labelled Tier A genes", "trrust_centre_labelled")):
        add(panel="B", source_file="verify/centring_sensitivity.csv", set="tierA", embedding=r["model"],
            task="node_trrust", estimator="full", statistic="AUROC", value=r[col], centring=conv,
            note="point estimate (no interval)")
    add(panel="B", source_file="verify/centring_sensitivity.csv", set="tierA", embedding=r["model"],
        task="node_trrust", estimator="full", statistic="AUROC labelled-centred minus Tier-A-centred",
        value=r["diff_labelled_minus_tierA"])
add_result_rows(tB, "C", "tierB_results.csv")
for _, r in pdl[pdl.set.isin(["tierA"] + [f"tierB_{c}" for c in CTX])].iterrows():
    ppi = r["task"] == "ppi"
    add(panel="B" if r["set"] == "tierA" else "C", source_file="pairwise_deltas_vs_scgpt.csv", set=r["set"],
        embedding=r["model_a"], label=r["label_a"], comparator=r["model_b"], comparator_label=r["label_b"],
        task=r["task"], estimator="full", statistic="delta (embedding minus comparator)", value=r["delta"],
        value_embedding=r["value_a"], value_comparator=r["value_b"],
        ci_lo=r["ci_lo_deleted"] if ppi else r["ci_lo"], ci_hi=r["ci_hi_deleted"] if ppi else r["ci_hi"],
        ci_method=("paired vertex subsampling m=n/2 B=200 (identical node subsets), delete-d sqrt(m/(n-m)); "
                   "pre-registered sqrt(m/n) interval in ci_*_preregistered") if ppi else r["ci_method"],
        ci_lo_preregistered=r["ci_lo"] if ppi else np.nan, ci_hi_preregistered=r["ci_hi"] if ppi else np.nan,
        ci_excludes_0=bool(r["ci_excludes_0_deleted"]) if ppi else bool(r["ci_excludes_0"]), B=r["B"])
for _, r in cdel.iterrows():
    add(panel="B" if r["set"] == "tierA" else "C", source_file="verify/centring_delta.csv", set=r["set"],
        embedding=r["model_a"], comparator=r["model_b"], task="node_trrust", estimator="full",
        statistic="delta (embedding minus comparator)", value=r["delta"], value_embedding=r["value_a"],
        value_comparator=r["value_b"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"],
        ci_method="paired gene bootstrap over the TRRUST-labelled genes, B=1000, percentile re-centred",
        n_pos=r["n_tf"], n_genes=r["n_lab"],
        centring={"centre_set": "centred over the analysed gene set",
                  "centre_labelled": "centred over the TRRUST-labelled genes"}[r["convention"]])
SELQ = ["legacy_value", "null_mean", "null_sd", "null_q025", "null_q50", "null_q975", "legacy_percentile",
        "p_two_sided", "legacy_minus_null_mean"]
for _, r in sel[sel.row_type == "subset_null"].iterrows():
    for q in SELQ:
        add(panel="D", source_file="selection_bias.csv", set=r["variant"], embedding=r["model"], label=r["label"],
            task="node_trrust", estimator="full", statistic=q, value=r[q], n_genes=r["n_legacy_genes_used"],
            p_N=r["null_N"], note="matched-subset null; tierA_value (Tier A centring) = %.6f" % r["tierA_value"])
for _, r in sel[sel.row_type == "kendall_tau"].iterrows():
    for q in ["kendall_tau_legacy_vs_tierA", "kendall_tau_subsetmean_vs_legacy", "kendall_tau_subsetmean_vs_tierA"]:
        add(panel="D", source_file="selection_bias.csv", set=r["variant"], embedding="ALL", task="node_trrust",
            statistic=q, value=r[q],
            ci_lo=r["tau_ci_lo"] if q == "kendall_tau_legacy_vs_tierA" else np.nan,
            ci_hi=r["tau_ci_hi"] if q == "kendall_tau_legacy_vs_tierA" else np.nan,
            ci_method=r["tau_ci_method"] + f", B={int(r['tau_boot_B'])}" if q == "kendall_tau_legacy_vs_tierA" else "",
            p=r["kendall_p_asymptotic"] if q == "kendall_tau_legacy_vs_tierA" else np.nan,
            n_models=r["n_models"],
            note=f"order_legacy={r['order_legacy']}; order_tierA={r['order_tierA']}; order_subset_mean={r['order_subset_mean']}")
# tau of legacy178 against labelled-centred Tier A (point estimate; trivial derived)
s178 = sel[(sel.variant == "legacy178") & (sel.row_type == "subset_null")].set_index("model")
cs = csens.set_index("model")
tau_lab = kendalltau(s178.legacy_value, cs.loc[s178.index, "trrust_centre_labelled"])[0]
tau_set = kendalltau(s178.legacy_value, cs.loc[s178.index, "trrust_centre_tierA"])[0]
assert abs(tau_lab - 0.964) < 5e-4, tau_lab
tau_stored = sel[sel.variant.eq("legacy178") & sel.row_type.eq("kendall_tau")].kendall_tau_legacy_vs_tierA.iloc[0]
assert abs(tau_set - tau_stored) < 1e-9, (tau_set, tau_stored)  # same computation reproduces the stored tau
add(panel="D", source_file="selection_bias.csv + verify/centring_sensitivity.csv (computed by build_S10.py)",
    set="legacy178", embedding="ALL", task="node_trrust",
    statistic="kendall_tau_legacy_vs_tierA_centred_over_labelled", value=tau_lab, n_models=len(s178),
    note="point estimate only")
for k in ("legacy_mean_trrust_degree", "pool_mean_trrust_degree"):
    for lab, v in smeta[k].items():
        add(panel="D", source_file="work/selection_meta.json", set="legacy186" if k.startswith("legacy") else
            "pool (TRRUST-labelled Tier A)", statistic=f"{k} ({lab})", value=v)
add(panel="D", source_file="work/selection_meta.json", set="pool", statistic="pool_size", value=smeta["pool_size"])
add(panel="D", source_file="work/selection_meta.json", set="pool", statistic="n_TF_pool", value=smeta["n_TF_pool"])
add(panel="D", source_file="work/selection_meta.json", set="legacy186", statistic="n_legacy_genes_not_in_lung_V",
    value=len(smeta["legacy_not_in_lungV"]), note=",".join(smeta["legacy_not_in_lungV"]))
for _, r in proj.iterrows():
    full = A_row(r["base_model"], r["task"])
    v, lo, hi = ival(full, r["task"])
    assert abs(v - r["full_value"]) < 1e-12
    for q in ["proj_mean", "proj_sd", "proj_min", "proj_max", "proj_minus_full"]:
        add(panel="D", source_file="projection_control.csv", set="tierA", embedding=r["family"],
            label=r["base_model"], d=r["d_projected"], task=r["task"], estimator="full (random Gaussian projection)",
            statistic=q, value=r[q], n_seeds=r["n_seeds"], note=f"seeds_values={r['seeds_values']}")
    add(panel="D", source_file="projection_control.csv + tierA_results.csv", set="tierA", embedding=r["base_model"],
        d=r["d_full"], task=r["task"], estimator="full", statistic="full_value", value=v, ci_lo=lo, ci_hi=hi,
        ci_method="as tierA_results.csv (PPI: delete-d)",
        ci_lo_preregistered=r["full_ci_lo"] if r["task"] == "ppi" else np.nan,
        ci_hi_preregistered=r["full_ci_hi"] if r["task"] == "ppi" else np.nan)

ORDER = ["panel", "source_file", "block", "set", "embedding", "model", "label", "kind", "ctx", "layer", "d",
         "comparator", "comparator_label", "task", "estimator", "centring", "metric", "statistic", "value",
         "value_embedding", "value_comparator", "ci_lo", "ci_hi", "ci_method", "ci_lo_preregistered",
         "ci_hi_preregistered", "ci_excludes_0", "B", "n_genes", "n_pos", "n_neg", "n_tf", "n_models", "n_seeds",
         "p", "p_N", "p_method", "best_band", "value_vocab_centred", "value_detection_stratified", "note"]
out = pd.DataFrame(rows)
for c in ORDER:
    if c not in out:
        out[c] = np.nan
out = out[ORDER]
# neutral names for the public CSV (values unchanged): the screen's 186- and 178-gene sets are the "curated" sets, as in
# the caption
for c in out.columns:
    if out[c].dtype == object:
        out[c] = out[c].map(lambda x: x.replace("(legacy and Tier A", "(curated set and Tier A").replace("legacy", "curated")
                            if isinstance(x, str) else x)
_txt = " ".join(out.astype(str).values.ravel()).lower() + " " + " ".join(out.columns).lower()
for bad in ("legacy", "submitted", "revision", "verifier", "round-"):
    assert bad not in _txt, bad
# the comparison class (other genes, target-only genes, rewired non-pairs or candidate non-edges, by task) is named as in
# the table header, not as negatives
assert "n_comparison" not in out.columns
out = out.rename(columns={"n_neg": "n_comparison"})
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_CSV, index=False)

# ---------------------------------------------------------------- LaTeX
L = []
W = "16.9cm"


def head(title):
    L.append(f"\\multicolumn{{7}}{{@{{}}p{{{W}}}@{{}}}}{{\\footnotesize\\textbf{{{title}}}}}\\\\*")


def colhead(*cells):
    L.append(" & ".join(cells) + "\\\\*")
    L.append("\\cmidrule(lr){1-7}")


def line(*cells):
    assert len(cells) == 7, cells
    L.append(" & ".join(cells) + "\\\\")


def span(first, second, text):
    L.append(f"{first} & {second} & \\multicolumn{{5}}{{>{{\\scriptsize\\raggedright\\arraybackslash}}p{{11.6cm}}@{{}}}}{{{text}}}\\\\")


# ---- Panel A
head("A. Gene sets: vocabularies, intersections, the two tiers and the curated set")
L.append("Gene set & Genes & \\multicolumn{5}{>{\\scriptsize\\raggedright\\arraybackslash}p{11.6cm}@{}}"
         "{Composition, or genes lost at this step}\\\\*")
L.append("\\cmidrule(lr){1-7}")
voc = att[att.stage == "vocabulary"].set_index("set").n_genes
vtxt = "; ".join(f"{NAME[m]} {intc(voc[m])}" + (" (181 all-zero rows dropped)" if m == "seq_features" else "")
                 for m in MODELS)
span("Vocabulary of each static representation (genes with a non-zero vector)",
     f"{intc(voc.min())}--{intc(voc.max())}",
     vtxt)
n3 = att_n("intersection", "Geneformer V1-10M & V2-104M & ESM2-UCE (submitted models)")
span("Geneformer V1-10M $\\cap$ V2-104M $\\cap$ ESM2-15B", intc(n3),
     "the vocabularies used to define the 186-gene curated set (below)")
a_l = A_row("scgpt_static", "node_lambert"); a_t = A_row("scgpt_static", "node_trrust")
a_p = A_row("scgpt_static", "ppi"); a_e = A_row("scgpt_static", "edge")
nA = att_n("Tier A", "intersection of all 11 static vocabularies")
assert nA == 17693 and a_l.n_pos == att_n("Tier A", "Lambert TFs in Tier A")
span("Tier A: all 11 vocabularies", intc(nA),
     f"Lambert TFs {intc(a_l.n_pos)} vs {intc(a_l.n_neg)} other genes; TRRUST TFs {intc(a_t.n_pos)} vs "
     f"{intc(a_t.n_neg)} target-only genes; STRING physical pairs {intc(a_p.n_pos)} vs {intc(a_p.n_neg)} rewired "
     f"non-pairs; TRRUST edges {intc(a_e.n_pos)} ({intc(a_e.n_tf)} TFs) vs {intc(a_e.n_neg)} candidate non-edges; "
     f"{intc(bmeta['tierA']['n_imputed'])} genes not observed in lung (assigned the median lung detection decile)")
for c in CTX:
    st = f"context {c}"
    span(f"Tier B, {CTXNAME[c]}", intc(att_n(st, "Tier B = V & Tier A")),
         f"gene axis $G$ {intc(att_n(st, 'gene axis (h5ad var)'))}; observed $V$ "
         f"{intc(att_n(st, 'V (observed, non-zero embedding)'))}; Lambert TFs in Tier B "
         f"{intc(att_n(st, 'Lambert TFs in Tier B'))}; genes in a TRRUST edge within Tier B "
         f"{intc(att_n(st, 'TRRUST-participating within Tier B'))}; {CELLS[c]} cells")
leg = att[att.stage == "legacy"].reset_index(drop=True)
assert list(leg.n_genes) == [209, 209, 187, 186, 178, 172]
hla = leg.note[4].split(": ")[1].replace(",", ", ")
span("Genes named in the curated kidney~A edge list", intc(leg.n_genes[0]),
     "1{,}152 TRRUST edge and sampled non-edge rows, so the list includes genes that occur only in non-edges; "
     "all 209 are on the kidney~A gene axis")
span("\\quad observed in kidney~A", intc(leg.n_genes[2]), "22 never-observed genes dropped")
m186 = bmeta["legacy186"]; m178 = bmeta["legacy178"]
assert m186["n_trrust_tf"] + m186["n_trrust_tgt"] == 186 and m178["n_trrust_tf"] + m178["n_trrust_tgt"] == 178
span("\\quad $\\cap$ Geneformer V1-10M, V2-104M, ESM2-15B: the 186-gene set", intc(leg.n_genes[3]),
     f"lost: NEAT1 (lncRNA); {m186['n_trrust_tf']} TRRUST TFs, {m186['n_trrust_tgt']} target-only genes; "
     f"{len(smeta['legacy_not_in_lungV'])} genes not observed in lung")
span("\\quad $\\cap$ Tier A: the 178-gene set", intc(leg.n_genes[4]),
     f"lost (absent from tGPT): {hla}; {m178['n_trrust_tf']} TRRUST TFs, {m178['n_trrust_tgt']} target-only genes")
span("Kidney~A $V \\cap {}$Tier A, in a TRRUST edge", intc(leg.n_genes[5]),
     f"the most a kidney~A-only analysis of the TRRUST tasks could use (kidney~A observes "
     f"{intc(att_n('context kidney_A', 'V (observed, non-zero embedding)'))} of its "
     f"{intc(att_n('context kidney_A', 'gene axis (h5ad var)'))} genes)")

# ---- Panel B
L.append("\\midrule")
head("B. Tier A (17{,}693 genes): every static representation, full dimension, centred over the analysed genes")
colhead("Representation", "Type; $d$", "Lambert: stratified (unstratified) [95\\% CI]",
        "TRRUST, centred over Tier A [95\\% CI]", "TRRUST, centred over labelled genes",
        "PPI [95\\% CI, delete-$d$]", "Edge, within TF [95\\% CI]")
for gi, (gname, ms) in enumerate(GROUPS):
    if gi:
        L.append("\\cmidrule(lr){1-7}")
    for m in ms:
        dr = det_row("tierA", m)
        lam = A_row(m, "node_lambert"); assert abs(dr.unstratified_value - lam.value) < 1e-12
        tr = A_row(m, "node_trrust"); pp = A_row(m, "ppi"); ed = A_row(m, "edge")
        assert abs(cs.loc[m, "trrust_centre_tierA"] - tr.value) < 1e-9
        pv, plo, phi = ival(pp, "ppi")
        line(NAME[m], f"{TYPE[m]}; {intc(lam.d)}",
             nl(f"{num(dr.value)} ({num(lam.value)})", ci(dr.ci_lo, dr.ci_hi)),
             nl(num(tr.value), ci(tr.ci_lo, tr.ci_hi)),
             num(cs.loc[m, "trrust_centre_labelled"]),
             nl(num(pv), ci(plo, phi)),
             nl(num(ed.value), ci(ed.ci_lo, ed.ci_hi)))
L.append("\\cmidrule(lr){1-7}")


def pd_row(s, a, b, task):
    r = pdl[(pdl.set == s) & (pdl.model_a == a) & (pdl.model_b == b) & (pdl.task == task)]
    assert len(r) == 1, (s, a, b, task)
    return r.iloc[0]


def cd_row(s, conv, a, b):
    r = cdel[(cdel.set == s) & (cdel.convention == conv) & (cdel.model_a == a) & (cdel.model_b == b)]
    assert len(r) == 1
    return r.iloc[0]


for m in ["esm2_uce", "gf_v2_104m"]:
    cells = []
    for task in ["node_lambert_detstrat", "node_trrust"]:
        r = pd_row("tierA", m, "scgpt_static", task)
        cells.append(nl(num(r.delta, sign=True), ci(r.ci_lo, r.ci_hi)))
    r = cd_row("tierA", "centre_labelled", m, "scgpt_static")
    r_set = cd_row("tierA", "centre_set", m, "scgpt_static")
    assert abs(r_set.delta - pd_row("tierA", m, "scgpt_static", "node_trrust").delta) < 1e-9
    cells.append(nl(num(r.delta, sign=True), ci(r.ci_lo, r.ci_hi)))
    r = pd_row("tierA", m, "scgpt_static", "ppi")
    cells.append(nl(num(r.delta, sign=True), ci(r.ci_lo_deleted, r.ci_hi_deleted)))
    r = pd_row("tierA", m, "scgpt_static", "edge")
    cells.append(nl(num(r.delta, sign=True), ci(r.ci_lo, r.ci_hi)))
    line(f"{NAME[m]} minus scGPT (paired)", "", cells[0], cells[1], cells[2], cells[3], cells[4])

# ---- Panel C
L.append("\\midrule")
head("C. Tier B: scGPT's residual stream ($L_0$--$L_{11}$) against three static representations on each "
     "extraction's observed genes")
colhead("Extraction: task", "Positives / comparison set", "Residual stream: range over 12 layers",
        "Best layer: AUROC [95\\% CI]", "ESM2-15B: AUROC; minus layer; $k_{+}/k_{-}$",
        "scGPT token table: AUROC; minus layer; $k_{+}/k_{-}$",
        "Geneformer V2-104M: AUROC; minus layer; $k_{+}/k_{-}$")
CTASK = [("node_lambert_detstrat", "Lambert (stratified)"), ("node_trrust", "TRRUST"), ("ppi", "PPI"),
         ("edge", "edge")]
for ci_, c in enumerate(CTX):
    s = f"tierB_{c}"
    if ci_:
        L.append("\\cmidrule(lr){1-7}")
    for task, tname in CTASK:
        layers = [f"{c}@L{k}" for k in range(12)]
        vals = []
        for emb in layers:
            if task == "node_lambert_detstrat":
                r = det_row(s, emb); vals.append((r.value, r.ci_lo, r.ci_hi))
            else:
                vals.append(ival(B_row(s, emb, task), task))
        v = np.array([x[0] for x in vals]); kb = int(np.argmax(v))
        if task == "node_lambert_detstrat":
            r0 = det_row(s, layers[0]); npos, nneg, ntf = r0.n_pos, r0.n_neg, np.nan
        else:
            r0 = B_row(s, layers[0], task); npos, nneg, ntf = r0.n_pos, r0.n_neg, r0.n_tf
        cnt = f"{intc(npos)} / {intc(nneg)}" + (f" ({intc(ntf)} TFs)" if task == "edge" else "")
        cells = [f"{CTXNAME[c]}: {tname}", cnt, rng(v.min(), v.max()),
                 nl(f"$L_{{{kb}}}$: {num(v[kb])}", ci(vals[kb][1], vals[kb][2]))]
        for m in ["esm2_uce", "scgpt_static", "gf_v2_104m"]:
            d = pdl[(pdl.set == s) & (pdl.model_a == m) & (pdl.task == task) & (pdl.model_b.isin(layers))]
            assert len(d) == 12 and set(d.model_b) == set(layers)
            lo = d.ci_lo_deleted if task == "ppi" else d.ci_lo
            hi = d.ci_hi_deleted if task == "ppi" else d.ci_hi
            assert lo.notna().all()
            kp, km = int((lo > 0).sum()), int((hi < 0).sum())
            # comparator's own value on the same genes
            if task == "node_lambert_detstrat":
                cv = det_row(s, m).value
            else:
                cv = B_row(s, m, task).value
            assert np.allclose(d.value_a, cv, atol=1e-9)
            # contextual values in the delta file agree with the per-layer values used above
            dv = d.set_index("model_b").loc[layers, "value_b"].to_numpy()
            assert np.allclose(dv, v, atol=1e-9), (s, task, m)
            cells.append(nl(num(cv), rng(d.delta.min(), d.delta.max()), f"{kp}/{km}"))
        line(*cells)

# ---- Panel D
L.append("\\midrule")
head("D. The 186- and 178-gene curated sets against matched random subsets of Tier A (TRRUST task), "
     "and a projection control")
colhead("Representation", "186 genes: AUROC", "186: null mean [2.5\\%, 97.5\\%]", "186: percentile; two-sided $p$",
        "186: excess over null mean", "178 genes: AUROC; percentile; excess",
        "186, degree-matched: null mean; percentile; excess")
sub = sel[sel.row_type == "subset_null"]


def srow(v, m):
    r = sub[(sub.variant == v) & (sub.model == m)]
    assert len(r) == 1
    return r.iloc[0]


for gi, (gname, ms) in enumerate(GROUPS):
    if gi:
        L.append("\\cmidrule(lr){1-7}")
    for m in ms:
        a, b, g = srow("legacy186", m), srow("legacy178", m), srow("legacy186_degmatched", m)
        mark = "$^{c}$" if int(a.n_legacy_genes_used) != 186 else ""
        line(NAME[m], num(a.legacy_value) + mark, nl(num(a.null_mean), ci(a.null_q025, a.null_q975)),
             nl(num(a.legacy_percentile, 4), f"$p = {pval(a.p_two_sided).strip('$')}$"),
             num(a.legacy_minus_null_mean, sign=True),
             nl(f"{num(b.legacy_value)}; {num(b.legacy_percentile, 4)}", num(b.legacy_minus_null_mean, sign=True)),
             nl(f"{num(g.null_mean)}; {num(g.legacy_percentile, 4)}", num(g.legacy_minus_null_mean, sign=True)))
L.append("\\cmidrule(lr){1-7}")
kt = sel[sel.row_type == "kendall_tau"].set_index("variant")
k186, k178 = kt.loc["legacy186"], kt.loc["legacy178"]
span("Kendall $\\tau$, order of the representations, curated set vs Tier A", "",
     f"186 genes ({int(k186.n_models)} representations, without tGPT): {num(k186.kendall_tau_legacy_vs_tierA)} "
     f"{ci(k186.tau_ci_lo, k186.tau_ci_hi)}; 178 genes ({int(k178.n_models)}): "
     f"{num(k178.kendall_tau_legacy_vs_tierA)} {ci(k178.tau_ci_lo, k178.tau_ci_hi)}; 178 genes vs Tier A centred "
     f"over its TRRUST-labelled genes: {num(tau_lab)} (point estimate)")
dl, dp = smeta["legacy_mean_trrust_degree"], smeta["pool_mean_trrust_degree"]
span("Mean TRRUST degree, 186-gene set vs pool", "",
     f"TFs {dl['TF']:.1f} vs {dp['TF']:.1f}; target-only genes {dl['TG']:.1f} vs {dp['TG']:.1f}")
L.append("\\cmidrule(lr){1-7}")
colhead("Projection control (Tier A): representation, task", "$d$: full $\\to$ projected",
        "Full: AUROC [95\\% CI]", "Projected: mean (SD) over 5 seeds", "Projected: min--max",
        "Projected minus full", "scGPT token table, full 512-$d$")
for _, r in proj.iterrows():
    full = A_row(r["base_model"], r["task"]); fv, flo, fhi = ival(full, r["task"])
    tn = {"node_lambert": "Lambert (unstratified)", "node_trrust": "TRRUST", "ppi": "PPI",
          "edge": "edge"}[r["task"]]
    ref = num(r["scgpt_static_full_value"]) if r["base_model"] != "scgpt_static" else "--"
    line(f"{NAME[r['base_model']]}: {tn}", f"{intc(r['d_full'])} $\\to$ {intc(r['d_projected'])}",
         nl(num(fv), ci(flo, fhi)), f"{num(r['proj_mean'])} ({num(r['proj_sd'])})",
         rng(r["proj_min"], r["proj_max"]), num(r["proj_minus_full"], sign=True), ref)

body = "\n".join(L)

# ---------------------------------------------------------------- caption
labdiff = (cs.trrust_centre_labelled - cs.trrust_centre_tierA)
expr = ["scgpt_static", "gf_v1_10m", "gf_v2_104m", "gf_v2_316m", "tgpt", "aido_100m"]
other = ["scfoundation", "esm2_uce", "esm2_650m", "seq_features", "hpa_coexp"]
lo_e, hi_e = labdiff[expr].min(), labdiff[expr].max()
hi_o = labdiff[other].abs().max()
assert hi_o <= 0.010
pA = det[det.set == "tierA"].p
assert np.allclose(pA, 1 / 2001)
cap = (
    "\\textbf{Cross-model comparison: gene sets, eleven static representations on 17{,}693 genes, scGPT's residual "
    "stream against static representations, and the bias of a 186-gene curated set.} "
    "Unless stated otherwise, values in Panels B--D are AUROCs of full-dimensional vectors centred over the analysed gene "
    "set and compared by cosine. "
    "\\emph{Lambert}: leave-one-out (LOO) mean cosine to the census transcription factors (TFs) of Lambert et al., TFs "
    "against all other genes; \\emph{TRRUST}: the same score, TRRUST TFs against genes that are only targets; "
    "\\emph{PPI}: STRING physical pairs (score $\\geq 700$) against non-pairs from one degree-preserving rewiring; "
    "\\emph{edge}: within-TF AUROC of the TF--target cosine, TRRUST edges against candidate non-edges (targets of other "
    "TFs, absent from TRRUST and DoRothEA). 95\\% intervals: node tasks, gene bootstrap ($B = 1{,}000$); PPI, vertex "
    "subsampling ($m = n/2$, $B = 200$) with delete-$d$ scaling; edge, two-way TF $\\times$ target bootstrap "
    "($B = 2{,}000$); paired differences reuse the same genes and draws for both representations. $d$, dimension; LM, "
    "language model; HPA, Human Protein Atlas; $G$, gene axis; $V$, observed genes. "
    "\\textbf{A}: Tier A, genes in all 11 vocabularies; Tier B, $V \\cap {}$Tier A; the curated set starts from the "
    "automated screen's kidney~A TRRUST edge list. "
    "\\textbf{B}: token table, the model's input gene-embedding table ($^{a}$with its LayerNorm applied); $^{b}$gene outputs of one forward pass on the HPA "
    "mean-expression profile (no gene-token table, so not like-for-like with the token tables); ESM2-15B, the "
    "5{,}120-dimensional per-protein embedding distributed with UCE; ESM2-650M, mean-pooled over the canonical "
    "sequence; HPA, single-cell expression over 154 cell types, 128 principal components. Lambert is pooled within "
    "lung detection deciles (unstratified in parentheses); the "
    f"{intc(bmeta['tierA']['n_imputed'])} genes not observed in lung all sit in the median decile, so this "
    "stratification is largely nominal (Panel C controls detection); within-decile permutation $p = 1/2{,}001$ "
    "throughout. TRRUST depends on the centring set: all Tier A genes, or only the "
    f"{intc(smeta['pool_size'])} labelled genes ({intc(a_t.n_pos)} TFs, {intc(a_t.n_neg)} target-only; point "
    "estimates; the convention of Panel D). Labelled-gene centring raises the single-cell models other than "
    f"scFoundation by {num(lo_e)}--{num(hi_e)} and moves the others by at most {num(hi_o)}, roughly halving the gap "
    "between protein and expression-trained representations (last two rows: paired differences from scGPT "
    "[95\\% CI], for Lambert of the stratified values). "
    "\\textbf{C}: scGPT whole-human residual stream, one extraction per context (192--256 cells), centred over the "
    "Tier B set; Lambert is stratified within the extraction's deciles of $\\log n_{\\mathrm{obs}}$ (the H2 "
    "statistic), the other tasks are unstratified. Positives / comparison set: TFs / other genes (Lambert), TFs / "
    "target-only genes (TRRUST), pairs / rewired non-pairs (PPI), documented edges / candidate non-edges (edge). "
    "Best layer: highest AUROC, chosen post hoc. Static columns: the "
    "representation's AUROC on the same genes; range over the 12 layers of the paired difference, static minus layer; "
    "$k_{+}/k_{-}$, layers whose interval lies entirely above / below 0. "
    "\\textbf{D}: TRRUST AUROC (centred over the subset) on the curated set and in 10{,}000 random subsets of the "
    f"{intc(smeta['pool_size'])} labelled Tier A genes matched on size, TF:target-only ratio "
    f"({m186['n_trrust_tf']}:{m186['n_trrust_tgt']} for 186) and lung mean-expression decile (degree-matched: "
    "TRRUST-degree quartiles and lung mean-expression quintiles; degree, TRRUST out-degree of a TF or in-degree of a "
    "target-only gene; pool, the labelled genes). Percentile: share of null subsets below the curated value (ties half); two-sided $p$ "
    "with $+1$ correction (floor $2.0\\times10^{-4}$); excess, curated value minus null mean; $^{c}$scored on the 178 "
    "genes in tGPT's vocabulary. Kendall $\\tau$: agreement of the representation orders, with a paired "
    "gene-bootstrap interval ($B = 1{,}000$). Projection: Gaussian random projection (5 seeds) at a 5:1 ratio, scored "
    "as in Panel B (Lambert unstratified). All comparisons are sensitivity analyses outside the pre-registered "
    "hypotheses. "
    "For PPI, the vertex-subsampling intervals with the pre-registered scaling factor $\\sqrt{m/n}$, which omits the "
    "finite-population term and so gives half-widths smaller by a factor $\\sqrt{2}$, are in the CSV (columns "
    "\\texttt{ci\\_lo\\_preregistered}/\\texttt{ci\\_hi\\_preregistered}). "
    "Full values: \\texttt{Table\\_S10\\_cross\\_model.csv}."
)

tex = f"""\\setcounter{{table}}{{9}}
\\begingroup
\\footnotesize
\\setlength{{\\tabcolsep}}{{3pt}}
\\renewcommand{{\\arraystretch}}{{1.05}}
\\begin{{longtable}}{{@{{}}>{{\\scriptsize\\raggedright\\arraybackslash}}p{{2.9cm}}>{{\\scriptsize\\raggedright\\arraybackslash}}p{{2.0cm}}*{{2}}{{>{{\\scriptsize\\raggedright\\arraybackslash}}p{{2.1cm}}}}*{{3}}{{>{{\\scriptsize\\raggedright\\arraybackslash}}p{{2.2cm}}}}@{{}}}}
\\caption{{{cap}}}\\label{{tab:S10}}\\\\
\\toprule
\\endfirsthead
\\multicolumn{{7}}{{@{{}}l}}{{\\footnotesize Table~\\thetable\\ (continued)}}\\\\
\\toprule
\\endhead
\\midrule
\\multicolumn{{7}}{{r@{{}}}}{{\\footnotesize continued on next page}}\\\\
\\endfoot
\\bottomrule
\\endlastfoot
{body}
\\end{{longtable}}
\\endgroup
"""
OUT_TEX.write_text(tex)
print("wrote", OUT_TEX, "and", OUT_CSV, f"({len(out)} rows)")
