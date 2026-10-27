"""E9 / A6 step 7 -- multiplicity audit over ALL screen records, and re-derivation of the 183-record audit.

Inputs (read-only):
  iterations/iter_XXXX/executor_hypothesis_screen.json          (every iteration that has one)
  revision_experiments/e5_multiple_testing/master_p.csv           (round-1, iters 1-63)
  revision2_experiments/master_p_round2.csv, audit_log_round2.json, headline_audit_round2.csv,
  revision2_experiments/headline_rebuilt.csv                      (round-2)

Procedure (identical to round 1 / round 2 so that the 183-record numbers can be re-derived):
  * primary p = round-1 regex extraction from the free-text `result_value` (empirical/permutation p first,
    else the minimum of any p=/p</p<= value), code copied verbatim from e5_multiple_testing/aggregate.py;
  * N = max integer before 'perm'/'permut'/'shuffle' in the record's `method` string (round-2 regex), else 1,000;
  * p == 0 -> resolution floor 1/(N+1);
  * Bonferroni at 0.05/n (n = 183 for the round-2 re-derivation, n = number of records with extractable p for
    the all-records audit) and BH at q = 0.05.

Outputs (this folder): audit_409_records.csv, audit_409.json
"""
from __future__ import annotations
import json, re
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
P41 = HERE.parents[1]
ITER = P41 / "iterations"
R1 = P41 / "revision_experiments/e5_multiple_testing"
R2 = P41 / "revision2_experiments"

# ---------------- round-1 extraction (verbatim logic from aggregate.py) ----------------
P_PATTERNS = [
    ("emp_p_eq",      re.compile(r"\bemp(?:irical)?_p\s*=\s*([0-9eE.+\-]+)")),
    ("emp_p_lt",      re.compile(r"\bemp(?:irical)?_p\s*<\s*([0-9eE.+\-]+)")),
    ("emp_p_le",      re.compile(r"\bemp(?:irical)?_p\s*≤\s*([0-9eE.+\-]+)")),
    ("perm_p_eq",     re.compile(r"\b(?:permutation|perm)_p\s*=\s*([0-9eE.+\-]+)")),
    ("perm_p_lt",     re.compile(r"\b(?:permutation|perm)_p\s*<\s*([0-9eE.+\-]+)")),
    ("perm_p_le",     re.compile(r"\b(?:permutation|perm)_p\s*≤\s*([0-9eE.+\-]+)")),
    ("p_perm_eq",     re.compile(r"\bp_perm\s*=\s*([0-9eE.+\-]+)")),
    ("p_perm_lt",     re.compile(r"\bp_perm\s*<\s*([0-9eE.+\-]+)")),
    ("p_perm_le",     re.compile(r"\bp_perm\s*≤\s*([0-9eE.+\-]+)")),
    ("p_emp_paren",   re.compile(r"\bp\s*=\s*([0-9eE.+\-]+)\s*\(\s*emp", re.IGNORECASE)),
    ("p_eq",          re.compile(r"(?<!_)\bp\s*=\s*([0-9eE.+\-]+)")),
    ("p_lt",          re.compile(r"(?<!_)\bp\s*<\s*([0-9eE.+\-]+)")),
    ("p_le",          re.compile(r"(?<!_)\bp\s*≤\s*([0-9eE.+\-]+)")),
    ("p_value_eq",    re.compile(r"\bp[-_ ]?value\s*=\s*([0-9eE.+\-]+)")),
]

def _coerce_p(s):
    s = s.rstrip(".,;)")
    try: v = float(s)
    except ValueError: return None
    return v if 0 <= v <= 1 else None

def extract_primary_p(rv):
    if not isinstance(rv, str): return None, None, []
    allp = []
    for label, pat in P_PATTERNS:
        ms = pat.findall(rv)
        if ms:
            ps = [p for p in (_coerce_p(m) for m in ms) if p is not None]
            if ps:
                allp.extend(ps)
                if "emp" in label or "perm" in label or "p_perm" in label:
                    return float(np.min(ps)), label, allp
    if allp: return float(np.min(allp)), "fallback_min_p", allp
    return None, None, []

NPERM_RE = re.compile(r'(\d[\d,]*)\s*[- ]?(?:perm|permut|shuffle)', re.I)     # round-2 regex, verbatim
def nperm_from_method(method):
    m = NPERM_RE.findall(method or "")
    return max(int(x.replace(",", "")) for x in m) if m else None

def bh_fdr(p):                      # same algorithm as stats3.bh_fdr / round-2 common.bh_fdr
    p = np.asarray(p, dtype=np.float64); n = p.size
    order = p.argsort(); q = np.empty(n); prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]; prev = min(prev, p[i] * n / (rank + 1)); q[i] = min(prev, 1.0)
    return q

# round-1 headline list (verbatim from aggregate.py HEADLINE_FINDINGS)
HEADLINES = {
    "SV1_extracellular_enrichment": ("iter_0006", "H02"),
    "SV1_secretory_pathway_layered": ("iter_0008", "H02"),
    "SV2_PPI_copole_string07": ("iter_0010", "H02"),
    "SV2_TRRUST_copole": ("iter_0011", "H01"),
    "TF_vs_target_AUROC_joint": ("iter_0056", "H01"),
    "celltype_marker_AUROC_851": ("iter_0023", "H01"),
    "edge_AUROC_decay": ("iter_0062", "H01"),
    "string_confidence_gradient": ("iter_0015", "H02"),
    "bcell_attractor_convergence": ("iter_0042", "H03"),
    "celltype_marker_AUROC_853_initial": ("iter_0022", "H02"),
    "TF_target_AUROC_cross_seed": ("iter_0057", "H01"),
}

def load_records():
    rows = []
    for p in sorted(ITER.glob("iter_*/executor_hypothesis_screen.json")):
        if p.name.startswith("._"): continue
        d = json.loads(p.read_text())
        it = p.parent.name
        for k, h in enumerate(d.get("hypotheses", [])):
            pp, src, allp = extract_primary_p(str(h.get("result_value", "")))
            rows.append(dict(iteration=it, iter_num=int(it[5:9]), aborted_duplicate_run=("incomplete" in it),
                             record_index=k, hypothesis_id=h.get("id", ""),
                             name=h.get("name", ""), family=h.get("family", "unknown"),
                             decision=h.get("decision", ""), result_direction=h.get("result_direction", ""),
                             method=h.get("method", ""), result_value=str(h.get("result_value", "")),
                             primary_p_raw=pp, p_source=src, n_p_in_string=len(allp),
                             n_perm_method=nperm_from_method(h.get("method", ""))))
    return pd.DataFrame(rows)

def audit(df, denom=None, alpha=0.05, default_n=1000):
    df = df.copy()
    df["n_perm_used"] = df.n_perm_method.fillna(default_n).astype(int)
    df["resolution_limited"] = df.primary_p_raw.notna() & (df.primary_p_raw == 0.0)
    df["p_floored"] = np.where(df.resolution_limited, 1.0 / (df.n_perm_used + 1), df.primary_p_raw)
    ok = df.p_floored.notna()
    n_ext = int(ok.sum()); m = denom if denom is not None else n_ext
    q = np.full(len(df), np.nan); q[ok.values] = bh_fdr(df.loc[ok, "p_floored"].to_numpy())
    df["bh_q"] = q
    df["bh_pass"] = df.bh_q < alpha
    df["bonferroni_alpha"] = alpha / m
    df["bonferroni_pass"] = df.p_floored < alpha / m
    return df, n_ext

def main():
    df = load_records()
    out = {"n_records_all": int(len(df)), "n_iterations_with_json": int(df.iteration.nunique()),
           "iterations_with_json_first_last": [df.iteration.min(), df.iteration.max()]}
    dup = df.groupby(["iteration", "hypothesis_id"]).size()
    out["duplicate_iteration_id_pairs"] = {f"{a}/{b}": int(c) for (a, b), c in dup[dup > 1].items()}

    # ---------------- (a) re-derive the 183-record round-1 / round-2 numbers ----------------
    d183 = df[df.iter_num <= 63].copy().reset_index(drop=True)
    r1 = pd.read_csv(R1 / "master_p.csv")
    # row-by-row identity with round-1 master_p.csv (same order: iteration, record order)
    same_n = len(r1) == len(d183)
    pr = d183.primary_p_raw.to_numpy(dtype=float); pm = r1.primary_p.to_numpy(dtype=float)
    p_identical = bool(same_n and np.array_equal(np.isnan(pr), np.isnan(pm)) and
                       bool(np.allclose(np.nan_to_num(pr, nan=-1), np.nan_to_num(pm, nan=-1), rtol=1e-12, atol=0)) and
                       list(d183.iteration) == list(r1.iter) and list(d183.hypothesis_id) == list(r1.hypothesis_id))
    ext = d183.primary_p_raw.notna()
    thr183 = 0.05 / 183
    r1_bonf = int((d183.primary_p_raw[ext] <= thr183).sum())             # round 1: '<=' and raw p (zeros kept)
    r1_bh = int((bh_fdr(d183.primary_p_raw[ext].to_numpy()) <= 0.05).sum())
    a183, n_ext183 = audit(d183, denom=183)
    r2 = pd.read_csv(R2 / "master_p_round2.csv")
    p2_identical = bool(np.allclose(np.nan_to_num(a183.p_floored.to_numpy(float), nan=-1),
                                    np.nan_to_num(r2.p_round2.to_numpy(float), nan=-1), rtol=0, atol=1e-15))
    n_identical = bool((a183.n_perm_used.to_numpy() == r2.n_perm.to_numpy()).all())
    # headline 8 -> 5 (round-2 substitutes re-measured p for five claims from headline_rebuilt.csv)
    reb = pd.read_csv(R2 / "headline_rebuilt.csv")
    best = lambda pre: (None if reb[reb.claim.str.startswith(pre)].empty else float(reb[reb.claim.str.startswith(pre)].p.min()))
    REMEAS = {"SV1_extracellular_enrichment": best("SV1_extracellular_L11"),
              "SV1_secretory_pathway_layered": best("SV1_extracellular_best_layer"),
              "SV2_PPI_copole_string07": best("SV2_STRING_copole"),
              "SV2_TRRUST_copole": best("SV2_TRRUST_copole"),
              "TF_vs_target_AUROC_joint": best("TF_target_AUROC_SV2_7_stouffer")}
    hl_r1 = hl_r2 = 0; hl_rows = []
    for tag, (it, hid) in HEADLINES.items():
        rr = d183[(d183.iteration == it) & (d183.hypothesis_id == hid)].iloc[0]
        p1 = rr.primary_p_raw
        p2 = REMEAS.get(tag) if REMEAS.get(tag) is not None else (1 / 1001 if p1 == 0 else p1)
        hl_r1 += int(p1 <= thr183); hl_r2 += int(p2 < thr183)
        hl_rows.append((tag, p1, p2))
    out["rederive_183"] = dict(
        n_total=int(len(d183)), n_extractable=int(ext.sum()), n_p_equal_zero=int((d183.primary_p_raw == 0).sum()),
        round1_bonferroni_pass_0p05_over_183=r1_bonf, round1_bh_pass=r1_bh,
        round2_bonferroni_pass=int(a183.bonferroni_pass.sum()), round2_bh_pass=int(a183.bh_pass.sum()),
        round2_bonferroni_pass_resolution_limited=int((a183.bonferroni_pass & a183.resolution_limited).sum()),
        headline_bonferroni_round1=hl_r1, headline_bonferroni_round2=hl_r2,
        primary_p_identical_to_round1_master=p_identical, floored_p_identical_to_round2_master=p2_identical,
        n_perm_identical_to_round2_master=n_identical,
        expected=dict(n_extractable=112, n_p0=18, bonferroni_r1=43, bonferroni_r2=25, bh=70, headline="8->5"),
    )
    o = out["rederive_183"]
    o["reproduces"] = bool(o["n_extractable"] == 112 and o["n_p_equal_zero"] == 18 and r1_bonf == 43 and
                           o["round2_bonferroni_pass"] == 25 and o["round2_bh_pass"] == 70 and hl_r1 == 8 and hl_r2 == 5)
    # sensitivity: (i) the confidence-gradient record p = 1.4e-24 is scipy's t-approximation at rho = 1 with n = 5;
    # the exact two-sided Spearman p is 2/5! = 0.0167; (ii) the Stouffer-combined TF/target p is invalid (dependent layers)
    sens = {}
    for tag, p1, p2 in hl_rows:
        q = p2
        if tag == "string_confidence_gradient": q = 2 / 120
        sens[tag] = dict(p_round2=p2, p_sensitivity=q, survives_round2=bool(p2 < thr183), survives_exact_gradient_p=bool(q < thr183))
    n_exact = sum(v["survives_exact_gradient_p"] for v in sens.values())
    n_exact_nostouffer = sum(v["survives_exact_gradient_p"] for k, v in sens.items() if k != "TF_vs_target_AUROC_joint")
    o["headline_round2_sensitivity"] = dict(per_headline=sens, survivors_with_exact_gradient_p=n_exact,
                                            survivors_with_exact_gradient_p_and_no_stouffer=n_exact_nostouffer)
    o["note_headline_key_mismatch"] = ("round-1 headline 'SV2_PPI_copole_string07' is keyed to iter_0010/H02, which is the "
                                       "TRRUST co-pole test (STRING co-pole records are iter_0012/H03 and iter_0013/H01)")
    o["note_round2_headline_remeasure"] = ("round 2 replaced 5 headline p by re-measured values; the TF/target one is a "
                                           "Stouffer combination across layers, which stats3.stouffer() now rejects "
                                           "because layers are dependent")

    # ---------------- (b) all-records audit ----------------
    aall, n_ext_all = audit(df, denom=None)
    out["all_records"] = dict(
        n_total=int(len(aall)), n_extractable=n_ext_all, n_p_equal_zero=int(aall.resolution_limited.sum()),
        bonferroni_alpha=0.05 / n_ext_all, n_bonferroni_pass=int(aall.bonferroni_pass.sum()),
        n_bonferroni_pass_resolution_limited=int((aall.bonferroni_pass & aall.resolution_limited).sum()),
        n_bh_pass=int(aall.bh_pass.sum()),
        n_resolution_limited_that_cannot_pass_bonferroni=int((aall.resolution_limited &
                                                               (1 / (aall.n_perm_used + 1) >= 0.05 / n_ext_all)).sum()),
        by_phase={k: dict(n=int(len(g)), n_extractable=int(g.p_floored.notna().sum()),
                          n_bonferroni=int(g.bonferroni_pass.sum()), n_bh=int(g.bh_pass.sum()))
                  for k, g in aall.groupby(np.where(aall.iter_num <= 63, "iter_1_63", "iter_64_plus"))},
        n_perm_source=dict(from_method=int(aall.n_perm_method.notna().sum()),
                           default_1000=int(aall.n_perm_method.isna().sum())),
    )
    # sensitivity: drop the 3 records of the aborted duplicate run iter_0064_incomplete_* (same IDs re-run in iter_0064)
    a406, n_ext406 = audit(df[~df.aborted_duplicate_run], denom=None)
    out["all_records_excluding_aborted_duplicate_run"] = dict(
        n_total=int(len(a406)), n_extractable=n_ext406, n_p_equal_zero=int(a406.resolution_limited.sum()),
        bonferroni_alpha=0.05 / n_ext406, n_bonferroni_pass=int(a406.bonferroni_pass.sum()), n_bh_pass=int(a406.bh_pass.sum()))
    # sensitivity: Bonferroni over every record (denominator = all records, as round 1 did with 183)
    out["all_records"]["n_bonferroni_pass_denominator_all_records"] = int((aall.p_floored < 0.05 / len(aall)).sum())
    hl = []
    for tag, (it, hid) in HEADLINES.items():
        r = aall[(aall.iteration == it) & (aall.hypothesis_id == hid)].iloc[0]
        hl.append(dict(headline=tag, iteration=it, hypothesis_id=hid, p_raw=r.primary_p_raw, n_perm=int(r.n_perm_used),
                       p_floored=float(r.p_floored), resolution_limited=bool(r.resolution_limited),
                       bh_q_all=float(r.bh_q), bh_pass_all=bool(r.bh_pass),
                       bonferroni_pass_all=bool(r.bonferroni_pass),
                       bonferroni_pass_denominator_409=bool(r.p_floored < 0.05 / len(aall)),
                       bonferroni_pass_183_round2=bool(a183[(a183.iteration == it) & (a183.hypothesis_id == hid)].bonferroni_pass.iloc[0])))
    out["headline_findings_all_records"] = hl
    out["headline_survivors_all_records_bonferroni"] = int(sum(h["bonferroni_pass_all"] for h in hl))
    alpha_all = 0.05 / n_ext_all
    out["headline_survivors_all_records_bonferroni_exact_gradient_p"] = int(sum(
        ((2 / 120) if h["headline"] == "string_confidence_gradient" else h["p_floored"]) < alpha_all for h in hl))
    out["headline_survivors_denominator_409_bonferroni"] = int(sum(h["bonferroni_pass_denominator_409"] for h in hl))
    # verifier addition: denominator = all records AND the exact gradient p (the 409-denominator count above still
    # uses the record's artefactual p = 1.4e-24 for string_confidence_gradient)
    out["headline_survivors_denominator_409_bonferroni_exact_gradient_p"] = int(sum(
        ((2 / 120) if h["headline"] == "string_confidence_gradient" else h["p_floored"]) < 0.05 / len(aall) for h in hl))
    out["headline_survivors_all_records_bh"] = int(sum(h["bh_pass_all"] for h in hl))
    aall.drop(columns=["method", "result_value"]).to_csv(HERE / "audit_409_records.csv", index=False)
    (HERE / "audit_409.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))

if __name__ == "__main__":
    main()
