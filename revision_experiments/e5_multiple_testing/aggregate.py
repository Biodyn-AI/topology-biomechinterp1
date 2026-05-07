#!/usr/bin/env python3
"""E5 — Family-wise multiple-testing audit across all 183 hypotheses.

Aggregates p-values from iter_0001..iter_0063 executor_hypothesis_screen.json,
applies BH-FDR (q=0.05) and Bonferroni (alpha=0.05/N), and reports per-family
and global outcomes. Headline findings are explicitly tracked.

Outputs:
  master_p.csv         per-hypothesis primary p, BH q, Bonferroni decision
  family_summary.csv   per-family pass rates
  headline_audit.md    headline findings + Bonferroni status
  audit_log.json       full extraction log (which p-values came from which strings)
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

def _df_to_md(df: "pd.DataFrame") -> str:
    """Tiny markdown table formatter (no tabulate dependency)."""
    if df.empty:
        return "_(empty)_"
    cols = list(df.columns)
    header = "| " + " | ".join(cols) + " |"
    sep = "|" + "|".join(["---"] * len(cols)) + "|"
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if v is None or (isinstance(v, float) and np.isnan(v)):
                cells.append("")
            elif isinstance(v, float):
                cells.append(f"{v:.3g}")
            else:
                cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([header, sep] + rows)


PROJECT = Path(__file__).resolve().parents[2]
ITER_DIR = PROJECT / "iterations"
OUT_DIR = Path(__file__).resolve().parent

PAPER_ITER_RANGE = range(1, 64)  # paper claims 63 iterations / 183 hypotheses

# Headline findings cited in the paper abstract or main-text claims.
# These are the ones whose Bonferroni status we explicitly track.
HEADLINE_FINDINGS = {
    "SV1_extracellular_enrichment": {
        "iter_id": ("iter_0006", "H02"),  # GO:0005615 OR=6.37 p=0.0003
        "claim": "SV1 separates extracellular vs cytosolic; extracellular OR=6.37, p=2.6e-4",
    },
    "SV1_secretory_pathway_layered": {
        "iter_id": ("iter_0008", "H02"),  # mito_L3 OR=23.25 p=1e-7
        "claim": "Mitochondria/ER/extracellular sequence across layers (paper §2.2)",
    },
    "SV2_PPI_copole_string07": {
        "iter_id": ("iter_0010", "H02"),  # SV2 obs=0.206 vs null=0.122, emp_p=0.000
        "claim": "STRING≥0.7 PPI co-pole on SV2 at 12/12 layers",
    },
    "SV2_TRRUST_copole": {
        "iter_id": ("iter_0011", "H01"),  # 12/12 layers p<0.05
        "claim": "TRRUST TF–target co-pole on SV2; activation 12/12, repression 1/12",
    },
    "TF_vs_target_AUROC_joint": {
        "iter_id": ("iter_0056", "H01"),
        "claim": "Joint SV2-SV7 TF/target AUROC=0.744 mean (max 0.789 at L3)",
    },
    "celltype_marker_AUROC_851": {
        "iter_id": ("iter_0023", "H01"),  # AUROC=0.851 with contamination control
        "claim": "Cell-type marker AUROC=0.851 (12/12 layers p<1e-6)",
    },
    "edge_AUROC_decay": {
        "iter_id": ("iter_0062", "H01"),
        "claim": "Edge-level AUROC peak 0.602 (L0); L0–L8 perm_p≤0.045",
    },
    "string_confidence_gradient": {
        "iter_id": ("iter_0015", "H02"),  # ρ=1.000, p=1.4e-24
        "claim": "STRING confidence gradient on SV2: rho=1.000, p=1.4e-24",
    },
    "bcell_attractor_convergence": {
        "iter_id": ("iter_0042", "H03"),  # BATF & BACH2 convergence
        "claim": "BATF/BACH2 rank decay across L0→L11 (ρ=-0.972, p<1e-4)",
    },
    "celltype_marker_AUROC_853_initial": {
        "iter_id": ("iter_0022", "H02"),  # AUROC=0.853 initial finding
        "claim": "Cell-type marker AUROC=0.853 (12/12 layers p<1e-6)",
    },
    "TF_target_AUROC_cross_seed": {
        "iter_id": ("iter_0057", "H01"),  # cross-seed validation
        "claim": "Joint SV2-7 cross-seed: 0.744/0.753/0.757; p_perm=0.000",
    },
}


# ---------- p-value extraction --------------------------------------------------

# Regexes (tried in priority order). Each captures a numeric p-value.
P_PATTERNS = [
    # Highest priority: empirical / permutation p-values (these are nulls done right)
    ("emp_p_eq",      re.compile(r"\bemp(?:irical)?_p\s*=\s*([0-9eE.+\-]+)")),
    ("emp_p_lt",      re.compile(r"\bemp(?:irical)?_p\s*<\s*([0-9eE.+\-]+)")),
    ("emp_p_le",      re.compile(r"\bemp(?:irical)?_p\s*≤\s*([0-9eE.+\-]+)")),
    ("perm_p_eq",     re.compile(r"\b(?:permutation|perm)_p\s*=\s*([0-9eE.+\-]+)")),
    ("perm_p_lt",     re.compile(r"\b(?:permutation|perm)_p\s*<\s*([0-9eE.+\-]+)")),
    ("perm_p_le",     re.compile(r"\b(?:permutation|perm)_p\s*≤\s*([0-9eE.+\-]+)")),
    # `p_perm=` and `p_perm<` (alternate ordering)
    ("p_perm_eq",     re.compile(r"\bp_perm\s*=\s*([0-9eE.+\-]+)")),
    ("p_perm_lt",     re.compile(r"\bp_perm\s*<\s*([0-9eE.+\-]+)")),
    ("p_perm_le",     re.compile(r"\bp_perm\s*≤\s*([0-9eE.+\-]+)")),
    ("p_emp_paren",   re.compile(r"\bp\s*=\s*([0-9eE.+\-]+)\s*\(\s*emp", re.IGNORECASE)),
    # Fall-back: explicit p=, p<, p≤
    ("p_eq",          re.compile(r"(?<!_)\bp\s*=\s*([0-9eE.+\-]+)")),
    ("p_lt",          re.compile(r"(?<!_)\bp\s*<\s*([0-9eE.+\-]+)")),
    ("p_le",          re.compile(r"(?<!_)\bp\s*≤\s*([0-9eE.+\-]+)")),
    # Fisher / parametric
    ("p_value_eq",    re.compile(r"\bp[-_ ]?value\s*=\s*([0-9eE.+\-]+)")),
]


def _coerce_p(s: str) -> Optional[float]:
    s = s.rstrip(".,;)")
    try:
        v = float(s)
    except ValueError:
        return None
    if not (0 <= v <= 1):
        return None
    return v


def extract_primary_p(result_value: str) -> tuple[Optional[float], Optional[str], list[float]]:
    """Return (primary_p, source_label, all_p_found).

    Strategy: collect all p's from highest-priority pattern that matches; pick the
    minimum. If none found, return (None, None, []).
    """
    if not isinstance(result_value, str):
        return None, None, []

    all_p_found: list[float] = []
    for label, pat in P_PATTERNS:
        matches = pat.findall(result_value)
        if matches:
            ps = [p for p in (_coerce_p(m) for m in matches) if p is not None]
            if ps:
                all_p_found.extend(ps)
                # If this is one of the empirical/permutation patterns, take it as primary.
                if "emp" in label or "perm" in label or "p_perm" in label:
                    return float(np.min(ps)), label, all_p_found

    # No empirical p found; fall back to any extracted p.
    if all_p_found:
        return float(np.min(all_p_found)), "fallback_min_p", all_p_found
    return None, None, []


# ---------- BH-FDR --------------------------------------------------------------

def bh_fdr(p_values: np.ndarray, q: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Return (adjusted q-values, reject_at_q array). Standard BH."""
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([]), np.array([], dtype=bool)
    order = np.argsort(p)
    ranked = p[order]
    # BH adjusted: q_i = min_{j>=i} ( ranked_j * n / j )
    raw = ranked * n / (np.arange(1, n + 1))
    cummin_rev = np.minimum.accumulate(raw[::-1])[::-1]
    adj_sorted = np.clip(cummin_rev, 0, 1)
    adj = np.empty_like(adj_sorted)
    adj[order] = adj_sorted
    reject = adj <= q
    return adj, reject


# ---------- main ----------------------------------------------------------------

def main() -> None:
    rows = []
    extraction_log = []

    for i in PAPER_ITER_RANGE:
        p = ITER_DIR / f"iter_{i:04d}" / "executor_hypothesis_screen.json"
        if not p.exists():
            continue
        try:
            d = json.loads(p.read_text())
        except Exception as e:
            extraction_log.append({"iter": i, "error": str(e)})
            continue
        for h in d.get("hypotheses", []):
            rv = h.get("result_value", "")
            primary_p, src, all_p = extract_primary_p(str(rv))
            rows.append(
                dict(
                    iter=f"iter_{i:04d}",
                    hypothesis_id=h.get("id", ""),
                    name=h.get("name", "")[:120],
                    family=h.get("family", "unknown"),
                    direction=h.get("result_direction", "unknown"),
                    primary_p=primary_p,
                    p_source=src,
                    n_p_in_string=len(all_p),
                    result_value=str(rv)[:280],
                )
            )

    df = pd.DataFrame(rows)
    print(f"Loaded {len(df)} hypothesis records from iters 1..63.")

    # Bonferroni and BH applied across only those with extractable p.
    extractable = df[df["primary_p"].notna()].copy()
    n_extract = len(extractable)
    n_total = len(df)
    print(f"  Extractable p: {n_extract}/{n_total}")
    print(f"  Direction breakdown: ")
    print(df["direction"].value_counts().to_string())

    # Global Bonferroni at alpha = 0.05 / 183 (paper-cited count, used as denominator
    # even when extraction yields fewer to be transparent; we report both).
    alpha = 0.05
    bonf_denom_paper = n_total  # 183 by construction
    bonf_thresh_paper = alpha / bonf_denom_paper
    bonf_denom_strict = n_extract
    bonf_thresh_strict = alpha / max(bonf_denom_strict, 1)

    # BH on extractable subset only (BH undefined for missing p)
    if n_extract:
        bh_q, bh_reject = bh_fdr(extractable["primary_p"].to_numpy(), q=alpha)
        extractable["bh_q"] = bh_q
        extractable["bh_reject_q05"] = bh_reject
        extractable["bonferroni_paper"] = (
            extractable["primary_p"] <= bonf_thresh_paper
        )
        extractable["bonferroni_strict"] = (
            extractable["primary_p"] <= bonf_thresh_strict
        )
    else:
        for col in ("bh_q", "bh_reject_q05", "bonferroni_paper", "bonferroni_strict"):
            extractable[col] = np.nan

    # Merge back
    df_out = df.merge(
        extractable[
            [
                "iter",
                "hypothesis_id",
                "bh_q",
                "bh_reject_q05",
                "bonferroni_paper",
                "bonferroni_strict",
            ]
        ],
        on=["iter", "hypothesis_id"],
        how="left",
    )
    df_out.to_csv(OUT_DIR / "master_p.csv", index=False)

    # Per-family BH (within-family correction)
    fam_rows = []
    for fam, grp in extractable.groupby("family"):
        if grp.empty:
            continue
        bh_q_fam, bh_rej_fam = bh_fdr(grp["primary_p"].to_numpy(), q=alpha)
        bonf_fam = grp["primary_p"].to_numpy() <= alpha / len(grp)
        fam_rows.append(
            dict(
                family=fam,
                n_hypotheses=len(grp),
                n_extractable_p=int(grp["primary_p"].notna().sum()),
                n_bh_pass_within=int(bh_rej_fam.sum()),
                n_bonferroni_pass_within=int(bonf_fam.sum()),
                n_bh_pass_global=int(grp["bh_reject_q05"].sum()) if "bh_reject_q05" in grp.columns else 0,
                n_bonferroni_pass_paperdenom=int((grp["primary_p"] <= bonf_thresh_paper).sum()),
            )
        )
    fam_df = pd.DataFrame(fam_rows).sort_values("n_hypotheses", ascending=False)
    fam_df.to_csv(OUT_DIR / "family_summary.csv", index=False)

    # Headline-finding audit
    headline_rows = []
    for tag, info in HEADLINE_FINDINGS.items():
        it, hid = info["iter_id"]
        hit = df_out[(df_out["iter"] == it) & (df_out["hypothesis_id"] == hid)]
        if hit.empty:
            headline_rows.append(
                dict(
                    headline=tag,
                    claim=info["claim"],
                    iter=it,
                    hypothesis_id=hid,
                    primary_p=None,
                    bh_q=None,
                    bonferroni_pass=None,
                    found=False,
                )
            )
            continue
        r = hit.iloc[0]
        headline_rows.append(
            dict(
                headline=tag,
                claim=info["claim"],
                iter=it,
                hypothesis_id=hid,
                primary_p=r["primary_p"],
                bh_q=r["bh_q"],
                bonferroni_pass=bool(r["bonferroni_paper"]) if pd.notna(r["bonferroni_paper"]) else None,
                found=True,
            )
        )
    head_df = pd.DataFrame(headline_rows)
    head_df.to_csv(OUT_DIR / "headline_audit.csv", index=False)

    # Markdown summary
    lines = []
    lines.append("# E5 — Multiple-testing audit across the 183 paper hypotheses\n")
    lines.append("## Global counts\n")
    lines.append(f"- Total hypotheses (iter 1..63 inclusive): **{n_total}** (paper claim: 183 — match: {'yes' if n_total == 183 else f'no, {n_total}'})")
    lines.append(f"- Hypotheses with extractable primary $p$: **{n_extract}** ({100*n_extract/max(n_total,1):.1f}%)")
    lines.append("")
    lines.append("## Direction breakdown (paper's own classification)")
    lines.append("```")
    lines.append(df["direction"].value_counts().to_string())
    lines.append("```\n")
    lines.append("## Multiple-testing thresholds")
    lines.append(f"- Bonferroni (paper denom $N=183$): $\\alpha = 0.05/183 = {bonf_thresh_paper:.3e}$")
    lines.append(f"- Bonferroni (extractable denom $N={n_extract}$): $\\alpha = 0.05/{n_extract} = {bonf_thresh_strict:.3e}$")
    lines.append("- BH-FDR at $q=0.05$ across the extractable set.\n")
    n_bh_pass = int(extractable["bh_reject_q05"].sum()) if n_extract else 0
    n_bonf_paper = int((extractable["primary_p"] <= bonf_thresh_paper).sum()) if n_extract else 0
    n_bonf_strict = int((extractable["primary_p"] <= bonf_thresh_strict).sum()) if n_extract else 0
    lines.append("## Global pass counts")
    lines.append(f"- BH-FDR $q\\le0.05$: **{n_bh_pass} / {n_extract}**")
    lines.append(f"- Bonferroni (paper denom): **{n_bonf_paper} / {n_extract}**")
    lines.append(f"- Bonferroni (strict denom): **{n_bonf_strict} / {n_extract}**\n")
    lines.append("## Per-family\n")
    lines.append(fam_df.pipe(_df_to_md))
    lines.append("\n## Headline-finding audit\n")
    lines.append(head_df.pipe(_df_to_md))
    lines.append("\n## Notes / caveats")
    lines.append("- p-values were extracted by regex from the executor's free-form `result_value`. "
                 "Where the executor reported empirical/permutation p-values, those were preferred over GO Fisher p-values.")
    lines.append("- 'Extractable' means at least one numeric p-value was parseable from the result_value string. "
                 "Hypotheses with `direction='negative'` or `'mixed'` may have no headline p-value, so missingness is informative.")
    lines.append("- Bonferroni with the paper's denominator (N=183) is the strictest defensible global correction, "
                 "and is the threshold against which abstract claims are audited.")
    lines.append("- Within-family BH/Bonferroni is reported because hypotheses within the same family are conceptually related "
                 "and external readers may prefer family-wise correction over global.")
    (OUT_DIR / "headline_audit.md").write_text("\n".join(lines))

    (OUT_DIR / "audit_log.json").write_text(
        json.dumps(
            dict(
                n_iterations_scanned=len(list(PAPER_ITER_RANGE)),
                n_total=n_total,
                n_extractable=n_extract,
                bonferroni_paper_threshold=bonf_thresh_paper,
                bonferroni_strict_threshold=bonf_thresh_strict,
                bh_q=alpha,
                n_bh_pass=n_bh_pass,
                n_bonferroni_paper_pass=n_bonf_paper,
                n_bonferroni_strict_pass=n_bonf_strict,
                extraction_errors=extraction_log,
            ),
            indent=2,
        )
    )
    print("\nWrote:")
    print(f"  {OUT_DIR/'master_p.csv'}")
    print(f"  {OUT_DIR/'family_summary.csv'}")
    print(f"  {OUT_DIR/'headline_audit.csv'}")
    print(f"  {OUT_DIR/'headline_audit.md'}")
    print(f"  {OUT_DIR/'audit_log.json'}")


if __name__ == "__main__":
    main()
