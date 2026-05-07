#!/usr/bin/env python3
"""E6 — Family-level retire/validate ratios with binomial selection-bias test.

Reviewer 1 (#4): "the autonomous two-agent hypothesis-screening loop lacks
sufficient methodological detail and may introduce uncharacterised selection
bias."

This script formalises the selection-bias check by computing the binomial
probability of seeing the observed positive-rate per family under a 50/50 null,
and reporting the family-level retire/validate decisions made by the
brainstormer agent. High positive rates indicate the brainstormer was effective
at retiring unproductive branches early — but also that the family-level
positive rate over-reports the rate at which a *random* hypothesis would be
positive.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
ITER_DIR = PROJECT / "iterations"
OUT = Path(__file__).resolve().parent


def binom_two_tailed(k: int, n: int, p: float = 0.5) -> float:
    """Two-tailed binomial p-value (exact)."""
    if n == 0:
        return 1.0
    obs = abs(k - n * p)
    total = 0.0
    for i in range(n + 1):
        prob = math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
        if abs(i - n * p) >= obs - 1e-9:
            total += prob
    return min(1.0, total)


def main() -> None:
    families: dict[str, dict] = {}
    decisions: dict[str, dict] = {}
    extends: dict[str, list[tuple[str, str, str]]] = {}  # finding -> [(iter,id,direction)]

    for i in range(1, 64):
        p = ITER_DIR / f"iter_{i:04d}" / "executor_hypothesis_screen.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        for h in d.get("hypotheses", []):
            fam = h.get("family", "unknown")
            direction = h.get("result_direction", "unknown")
            decision = h.get("decision", "unknown")
            ext = h.get("extends_finding", "")
            slot = families.setdefault(
                fam, {"positive": 0, "negative": 0, "mixed": 0, "inconclusive": 0, "total": 0}
            )
            slot[direction] = slot.get(direction, 0) + 1
            slot["total"] += 1
            decisions.setdefault(fam, {})[decision] = decisions.setdefault(fam, {}).get(decision, 0) + 1
            if ext:
                extends.setdefault(ext, []).append((f"iter_{i:04d}", h.get("id", ""), direction))

    rows = []
    for fam, d in sorted(families.items(), key=lambda x: -x[1]["total"]):
        n = d["total"]
        pos = d["positive"]
        neg = d["negative"]
        binom_p = binom_two_tailed(pos, n, 0.5)
        rows.append(
            dict(
                family=fam,
                n=n,
                positive=pos,
                negative=neg,
                mixed=d["mixed"],
                inconclusive=d["inconclusive"],
                positive_rate=round(pos / max(n, 1), 3),
                binomial_p_50null=round(binom_p, 4),
                p_pos_half_or_higher=round(pos / max(n, 1), 3) >= 0.5,
            )
        )
    fam_df = pd.DataFrame(rows)
    fam_df.to_csv(OUT / "family_retire_validate.csv", index=False)

    # extends_finding distribution
    ext_rows = []
    for finding, items in extends.items():
        n = len(items)
        pos = sum(1 for _, _, dir_ in items if dir_ == "positive")
        ext_rows.append(dict(finding=finding, n=n, positive=pos, positive_rate=round(pos / max(n, 1), 3)))
    ext_df = pd.DataFrame(ext_rows)
    if not ext_df.empty:
        ext_df = ext_df.sort_values("n", ascending=False)
    ext_df.to_csv(OUT / "extends_finding_distribution.csv", index=False)

    print("Family retire/validate:")
    print(fam_df.to_string(index=False))
    print("\nExtends-finding distribution:")
    print(ext_df.to_string(index=False))


if __name__ == "__main__":
    main()
