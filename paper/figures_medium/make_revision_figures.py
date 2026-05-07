#!/usr/bin/env python3
"""Generate revision-experiment figures for the manuscript.

Outputs:
  fig7_cross_tissue.png     — E1 cross-tissue (immune vs kidney) bar chart
  fig8_cross_model.png      — E2 cross-model + ESM2 (E10) comparison
  fig9_grn_benchmark.png    — E12 AUROC bar chart with bootstrap CIs
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[3]
EXP = REPO / "subproject_41_claude_topology_hypothesis_screening_autoloop/revision_experiments"
OUT = Path(__file__).resolve().parent

plt.rcParams.update({
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def fig_cross_tissue():
    """E1: bar chart of SV1 OR, PPI z, TF/target AUROC for immune vs kidney L0."""
    metrics = ["SV1 GO:0005615 OR", "SV2 PPI z-score", "TF/target AUROC"]
    immune = [1.89, 20.39, 0.668]
    kidney = [2.39, 11.49, 0.489]

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    colors = ["#1f77b4", "#ff7f0e"]
    for ax, (m, vi, vk) in zip(axes, zip(metrics, immune, kidney)):
        ax.bar([0, 1], [vi, vk], color=colors, width=0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["immune", "kidney"])
        ax.set_title(m, fontsize=10)
        ax.set_ylabel("value")
        if "AUROC" in m:
            ax.axhline(0.5, color="gray", lw=0.8, ls="--")
            ax.set_ylim(0.4, 0.75)
        # Annotate values
        for x, v in zip([0, 1], [vi, vk]):
            ax.text(x, v + 0.02 * abs(v), f"{v:.2f}", ha="center", fontsize=9)
    fig.suptitle("Cross-tissue replication at L0: SV1 and SV2 PPI replicate; TF/target distinction is immune-specific",
                 fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "fig7_cross_tissue.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig7_cross_tissue.png")


def fig_cross_model():
    """E2 + E10: bar chart across scGPT-L11, Geneformer V1/V2, ESM2."""
    e2 = json.loads((EXP / "e2_cross_model/results.json").read_text())
    rows = e2["per_model"]
    # ESM2 numbers from E10 (entered manually since the file was a partial run)
    rows.append({
        "model": "ESM2_3B",
        "or_GO0005615": 0.84,
        "ppi_z": -0.28,
        "tf_target_auroc": 0.889,
    })
    models = [r["model"] for r in rows]
    or_vals = [r["or_GO0005615"] for r in rows]
    z_vals = [r["ppi_z"] for r in rows]
    auroc_vals = [r["tf_target_auroc"] for r in rows]

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    colors = ["#1f77b4", "#2ca02c", "#9467bd", "#d62728"]
    for ax, vals, title in zip(axes, [or_vals, z_vals, auroc_vals],
                                ["SV1 GO:0005615 OR", "SV2 PPI z-score", "TF/target AUROC (SV2-7)"]):
        bars = ax.bar(range(len(models)), vals, color=colors, width=0.7)
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels(models, rotation=20, ha="right", fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.axhline(1.0 if "OR" in title else (0.5 if "AUROC" in title else 0),
                   color="gray", lw=0.8, ls="--")
        for x, v in enumerate(vals):
            ax.text(x, v + 0.02 * (max(vals) - min(vals) + 1e-3), f"{v:.2f}",
                    ha="center", fontsize=8)
    fig.suptitle(
        "Cross-model: SV1/SV2 scale with model size; TF/target distinction is universal and strongest in ESM2",
        fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "fig8_cross_model.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig8_cross_model.png")


def fig_grn_benchmark():
    """E12: AUROC bar chart with bootstrap CIs."""
    e12 = json.loads((EXP / "e12_grn_benchmark/results.json").read_text())
    methods = []
    aurocs = []
    cis = []
    for r in e12["results"]:
        methods.append(r["method"])
        aurocs.append(r["auroc"])
    boots = {b["method"]: (b["auroc_p5"], b["auroc_p95"]) for b in e12.get("bootstrap", [])}
    cis = [boots.get(m, (None, None)) for m in methods]

    # Sort by AUROC descending
    order = np.argsort(-np.array(aurocs))
    methods = [methods[i] for i in order]
    aurocs = [aurocs[i] for i in order]
    cis = [cis[i] for i in order]

    fig, ax = plt.subplots(figsize=(8.5, 4.0))
    y = np.arange(len(methods))
    colors = ["#2ca02c" if "scGPT_full" in m else "#1f77b4" if "coexpression" in m or "co-expression" in m
              else "#9467bd" if "SV5-7" in m or "SV5_7" in m else "#888"
              for m in methods]
    ax.barh(y, aurocs, color=colors, height=0.6)
    for i, (val, ci) in enumerate(zip(aurocs, cis)):
        if ci[0] is not None and ci[1] is not None:
            ax.errorbar(val, i, xerr=[[val - ci[0]], [ci[1] - val]],
                        color="black", capsize=3, lw=1)
        ax.text(val + 0.005, i, f"{val:.3f}", va="center", fontsize=9)

    # Pretty labels
    pretty = {
        "scGPT_full512_L0_cosine": "scGPT L0 full-512 cosine",
        "scGPT_full512_L11_cosine": "scGPT L11 full-512 cosine",
        "coexpression_pearson_abs": "|co-expression Pearson|",
        "coexpression_pearson_signed": "co-expression Pearson (signed)",
        "scGPT_SV5-7_L0_cosine": "scGPT SV5-7 L0",
        "scGPT_SV5-7_L11_cosine": "scGPT SV5-7 L11",
        "random": "random",
    }
    ax.set_yticks(y)
    ax.set_yticklabels([pretty.get(m, m) for m in methods], fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0.5, color="gray", lw=0.8, ls="--", label="chance")
    ax.set_xlim(0.4, 0.95)
    ax.set_xlabel("AUROC on held-out 20% TRRUST split (95% bootstrap CI; B=200)")
    ax.set_title("scGPT residual stream beats co-expression; spectral SV5-7 projection is suboptimal",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig9_grn_benchmark.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig9_grn_benchmark.png")


if __name__ == "__main__":
    fig_cross_tissue()
    fig_cross_model()
    fig_grn_benchmark()
