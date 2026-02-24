# Brainstormer Structured Feedback — iter_0025

**Date:** 2026-02-22
**Research gate passed:** YES

---

## Iteration Assessment

### What was strong

**H01 (multi-predictor joint model)** is the most important result to date. Three independent biological databases (STRING PPI, Dorothea regulatory confidence, GO CC compartment) each contribute non-redundant geometric signal with VIF < 2. Total partial R² = 1.7% is meaningful given 512-dimensional noise. This demonstrates convergent validity: functionally distinct annotation sources converge on the same embedding geometry.

**H02 (layer timeline)** is a mechanistic finding with publication-grade narrative: structural information (compartment, L5) encodes before functional process information (GO BP, L7), suggesting a biologically plausible representational hierarchy within the transformer. The 7-layer span across anchors is strong.

**H03 (chromosomal control)** did its job perfectly. AUROC 0.515 vs functional signals at 0.60–0.68 is a 10× gap. This is the cleanest negative control in the project so far.

### Critical gaps now exposed

The project has characterized *what* is encoded in scGPT's manifold. The next scientific questions are:
1. **Generalizability**: Is this scGPT-specific or common to gene LLMs? (Geneformer replication)
2. **Mechanism**: Which transformer components drive the layer-resolved biology? (Attention attribution)
3. **Completeness**: What fraction of total biological structure is captured? Is 1.7% R² a ceiling or a floor?
4. **Topology**: Do functionally coherent gene sets form topologically stable clusters (persistent homology), not just pairwise distance signals?
5. **Causal direction**: Do regulatory interactions (Dorothea) encode direction (activator vs. repressor), not just confidence magnitude?

### What to retire

See roadmap file. Short version: no direction needs hard retirement — the project is in a productive phase. Some directions (chromosomal confound, GO BP vs GO CC redundancy) are now resolved and need no further work.

### Iteration quality

Clean execution, no blockers, three meaningful experiments. Iteration was cost-efficient and built directly on prior positive results. The next iteration should be more ambitious.

---

## Numerical Anchors for Next Hypotheses

| Signal | Metric | Value |
|--------|--------|-------|
| STRING PPI | Spearman ρ (L8) | +0.152 (univariate) / +0.098 (partial) |
| Dorothea confidence | AUROC (L7) | 0.679 |
| GO CC | Spearman ρ (L5) | +0.124 |
| GO BP | Spearman ρ (L7) | +0.083 |
| TRRUST activation | AUROC (L5) | 0.659 |
| Joint model | Total partial R² (L8) | 0.01695 |
| Chromosomal proximity | AUROC (L8) | 0.515 (null validated) |

Current explained variance ceiling from known anchors: ~1.7% partial R². Substantial unexplained structure remains.
