# Brainstormer Structured Feedback — iter_0061

**Date**: 2026-02-23
**Research gate**: PASSED (3 hypotheses tested, all with permutation-corrected nulls)

---

## Iteration Quality Assessment

**Strengths:**
- H03 FFL retirement is correct and clean. The null t_mean≈0.50 insight is a genuine methodological contribution — all future geometric "betweenness" or "ordering" tests must use this corrected null, not t_mean>0.
- H01 CKA is well-executed: three null types (consecutive, cross-seed, vs-L0) give a complete picture.
- H02 correctly generalizes the HIF1A finding downward — closes the branch with proper N.

**Key emergent finding (not hypothesis-driven):**
Cross-seed CKA drops sharply at L10-L11 (0.870→0.779), an unusually large 9-point drop after a smooth monotonic decline. This suggests deep layers are either (a) encoding sample-specific features, or (b) operating in a less-constrained regime where random initialization diverges. This is the most actionable side-finding from iter_0061.

**Methodological lesson solidified:**
For any geometric "ordering" or "betweenness" test in high-dimensional embedding space, the permutation null is not 0 — it's ~0.50 (by symmetry). All prior FFL results using parametric t>0 are invalidated. This is now locked.

---

## Status of Active Research Threads

| Thread | Status | Evidence |
|--------|--------|----------|
| TF-target SV5-7 separation (AUROC~0.57-0.74) | **Active — underpowered** | Cycle1: 288 pairs; cycle4 has 735 — test at full scale |
| Deep-layer (L10-L11) geometric divergence | **New — unexplored** | CKA drop 0.87→0.78; no topology test yet |
| HIF1A L8 flip | **Active — n=1 outlier** | Biologically interpretable but needs neighborhood test |
| FFL geometric ordering | **Retired** | Permutation null kills this |
| L8 CKA boundary | **Retired** | Consecutive CKA uniform 0.977-0.983 |
| Dual-role gene ambiguity | **Retired** | All consistently target-like |
| Stable low-CV TF anchors (STAT4, BACH2, ZEB1, NFATC2) | **Active — unexplored geometry** | Found in iter_0059, not followed up |
| eff_rank → AUROC confound | **Closed** | Layer-depth drives both |
| Cross-cycle transfer | **Unexplored** | cycle4 vs cycle1 classifier transfer |

---

## Critical Decision: Next Priority

The most important unresolved question: **Is the TF-target AUROC signal real biological structure or a layer-depth artifact?**

Evidence so far: partial Spearman(AUROC, eff_rank | layer) = -0.045 (iter_0059/H03). The confound is layer-depth itself, not eff_rank. But we haven't tested sparsity or gene-count as confounds at layer level for cycle4. This should be done before claiming any AUROC as biologically meaningful.

Second priority: Scale up to cycle4_immune 735 pairs. The 288-pair cycle1 dataset underpowers all AUROC estimates.

Third priority: Topological characterization of L10-L11 divergence.
