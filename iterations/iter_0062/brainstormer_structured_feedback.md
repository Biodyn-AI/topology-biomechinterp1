# Brainstormer Structured Feedback — iter_0062

**Date**: 2026-02-23
**Gate status**: PASSED (passed_min_research_gate = true)

---

## Assessment of iter_0062 Results

### M1: cycle4_immune edge-level AUROC (SV5-7)
**Signal quality**: Genuine, modest. AUROC peaks at 0.602 (L0), declines monotonically to chance at L9–L11. Spearman(AUROC, layer) = −0.958, p < 1e-6. Cross-seed agreement tight at L3 (range < 0.01), acceptable at L2. **This is a publication-quality finding**: edge-level TF-target co-embedding is geometrically encoded in early-to-mid layers of scGPT, not in late layers.

SV2-4 sub-null AUROC (<0.50) is an unexpected and interpretable phenomenon: TF-target pairs may be actively repelled in SV2-4 space (the PPI/subcellular compartment axes from iter_0014-0017). This is the most interesting unexplained anomaly in the current data.

### M2: AUROC mechanism diagnostic
Sparsity eliminated as confound. Nz_count = 2039 at all 12 layers. Depth is the real driver. This is a clean negative control that upgrades M1 from "interesting" to "defensible." Worth a paragraph in the paper.

### C1: Cross-cycle Procrustes (design flaw, negative)
The flaw is fundamental: cosine similarity is rotation-invariant, so Procrustes rotation cannot affect AUROC. The secondary observation (c1_own AUROC = 0.54–0.56 vs c4 = 0.50–0.60) is suggestive of cycle-specific geometry but is confounded by small shared-gene count (80). This hypothesis needs complete redesign with a non-rotation-invariant metric.

---

## Key Open Questions Surfaced

1. **SV2-4 sub-null AUROC**: Why are TF-target pairs repelled in the PPI/compartment axes? This contradicts prior iter results where TF-target activation pairs showed positive STRING co-pole enrichment. The two subspaces may encode orthogonal biological dimensions, and TF-target pairs cluster in SV5-7 precisely because they do NOT cluster by subcellular compartment.

2. **L0 peak vs L3 classification peak**: Edge-level similarity peaks at L0; TF/target class-level AUROC (prior iters) peaked at L2–L3. These measure different structure: pairwise co-embedding vs class geometry. The L0 peak is consistent with raw transcriptional co-expression being most preserved at input layers.

3. **Cross-cycle geometry**: The 80 shared gene finding means cycle1 and cycle4 share only a small regulatory vocabulary in embedding space. This suggests training cycle induces substantial gene-position rearrangement, not just metric scaling.

---

## Trajectory Assessment

The research is now in a late-stage characterization phase. The core claims (PPI/regulatory manifold organization, multi-axis structure, depth encoding timeline) are established. The current thread (edge-level AUROC in cycle4) is adding quantitative detail to these claims. Two productive directions remain:
1. **Mechanistic depth**: Why does the signal decay by L9? Is it cell-type specialization outcompeting regulatory co-expression?
2. **Generalization**: Do these signals transfer across cycles, conditions, or model families?

The SV2-4 anti-correlation is underexplored and could be the highest-novelty finding in the portfolio.
