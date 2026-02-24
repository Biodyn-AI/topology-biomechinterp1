# Brainstormer Structured Feedback — iter_0019

## Gate Status
`passed_min_research_gate: true` — all three hypotheses produced quantitative results with positive outcomes.

---

## Iteration Quality Assessment

### What worked
- **H01** is the strongest result this project has produced: Geneformer static token embeddings independently encode PPI geometry at MW_p=7.8e-127, and cross-model Spearman_abs=0.446 with scGPT SV2 is a genuine alignment signal. The sign-flip between layers (L00: -0.842, L11: +0.702) is a real geometric phenomenon worth exploiting.
- **H02** delivers the clearest mechanistic insight to date: the gene manifold is genuinely multi-dimensional. Monotonic 0.86x→2.18x enrichment gradient across 0–3 co-polar axes confirms that PPI proximity is distributed across multiple orthogonal SVD components, not concentrated in a single axis.
- **H03** delivers the cleanest dissociation result: SVD geometry → PPI, attention → TF regulation (direction-agnostic). The fact that STRING pairs are NOT enriched in attention (p=0.084 NS) while TRRUST pairs are (~2x, p<0.01) is a publishable mechanistic claim.

### What needs immediate follow-up
- H01: Static embedding result is strong but incomplete. Contextualized GF inference on the same lung dataset is required to close the cross-model story. Without it, the alignment is between pre-training token distributions, not true representation geometry.
- H02: Shuffle null is missing. The 2.18x enrichment at 3-axis co-polarity has no validated control. Must add axis-permutation null before this can be a publication-grade claim.
- H03: The attention-SVD dissociation has been established but not exploited. A joint predictor combining both scores has not been tested.

### What is stale
- All SVD signed-direction tests: definitively negative across iter_0013–0017. No rescue value.
- GO term sweeps: negative in iter_0011. No topological signal found in GO membership.
- 1D SV2 distance ranking: 1.2x enrichment (iter_0018) — superseded entirely by multi-axis composite (2.18x). Retire.
- Random Gaussian null: confirmed negative in iter_0018. No need to retest.

---

## Strategic Assessment

The project now has three orthogonal but potentially integrable results:
1. **Cross-model manifold alignment** (SVD geometry is model-invariant)
2. **Multi-dimensional PPI manifold** (3-axis composite outperforms 1D)
3. **Attention/SVD dissociation** (two complementary geometric encodings)

The natural next step is to combine and validate these into a unified geometric model of the scGPT gene representation space. The cross-model validation (GF) would provide external support for the SVD geometry claims.

The weakest gap: no topology tests (persistent homology, Mapper) have been attempted despite this being the project's original mandate. This iteration should re-prioritize topological methods, which are now better motivated by the confirmed multi-dimensional structure.
