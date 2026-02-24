# Brainstormer Structured Feedback — iter_0057

**Date**: 2026-02-23
**Research gate**: PASSED

---

## Outcome Assessment

### H01 — Joint 6D Cross-Seed Validation: STRONG POSITIVE
Mean AUROC=0.751 across 3 seeds, p_perm=0.000 at all 12 layers. This converts Claim 50 from "promising" to fully validated. The key architectural fact established: the *combined span* of SV2-7 encodes TF/target identity robustly, even though individual subspace orientations (SV5-7 directions) drift 13–41° across seeds.

Actionable follow-ups unlocked by this result:
- Which TF genes are consistently highest-ranked across seeds? (boundary consistency test)
- Can we explain *why* early layers are better discriminators geometrically?
- Does this generalize beyond immune cell type?

### H02 — Basis Permutation: DEFINITIVELY NEGATIVE
Ruled out cleanly. SV5-7 directions drift substantially. Main SV2-4 is near-orthogonal to seed43 SV5-7 (mean PA=74.8°). This is a clean negative — retire the specific mechanism, but the observation that "combined span is stable while axes are not" is scientifically interesting and itself worth reporting.

### H03 — SV Energy Fraction vs AUROC: NOVEL STRUCTURAL FINDING
rho=-0.93 is unexpectedly strong. Early layers: low fraction of variance in SV2-7 but high AUROC. This means early-layer discriminative information is *distributed across many singular components*, not concentrated. This connects to effective rank / participation ratio of the embedding — a direction that has not been directly tested yet.

The negative sign is the key: concentration of variance in top SVs is *inversely* related to regulatory discriminability. This pattern should be tested in cross-seed to confirm it's not a single-run artifact.

---

## Strategic Observations

1. The core validated finding (6D SV2-7 TF/target discrimination, L0-L3 peak, L8 minimum) is now seed-stable. The axis of inquiry should now pivot toward:
   - **Mechanistic explanations** of the layer profile (why peak at L2-L3? why minimum at L8?)
   - **Biological interpretation** (which TF families drive the boundary? does GO enrichment separate the SV2-7 projections?)
   - **Generalization probes** (cell types, models, gene sets)

2. The rho=-0.93 finding (H03) suggests effective rank is a better layer-descriptor than energy concentration. This is a clean hypothesis to extend.

3. Untested high-value hypotheses from iter_0056 portfolio that should be carried forward:
   - TRRUST Graph Laplacian spectral alignment (B3)
   - PH on TF vs target gene sets (B1)
   - OOD generalization train-on-immune test-on-lung (C3)
   - TwoNN intrinsic dimension TF vs target (B2)

4. No hypothesis in the current portfolio has addressed *which genes* sit at the TF/target decision boundary in 6D space. This is the most direct path to biological interpretation.
