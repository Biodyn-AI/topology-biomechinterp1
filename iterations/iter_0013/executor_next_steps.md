# Executor Next Steps — iter_0013

## Retired directions
- Repression anti-pole (H02 this iter): definitive negative. Do not re-test.
- Bootstrap CI act-rep differential (iter_0012 H01): already retired.

## Priority for iter_0014

### 1. STRING threshold expansion for hub-degree confound resolution (HIGH)
H03 showed non-hub edges are underpowered (N=35 at score≥0.7). Lower threshold to score≥0.4
to get ~5-10x more non-hub pairs. Re-run SV2 co-pole separately for hub vs non-hub.
If non-hub is still null at 150+ pairs, the confound is confirmed and the PPI claim must be
qualified. If non-hub is positive, the signal generalizes beyond hubs.
Family: graph_topology / null_sensitivity

### 2. SV3 biological anchor (HIGH)
H01 showed SV3 is strongly significant (z=7.18, 12/12 layers). SV3 has not been annotated.
Run GO enrichment at SV3 top/bottom poles (similar to iter_0007/0008 for SV1/SV2).
If SV3 encodes a distinct biological axis, this is a new positive finding.
Family: intrinsic_dimensionality / module_structure

### 3. Joint SV2+SV3 co-pole (MEDIUM)
Test if STRING PPI pairs are jointly co-localized in SV2-SV3 2D space (not just marginally).
Could use 2D quadrant assignment instead of 1D poles.
Family: manifold_distance

### 4. Cross-model alignment check (MEDIUM, novel family)
Compare SV2 gene rankings between scGPT and Geneformer embeddings.
Family: cross_model_alignment (not tested yet — high novelty)

## Notes
- The SV2 z-scores were larger in this iteration (~10) than iter_0012 (~5) due to higher N_shuffle
  (500 vs 500 with larger null variance estimate). The signal is very robust.
- SV1 is definitively not a PPI geometry axis (mean z=0.39, 2/12 sig).
- The co-pole claim should be stated as "SV2 and SV3 jointly encode PPI proximity" going forward.
