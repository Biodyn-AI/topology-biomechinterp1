# Next Steps: iter_0045 → iter_0046

## Key finding this iteration
Full-manifold TwoNN on all 4941 gene vectors shows monotone ID compression: L0=32.57 → L11=18.05
(−44.6%), far below Gaussian null (~123D). Every layer compresses. This is the main carry-forward.

## Retirements
- **H01 (kNN lineage purity with n=1 per lineage)**: Retired. Underpowered by design. Do not retry
  without at least 4 in-vocab genes per lineage.

## Priority 1: Cross-run replication of H02 (ID compression)
- Compute TwoNN on cycle1 embeddings [12,4803,512] at same layer points
- Compute on cycle4_seed43 and cycle4_seed44 for multi-seed consistency
- Also test: does ID compression rate differ across seeds? (stability check)
- Command: conda run -n subproject40-topology python [replication script]

## Priority 2: Bootstrap stability of ID compression
- Resample n_samples in {500, 1000, 1500, 2000, 3000} at L0 and L11
- Confirm monotone compression is not subsample-size artifact
- 95% CI for ID at each layer from 10 bootstrap resamples

## Priority 3: Functional gene subset ID (within-context TwoNN)
- For B-cell markers, TFs, and random control sets: compute local ID (TwoNN on n=50 nearest neighbors
  in full embedding space) at each layer
- Does B-cell subset ID compress faster than random?

## Priority 4: CD28 trajectory replication + expansion
- Verify CD28 rank-to-T-cell-centroid trajectory (L0=595→L6=86→L11=149) in cycle1 and cycle4_seed43
- Add FOS, JUN, JUNB to circuit if in vocab (AP1 transcription factor family)
- Also check: does CD28 U-shape peak at L6 in other cycles?

## Recommended portfolio for iter_0046
1. **H1 (intrinsic_dimensionality, new_method)**: Bootstrap stability of ID compression (n_samples
   sensitivity) + cross-seed replication
2. **H2 (intrinsic_dimensionality, new_method)**: Local ID for functional gene subsets vs random
3. **H3 (manifold_distance, refinement)**: CD28 trajectory cross-run validation + AP1 expansion
