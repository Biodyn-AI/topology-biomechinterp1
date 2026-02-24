# Next Steps — iter_0061

## Retirements confirmed this iteration
- FFL geometric ordering branch (all subspaces, permutation-corrected) → RETIRED
- Dual-role gene L8 margin flip → closed
- L8 CKA boundary → falsified

## Key insight from this iteration
**Geometric baseline correction**: Random gene triplets in embedding space have t_mean ~0.50 (not 0), because by geometry ~50% of random triplets show the intermediate gene "between" the other two. All parametric tests against 0 in FFL/geometry experiments are misleading. Future geometric tests must use permutation-corrected baselines.

## Top 3 for iter_0062

### H01: Deep-layer (L10-L11) cross-seed divergence — persistent homology
**Rationale**: Cross-seed CKA drops from 0.93 at L7 to 0.78 at L11. This is the most significant pattern from this iteration. Are deep-layer representations topologically less consistent (different Betti curves)?
**Method**: For each layer (all 12): compute Vietoris-Rips persistence on the k-NN neighborhood of TRRUST-relevant genes. Compare Betti-0 and Betti-1 curves between seeds (Wasserstein distance or L2 of Betti curves). Test: does between-seed Betti-curve distance increase at L10-L11?
**Family**: persistent_homology
**Cost**: medium (PH is slower but N~300 genes is tractable)

### H02: Local intrinsic dimensionality by layer and by gene neighborhood
**Rationale**: Cross-seed CKA declining could reflect that deep layers expand the intrinsic dimensionality. Test: local ID via TwoNN or participation ratio within k-NN neighborhoods of TF genes.
**Method**: For k=10, 20, 50: build k-NN graph of TRRUST genes at each layer. Estimate local ID via participation ratio of pairwise distances. Compare L0-L7 vs L8-L11.
**Family**: intrinsic_dimensionality
**Cost**: low

### H03: Novel — kNN graph structure of TRRUST genes vs null (STRING co-expression as baseline)
**Rationale**: Entirely new direction. Build kNN graph (k=10) of TRRUST-relevant genes at each layer. Test if TRRUST edges are enriched in kNN graph edges vs random (network assortativity / edge overlap). Compare to STRING co-expression baseline.
**Family**: graph_topology
**Cost**: low-medium

## Methodological note for all future geometric tests
Always run permutation test (gene position shuffle) alongside parametric test. Use null-corrected effect size: (real − null_mean) / null_std as a Z-score.

## Retired directions (do not revisit without material change)
- FFL geometric ordering (all subspaces): perm-p > 0.29
- Signed displacement projection (AUROC): ceiling ~0.56
- Topology stability (bootstrap persistence): iter_0056/0057 negative
- bHLH/HIF1A margin generalization: closed
