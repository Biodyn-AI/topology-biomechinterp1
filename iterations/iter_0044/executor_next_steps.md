# Next Steps: iter_0044 → iter_0045

## Retirements this iteration
- **H01 (cross-lineage compression)**: Retired. 5-point TwoNN is methodologically unsuitable for detecting layer-wise compression patterns. The iter_0043 finding was likely computed differently (full-manifold TwoNN). Do not revisit without a full-manifold TwoNN implementation.
- **H03 (BCL6 metabolic divergence)**: Retired. BCL6 co-moves with both tracks (rho=0.5), confirming no selective metabolic attractor. Second negative result combined with iter_0043 H01/H02.

## Open directions

### Priority 1: Fix TwoNN methodology (foundational)
The correct TwoNN for layer-wise intrinsic dimensionality should:
1. Compute TwoNN on the FULL gene manifold (~4941 genes) at each layer
2. Report the global ID estimate, OR local ID in the vicinity of target gene subsets
This would allow proper replication of iter_0042/0043 positive results and meaningful comparison.

### Priority 2: TCR circuit with clean gene separation (H02 follow-up)
- Remove CD247 from both centroid and circuit definitions
- Use clean T-cell centroid: CD3E, CD3D, TRAC, TRBC1 (CD247 excluded)
- Use pure signaling circuit: CD28, LCK (2 genes only, or wait for more)
- Or try a different signaling axis: CD2/CD58 interaction, or PTPRC/CD45

### Priority 3: New family — graph topology (community detection)
The module_structure family has been positive in iter_0036 (B-cell kNN communities). Extend to:
- Build kNN graph (k=10) on all ~4941 genes at L2 and L11
- Detect communities (Louvain or spectral clustering)
- Test if T-cell/myeloid/B-cell markers sort into distinct communities
- Compare modularity L2 vs L11

### Priority 4: Persistent homology on gene embedding subsets
Family not tested in recent iterations. Try Rips complex on B-cell marker subset (n=20 nearest neighbors at L11) and compute Betti numbers. Compare to random subset of same size.

## Recommended portfolio for iter_0045
1. **H1 (new family, graph_topology)**: kNN community structure at L2 vs L11; lineage sorting test
2. **H2 (new method, intrinsic_dimensionality)**: Full-manifold TwoNN on all 4941 genes at each layer; measure overall ID trajectory
3. **H3 (new method, manifold_distance)**: TCR circuit with CD247 removed from centroid; clean test of signaling convergence
