# Next Steps — iter_0021

## Priority 1 (Follow-up H03 — high upside)
**H1 persistence loop analysis**: iter_0020 collected H0+H1 diagrams from ripser. Next: extract H1 (loop) lifetimes per layer, compare H1 Betti number between STRING-dense subgraphs vs STRING-sparse subgraphs. Test whether GO-annotated gene clusters have distinct H1 topology. This is the core "topological structure" claim the project was designed for.

## Priority 2 (Consolidate validated findings)
**Layer-stratified bootstrap for multi-axis composite**: H01 validated with permutation z=4.88. Now run bootstrap CI across all 12 layers to get uncertainty bands on the enrichment-vs-count curve. This produces a publication-quality figure.

## Priority 3 (New family — intrinsic dimensionality)
**Local intrinsic dimensionality by biological annotation**: Use TwoNN or PCA-based local ID estimator on the 209-gene neighborhood at each layer. Test if genes in the same GO/STRING community have lower local ID (tighter cluster) than genes in different communities. Relates to manifold_distance family which has not been tested in this precise way.

## Retired directions (do not revisit without rescue rationale)
- TRRUST signed SVD (consistently negative, >=3 tries)
- 1D single-SV TRRUST co-polarity
- GO sweep for module structure (negative/inconclusive >=2)
- Attention for STRING (AUROC 0.48, below chance — definitively negative)

## Brainstormer brief
Key results to communicate:
1. Multi-axis composite z=4.88 (perm null) + magnitude quintile monotone (r=0.900) — validated claim
2. H03: STRING pairs geometrically closer on unit sphere, d=−0.237, all 12 layers, perm control clean
3. H02: dissociation of attention (TRRUST) vs SVD (STRING) confirmed by AUROC profile
4. Next: H1 topological loop analysis, local ID by community, layer bootstrap CIs
