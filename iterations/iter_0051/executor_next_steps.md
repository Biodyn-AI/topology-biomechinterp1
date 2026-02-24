# Next Steps: iter_0052

## Immediate follow-ups

### Priority 1: SV1 housekeeping gene hypothesis (from H02)
- The strong depletion of TFs from high-SV1 genes (OR=0.108 at L0) suggests SV1 captures housekeeping vs regulatory dichotomy.
- Download or hardcode a housekeeping gene set (e.g., Eisenberg & Levanon 2013, ~3804 human housekeeping genes).
- Check enrichment of housekeeping genes in the top-20% SV1-high group vs bottom-20%.
- **Expected outcome**: high-SV1 genes enriched for housekeeping; if confirmed, this mechanistically grounds the SV1 axis.

### Priority 2: Spectral direction search for regulatory-specific content
- H01 showed SV2-4 is co-expression confounded; only 1/12 layers (L8) retains residual signal.
- Try: compute first 15 SVs and scan which SV group (5-8, 9-12) shows TRRUST proximity AFTER co-expression regression.
- Test whether any orthogonal spectral axis carries regulatory signal independent of co-expression.

### Priority 3: H1 Betti curves (loops) on circuit genes
- H03 found circuit genes are spatially compact in SV2-4 but don't form topologically distinct 0-dim islands.
- Try H1 persistent homology (loops) using ripser/gudhi — circuits may form ring-like topological structures.
- Focus on L8 (best rbc layer) and compare circuit vs non-circuit Betti-1 lifetimes.

## Paper update requirements
- Update SV2-4 regulatory claim: was "encodes regulatory proximity," now "encodes co-expression structure with TF enrichment; regulatory-specific content is largely a co-expression confound."
- Add SV1-TF depletion finding: mechanistic identity of the SV1 axis.
- Provide precise H02 Fisher statistics in paper.

## Retirement candidates
- `module_structure` (SV2-4 regulatory proximity): retire. Replaced by co-expression framing.
- `manifold_distance` family: >=2 negative/inconclusive without new method — retire unless new anchor.

## New directions to explore
1. Cross-model alignment (scGPT vs Geneformer spectral axes) — not yet tested.
2. GO enrichment of SV1-high vs SV1-low gene sets — biological grounding.
3. Layer-specific kNN graph rewiring null to test graph topology changes.
