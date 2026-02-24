# Executor Next Steps — iter_0015

## High-priority follow-ups (iter_0016)

### Priority 1: SV4 GO biology characterization
- Analogous to iter_0014 H01 (SV3 GO enrichment)
- For SV4: at each of 12 layers, find top-K=52 and bottom-K=52 genes, run Fisher exact vs GO BP+CC
- SV4 is prominent in early layers (L0-L6); expect a distinct biological theme
- Zero new infrastructure; purely extends proven pipeline
- Expected: positive (SV4 follows same pattern as SV3)

### Priority 2: STRING confidence gradient replication on SV3
- H02 showed Spearman rho=1.000 for SV2. Replicate for SV3 (which also has strong PPI signal)
- Test: does the gradient hold equally, more strongly, or differently for SV3?
- Comparison of SV2 vs SV3 gradient slopes = quantitative geometry comparison

### Priority 3: GO-PPI driver with extreme quintiles
- H03 used median split (GO Jaccard threshold=0.042, quite low mean=0.063)
- Repeat with top 20% GO Jaccard vs bottom 20% STRING score for "pure GO" group
- Expected: cleaner separation; test whether GO-only signal remains positive

### Priority 4: Cross-model alignment (SV2 of scGPT vs Geneformer)
- **New family** (cross_model_alignment) — not yet explored
- If Geneformer embeddings available: compute SV2 of Geneformer on same 209 genes
- Test: does STRING co-pole enrichment hold in Geneformer SV2? Does Procrustes align the two SV2 axes?
- High-risk/high-reward; requires checking data availability first

## Retired directions (do not revisit without rescue rationale)
- SV1 PPI co-pole (consistently z≈0; no evidence)
- Repression anti-pole (negative, iter_0013 H02)
- GO BP enrichment in SV2 poles (negative, iter_0011 H03)
- TRRUST bootstrap CI on act-rep differential (underpowered, iter_0012 H01)
- Layer-depth Spearman monotonicity (neutral, iter_0014 H02)

## State of confirmed claims
1. **SV2 PPI co-pole**: STRING score≥0.4, 12/12 layers sig, mean z~3-10; hub-confound cleared (iter_0014 H03)
2. **SV3 PPI co-pole**: 12/12 layers sig, mean z~7 (iter_0013 H01)
3. **SV4 PPI co-pole**: 7/12 layers sig, mean z=2.21 (iter_0015 H01) — NEW
4. **SV2 confidence gradient**: rho=1.000, Q1 z=1.0 → Q5 z=5.0 (iter_0015 H02) — NEW
5. **PPI-only geometry**: z=2.66, 12/12 sig, rules out GO co-annotation as exclusive driver (iter_0015 H03) — NEW
6. **SV3 GO biology**: kinase/signaling (early) → MHC-II/immune (late), 12/12 sig (iter_0014 H01)
7. **TRRUST activation co-pole**: 12/12 sig (iter_0011 H01)
