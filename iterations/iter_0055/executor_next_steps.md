# Next Steps — iter_0055

## Key findings to follow up

### 1. Complementary subspace encoding (HIGH PRIORITY)
H02 reveals a striking complementary pattern:
- SV5-7: high AUROC at L0-L3 and L9-L11, near-random at L6-L8
- SV2-4: high AUROC at L4-L8, drops at L9-L11

**Next**: Combine both subspaces (SV2-7 jointly) across all layers. Test if the combined 6D projection maintains >0.73 AUROC at all layers, implying full information preservation.

### 2. Layer-specific crossover biological interpretation (HIGH PRIORITY)
H03 shows SV5-7 directionality surpasses SV2-4 at L9. H02 shows SV5-7 AUROC recovers at L9.
These two findings together suggest L9 is a mechanistically distinct layer.

**Next**: Extract L9-specific gene displacement vectors (top movers in SV5-7 at L9 vs L0). Check if these genes are enriched for known late-stage regulatory targets or binding factors. Use GO enrichment.

### 3. Full SV1-7 joint classification (MEDIUM PRIORITY)
Current results use only 3-axis subspaces. SV1 has been separately shown to encode non-circuit vs circuit information.

**Next**: Joint LR classifier on SV1-7 (7D) across all 12 layers. Expected to show monotonically high AUROC if subspaces encode complementary but consistent signals.

### 4. Cross-seed directionality trajectory replication (MEDIUM PRIORITY)
H01 cross-seed AUROC replicates well. H03 directionality uses only main seed.

**Next**: Run H03 on seed43 and seed44. Verify that the SV2-4 decrease and SV5-7 increase trajectories are reproducible across seeds.

### 5. Null sensitivity for directionality (LOW PRIORITY — if time allows)
H03 shows significant directionality in BOTH subspaces. A label-shuffle null would confirm these are not artifacts of SVD orientation.

**Next**: Shuffle TRRUST source/target labels (keeping gene identities), recompute displacement magnitude. Expected ~0 if signal is real.

## Retired directions (no further investment)
- persistent_homology (>=2 negative)
- intrinsic_dimensionality basic TwoNN (negative)
- module_structure co-expression proxy (negative)
