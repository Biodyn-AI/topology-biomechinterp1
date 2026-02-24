# Executor Next Steps — iter_0039

## Key findings to build on

1. **Drift endpoint resolved**: B-cell centroid drifts to GC-TF neighborhood (BATF rank-2, SPIB rank-5 at L11). This is the first direct evidence that scGPT layer processing moves B-cell identity representations toward the GC-TF cluster.

2. **Extended GC cluster validated**: BCL6 and PAX5 are in cycle4_immune embeddings and cluster with BATF/BACH2 (all showing strong negative z). The GC master regulator panel is now: BATF, BACH2, BCL6, PAX5 (4 genes).

3. **PRDM1 artifact identified**: The cycle1 result where plasma z went positive was an artifact of PRDM1 appearing in both the B-cell anchor AND the plasma panel. With clean anchor (no PRDM1), plasma divergence is RELATIVE (z increases monotonically: −1.35 → −0.81) but never positive. **Correct the iter_0038 claim in the paper.**

4. **CD19 + BLK as minimal anchor**: LOO shows CD19 and BLK are the critical anchor genes for GC-TF proximity signal. MS4A1, CD79A, PRDM1 are redundant for precision@10. Future experiments should use CD19+BLK as minimal 2-gene B-cell anchor.

## Priorities for iter_0040

### Experiment 1 (HIGH-VALUE): Cross-model validation in Geneformer token space
- Locate Geneformer token embeddings (check subproject_38 outputs; or use the `geneformer_embedding_geometry_audit.py` script)
- Test: do BATF, BACH2, BCL6, PAX5 cluster near B-cell identity tokens in Geneformer?
- If unavailable: use cycle4 seed replicate (cycle4_immune_seed43, cycle4_immune_seed44) to check seed robustness

### Experiment 2 (EXTENSION): Minimal 2-gene anchor (CD19+BLK) precision screen
- Verify CD19+BLK yields strong GC-TF precision even with just 2 anchor genes
- Test across all 12 layers (not just probed layers) to identify the transition point

### Experiment 3 (NEW FAMILY): Principal angle between GC-TF and plasma-TF subspaces
- At L11 (cycle4_immune): compute SVD of GC-TF embedding matrix [4, 512] and plasma-TF [2, 512]
- Compute principal angles between the subspaces
- Compare to 500 random same-size gene-pair subspaces
- If GC and plasma subtend near-orthogonal subspaces, this corroborates the 94° angle finding from iter_0038

### Retired directions (do not revisit)
- "Plasma z goes positive" claim (PRDM1 artifact; replace with relative divergence framing)
- NK/myeloid cell type screens (vocabulary limitations in cycle1)
- directional drift cosine alignment (cosines near 0, confirmed not directional)
- STRING PPI Euclidean proximity
- intrinsic_dimensionality SVD AUROC

## Paper updates needed for iter_0040
1. Correct the PRDM1 artifact: plasma z goes relatively less negative (not absolutely positive)
2. Add drift endpoint finding: B-cell centroid converges to GC-TF cluster
3. Add extended GC cluster result: BCL6 and PAX5 confirmed
4. Add LOO ablation result: CD19 and BLK are the critical anchor genes
