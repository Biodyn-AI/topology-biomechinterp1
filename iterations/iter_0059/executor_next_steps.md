# Executor Next Steps — iter_0059

## Completed This Iteration
- H01: Cross-seed TF boundary anchor stability confirmed (CV < 0.15 for top TFs). Out-degree hypothesis falsified (r=-0.012).
- H02: Pairwise TF→target 6D distance AUROC = 0.565 at best — real but below 0.62 threshold.
- H03: **Critical**: eff_rank vs AUROC correlation is a layer-depth confound (partial_r=-0.045 after controlling layer). The raw r=0.855 is spurious.

## Retire Directions
- **eff_rank as AUROC predictor**: fully confounded by layer. Retired.
- **Out-degree as boundary anchor explanation**: r=-0.012. Retire.
- **Pairwise 6D dist AUROC > 0.62**: signal too weak (0.565). Retire unless method changes substantially.

## Priority Hypotheses for iter_0060

### High Priority (novel + high upside)
1. **Within-layer seed variance vs eff_rank** (new_method, intrinsic_dimensionality):
   - Partial r= -0.045 showed eff_rank vs AUROC is confounded *across* layers.
   - But *within* each layer, do seed-to-seed AUROC differences (3 seeds) correlate with seed-to-seed eff_rank differences?
   - N=12 layer×3-seed = 36 obs; within-layer residuals free from layer confound.

2. **GO/functional annotation overlay on stable boundary anchors** (new_family, module_structure→biological_anchor):
   - STAT4, BACH2, ZEB1, RUNX1 are consistently high-margin across seeds.
   - Test: do these TFs share GO terms (immune response, cell migration) with their top targets?
   - Use goatools or manual lookup vs gene sets.

3. **SV2-4 vs SV5-7 unified AUROC comparison** (refinement, module_structure):
   - Prior H02 used SV5-7. Test if SV2-4 or joint SV2-7 gives better pairwise AUROC.
   - H02 result (AUROC=0.565 in SV5-7) is the baseline.

4. **TwoNN intrinsic dimension per layer** (new_method, intrinsic_dimensionality):
   - Replace effective rank with TwoNN estimator (model-free, local).
   - Test whether TwoNN ID differs between TF and target subsets at each layer.

### Medium Priority
5. **Community structure vs GO modules**: build kNN graph on SV5-7 embedding, detect communities (Leiden/Louvain), overlap with GO biological process.
6. **Cross-model validation (Geneformer)**: if Geneformer embeddings available, compute same SV5-7 AUROC for comparison.

## Critical Conceptual Update
The iter_0059 H03 finding changes the narrative: the layer-depth AUROC trend is NOT explained by dimensionality reduction (effective rank). The question is what actually drives the geometric separation of TF vs target genes across layers. Candidates:
- Layer-specific attention patterns concentrating on regulatory neighbors
- Progressive nonlinearity increasing class-conditional covariance
- Embedding normalization artifacts by layer

Next iteration should test at least one mechanistic explanation for the layer-depth AUROC trend.
