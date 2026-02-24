# Next Steps — iter_0036

## Priority 1 (High): Cross-model Geneformer validation
The B-cell geometry finding is now complete within scGPT. The highest-value next step is cross-model validation:
- Obtain Geneformer residual embeddings for the same gene set
- Compute kNN precision@10 for B-cell markers in Geneformer representations
- Compare z-scores across layers: does B-cell geometry emerge at early layers in Geneformer too?

## Priority 2 (High): Expand marker set + confirm T-cell specificity
- Current B-cell set: 7 in-vocab. Explore adding more (EBF1, XBP1, BLNK were in list but not in vocab — check full vocab for additional B-cell genes).
- Confirm T-cell specificity: run permutation null for T-cell markers to show their low z=0.29 is not artifact-driven (complementary negative control).

## Priority 3 (Medium): SPIB/BACH2/BATF biological coherence
- Check whether SPIB and BACH2 appear as edges in the TRRUST/Dorothea dataset used in this project.
- Both are B-cell identity TFs regulated by the BCR signaling pathway. Confirm STRING connections.
- This would strengthen the H03 finding from "top neighbors happen to be B-cell TFs" to "top neighbors are known B-cell TF regulators."

## Priority 4 (Lower): Persistent homology on B-cell neighborhood
- Retired directions (intrinsic_dimensionality, module_structure with TRRUST/GO annotation) remain retired.
- Consider: Betti-0 curve of B-cell marker subgraph as filtration radius varies — does B-cell gene set form a persistent connected component before non-B-cell genes merge?
- This would be a new_method in persistent_homology family (not previously tested on cell-type marker subgraphs).

## Key claims now established (paper-ready)
1. B-cell geometric clustering in scGPT L2: precision@10 z=4.55, survives permutation null (p=0.005)
2. Cell-type specificity: z=4.55 (B-cell, n=7) vs z=0.29 (T-cell, n=12) — 15x difference
3. Biological neighborhood content: SPIB, BACH2, BATF in top neighbors (known B-cell TFs)
4. Layer trajectory: clustering peaks at early layers (L2) and attenuates toward L11

## Retired directions (do not revisit without rescue rationale)
- TRRUST/GO co-annotation → embedding proximity: 2 failures (iter_0033, iter_0035)
- intrinsic_dimensionality (SVD AUROC for TFs, etc.): repeated negatives
