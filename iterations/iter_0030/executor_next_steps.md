# Next Steps — iter_0030

## Critical pivot: Remove OOV genes from all future analyses

The 14 "diverging" genes (FOS, JUNB, TNF, PTGS2, HLA-A, HLA-DPB1, KLF6, LDHA, LGALS1, NCAM1, NCOA3, NR4A3, PAX5, TBXAS1) are OOV tokens in scGPT with zero L0 embeddings. All prior bifurcation analyses are OOV artifacts. Future experiments must use the **195 in-vocabulary named genes only**.

---

## Priority hypotheses for iter_0031

### H01 (HIGH): kNN graph topology on 195 in-vocab genes
- **Family:** graph_topology
- **Rationale:** Re-run the iter_0028 H03 / iter_0029 H01 analyses (spectral gap, connected components, Fiedler vector) on the 195 in-vocab genes. Does meaningful topology survive without OOV contamination?
- **Method:** k=10 kNN graph on 195 in-vocab genes per layer. Spectral gap, n_components, Fiedler partition. Compare to random Gaussian null.
- **Novelty:** New method (corrected subset)

### H02 (HIGH): STRING distance-proximity correlation on 195 in-vocab genes
- **Family:** manifold_distance
- **Rationale:** The iter_0022 STRING-distance Spearman result (+0.3–0.4 AUROC) was computed on all 209 genes. Retest on 195 in-vocab genes to check if it holds without OOV contamination.
- **Method:** Pairwise Spearman(STRING_score, L2_distance) at each of 12 layers for all STRING pairs among 195 in-vocab genes.
- **Novelty:** Refinement (corrected subset)

### H03 (MEDIUM): Intrinsic dimensionality of 195 in-vocab genes by layer
- **Family:** intrinsic_dimensionality
- **Rationale:** The iter_0026 H03 PCA/PRF analysis was also contaminated by OOV genes. Retest on the corrected 195-gene set.
- **Method:** Participation ratio (sum(sv)^2 / sum(sv^2)) on PCA of 195 in-vocab gene embeddings per layer. Also compute PCA variance explained (PC1/PC2).
- **Novelty:** Refinement (corrected subset)

---

## Retired directions

- "Diverging cluster" topology (bifurcation, trajectory clustering) — OOV artifact, retired
- Slope-vs-STRING-degree — no correlation found, retire
- PC1 cell-type axis — tested negative iter_0027
- TRRUST polarity — tested negative iter_0027
- Spectral gap k-robustness — was positive but contaminated by OOV; re-test on corrected set

---

## Long-term agenda

After correcting for OOV, the 195-gene in-vocab analyses should be the main focus:
1. Verify STRING-distance correlation on corrected set
2. Verify spectral structure on corrected set
3. Cross-model alignment (scGPT vs Geneformer) on in-vocab genes — untested
4. Persistent homology on in-vocab set — high novelty
