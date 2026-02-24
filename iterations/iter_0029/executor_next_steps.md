# Executor Next Steps — iter_0029

## Top Priority for iter_0030

### 1. Biological characterization of the 14 diverging genes (HIGH VALUE)
The 14 outlier genes (FOS, JUNB, TNF, PTGS2, HLA-A, HLA-DPB1, KLF6, LDHA, LGALS1, NCAM1, NCOA3, NR4A3, PAX5, TBXAS1, TNF) diverge from centroid at deep layers. Next:
- Fetch GO enrichment for these 14 genes (mygene API) and compare to the 195 converging genes
- Test if they are enriched for "immediate early response" / "stress response" GO terms
- Check if these 14 genes appear as cell-type marker genes more frequently in published datasets

### 2. Bootstrap confidence intervals for trajectory clusters (VALIDATION)
- Bootstrap the trajectory cluster assignments (500x, resample layers) to confirm the 2-cluster solution is stable
- Compute Jaccard similarity between bootstrap cluster assignments and original

### 3. Hub gene centrality at deep layers — cross with trajectory clusters
- From iter_0028 H01: STRING hub genes migrate centrally at deep layers (rho=+0.225 at L10)
- Test: Are the 14 diverging genes systematically *non-hubs* vs the 195 converging genes?
- MW test STRING degree: diverging vs converging (already done, p=0.35 NS — document as control)

### 4. Cross-model validation (if Geneformer embeddings available)
- Do Geneformer embeddings also produce a 2-trajectory split at deep layers?
- If yes, the 14-gene outlier group has cross-model validity

## Retirement Status Updates
- PC1 T-cell/APC axis: RETIRED (iter_0028)
- TRRUST polarity: RETIRED (iter_0027)
- Dorothea confidence-tier distance (raw): RETIRED

## Open Questions
- Are the 14 diverging genes the "most context-dependent" genes in the model?
- Does the spectral gap collapse at L11 reflect the convergence of 195 genes to a near-degenerate manifold?
- Can we characterize the "converging" manifold at L11 as approximately 1D (consistent with iter_0026 H03 PC1 dominance finding)?
