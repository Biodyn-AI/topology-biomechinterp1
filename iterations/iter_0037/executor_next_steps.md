# Executor Next Steps — iter_0037

## Key findings to build on

1. **B-cell z-score = 7.55** (strongest seen; 5-gene panel, L2) — the B-cell geometric signal is robust.
2. **Germinal-center TFs (BATF, SPIB, BACH2) cluster near B-cell centroid; plasma TFs (IRF4) do not** — sub-lineage distinction within B-cell program.
3. **T-cell permutation null confirms specificity** — T-cell emp_p=0.96 vs B-cell emp_p=0.005.
4. **Centroid distance matrix**: B-cell↔T-cell=4.82, B-cell↔DC=6.79, B-cell↔Plasma=6.18.

## Top priorities for iter_0038

### H-A (HIGH PRIORITY): Geneformer cross-model validation
- Load Geneformer embeddings (if available) and test B-cell precision@10 at equivalent layers.
- If B-cell clustering survives in Geneformer, the finding generalizes from scGPT-specific to universal across transformer gene models.
- This is the most important remaining gap for the paper.
- Method: same precision@k=10 protocol; compare z-scores across models.

### H-B (MEDIUM PRIORITY): Germinal-center vs plasma geometric split
- The iter_0037 H02 result shows BATF/SPIB/BACH2 cluster near B-cell centroid while IRF4/IRF8 are far.
- Next: compute a GC-TF subcluster centroid vs plasma-TF subcluster centroid distance.
- Test whether this separation is layer-dependent (grows/shrinks from L0→L11).
- Biological anchor: this would map the B-cell differentiation axis geometrically.

### H-C (MEDIUM PRIORITY): STRING TF connectivity with actual STRING data
- Download STRING human PPI file (~4M edges) and compute full scoring for B-cell neighborhood.
- Compute odds ratio: STRING-connected neighbors vs random using Fisher's exact test.
- This provides a quantitative, reproducible biological coherence metric.

### H-D (LOW PRIORITY): PAX5/EBF1 expanded search
- Search full scGPT vocabulary (not just named genes) for PAX5, EBF1, BCL6, IKZF1, IKZF3.
- If in full vocabulary (4803 genes), test their distance to B-cell centroid.
- PAX5 and EBF1 are master B-cell lineage TFs — strong prior expectation they cluster near B-cell markers.

## Retired directions
- intrinsic_dimensionality SVD AUROC (confirmed negative in multiple iterations)
- TRRUST/GO proximity (confirmed negative)
- STRING PPI Euclidean proximity (confirmed negative in iter_0031 H02)
- graph_topology Laplacian spectral gap (tested, not biologically anchored)
- module_structure GO BP (negative)

## Paper status
- The B-cell geometric story has 15+ supporting findings across 7 iterations.
- Main gap: cross-model validation (Geneformer) and a broader cell-type panel.
- Current paper claim should be: "scGPT embeds B-cell identity markers in a geometrically coherent subspace distinguishable from T-cell, myeloid, and DC programs."
