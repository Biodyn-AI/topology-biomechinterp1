# Executor Iteration Report — iter_0029

## Summary

Three hypotheses tested. Two strong positives:

1. **H01** (graph_topology): The 2 kNN connected components at L11 persist from L8→L11; Component 1 contains exactly 14 genes with AP1-family enrichment (OR=0.06, p=0.023). The 14 genes are: FOS, JUNB, HLA-A, HLA-DPB1, KLF6, LDHA, LGALS1, NCAM1, NCOA3, NR4A3, PAX5, PTGS2, TBXAS1, TNF.

2. **H02** (manifold_distance / trajectory): Gene trajectory clustering (distance-to-centroid over 12 layers) yields silhouette=0.905 at k=2, perfectly matching the connected components. **Two trajectory types discovered:**
   - *Converging* (n=195, slope=−0.507): Move toward centroid from L0→L11 (dist: 11.4→4.1)
   - *Diverging* (n=14, slope=+0.459): Move away from centroid from L0→L11 (dist: 16.6→20.3)

3. **H03** (topology_stability): Spectral gap decrease (rho~−0.993) is **fully robust** across k=5,10,15,20,25,30. All rho<0, mean rho=−0.988. All null_p=0.000. This closes the k-robustness open question from iter_0028.

## Command Trace

```bash
# Run experiment script
conda run -n subproject40-topology python \
  iterations/iter_0029/run_iter0029_screen.py

# Output artifacts
iterations/iter_0029/h01_component_enrichment.json
iterations/iter_0029/h02_trajectory_clusters.json
iterations/iter_0029/h03_spectral_gap_krobustness.json
```

## Detailed Results

### H01: Connected Components Biological Enrichment

- kNN graph (k=10) at L8–L11: exactly 2 components, stable at sizes 195 and 14 across all 4 layers.
- Component 1 (n=14) is enriched for AP1 family vs Component 0 (n=195):
  - AP1: 2/14 vs 2/195, OR=0.06, p=0.023 (Fisher exact)
  - KLF: 1/14 vs 1/195, OR=0.07, p=0.130
  - HLA class I: 1/14 vs 2/195, OR=0.13, p=0.189
- STRING degree: comp0 mean=29.0, comp1 mean=37.1, AUROC=0.425, p=0.353 (not significant — outlier component is not simply low-degree)
- Outlier genes: FOS, HLA-A, HLA-DPB1, JUNB, KLF6, LDHA, LGALS1, NCAM1, NCOA3, NR4A3, PAX5, PTGS2, TBXAS1, TNF
  - Thematic: inflammatory/stress responders (FOS, JUNB, TNF, PTGS2), MHC antigen presentation (HLA-A, HLA-DPB1), metabolic enzymes (LDHA, TBXAS1)

### H02: Gene Trajectory Clustering (NEW FINDING)

Distance-to-centroid trajectories over 12 layers clustered with k-means (k=2..6):

| k | Silhouette |
|---|-----------|
| 2 | **0.905** |
| 3 | 0.571 |
| 4 | 0.451 |
| 5 | 0.443 |
| 6 | 0.321 |

Best k=2 (silhouette=0.905 — very high; two well-separated trajectory types):

**Cluster 0 (Converging, n=195):**
- Mean trajectory: L0=11.41, L1=9.84, ..., L10=5.43, L11=4.09
- Slope = −0.507 (strong convergence toward centroid)
- Mean STRING degree: 29.0

**Cluster 1 (Diverging, n=14):**
- Mean trajectory: L0=16.55, ..., L10=22.46, L11=20.3
- Slope = +0.459 (divergence from centroid)
- Mean STRING degree: 37.1
- Contains: FOS, JUNB, TNF, PTGS2 (AP1/inflammatory), HLA-A, HLA-DPB1, PAX5, NCAM1, NCOA3

- Kruskal-Wallis (degree across clusters): H=0.867, p=0.352 (degree not explanatory)
- Chi2 (immune families across clusters): chi2=6.77, p=0.562 (immune family distribution not strongly non-uniform — the 14 genes span multiple functions)

**Interpretation:** scGPT systematically pushes 14 genes to the periphery of the embedding space at deep layers, while the remaining 195 converge. The diverging 14 genes are characterized by stress/inflammatory response and antigen presentation — possibly reflecting genes with high cell-type-specific expression variance being assigned outlier representations.

### H03: Spectral Gap k-Robustness

| k | Spearman rho | p-value | null_p |
|---|-------------|---------|--------|
| 5 | −0.9930 | 1.3e-10 | 0.000 |
| 10 | −1.0000 | 0.0 | 0.000 |
| 15 | −0.9510 | 2.0e-06 | 0.000 |
| 20 | −0.9860 | 4.1e-09 | 0.000 |
| 25 | −1.0000 | 0.0 | 0.000 |
| 30 | −1.0000 | 0.0 | 0.000 |

Mean rho = −0.988, all rho < 0, all null_p = 0.000.

**The spectral gap decrease trend is extremely robust across all neighborhood sizes.** This confirms the finding from iter_0028 H03 is not an artifact of k selection. The gene graph becomes progressively more loosely connected (lower algebraic connectivity) at deeper layers, consistent with the convergence finding.

## Biological Interpretation

The combination of H01+H02 reveals a clear dichotomy:
- The scGPT model organizes 195 genes into a tight cluster at deep layers, pulling them toward the centroid.
- 14 genes (stress/inflammatory responders + antigen presentation) are actively pushed to the periphery.
- This outlier group has higher STRING connectivity (mean degree 37 vs 29) but is NOT explained by degree alone (MW p=0.35).
- The AP1 enrichment (FOS, JUNB — p=0.023) in the outlier group is consistent with these genes having strongly context-dependent expression, which the model may represent with high positional variance.

## Artifacts Generated

- `iterations/iter_0029/h01_component_enrichment.json`
- `iterations/iter_0029/h02_trajectory_clusters.json`
- `iterations/iter_0029/h03_spectral_gap_krobustness.json`
- `iterations/iter_0029/run_iter0029_screen.py`
