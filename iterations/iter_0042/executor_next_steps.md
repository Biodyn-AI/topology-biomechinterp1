# Executor Next Steps — iter_0042

## Immediate priorities for iter_0043

### Priority 1 (H03 follow-up): Full GC-TF triangle visualization
BATF, BACH2, PAX5 all converge near each other by L11 (distances ~4–5 each).
- Compute full pairwise distance matrix among {PAX5, BATF, BACH2, BCL6} at each layer.
- Show: PAX5–BATF, PAX5–BACH2, BATF–BACH2 all decrease; BCL6–PAX5, BCL6–BATF, BCL6–BACH2 remain large or decrease more slowly.
- This would produce a "geometric triangle" evolving across layers — highly visual and publishable.
- Family: manifold_distance, novelty: new_method.

### Priority 2 (H02 follow-up): GO annotation of BCL6 metabolic cluster
Need to annotate what the BCL6 metabolic cluster means biologically:
- Use MyGene.info or manual annotation for: NAMPT, GLUL, PFKFB3, ACSL1, NIBAN1, FNDC3B, VMP1, STAT3, CEBPD, TRIB1.
- Common themes: glycolysis, fatty acid synthesis, metabolic stress response.
- Confirm BCL6's known role in Warburg effect / aerobic glycolysis regulation.
- Family: module_structure, novelty: new_method.

### Priority 3 (New family — null sensitivity): Bootstrap stability of BATF convergence
BATF has rho=-0.972 — exceptionally strong. Bootstrap test:
- Sample 80% of all genes (keeping BATF + B-cell markers), recompute BATF rank near centroid.
- 200 bootstrap replicates. Check stability of convergence trajectory.
- Family: null_sensitivity, novelty: new_method.

## Retired directions (do not re-test without major change)
- Cross-model alignment (scGPT vs Geneformer input embeddings) — blocked by vocab mismatch
- GC-plasma subspace angles — inconclusive in iter_0040
- Chromosomal proximity — negative
- PC2/PC3 axes — retired

## Paper update priorities
- Add H03 BATF/BACH2 convergence as new subsection: "Dynamic recruitment of GC transcription factors"
- Update H01 (TwoNN changepoint L3) in intrinsic dimensionality section with slope-change table
- Add H02 result (BCL6 specificity) to BCL6 subsection: "BCL6 metabolic isolation is gene-specific"
