# iter_0041 Next Steps

## Confirmed findings to promote to paper

1. **GC attractor uniqueness (cross-lineage)**: B-cell/GC-TF attractor onset at L3 is unique vs T-cell/myeloid TFs which are pre-wired from L0. This elevates the finding from B-cell-specific to a general principle about *delayed vs. immediate* lineage encoding in scGPT.

2. **PAX5 L0 pre-wiring**: PAX5 neighbors are B-cell surface receptor genes (FCRL1, CD22, VPREB3) from L0 — biologically coherent. PAX5 is pre-encoded as B-cell master regulator.

3. **BCL6 metabolic isolation**: BCL6 lives in a metabolic stress neighborhood (NAMPT, GLUL, PFKFB3) not B-cell/GC cluster throughout all layers. Remove BCL6 from all GC attractor claims.

4. **B-cell manifold dimensionality reduction**: TwoNN ID decreases 8.16→5.09 across layers (rho=-0.951, p<0.0001). Unique to B-cell, not seen in T-cell or myeloid. First geometric evidence of manifold compression coinciding with attractor formation.

## High-priority next experiments

### H-A: B-cell ID reduction onset alignment
- Does the TwoNN ID decrease show an inflection near L3 (attractor onset)?
- Compute dID/dL per layer, find elbow or change point.
- Low cost, directly links ID reduction to attractor finding.

### H-B: T-cell vocabulary expansion
- Cycle1 vocab has different gene coverage. Check if FOXP3/GATA3/TBX21 are in cycle1.
- If yes, run T-cell attractor test in cycle1 to get fuller TF coverage.
- Important: currently only RUNX3 of 5 T-cell TFs available in cycle4.

### H-C: PAX5 L0 pre-wiring null test
- Shuffle gene-to-embedding assignment 200 times.
- Recompute PAX5 neighbor B-cell content at L0.
- Tests whether B-cell pre-wiring is statistically significant vs random.

### H-D: BCL6 GO annotation of metabolic neighborhood
- Map BCL6's top-20 neighbors to GO terms.
- Do they cluster in "cellular stress response", "metabolic process" GO terms?
- Provides biological interpretation of BCL6's non-GC placement.

### H-E: Multi-lineage ID comparison (bootstrap)
- Bootstrap resample each lineage marker set 100x.
- Compute uncertainty bounds on TwoNN ID per layer.
- Confirm B-cell ID decrease is robust vs. within-set resampling.

## Retirement candidates
- Cross-model alignment with Geneformer input token embeddings: retired (H02 iter_0040 negative). Need layer-wise Geneformer embeddings.
- GC-plasma subspace angle test: retired (H03 iter_0040 negative, near-orthogonal throughout).

## Paper structure update needed
- Add "Multi-lineage Comparison" section with H01 results.
- Add "Geometric Compression" section with H03 TwoNN ID results.
- Update B-cell attractor claim to distinguish BATF/BACH2/PAX5 (true attractor) from BCL6 (metabolic context).
