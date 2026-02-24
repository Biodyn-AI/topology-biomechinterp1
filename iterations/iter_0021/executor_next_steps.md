# Executor Next Steps — iter_0021

## High Priority (iter_0022)

### 1. TRRUST vs STRING Overlap-Corrected Test (new_method)
- H02 found TRRUST pairs are also geometrically closer. Key question: Is this independent of STRING overlap?
- Many TRRUST TF-target pairs may also be STRING edges.
- Compute TRRUST-only pairs (remove any pair also in STRING). Test: is proximity still significant?
- Expected: If still significant → regulation adds independent signal. If not → TRRUST proximity is mediated by shared STRING membership.
- Family: manifold_distance, novelty: new_method

### 2. STRING Score as Continuous AUROC Predictor (refinement of H01)
- H01 Spearman with 5 quintile means was underpowered.
- Use raw STRING score as continuous predictor. Compute Spearman(score, distance) for all 3092 pairs at each layer. Also test AUROC for thresholds 0.5, 0.6, 0.7, 0.8 vs 0.4 baseline.
- Family: manifold_distance, novelty: refinement

### 3. Layer-Stratified Bootstrap CIs for Co-Polarity Enrichment (validation)
- H03 computed bootstrap CI at layer 8 only.
- Run bootstrap CIs at all 12 layers to confirm early-layer enrichment (1.55x) is significantly higher than late (1.33x).
- Report CI for each layer. Test: CI at layers 1-4 vs CI at layers 9-11 — do they overlap?
- Family: persistent_homology, novelty: refinement

## Medium Priority

### 4. Cross-model alignment: Geneformer vs scGPT proximity (cross_model_alignment)
- iter_0019 H01 found Geneformer word embeddings have high Procrustes alignment with scGPT.
- Compute pairwise distances in Geneformer space for STRING pairs. Do STRING pairs show the same proximity in Geneformer as in scGPT?
- If YES: geometry is model-agnostic → strong evidence for biological signal.
- Family: cross_model_alignment, novelty: new_method

### 5. Intrinsic dimensionality by layer (intrinsic_dimensionality)
- H0 lifetime declines from 0.667 to 0.242 across layers (compaction signature).
- Compute TwoNN or PCA effective dimension for named gene embeddings at each layer.
- Test: does intrinsic dimension correlate with H0/H1 lifetime (confirming compaction interpretation)?
- Family: intrinsic_dimensionality, novelty: new_family

## Retirement Status

- `null_sensitivity` (permutation shuffles): Sufficient coverage. H01 iter_0020 established permutation null z=4.88. No need for more iterations here.
- `graph_topology` (SVD co-polarity direction): Well established across 6 iterations. Focus now on refinements and cross-model.

## Key Claims for Paper (ready to write)

1. STRING pairs are geometrically closer in scGPT unit-sphere space across all 12 layers (d=−0.237, bootstrap CI confirmed).
2. 3-axis co-polarity enrichment: 1.33x at layer 8, 1.55x at layers 1-4. Bootstrap 95% CI=[1.23, 1.41] excludes null.
3. H1 persistence lifetime declines monotonically with layer (rho=−0.916, p<0.0001): progressive topological compaction.
4. BOTH STRING (PPI) and TRRUST (TF regulation) pairs show significant proximity (effects −0.048 and −0.043 respectively): geometry encodes general biological interaction proximity.
