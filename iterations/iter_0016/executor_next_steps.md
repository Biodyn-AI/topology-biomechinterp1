# Executor Next Steps — iter_0016

## Top Priorities for iter_0017

### Priority 1: Cross-model alignment (Geneformer)
**H-G (Geneformer SV2 co-pole)**: Check if Geneformer residual embeddings are available.
If present, compute SV2 co-pole z-score for STRING pairs (K=52, N=3092, score≥0.4) and compare
to scGPT SV2 z=5+. Conservation of PPI geometry across two distinct biological LLMs would be the
highest-upside single experiment in the project.
- Action: `find /Volumes/Crucial\ X6/MacBook/biomechinterp/biodyn-work -name "*.npy" | grep -i geneformer`

### Priority 2: Multi-axis PPI gradient (SV4 quintile sweep)
SV4 shows mean_z=2.21 (iter_0015). Run the full quintile gradient test for SV4 to confirm
whether SV4 also encodes a monotonic confidence gradient (Q1→Q5). Expected: weaker but monotonic.
- Minimal cost: extend H01 script to SV4_IDX=3.

### Priority 3: Axis independence with biological labels
Use the axis-dominant pair labels from H02 to test whether the 3 axis clusters differ in GO
or TRRUST category composition. E.g., SV2-dominant pairs → more signal transduction; SV4-dominant
→ more apoptosis/endosome. This would make the orthogonality story mechanistically interpretable.
- Method: compute GO Jaccard within each axis-dominant subset vs background.

### Priority 4: SV3 GO biology profile
Analogous to H03 but for SV3. This would complete the biological characterization of all 4 axes
(SV2 already done in iter_0010/011, SV3 partial in iter_0013/014, SV4 done here).
- Zero new infrastructure: copy H03 script, swap SV4_IDX=2 → SV3_IDX=2.

## Evidence Consolidation
The project now has a coherent, multi-axis story ready for paper:
1. scGPT residual stream encodes STRING PPI proximity on at least 4 orthogonal axes (SV2–SV5)
2. SV2 and SV3 show monotonic confidence gradient (rho=1.00, 0.90 respectively)
3. SV2/SV3/SV4 are largely independent (max r=0.247 inter-axis copole correlation)
4. PPI geometry is driven by physical interaction independent of GO co-annotation (iter_0015 H03)
5. SV4 GO profile: endosome/apoptosis vs antiviral/T-cell biology — distinct from SV2/SV3

## Retirement Status
- module_structure single-gene permutation tests: retired (2+ negatives; focus on gene-pair level)
- TRRUST TF→target direction (SV2/SV3): retired (iter_0013 H02 negative)
