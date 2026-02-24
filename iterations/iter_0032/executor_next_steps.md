# Executor Next Steps: iter_0033

## Priority 1 (HIGH): Regulatory Geometry Decay — Activation vs Repression Split
**Rationale:** H02 showed high-conf regulatory AUROC decays from 0.57 at L0 to 0.50 at L11. Now split by regulation type (activation vs repression) to test if decay differs. Hypothesis: activation pairs might preserve geometry longer (positive feedback loops reinforce co-expression).
- Source: Dorothea with sign annotation (mor > 0 = activation, mor < 0 = repression)
- Test: AUROC per layer separately for Act vs Rep; compare rho values

## Priority 2 (HIGH): Community 0 vs Community 1 — Differential Analysis
**Rationale:** H01 found very strong 2-community structure at L11 (z=34) but biology not clear. Now characterize communities functionally:
- What is the mean L2 norm of each community (norm CV was rising per iter_0031)?
- Are there STRING degree differences between communities?
- Do community members differ in: (a) TF vs non-TF, (b) cell-surface vs intracellular, (c) mean pairwise distance within community
- Compute PC1 scores per community (community 1 = negative/B-cell pole?)

## Priority 3 (MEDIUM): Dorothea OOV-Corrected Full Cross-Layer Profile
**Rationale:** The decay from AUROC=0.57 to 0.50 was measured with 5000 random pairs. Verify:
- Check if decay is robust to different random seeds (3 seeds)
- Run on Geneformer embeddings (if available) for cross-model comparison
- Test whether the "transition layer" (where AUROC crosses 0.5) correlates with spectral gap changes

## Priority 4 (LOWER): PC1 Polarity with Expanded Gene Sets
**Rationale:** H03 underpowered (only 3 B-cell markers in 195-gene vocab). Expand:
- Take all B-lineage genes in full vocab (CD19, MS4A1, CD79A, CD79B, IGHM, IGHD, IGKC, EBF1, PAX5, MZB1, etc. — check which are in gene_list.txt)
- Compute PC1 projection for ALL ~4800 vocab genes at L11 and test larger biological groups

## Retirement Candidates
- intrinsic_dimensionality: PR collapse and PC1 structure established. No new hypothesis.
- null_sensitivity: chromosome co-location null was negative (iter_0025). Don't retry.
- STRING→L2 distance: definitively retired (iter_0031 H02, multiple prior nulls).

## Artifacts to Produce in iter_0033
- `h01_community_differential_analysis.json` — characterize 2 communities
- `h02_dorothea_activation_repression_decay.json` — activation vs repression decay
- `h03_pc1_expanded_vocab.json` — expanded B/T cell gene sets on PC1
