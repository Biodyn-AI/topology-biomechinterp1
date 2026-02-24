# Executor Next Steps — iter_0022

## Top Priority (iter_0023)

### 1. Cell-type marker separation — power extension (HIGH VALUE)
- Current: 3 cell types, 30 within pairs, 61 cross pairs
- Next: Add macrophage (ALOX5+PTPRC), myeloid (CD14, LYZ), NK (NCAM1+PRF1) markers if present
- Also test: within-cell-type gene cluster UMAP visualization to document geometry
- Check against Geneformer word embeddings (same cell-type clustering?)
- Consider: PBMC or Tabula Sapiens data for broader cell-type coverage

### 2. TRRUST-exclusive — directionality test (PENDING)
- Split 141 TRRUST-exclusive pairs by direction (activation vs repression)
- Does repression vs activation have different proximity? (iter_0019 found no attention difference, but embedding may differ)
- Also test: TF hub gene effect (do TFs with >5 targets cluster together?)

### 3. Persistence entropy — proper computation (BLOCKED on ripser data)
- Need actual H1 lifetime arrays per gene pair (not stored from iter_0020)
- Rerun ripser at each layer and store full dgms[1] lifetime array
- Compute true Shannon entropy: -sum(p_i * log p_i) where p_i = lifetime_i / total

## Secondary (iter_0024+)

### 4. Cross-model alignment — cell-type geometry replication
- Do Geneformer word embeddings show the same cell-type cluster structure?
- CKA / Procrustes alignment between scGPT and Geneformer T-cell marker cluster
- If cell-type geometry is shared across models → architecture-agnostic biological encoding

### 5. STRING continuous — fine-grained analysis
- Current: Spearman rho=−0.093 (uniform across layers, weak but significant)
- Next: Is there a nonlinear threshold? Compute mean distance by score decile (10 points vs 5)
- Compare high-confidence (>0.7) vs medium (0.4-0.6) vs background

### 6. Intrinsic dimensionality by cell type
- For T cell marker genes: does their local neighborhood have lower ID than background?
- Two-NN estimator or PCA explained variance for within-type vs cross-type subspaces

## Retired Directions
- NULL sensitivity / permutation tests: covered, no further value
- Quintile-binned STRING gradient: replaced by continuous AUROC/Spearman
- Standalone SVD direction counting: covered in co-polarity

## Paper Updates Needed
- Add cell-type marker separation as a new section ("Geometric Encoding of Cell-Type Identity")
- Update TRRUST section to note independence from STRING overlap (iter_0022 H01)
- Add per-layer co-polarity bootstrap table (all 12 layers, CI>1)
