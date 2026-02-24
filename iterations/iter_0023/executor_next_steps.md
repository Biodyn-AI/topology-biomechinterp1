# Executor Next Steps — iter_0023

## Top Priority (iter_0024)

### 1. TRRUST activation asymmetry — consolidation and extension (HIGH VALUE)
- Current: Activation AUROC=0.640 (12/12 sig), Repression AUROC=0.459 (0/12 sig)
- Next: Test with Dorothea (dorothea_human.tsv) — larger TF-target dataset with confidence levels
- Check: is the asymmetry consistent with the Dorothea A-E confidence tiers?
- Also: TF hub centrality — do TFs with many activation targets occupy central embedding positions?
- Geneformer replication: does word embedding space show same activation/repression asymmetry?

### 2. GO Biological Process proximity — deeper analysis (PROMISING)
- Current: rho=-0.077, AUROC=0.557 (3rd biological anchor confirmed)
- Next: Test GO BP with pathway-level grouping — are specific biological processes (immune signaling, ECM, metabolism) geometrically clustered?
- Compare GO BP vs GO MF vs GO CC — which ontology predicts embedding proximity best?
- Effect at each layer: rho deepens from L0=-0.063 to L7=-0.090 then returns (interesting shape — peak at L7)

### 3. Cross-model validation of cell-type geometry (HIGH VALUE)
- Do Geneformer word embeddings cluster same cell-type markers?
- Use iter_0019 Geneformer word embeddings (already loaded in that iteration)
- CKA alignment between scGPT cell-type marker clusters and Geneformer
- If yes → cell-type geometry is architecture-agnostic

## Secondary (iter_0025+)

### 4. Cell-type pair geometry — why T_cell vs B_cell only partially separated?
- T_cell vs B_cell AUROC=0.678 (within-immune cluster tight)
- T_cell vs fibroblast AUROC=0.195 (fibroblasts far from T cells, but AUROC<0.5 means CROSS is CLOSER than random!)
- Explanation: fibroblast marker genes (DCN, COL1A1, VIM) may be highly specific and far from all immune genes
- Test: compute mean distance from fibroblast markers to all named genes vs T-cell markers

### 5. Intrinsic dimensionality by cell-type
- Do T-cell marker gene neighborhoods occupy lower-dimensional local subspaces?
- Two-NN estimator or PCA for within-type (N=21 T-cell pairs) vs cross-type subspace

## Retired Directions
- Proxy persistence entropy (entropy of distance histograms) — invalid proxy, confirmed negative
- Quintile-binned STRING: superseded by continuous Spearman
- Standalone SVD direction counting: covered in co-polarity

## Paper Updates Needed
- Add TRRUST directional asymmetry as new section ("Regulatory Direction Asymmetry")
- Update GO section with H02 results (3rd biological anchor)
- Update cell-type section with contamination control confirmation
- Update Abstract/Conclusions to include all 4 confirmed anchors: STRING, TRRUST-excl, GO-BP, cell-type markers
