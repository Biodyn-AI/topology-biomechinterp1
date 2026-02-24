# Executor Next Steps — iter_0005

## Immediate priorities (iter_0006)

### 1. Recover full gene vocabulary for GO clustering (highest priority)
The 4803-gene embedding includes the full scGPT gene vocabulary but we only have 209 named genes. To unlock GO clustering and drift enrichment at full power, need to recover gene names for all 4803 positions.

**Action:** Search subproject_38 code/data for a gene vocabulary file or the scGPT model's vocab. Expected path: `subproject_38/.../scgpt_vocab.json` or similar. If unavailable, use the scGPT model's built-in `vocab.json`.

### 2. MADA cross-check for intrinsic dimensionality
iter_0004 recommended MADA as a cross-check for TwoNN. Now H03 (effective rank) confirms the trend independently. Add MADA to complete the geometric dimensionality picture.

### 3. Persistent homology (Betti curves) on GO-term subgraphs
With the biological anchoring signal from H02 (immune/transcription genes show high drift), test whether the topological structure (persistent homology H0/H1) of these gene subsets differs from matched random sets.

**Family:** persistent_homology (not yet tested this project)

### 4. Validate H02 drift enrichment with larger gene pool
Once full gene vocabulary is recovered, rerun residual drift + GO enrichment with all 4803 genes (split top/bottom 500). The 7 nominally significant terms are biologically coherent; this needs FDR-level replication.

### 5. Cross-layer CKA derivative signal
H03 shows ER drops monotonically but TwoNN is more irregular. Run CKA similarity between consecutive layers (already have the matrix from iter_0004) to see if layer pairs with the steepest ER drop also show lowest consecutive CKA.

## Retired directions
- TRRUST co-target distance clustering (small gene pool, inconclusive x2)
- Cross-layer linear CKA replication (fully answered iter_0004)
- Cross-model feature-vector alignment (no Geneformer embeddings available)

## Hypothesis candidates for iter_0006

| ID | Name | Family | Novelty |
|----|------|--------|---------|
| H01 | Persistent homology (H0/H1) on high-drift gene subgraph | persistent_homology | new_family |
| H02 | Full-vocabulary GO clustering (if vocab recovered) | module_structure | new_method |
| H03 | MADA intrinsic dimensionality per layer | intrinsic_dimensionality | new_method |
