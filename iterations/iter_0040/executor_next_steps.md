# iter_0040 Next Steps

## High-priority follow-ups

### 1. BCL6 correction + revised GC attractor membership (IMMEDIATE)
- BCL6 ranks 750-1500 throughout all 12 layers — NOT part of GC attractor
- Revised GC attractor = {BATF, BACH2, PAX5} (3 genes confirmed)
- Update all paper claims that list BCL6 in the GC cluster
- Verify in cycle1 embeddings: does BCL6 also fail there?

### 2. PAX5 early proximity characterization (HIGH VALUE)
- PAX5 rank 36 at L0 already — anomalously close to B-cell centroid from the start
- Is this pre-existing or does it strengthen through layers?
- Test: what are the top-20 nearest genes to B-cell centroid at L0 that include PAX5? Any other TFs?
- Biological anchor: PAX5 is master B-cell TF (EBF1 target), explains L0 proximity

### 3. Geneformer layer-wise gene embeddings (CROSS-MODEL)
- Input token embeddings gave negative result (H02 this iteration)
- For a valid cross-model test: extract Geneformer layer-wise representations by running forward pass on immune scRNA-seq data with per-cell per-gene token tracking
- This requires engineering: run Geneformer on tabula_sapiens_immune, aggregate gene embeddings per layer across cells
- Mark as future engineering task; include in paper as limitation

### 4. IRF4 plasma-bridge test (REQUIRES CYCLE1)
- IRF4 is OOV in cycle4 but may be in cycle1 (4803 genes)
- If in cycle1: check IRF4 rank trajectory across layers — does it diverge from GC-TFs?
- Expected: if IRF4 bridges GC→plasma program, it should not converge as strongly as PAX5/BATF/BACH2

### 5. Onset layer L3 mechanistic probe
- Both H01 and H03 independently confirm onset at L3
- What changes at L3? Check attention head patterns, embedding norm, or kNN graph topology at L3 vs L2
- Is L3 the first attention layer that processes regulatory context beyond surface co-expression?

## Retired directions
- GC-plasma subspace convergence: null (rho=-0.098) — retire for at least 2 iterations
- BCL6 as GC attractor gene: corrected, remove from claims

## Open questions
- Does the GC attractor onset at L3 replicate in cycle1 embeddings?
- Is PAX5's L0 proximity explained by training data co-expression with B-cell markers, or by regulatory network position?
- What is the within-cell-type composition of the immune dataset that drives these patterns?
