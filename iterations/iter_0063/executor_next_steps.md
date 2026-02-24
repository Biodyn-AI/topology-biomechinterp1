# Next Steps — iter_0063

## High Priority

### 1. H-G follow-up: Signed regulation with permutation control
- Repression AUROC (0.620) > Activation AUROC (0.599) in SV5-7 is an interesting finding.
- Need: proper permutation null (shuffle regulation labels, not gene positions) to confirm this difference is not due to sample size (n=141 vs 270).
- Also: cross-seed replication of the activation/repression AUROC gap.

### 2. H-A refinement: TF-class hypothesis across layers
- TF-TF AUROC in SV2-4 (0.539) vs TF-target AUROC (0.485) at L0 supports the class-separation hypothesis.
- Question: does TF clustering in SV2-4 also follow the layer-depth decay observed for TF-target cosine AUROC?
- Test: Run the TF-TF AUROC in SV2-4 across all 12 layers to see if the class-separation signal is consistent.

### 3. H-B redo with higher N permutations (N=500)
- With N=100, minimum detectable p=0.01; Bonferroni threshold is 0.005 for 10 windows.
- Redo with N=500 (minimum detectable p=0.002) to give Bonferroni-corrected test real power.
- Expected runtime: ~5 min.

## Medium Priority

### 4. Split-regime validation of H-G finding
- Repression vs activation AUROC difference — does it replicate in cycle1 (different gene set)?
- Quick check: run same comparison on cycle1_main with trrust_human.

### 5. SV2-4 repulsion: non-TF target clustering
- If SV2-4 separates TF class from non-TF, test: AUROC for random TF-nonTF pairs (should also be < 0.50).
- Compare with TF-target edge AUROC in SV2-4.

## Retired / Deprioritized
- Cross-model alignment (scGPT vs Geneformer): 2 consecutive negatives.
- topology_stability: negative.
- Procrustes transfer (C1): design flaw, retired.
