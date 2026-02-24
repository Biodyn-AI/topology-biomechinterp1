# Next Steps — iter_0058

## Confirmed findings this iteration

1. **TF boundary anchors** (H01): BCL11A, NFKB1, FOXO3, ZEB1, STAT3 are top boundary TFs at peak layers. Forkhead/STAT families well-separated; bZIP/C2H2-ZF near target centroid. Biologically coherent for immune dataset.
2. **Laplacian alignment** (H02): Statistically significant but practically weak alignment (Z=7.81, Δ≈3° from random). TRRUST graph structure does not strongly imprint on SV subspace geometry.
3. **Effective rank → AUROC** (H03): rho=0.855 across 36 (seed×layer) observations. Higher full-spectrum dimensionality → better TF/target discrimination. Strong cross-seed reproducible result.

---

## Priority hypotheses for iter_0059

### Top priority (P1): Cross-seed boundary anchor consistency
- For each top-10 TF anchor (BCL11A, NFKB1, RB1, FOXO3, ZEB1), compute margin in seed43 and seed44 at L2/L3.
- Test: are the same TFs consistently top boundary anchors across seeds?
- Family: module_structure | novelty: refinement of H01

### Top priority (P2): Effective rank within-seed partial analysis
- Partial out seed variance: test rho(eff_rank, AUROC) within each seed separately (12 layer obs).
- Also: test whether layer-normed effective rank (relative to L0) predicts AUROC drop-off.
- Family: intrinsic_dimensionality | novelty: refinement of H03

### Novel hypothesis (P3): Circuit-restricted Laplacian alignment
- Restrict to the 295 circuit genes (TF+target) only — denser subgraph, potentially better aligned.
- Recompute Laplacian on 295×295 adjacency and compare to SV5-7 projections of those 295 genes.
- Z-score for H02 was 7.81 but absolute alignment weak — circuit-restricted may be different.
- Family: module_structure | novelty: new_method

### Novel hypothesis (P4): bZIP/C2H2-ZF family target co-expression test
- H01 found bZIP and C2H2-ZF TFs have negative margin (near target centroid).
- Hypothesis: these TF families regulate many of their own targets in the dataset, explaining proximity.
- Test: TRRUST degree distribution for bZIP vs Forkhead/STAT families; co-expression proxy.
- Family: null_sensitivity | novelty: new_family

---

## Retired directions (do not revisit)
- SV basis permutation (H02 iter_0057: definitively negative)
- Persistent homology on circuit genes (iter_0052: negative)
- SV8-14 subspace structure (deprioritized)
- Raw graph rewiring null for topology (iter_0050: retired)

---

## Key open question
**Why does effective rank correlate with AUROC?** Two possibilities:
1. High-dimensional diffuse representations preserve fine-grained gene identity signals (general capacity).
2. Specific TF/target-separating subspaces (SV2-7) are more orthogonal in higher-rank layers, enabling better LR discrimination.
The SV2-7 effective rank has weaker correlation (rho=0.539), suggesting (1) — global dimensionality, not just the TF subspace. Next iteration should disentangle these.
