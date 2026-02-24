# Brainstormer Hypothesis Roadmap — iter_0006 → iter_0007

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| Rewiring-null survival for H1 persistent homology | **RETIRE** | 0/24 significant across 4 rewiring variants (historical iters 0006–0008). No rescue potential. |
| Geneformer cross-model alignment | **DEPRIORITIZE** | Requires residual-level Geneformer embeddings not in workspace. Low immediate yield. |
| TRRUST co-target clustering (n=209) | **RETIRED** (prior iter) | Too small gene pool. |
| GO cosine distance clustering (n=209) | **RETIRE** | 0/12 significant; retire in 209-gene form. |

---

## New Hypothesis Portfolio

### A. Direct follow-ups from iter_0006 (high confidence, clear test)

**A1. SV1 feature-shuffle null — biological specificity confirmation**
- *Hypothesis*: The extracellular-vs-nuclear GO enrichment on SV1 loadings is not an artifact of gene-label assignment; it reflects genuine spectral organization.
- *Test*: Shuffle which gene name maps to which embedding index; recompute SVD and SV1 GO enrichment. Run 100 permutations; compute empirical p.
- *Signal if true*: Empirical p < 0.05 (true enrichment disappears in shuffle, or true OR >> shuffle mean OR).
- *Null*: Identical OR in shuffled null → enrichment is positional artifact.
- *Value*: high | *Cost*: low

**A2. Layer-wise SVD trajectory — when does SV1 dominance emerge?**
- *Hypothesis*: SV1 dominance accumulates rapidly in layers 1–4 (matching ER curvature breakpoint), then plateaus; the extracellular-vs-nuclear axis is already resolved by layer 4.
- *Test*: Compute SVD for all 12 layers. Track SV1/(SV1+SV2) ratio and top-5-SV cumvar. Compute GO enrichment of SV1 loadings per layer or at key breakpoint layers (0, 1, 4, 11).
- *Signal if true*: SV1 dominance ratio rises sharply layers 1–4, then saturates. GO enrichment for extracellular axis appears by layer 4.
- *Null*: Linear SV1 dominance accumulation across all 12 layers (no breakpoint).
- *Value*: high | *Cost*: medium

**A3. SV1 bottom-quartile enrichment — nuclear/regulatory complement**
- *Hypothesis*: Genes at the negative pole of SV1 are enriched in nuclear/DNA-binding/transcription-factor GO terms, completing the full extracellular→nuclear biological axis interpretation.
- *Test*: From existing h02_svd_gene_projections.csv, take bottom-25% SV1 loading genes. Run Fisher exact GO enrichment vs top-25% (already done for extracellular). Focus on GO:CC nuclear terms and GO:MF DNA-binding/TF terms.
- *Signal if true*: nuclear/TF GO terms significantly enriched at negative SV1 pole (OR > 2, p < 0.01).
- *Null*: No GO enrichment on negative SV1 pole; axis is one-sided.
- *Value*: high | *Cost*: low (data already on disk)

**A4. Drift enrichment feature-shuffle null**
- *Hypothesis*: GO enrichment of high-drift genes (TF activity, OR=3.30) is not an artifact of the named-gene selection from the 4803-gene vocabulary.
- *Test*: Randomly permute which gene names map to which embedding indices. Recompute L2 drift for 209 named genes in permuted context. Rerun GO enrichment top-50 vs bottom-50. Repeat 100×. Compare true OR to null distribution.
- *Signal if true*: True OR = 3.30 >> mean null OR (~1.0); empirical p < 0.05.
- *Null*: True OR ≤ 95th percentile of null ORs → enrichment is selection artifact.
- *Value*: high | *Cost*: low

### B. New geometric/manifold directions

**B1. Spectral gap dynamics — SV1/SV2 ratio across layers**
- *Hypothesis*: The 7.7× SV1/SV2 gap at layer 11 is a late-emerging phenomenon that correlates with ER compression; early layers have a small spectral gap that grows monotonically.
- *Test*: Compute SV1/SV2 ratio from full-vocab SVD at each of 12 layers. Compute Pearson/Spearman correlation with ER trajectory. Test whether ratio growth is front-loaded (matching ER) or back-loaded.
- *Signal if true*: Spectral gap ratio monotonically increases; Pearson r > 0.85 with ER; front-loaded growth pattern.
- *Null*: Spectral gap is already large at layer 0 or grows uniformly.
- *Value*: high | *Cost*: low (needs per-layer SVD, already planned in A2)

**B2. Persistent homology on SV1 projection — does topology survive spectral collapse?**
- *Hypothesis*: After projecting onto the top-k SV subspace (e.g., top-2 or top-5 SVs), persistent H1 loops are absent or suppressed relative to original embedding space, confirming that the topological structure lives in the complement of the dominant subspace.
- *Test*: Project 209 named genes onto (a) SV1-only, (b) SV2-SV5, (c) residual (remaining 507 dims). Compute PH H1 Betti persistence for each. Compare to full-embedding H1. Use feature-shuffle null per subspace.
- *Signal if true*: H1 signal is concentrated in the residual subspace, not in SV1. Top-SV projection kills topology.
- *Null*: H1 signal is invariant to projection choice.
- *Value*: medium | *Cost*: medium

**B3. Gene neighborhood rank correlation across layers — manifold deformation tracking**
- *Hypothesis*: Gene-neighbor relationships are significantly reorganized between layers 0 and 4 (the breakpoint regime) but minimally reorganized from layers 4 to 11.
- *Test*: For each gene, compute rank correlation of its k-nearest-neighbor distances at each pair of consecutive layers. Aggregate Spearman r by layer-pair. Test whether layer-pair (0,4) has lower mean Spearman than (4,11).
- *Signal if true*: Mean rank correlation (0→4) << mean rank correlation (4→11), consistent with front-loaded neighborhood reorganization.
- *Null*: Uniform rank correlation across all consecutive layer-pairs (CKA ≈ 1.0 is consistent with small updates but not full rank-preservation).
- *Value*: medium | *Cost*: low

**B4. Local intrinsic dimension variation across genes — which genes occupy locally higher-D neighborhoods?**
- *Hypothesis*: Regulatory/TF genes (high-drift cluster) occupy locally higher-dimensional neighborhoods than secreted/cytokine genes at layer 11, reflecting encoding complexity differences.
- *Test*: Compute per-gene local TwoNN ID using k=20 neighbors at layer 11. Compare median local-ID for (a) high-drift genes, (b) top SV1 genes (extracellular), (c) bottom SV1 genes (nuclear/TF). Mann-Whitney U test.
- *Signal if true*: Nuclear/TF genes have significantly higher local ID than secreted genes (p < 0.05), reflecting multi-dimensional TF interaction space vs simpler extracellular localization signals.
- *Null*: No local ID difference between biological subgroups.
- *Value*: medium | *Cost*: low

**B5. Gene covariance matrix spectral evolution — Marchenko-Pastur comparison**
- *Hypothesis*: The eigenvalue distribution of the gene-gene covariance matrix transitions from Marchenko-Pastur (random noise) at layer 0 toward a few-outlier-dominant spectrum at layer 11, quantifying the emergence of structured co-representation.
- *Test*: Compute covariance matrix of the 209 named gene embeddings at each layer. Compare eigenvalue distribution to MP distribution with matched parameters. Compute number of eigenvalues exceeding MP bulk upper edge. Track across 12 layers.
- *Signal if true*: Number of bulk-exceeding eigenvalues increases across layers, especially in layers 1–4; near-rank-1 at layer 11 shows single outlier eigenvalue >> all others.
- *Null*: Eigenvalue distribution remains MP-like across all layers.
- *Value*: high | *Cost*: medium

**B6. Betti-0 (connected components) evolution under adaptive filtration**
- *Hypothesis*: scGPT layer embeddings show progressive connected-component consolidation (Betti-0 decreases across layers) under Rips filtration, reflecting increasing representational coherence.
- *Test*: Compute Betti-0 persistence (number of components at representative filtration radii) for 209 gene embeddings at each of 12 layers. Test whether Betti-0 decreases monotonically. Use feature-shuffle as null.
- *Signal if true*: Betti-0 decreases from layer 0 to layer 11; observed Betti-0 at layer 11 is significantly below feature-shuffle null.
- *Null*: No layer-monotone Betti-0 trend.
- *Value*: medium | *Cost*: medium

### C. Biological anchoring / mechanistic directions

**C1. SV2 biological characterization — what does the secondary axis encode?**
- *Hypothesis*: SV2 (singular value 96, 8× smaller than SV1) encodes a biologically interpretable secondary functional axis distinct from the extracellular/nuclear SV1 axis.
- *Test*: From existing h02_svd_gene_projections.csv, compute GO enrichment of top/bottom quartile genes by SV2 loading. Test SV2 vs SV3 for known cell-type markers (B-cell, NK, monocyte) using Fisher exact.
- *Signal if true*: SV2 top genes enriched in a coherent GO functional cluster (e.g., immune effector function, cell-surface receptor activity) distinct from SV1.
- *Null*: SV2 has no significant GO enrichment (noise subspace).
- *Value*: medium | *Cost*: low (data already on disk)

**C2. Layer-wise gene cluster evolution — do discrete cell-type-specific gene groups emerge?**
- *Hypothesis*: k-means clustering (k=4–8) of gene embeddings at layers 0, 4, 8, 11 yields increasingly enriched cell-type-marker clusters; the clusters at layer 11 better separate B-cell, NK, monocyte, and secreted-factor gene groups than at layer 0.
- *Test*: Run k-means (k=6) on 209 named gene embeddings at layers 0, 4, 8, 11. For each cluster, test cell-type-marker gene enrichment (CD3D/CD3E/CD3G = T-cell; CD19/CD79A/MS4A1 = B-cell; CD14/LYZ = monocyte; etc.) using Fisher exact. Compare cluster purity across layers.
- *Signal if true*: Layer 11 clusters have higher cell-type marker enrichment (Fisher OR > 3, p < 0.01) than layer 0; improvement is front-loaded (layer 0→4 bigger gain than 4→11).
- *Null*: Cluster cell-type enrichment does not increase across layers.
- *Value*: high | *Cost*: medium

**C3. STRING network distance vs embedding distance correlation — biological network structure encoding**
- *Hypothesis*: At layer 11, genes that are connected in STRING protein-protein interaction network (high confidence ≥ 700) are significantly closer in embedding space than disconnected gene pairs.
- *Test*: Load STRING high-confidence interactions for the 209 named genes. Compute pairwise L2 distances in layer-11 embeddings. Test whether mean distance for STRING-connected pairs < STRING-disconnected pairs (Mann-Whitney U). Run permutation null (shuffle gene labels).
- *Signal if true*: STRING-connected pairs have significantly lower embedding distance (p < 0.01, effect size > 0.2 Cohen d). Permutation null p < 0.05.
- *Null*: No distance difference between connected and disconnected STRING pairs.
- *Value*: high | *Cost*: medium

**C4. scGPT embedding SV1 vs single-cell cell-type proportions — macro-level biological alignment**
- *Hypothesis*: The SV1 axis ordering of genes (extracellular→nuclear) correlates with the average expression rank of those genes across Tabula Sapiens immune cell types; secreted genes should be highest in secretory cells, nuclear regulators in progenitors.
- *Test*: Use mean cell-type expression profiles for the 209 genes (from Tabula Sapiens data if available, else from public single-cell references). Compute Spearman correlation between SV1 loading and mean expression in secretory cell types (plasma cells, mast cells) vs progenitors (HSC, common lymphoid progenitor).
- *Signal if true*: Spearman r > 0.3 (p < 0.01) between SV1 loading and secretory/progenitor expression ratio.
- *Null*: No correlation between SV1 loading and expression profile.
- *Value*: high | *Cost*: medium

### D. Algorithmic / model-mechanistic signatures

**D1. Layer-wise attention concentration — do early layers (1–4) show broad vs focused attention?**
- *Hypothesis*: Layers 1–4 (the compression regime) have higher attention entropy (broader attention distribution) than layers 5–11, consistent with global structure reorganization; later layers narrow attention to refine fine-grained features.
- *Test*: If scGPT attention weight tensors are accessible, compute per-head entropy across gene positions for each layer. Compare mean entropy layers 1–4 vs 5–11. Mann-Whitney U test.
- *Signal if true*: Mean attention entropy in layers 1–4 significantly higher than layers 5–11 (p < 0.05).
- *Null*: Uniform attention entropy across all layers.
- *Value*: medium | *Cost*: medium (requires attention weight access)

**D2. Gradient-weighted embedding saliency — which embedding dimensions carry the H1 signal?**
- *Hypothesis*: The H1 persistent homology signal is concentrated in a sparse subset of the 512 embedding dimensions, identifiable by gradient-based or PCA-based variance attribution.
- *Test*: For 209 named genes at layer 11, compute the top-20 embedding dimensions by variance. Test whether PH computed on these 20 dimensions reproduces the full-space H1 signal (Fisher test vs feature-shuffle null). Compare with bottom-20 variance dimensions.
- *Signal if true*: Top-20 variance dimensions reproduce H1 signal; bottom-20 do not.
- *Null*: H1 signal requires full 512 dimensions.
- *Value*: medium | *Cost*: low

---

## Top 3 for Immediate Execution

### Candidate 1 — High-Probability Discovery (SV1 null validation + layer trajectory)
**Execute H: A1 + A2 + A3 combined as one script**
- A1 (SV1 feature-shuffle null): Confirms biological specificity of the strongest result.
- A2 (layer-wise SVD trajectory): Fills the critical mechanistic gap — when does the extracellular axis emerge?
- A3 (SV1 bottom-quartile nuclear enrichment): Completes the biological axis interpretation.
- Together these three tests can produce a complete, null-validated, mechanistically resolved account of the dominant spectral structure.
- **All three use data already on disk; cost is low-medium.**

### Candidate 2 — High-Risk / High-Reward (Marchenko-Pastur spectral comparison)
**Execute B5 (covariance spectral evolution)**
- If the eigenvalue spectrum transitions from MP bulk to outlier-dominated, this gives a quantitative statistical mechanics framing for the spectral collapse — publishable and novel.
- Risk: the 209-gene n may produce noisy MP comparison; the test may need the full 4803-gene vocabulary for reliable MP fitting.
- **Reward**: unique framing connecting random matrix theory to representation learning.

### Candidate 3 — Cheap Broad Screen (SV2 characterization + drift shuffle null)
**Execute C1 + A4 combined**
- A4 (drift enrichment feature-shuffle null): 30 min to code, reuses existing GO enrichment code, decisively validates the H01 claim.
- C1 (SV2 GO enrichment): reuses h02_svd_gene_projections.csv, costs ~10 min to code.
- Combined, these two address remaining validation gaps without requiring new data collection.
