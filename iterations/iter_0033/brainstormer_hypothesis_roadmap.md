# Brainstormer Hypothesis Roadmap: iter_0033 → iter_0034+

**Date:** 2026-02-23

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| STRING PPI → L2 distance | Multiple definitive negatives across iterations | **retire_now** |
| Dorothea activation vs repression sign-split | Blocked: available files lack sign/mor column | **retire_now** (unless new data obtained) |
| Raw Betti number counting (H0/H1 PH) without biological anchoring | Signal was weak and hard to interpret without grounding | **deprioritize** (rescue only if combined with bio anchor) |
| Intrinsic dimensionality PR as standalone | Finding established (6.1× collapse L0→L11); no new hypothesis value | **deprioritize** (use as covariate only) |
| Cross-model CKA without matched gene vocabulary | Consistently inconclusive due to vocabulary mismatch | **retire_now** (revive only with confirmed 195-gene overlap) |

---

## New Hypothesis Portfolio

### H-A: Community 0 sub-identity — T-cell vs myeloid subcommunities
**Hypothesis:** The L11 k=10 Community 0 (n=107, non-B-cell) contains a T-cell subcommunity and a myeloid subcommunity separable at k=15 (3-way partition).
**Test:** Run k=15 partition at L11. Test T-cell markers (n=15 in vocab) and myeloid markers (n=6 in vocab) against each resulting community via Fisher exact.
**Expected signal:** T-cell markers enriched in one sub-community, myeloid in another, OR≥5, p<0.05.
**Null:** Random permutation of community labels.
**Value:** high | **Cost:** low

### H-B: Layer-wise B-cell community purity emergence
**Hypothesis:** The B-cell-dominant community becomes progressively purer (higher B-cell marker OR and lower p-value) as layers progress from L0 to L11, tracking the monotonic participation-ratio collapse.
**Test:** At each layer L0–L11, run k=10 greedy-modularity partition; identify the community with highest B-cell marker enrichment; compute OR and Fisher p per layer. Spearman rho(layer, OR).
**Expected signal:** Monotonic increase in OR with layer depth; Spearman rho>0.7, p<0.05 (n=12 layers).
**Null:** Shuffle gene-to-community assignments within each layer.
**Value:** high | **Cost:** medium

### H-C: PC2 axis — T-cell or myeloid identity
**Hypothesis:** PC2 (second SVD axis) of the 195-gene embedding encodes T-cell or myeloid identity, with the relevant markers enriched at one pole of PC2 at L11.
**Test:** Compute SVD top-3 axes per layer. Test T-cell markers (n=15), myeloid (n=6), NK markers at each pole of PC2 and PC3 via Mann-Whitney and AUROC. Test all layers.
**Expected signal:** AUROC<0.35 or >0.65 for at least one cell-type group on PC2/PC3, p<0.05.
**Null:** Background genes (all 195); random marker sets of same size.
**Value:** high | **Cost:** low

### H-D: Regulatory TF-target community co-membership
**Hypothesis:** TF-target gene pairs from Dorothea high-confidence network co-segregate into the same community more often than random gene pairs, especially at late layers.
**Test:** For 205 Dorothea high-conf pairs in 195-gene vocab: compute fraction in same community per layer. Compare to 1000 random same-size gene pairs (AUC over community co-membership). Spearman rho(layer, co-membership rate).
**Expected signal:** Co-membership rate > random for late layers (L8–L11), Spearman rho>0.5.
**Null:** Random pairs from same gene set.
**Value:** high | **Cost:** low

### H-E: Cross-model validation — Geneformer PC1 B-cell polarity
**Hypothesis:** Geneformer gene embeddings show the same B-cell negative PC1 axis as scGPT, indicating this is a universal property of transformer-based gene models.
**Test:** Load Geneformer embeddings for matched 195-gene vocabulary. Compute SVD PC1. Test B-cell markers vs background at negative pole. Report AUROC and p-value. Compare scGPT vs Geneformer PC1 direction cosine similarity.
**Expected signal:** B-cell markers at negative PC1 pole in Geneformer (AUROC<0.35, p<0.05); PC1 cosine similarity between models > 0.5.
**Null:** Permuted gene labels; random same-size marker sets.
**Value:** high | **Cost:** medium (contingent on Geneformer embedding availability)

### H-F: PC1 axis is interpretable as the first discriminant axis for cell-type classification
**Hypothesis:** PC1 directions across layers converge toward the direction that best separates B-cell from non-B-cell genes, measurable as increasing linear discriminant AUROC with layer depth.
**Test:** At each layer, compute 1D LDA (project gene embeddings onto direction maximizing B-cell vs non-B-cell separation). Compare LDA direction cosine similarity with PC1. Compute LDA AUROC per layer; test Spearman rho(layer, AUROC).
**Expected signal:** rho>0.6, p<0.05; LDA-PC1 cosine similarity increases with layer.
**Null:** Shuffle B-cell vs non-B-cell labels before LDA.
**Value:** high | **Cost:** low

### H-G: Community structure in higher-k regime — does a 5-community partition map to 5 immune cell types?
**Hypothesis:** At k=5 (5 communities at L11 per H02), the 5 communities map to B-cell, T-cell, myeloid, NK, and TF-heavy groups.
**Test:** Use the k=5 partition (already available from H02 data). Test each community for B-cell, T-cell, myeloid, NK, TF enrichment via Fisher exact.
**Expected signal:** At least 3 out of 5 communities show Fisher p<0.05 for one cell-type marker set.
**Null:** Random community assignment of same size.
**Value:** medium | **Cost:** low (k=5 partition already computed in H02)

### H-H: Filtration of kNN graph by edge weight — persistence of B-cell cluster
**Hypothesis:** When the kNN graph is filtrated by edge weight (cosine similarity threshold swept from 0 to 1), the B-cell community persists across a wider filtration range than expected under null, indicating a denser, more coherent cluster.
**Test:** Build weighted kNN graph at L11. Sweep edge-weight thresholds. At each threshold, compute fraction of B-cell markers still connected in one component vs fragmented. Compute "persistence length" (range of thresholds over which B-cell cluster remains connected). Compare to 1000 random same-size gene sets.
**Expected signal:** B-cell persistence length z-score > 2.
**Null:** Random gene sets matched for vocabulary membership.
**Value:** medium | **Cost:** medium

### H-I: Intrinsic dimension estimate per community — B-cell cluster is lower-dimensional
**Hypothesis:** The B-cell community genes lie in a lower-dimensional submanifold than Community 0 genes, reflecting a more stereotyped/constrained representation.
**Test:** Compute participation ratio (PR) separately for genes in Community 1 vs Community 0 at each layer. Test if PR(Community 1) < PR(Community 0) significantly.
**Expected signal:** PR(B-cell community) < PR(non-B-cell) at multiple layers; Mann-Whitney p<0.05.
**Null:** Random balanced gene splits matched for community size.
**Value:** medium | **Cost:** low

### H-J: GO biological process enrichment of community 1 vs community 0
**Hypothesis:** Community 1 genes (B-cell community, n=88) are significantly enriched for immune-related GO biological process terms relative to Community 0 (n=107), confirming biological coherence beyond known B-cell marker genes.
**Test:** Run Fisher exact test for GO BP term membership in Community 1 vs Community 0. Use 195-gene as background. Report top 10 enriched terms in each community.
**Expected signal:** Community 1 enriched for B-cell differentiation, lymphocyte activation, immunoglobulin production (FDR < 0.05).
**Null:** Random community assignment matched for sizes.
**Value:** medium | **Cost:** low (pure annotation lookup)

### H-K: Early-layer (L0) B-cell PC1 vs gene expression level correlation
**Hypothesis:** The B-cell PC1-negative signal at L0 is driven by gene expression level (highly expressed genes in B-cells cluster at one pole), not cell-type geometry. Later layers preserve the signal but for biological reasons.
**Test:** Correlate gene mean expression level in B-cell samples (from scRNA-seq metadata) with PC1 score at L0 and L11. If correlation is high at L0 but lower at L11, expression drives early geometry.
**Expected signal:** |rho(expr, PC1)| > 0.4 at L0 but < 0.2 at L11, confirming progressive decoupling.
**Null:** Shuffle gene labels; random gene subsets.
**Value:** medium | **Cost:** medium (requires expression data)

### H-L: Geodesic distance (shortest-path on kNN graph) vs Euclidean distance for B-cell markers
**Hypothesis:** B-cell marker genes are closer to each other in geodesic (graph shortest-path) distance than Euclidean distance relative to null, indicating they lie along a manifold corridor rather than being Euclidean near-neighbors.
**Test:** Compute pairwise geodesic distances on kNN graph at L11. Compare mean pairwise geodesic distance for B-cell gene pairs vs random pairs. Also compute Euclidean distance. Report ratio (geodesic/Euclidean) for B-cell vs random.
**Expected signal:** Geodesic/Euclidean ratio is smaller for B-cell pairs than random, p<0.05.
**Null:** Random gene pairs matched for Euclidean distance.
**Value:** medium | **Cost:** medium

---

## Top 3 for Immediate Execution

### #1 — High-probability discovery candidate
**H-B: Layer-wise B-cell community purity emergence**
- Directly extends the strongest current finding (2-community structure at L11).
- Uses already-computed infrastructure (kNN + greedy modularity per layer).
- If Spearman rho(layer, B-cell OR) is significant, this provides a mechanistic narrative: community consolidation = progressive cell-type disambiguation in the forward pass.
- Expected: 1–2 hrs compute, low overhead.

### #2 — High-risk/high-reward candidate
**H-E: Cross-model validation — Geneformer PC1 B-cell polarity**
- If Geneformer shows the same B-cell/PC1-negative axis, the finding becomes model-class universal, greatly elevating publishable significance.
- Contingent on Geneformer embedding availability in the reference root. Risk: embeddings may not be readily usable for 195-gene vocab without additional preprocessing.
- If succeeds: transforms paper from "scGPT observation" to "general transformer gene model property."

### #3 — Cheap broad-screen candidate
**H-C: PC2/PC3 axis — T-cell or myeloid identity**
- Re-uses already-computed SVD (from H03); only new work is testing PC2 and PC3 instead of PC1.
- Near-zero additional compute.
- If T-cell or myeloid markers appear on PC2/PC3 at any layer, this immediately extends the framework to a multi-axis cell-type geometry model.
- Expected: <30 min additional analysis on existing SVD results.
