# Brainstormer Hypothesis Roadmap — iter_0036

---

## Retire / Deprioritize

| Direction | Iterations with negatives | Verdict |
|-----------|--------------------------|---------|
| TRRUST/GO co-annotation → kNN proximity | iter_0033, iter_0035 | `retire_now` |
| Intrinsic dimensionality (SVD AUROC for TF sets) | iter_0033, iter_0035 | `retire_now` |
| GO BP enrichment in SV2 poles | iter_0011 | `retire_now` (already retired) |
| Repression anti-pole geometry | iter_0013 | `deprioritize` (signal absent, no strong rescue path) |
| Persistent homology H1 on general gene sets | iter_0003 | `deprioritize` (positive but not tied to biology yet) |

No `rescue_once_with_major_change` candidates at this time — retired directions had clean nulls at adequate sample sizes.

---

## New Hypothesis Portfolio

### H-A: Geneformer B-Cell kNN Precision@10 Cross-Model Validation
**Hypothesis**: The B-cell geometric clustering signal (kNN precision@10, z>4) is present in Geneformer residual embeddings at the same layers (early > late).
**Test**: Load Geneformer layer embeddings for the same 195 in-vocab genes. Compute precision@k=10 for the 7 B-cell markers across Geneformer layers. Bootstrap null: 500 random same-sized gene sets. Compare z-scores to scGPT z=4.55.
**Expected signal if true**: z > 2 at early Geneformer layers, declining toward final layer. Concordance with scGPT layer trajectory direction.
**Null/control**: Gene-name permutation null (same protocol as H02 iter_0036).
**Value**: `high` | **Cost**: `medium` (requires Geneformer embedding access)

### H-B: SPIB/BACH2/BATF STRING PPI Quantification
**Hypothesis**: The top B-cell neighborhood genes (SPIB, BACH2, BATF) have significantly higher STRING PPI confidence scores with B-cell markers than with randomly drawn in-vocab genes.
**Test**: For each B-cell marker (CD19, CD79A, BLK, etc.), query the existing STRING PPI dataset (from iter_0012, 9905 edges) for mean edge weight to {SPIB, BACH2, BATF} vs. random gene triples from in-vocab. Bootstrap null: 500 random gene triples. Compute z-score.
**Expected signal if true**: STRING mean score to B-cell TFs is significantly above random (z > 2), confirming that geometric neighborhood recapitulates known PPI biology.
**Null/control**: Same-sized random gene triples as "TF" set.
**Value**: `high` | **Cost**: `low` (STRING data already on disk)

### H-C: Expand Cell-Type Panel: NK, DC, Plasma Cells
**Hypothesis**: NK cell markers (NCAM1/CD56, NKG7, GNLY, PRF1, KLRB1) and dendritic cell markers (ITGAX/CD11c, FLT3, CLEC9A) show non-zero kNN clustering, while plasma cell markers (MZB1, JCHAIN, SDC1) cluster near B-cell markers.
**Test**: Find in-vocab genes for NK (target: 5+), DC (target: 3+), plasma cell (target: 3+) marker sets. Compute kNN precision@10 and z-score at L2/L5/L8/L11. Also measure centroid distance between plasma cell and B-cell centroids vs random pairs.
**Expected signal if true**: Plasma cells near B-cells (both are B-lineage). NK/DC clustering non-zero but weaker than B-cell z=4.55.
**Null/control**: Bootstrap random same-sized gene sets.
**Value**: `high` | **Cost**: `low`

### H-D: B-Cell Subspace Persistent Homology (Betti-0)
**Hypothesis**: The 7 B-cell markers form a persistent single connected component in the kNN graph at smaller filtration radii than the full 195-gene set, indicating geometric cohesion.
**Test**: Build kNN graph on 195 in-vocab genes at L2. Apply Betti-0 persistence (connected components as function of edge-weight threshold). Measure the radius at which the 7 B-cell markers merge into one component vs radius at which all 195 genes connect. Null: 500 random 7-gene subsets.
**Expected signal if true**: B-cell subgraph Betti-0 drops to 1 at significantly smaller radius than random sets (z > 2).
**Null/control**: Random 7-gene subsets from the 195 in-vocab genes.
**Value**: `medium` | **Cost**: `low`

### H-E: Layer-Specific Geometry — What Changes at L2→L11
**Hypothesis**: The attenuation of B-cell clustering from L2 (z=4.55) to L11 (z=1.65) reflects a progressive dispersion of B-cell markers into the broader embedding manifold, not random noise drift.
**Test**: Compute B-cell centroid movement (L2 norm of PC1-3 centroid displacement per layer step) and intra-cluster scatter (mean pairwise distance among B-cell markers) per layer. Test whether scatter growth is correlated with z-score decline.
**Expected signal if true**: Monotonic scatter increase (rho > 0.9) correlated with z-score decline; centroid itself stable or moving in one direction.
**Null/control**: Random same-sized gene sets tracked the same way.
**Value**: `medium` | **Cost**: `low`

### H-F: T-Cell Permutation Null Control
**Hypothesis**: T-cell kNN precision z=0.29 at L2 is a genuine null (not an artifact of small gene set or biased embeddings).
**Test**: Apply the same gene-name permutation null (200 permutations) used for B-cells to T-cell markers (n=12). Confirm perm z-scores mean ≈ 0 and real T-cell z=0.29 falls within permutation distribution (z_real_vs_perm < 1).
**Expected signal if true**: Real T-cell z=0.29 is within the permutation null distribution — confirms T-cell non-clustering is genuine, not masking an artifact.
**Null/control**: Permutation distribution itself.
**Value**: `medium` | **Cost**: `low` (code already written for H02 iter_0036)

### H-G: Cross-Model SV2 PPI Co-Pole in Geneformer
**Hypothesis**: The SV2 PPI co-pole signal (scGPT: 12/12 layers, z=3.3–6.5) is also present in Geneformer residual representations.
**Test**: Load Geneformer embeddings for the 209 named genes (same set as Track A). Compute SVD. Test STRING high-confidence PPI co-pole enrichment on Geneformer SV1/SV2/SV3 axes at each layer (same protocol as iter_0012/0013).
**Expected signal if true**: Geneformer has an SV axis (not necessarily SV2) with PPI co-pole z > 3 at multiple layers — same biology, possibly different axis ordering.
**Null/control**: Gene-label shuffles (N=500).
**Value**: `high` | **Cost**: `medium`

### H-H: B-Cell Geometry in Held-Out Marker Set
**Hypothesis**: The B-cell kNN clustering generalizes to a held-out set of B-cell markers not used in the original marker panel.
**Test**: Search full scGPT vocab for additional confirmed B-cell genes (EBF1, IRF4, PAX5, XBP1, BLNK, FCRL4, MS4A1/CD20, FCRLA). Test which are in-vocab. Split into discovery set (original 7) and held-out set (2-4 new genes). Test kNN precision@10 on held-out set alone.
**Expected signal if true**: Held-out B-cell genes also cluster above null (z > 2), confirming generalization beyond the original 7-marker panel.
**Null/control**: Bootstrap random same-sized gene sets.
**Value**: `high` | **Cost**: `low`

### H-I: Signed Regulatory Geometry in B-Cell Subspace
**Hypothesis**: TRRUST TF→target pairs involving B-cell TFs (SPIB, BACH2, BATF) show co-pole enrichment in the B-cell PC1 axis, analogous to the SV2 activation co-pole signal in Track A.
**Test**: Extract TRRUST activation edges where TF or target is in {SPIB, BACH2, BATF, CD19, CD79A, PAX5}. At L2, compute fraction of these pairs with both genes on the negative (B-cell) side of PC1. Null: 500 gene-label shuffles. Compare to all-edge TRRUST co-pole z.
**Expected signal if true**: B-cell TF→target pairs co-localize in PC1 negative pole at z > 2, extending the signed regulatory geometry finding to the cell-type-specific B-cell axis.
**Null/control**: TRRUST edges involving non-B-cell TFs.
**Value**: `medium` | **Cost**: `low`

### H-J: Geometric Distance Matrix Between Cell-Type Centroids
**Hypothesis**: The pairwise centroid distances between cell types (B, T, Myeloid, NK, DC, Plasma) in the embedding space reflect known immune lineage relationships: B-lineage (B, plasma) closer to each other than to myeloid, T closer to NK than to B.
**Test**: For all available cell-type marker sets (≥3 in-vocab genes), compute centroid (mean embedding) at L2 and L11. Build pairwise centroid distance matrix. Compare to hematopoietic lineage tree structure (B/Plasma closer, T/NK closer, Myeloid/DC closer).
**Expected signal if true**: Mantel test between embedding distance matrix and lineage distance matrix: Spearman rho > 0.5, p < 0.05.
**Null/control**: Random gene sets as pseudocentroids.
**Value**: `high` | **Cost**: `low` (if marker sets are found)

### H-K: Multi-Resolution kNN Stability (k=5, 10, 20, 50)
**Hypothesis**: The B-cell kNN precision signal is stable across scales (k=5 to k=50), indicating the clustering is not a k-specific artifact.
**Test**: Compute B-cell precision@k for k ∈ {5, 10, 20, 50} at L2. Bootstrap null at each k. Plot z-score vs k.
**Expected signal if true**: z > 2 at all k values tested, confirming scale-invariant geometric separation.
**Null/control**: Bootstrap random same-sized gene sets at each k.
**Value**: `medium` | **Cost**: `low`

### H-L: Geneformer vs scGPT B-Cell Geometry Rank Correlation
**Hypothesis**: Genes that are geometrically closer to B-cell markers in scGPT are also closer in Geneformer, i.e., both models learn a similar B-cell-anchored geometry.
**Test**: For each of 195 in-vocab genes, compute mean distance to B-cell centroid in scGPT (L2) and in Geneformer (best layer). Spearman correlation of these 195 ranks between models.
**Expected signal if true**: Spearman rho > 0.3, p < 0.001 — both models rank genes similarly by B-cell distance despite different architectures and training data.
**Null/control**: Permuted Geneformer embeddings; expected rho ≈ 0.
**Value**: `high` | **Cost**: `medium`

---

## Top 3 for Immediate Execution

### Candidate 1 — High-Probability Discovery: H-C (Expand Cell-Type Panel: NK, DC, Plasma)
**Rationale**: Low cost, in-vocab search is the bottleneck. If plasma cell markers cluster near B-cell markers, this immediately extends the biological story to a clinically relevant B-lineage differentiation axis. Uses code already written. Very likely to produce at least one positive result. Discovery probability: high.

### Candidate 2 — High-Risk/High-Reward: H-A (Geneformer B-Cell kNN Precision@10)
**Rationale**: If Geneformer also shows B-cell geometric clustering at early layers, the paper's central claim upgrades from "scGPT encodes B-cell geometry" to "genomic foundation models share a conserved B-cell geometric representation." This would be the highest-impact result in the entire subproject. Risk: Geneformer embeddings may not be readily accessible or may require non-trivial extraction.

### Candidate 3 — Cheap Broad Screen: H-B + H-F combined (STRING TF quantification + T-cell perm null)
**Rationale**: Both are near-zero marginal cost (STRING data on disk, permutation null code already written). H-B converts H03 from qualitative to quantitative. H-F completes the negative-control story for T-cells. Can be run as a single script in one iteration slot, leaving the other slot free for H-C or H-A.
