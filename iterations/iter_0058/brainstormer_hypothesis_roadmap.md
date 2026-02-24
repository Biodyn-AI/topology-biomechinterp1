# Brainstormer Hypothesis Roadmap — iter_0058

**Date**: 2026-02-23

---

## Retire / Deprioritize

| Hypothesis | Reason | Action |
|-----------|--------|--------|
| SV basis permutation (seed43) | Definitively ruled out by principal angle analysis (PA=74.8° cross-subspace) | RETIRED |
| SV rotation angle → AUROC | rho=-0.27, p=0.42 (iter_0056) | RETIRED |
| PH on full circuit genes | Negative (iter_0052) | RETIRED |
| Raw graph rewiring null for topology | iter_0050 retired | RETIRED |
| TRRUST Laplacian alignment (whole-graph) | Z=7.81 but Δ=3° absolute vs random; near-orthogonal geometry regardless | RESCUE-ONCE: circuit-only subgraph; retire if still weak |
| SV8-14 subspace structure | No positive signal found, high cost | DEPRIORITIZED |

---

## New Hypothesis Portfolio

### A. Refinement of Active Strong Threads

**A1. Cross-Seed Boundary Anchor Stability**
- *Hypothesis*: The top-10 TF boundary anchors (BCL11A, NFKB1, RB1, FOXO3, ZEB1) are stable across seeds in the joint 6D SV2-7 space at L2/L3.
- *Test*: Compute per-TF margin (dist_to_target_centroid − dist_to_TF_centroid) in seed43 and seed44 at L2 and L3. Spearman rank correlation of TF margin rankings across seeds.
- *Expected signal*: Spearman rho > 0.7 for main vs seed43 and main vs seed44 margin rankings.
- *Null*: Random permutation of TF-to-margin assignment.
- *Value*: high | *Cost*: low

**A2. Effective Rank Mechanistic Decoupling (Layer vs Dimensionality)**
- *Hypothesis*: The effective rank → AUROC correlation (rho=0.855) is driven by global dimensionality independent of layer depth; i.e., partial correlation controlling for layer index remains significant.
- *Test*: Partial Spearman correlation between eff_rank and AUROC with layer index partialled out across 36 (seed × layer) observations. Also: test within-seed rank correlation (12 observations per seed separately).
- *Expected signal*: Partial rho remains > 0.5 and significant.
- *Null*: Full correlation driven by layer index alone (both eff_rank and AUROC are monotone in layer depth).
- *Value*: high | *Cost*: low

**A3. bZIP/C2H2-ZF Proximity Explanation (TRRUST Degree)**
- *Hypothesis*: bZIP and C2H2-ZF TFs have negative boundary margin because they are high-degree TRRUST hubs (many targets in the dataset), pulling their embedding toward the target cloud centroid.
- *Test*: Correlate per-TF margin (L2) with TRRUST out-degree (number of TRRUST edges to targets in the nonzero gene set). Separately for bZIP and C2H2-ZF vs other families.
- *Expected signal*: Spearman rho < -0.4 between degree and margin; bZIP/C2H2-ZF families have significantly higher degree.
- *Null*: Margin is uncorrelated with TRRUST degree.
- *Value*: high | *Cost*: low

### B. New Geometric / Topological Directions

**B1. Individual Pair Distance as Edge Probability Predictor**
- *Hypothesis*: Pairwise Euclidean distance between individual TF and target gene embeddings in 6D SV2-7 space at L2 predicts held-out TRRUST edge probability (positive vs negative pairs).
- *Test*: For all TF-target gene pairs: compute 6D SV2-7 distance. Compare distance distributions for positive (589 TRRUST pairs) vs random negative pairs (same TF, random non-target). AUROC of edge prediction from distance alone. Cross-validate with 5-fold stratified on TF identity.
- *Expected signal*: AUROC > 0.65 for distance-based edge prediction; positive pairs closer than negatives.
- *Null*: Shuffled positive labels; distance AUROC ≈ 0.5.
- *Value*: high | *Cost*: medium

**B2. Geodesic Distance on Embedding Manifold (kNN Graph)**
- *Hypothesis*: Geodesic distances on the kNN graph of nonzero gene embeddings (L2) better separate TF classes than Euclidean distances in 6D SV2-7.
- *Test*: Build kNN graph (k=10, 15, 20) on [2039, 512] embeddings at L2. Dijkstra geodesic distances between all pairs. Project to 2D via MDS on geodesic distances. Compute TF vs target AUROC in 2D geodesic space vs 6D Euclidean.
- *Expected signal*: Geodesic AUROC > 6D Euclidean AUROC at L2 (improvement ≥ 0.02).
- *Null*: Geodesic AUROC ≤ Euclidean AUROC; shuffled labels.
- *Value*: medium | *Cost*: medium

**B3. Intrinsic Dimension Trajectory via TwoNN**
- *Hypothesis*: The intrinsic dimension (ID) of the nonzero gene embedding manifold decreases with layer depth, and this ID trajectory correlates with AUROC more tightly than effective rank.
- *Test*: Apply TwoNN estimator (or MLE-based) to [2039, 512] embeddings at each of 12 layers, 3 seeds. Spearman(ID, AUROC) across 36 observations. Compare rho to eff_rank rho=0.855.
- *Expected signal*: TwoNN ID and effective rank highly correlated; if TwoNN rho > 0.855, it's a stronger predictor (or the measures are distinct and complementary).
- *Null*: TwoNN ID uncorrelated with AUROC (rho ≈ 0).
- *Value*: medium | *Cost*: low

**B4. Persistent Homology on SV2-7 Subspace (Not Full Space)**
- *Hypothesis*: PH on the 6D SV2-7 projection of nonzero genes reveals a topological feature (H0 gap = cluster separation, or H1 loops) that correlates with AUROC across layers.
- *Test*: At each layer, compute PH (Vietoris-Rips) on the [2039, 6] SV2-7 projected gene cloud. Extract H0 gap (second − first death time) and H1 Betti number. Spearman(PH_feature, AUROC) across 12 layers.
- *Expected signal*: H0 gap (TF-cluster vs target-cluster separation) correlates with AUROC (rho > 0.6).
- *Null*: PH features uncorrelated with AUROC; layer-shuffle control.
- *Value*: medium | *Cost*: medium

**B5. Local Covariance Structure: Fisher Information Proxy**
- *Hypothesis*: At L2 (peak AUROC), TF gene embeddings have lower local covariance (more concentrated in the 6D SV2-7 subspace) than target-only genes, reflecting greater regularity of TF representations.
- *Test*: For each class (TF, target-only), compute covariance matrix in 6D SV2-7 at L2. Compare trace (total variance), condition number, and effective rank of each class covariance matrix. Test if TF class is more "concentrated" (lower condition number).
- *Expected signal*: TF class has lower condition number and lower covariance trace than target-only class.
- *Null*: No difference between class covariance structures.
- *Value*: medium | *Cost*: low

### C. Cross-Model and Biological Anchoring

**C1. scGPT vs Geneformer Effective Rank Alignment**
- *Hypothesis*: If Geneformer embeddings are available for the same immune gene set, the effective rank → AUROC relationship holds in Geneformer, making this a model-agnostic transformer property.
- *Test*: Extract Geneformer layer embeddings for the 2039 nonzero genes. Compute effective rank per layer. Compute TF/target AUROC in Geneformer 6D analogous subspace. Spearman correlation.
- *Expected signal*: Geneformer rho(eff_rank, AUROC) > 0.6.
- *Null*: scGPT-specific; Geneformer rho ≈ 0.
- *Value*: high | *Cost*: high (requires Geneformer inference)

**C2. TF Boundary Anchor Prediction from Gene Ontology**
- *Hypothesis*: The TF boundary margin (6D distance from target centroid) is predictable from GO biological process annotations — specifically, TFs annotated to "regulation of immune response" or "lymphocyte differentiation" have higher margin.
- *Test*: Fetch GO:BP terms for all 73 TF genes (via goatools or mygene). Binary feature: annotated to immune regulation terms. Mann-Whitney U test: margin of annotated vs unannotated TFs. Also: regression of margin on number of immune-related GO terms.
- *Expected signal*: Immune-annotated TFs have significantly higher margin (p < 0.01).
- *Null*: Margin uncorrelated with GO annotation category.
- *Value*: high | *Cost*: low

**C3. STRING Protein Interaction Network Alignment**
- *Hypothesis*: The 6D SV2-7 gene embedding distances are more correlated with STRING protein-protein interaction network distances than TRRUST distances, because STRING captures broader functional similarity.
- *Test*: Fetch STRING interaction scores for the 2039 nonzero genes. Compute graph distances (BFS) or similarity (co-occurrence) from STRING. Correlate with 6D Euclidean distances in SV2-7 at L2. Compare correlation strength to TRRUST Laplacian alignment.
- *Expected signal*: Mantel test or Spearman correlation between STRING and embedding distances is stronger than TRRUST correlation.
- *Null*: STRING distances uncorrelated with SV2-7 distances; no improvement over TRRUST.
- *Value*: medium | *Cost*: medium

**C4. Circuit-Restricted Laplacian Alignment (TRRUST Rescue)**
- *Hypothesis*: Restricting to the 295 circuit genes (TF + their TRRUST targets) and computing the denser TRRUST subgraph Laplacian shows stronger alignment with SV5-7 subspace projections of those 295 genes.
- *Test*: Build [295×295] TRRUST adjacency restricted to circuit genes. Compute Laplacian eigenvectors 1–9. Extract SV5-7 projections of only the 295 circuit genes at each layer. Compute principal angles. Compare to whole-graph baseline (H02, best PA=84.88°).
- *Expected signal*: Circuit-restricted PA < 80° (5+ degrees improvement over whole-graph best).
- *Null*: Circuit-restricted PA ≈ 84–88° (same as whole-graph); no improvement.
- *Value*: medium | *Cost*: low

### D. Mechanistic and Algorithmic Signatures

**D1. Attention Head Contribution to SV2-7 Subspace**
- *Hypothesis*: Specific attention heads in scGPT (at L2/L3) are responsible for writing the TF/target-separating SV2-7 geometry into the residual stream.
- *Test*: If scGPT attention weights are accessible: compute per-head residual contribution (attention × value matrix) to each gene embedding. Correlate per-head contribution magnitudes with TF boundary margin (H01). Identify top-contributing heads.
- *Expected signal*: A small subset of heads (< 4 out of all heads) account for > 60% of SV2-7 variance.
- *Null*: All heads contribute equally; no head-to-margin correlation.
- *Value*: high | *Cost*: high (requires model internals access)

**D2. Spectral Gap of Normalized Laplacian vs AUROC**
- *Hypothesis*: The spectral gap (second eigenvalue of normalized Laplacian) of the kNN graph built on gene embeddings correlates with TF/target AUROC across layers — capturing the clustering structure that drives classification.
- *Test*: Build kNN graph on [2039, 512] embeddings at each layer. Compute normalized Laplacian. Extract second eigenvalue (Fiedler value). Spearman(Fiedler, AUROC) across 12 layers × 3 seeds = 36 observations.
- *Expected signal*: Fiedler value correlates with AUROC (rho > 0.5); larger spectral gap = more separated clusters.
- *Null*: Fiedler value uncorrelated with AUROC.
- *Value*: medium | *Cost*: medium

**D3. Lipschitz Regularity of Layer Map (Embedding Smoothness)**
- *Hypothesis*: The Lipschitz constant of the embedding transformation from layer L to layer L+1 (estimated from pairwise distances) decreases with depth, and low-Lipschitz layers correspond to higher AUROC.
- *Test*: For consecutive layer pairs (L0→L1, ..., L10→L11): compute empirical Lipschitz estimate = max(|f(x_i)−f(x_j)|/|x_i−x_j|) for random gene pairs. Spearman(Lipschitz_L, AUROC_L) across 11 layer transitions.
- *Expected signal*: Low Lipschitz layers = smoother transformation = higher AUROC (rho < -0.4).
- *Null*: No correlation between Lipschitz constant and AUROC.
- *Value*: medium | *Cost*: low

---

## Top 3 for Immediate Execution

### Candidate 1 — High-Probability Discovery
**A1 + A3 combined: Cross-seed boundary anchor stability + TRRUST degree explanation**

Both can run from existing data in a single script.
- A1: Compute margin in seed43/seed44 at L2/L3 for top-10 TFs; Spearman rank correlation across seeds.
- A3: Load TRRUST edge file; count out-degree per TF; Spearman(degree, margin) for all 73 TFs; compare bZIP/C2H2-ZF vs other families.

*Why*: A1 converts H01 from a single-seed finding to a robust cross-seed claim (elevates to Claim 54). A3 provides mechanistic explanation for the family stratification. Both are low-cost and directly extend existing confirmed results.

### Candidate 2 — High-Risk / High-Reward
**B1: Individual pair distance as edge probability predictor**

This pivots from class-level TF/target separation to individual edge prediction — the step that converts a geometric finding into a functional regulatory inference tool.

*Why*: If pairwise 6D SV2-7 distance predicts individual TRRUST edges (AUROC > 0.65), this is a qualitative leap — from "TFs and targets are geometrically separable" to "embedding distance encodes regulatory relationships." This is publication-level strong if it holds, and the test is straightforward.

### Candidate 3 — Cheap Broad Screen
**A2 + D3 combined: Mechanistic decoupling of effective rank + Lipschitz regularity**

Both run from existing data with minimal compute.
- A2: Partial Spearman(eff_rank, AUROC | layer_index) across 36 observations; within-seed rho for each of 3 seeds separately.
- D3: Pairwise distance ratio for consecutive layer pairs; Spearman(Lipschitz_estimate, AUROC).

*Why*: A2 determines whether rho=0.855 is a causal-mechanistic property or a layer-depth confound — a critical distinction for paper narrative. D3 adds an orthogonal angle on layer-wise transformations. Both complete in one experiment run.
