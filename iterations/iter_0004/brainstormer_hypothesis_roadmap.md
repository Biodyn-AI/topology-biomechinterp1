# Brainstormer Hypothesis Roadmap — iter_0004

---

## Retire / Deprioritize

| Hypothesis | Reason | Status |
|-----------|--------|--------|
| kNN transitivity | Strongly negative across all layers/seeds; no rescue potential | **retire_now** |
| TRRUST co-target distance (209-gene pool) | Gene pool too small; same-pool null is inappropriate; 1.0% sig rate | **retire_now** |
| Cross-model feature vector alignment (summary vectors) | Spearman 0.83 but Fisher p = 0.41; low-dim summary features are too coarse | **rescue_once_with_major_change** (rescue = use full residual embedding tensors, CKA/Procrustes, not feature effect vectors) |
| Cross-layer linear CKA | Already fully resolved (CKA ≈ 1.0); further replication adds nothing | **retire_now** (result incorporated, move to derivative questions) |

---

## New Hypothesis Portfolio

### H-A: Per-Gene Layer Displacement (Residual Drift Outliers)
**Hypothesis:** Despite population-level CKA ≈ 1.0, individual genes vary substantially in how much their embedding changes across layers; the high-drift genes are functionally distinct (e.g., TF binding hubs, signaling intermediaries).
**Test:** For each gene, compute L2 norm of (layer_k_embedding − layer_0_embedding) for k=1..11; rank genes by cumulative drift; compare top-10% vs bottom-10% for GO enrichment and known regulatory role.
**Expected signal if true:** High-drift genes are enriched for known regulators/TFs; low-drift genes are housekeeping.
**Null/control:** Permute gene labels across drift ranks; check whether top/bottom 10% GO enrichment persists.
**Value:** high | **Cost:** low

### H-B: GO Term Embedding Clustering (Full Vocabulary)
**Hypothesis:** Genes sharing GO Biological Process annotations cluster more tightly in scGPT embedding space than randomly assembled gene groups of equal size.
**Test:** Use full scGPT gene vocabulary (~4000 genes). Select GO BP terms with 10–100 member genes each. For each term × layer, compute mean pairwise cosine distance within the GO group. Null: 50 random gene groups of matched size. Z-score per term-layer. Report fraction of significant tests and top clusters.
**Expected signal if true:** >5% of GO terms show z < -1.96 in ≥3 layers; biologically coherent terms (metabolic pathways, immune response) show strongest clustering.
**Null/control:** Random size-matched gene groups from the full vocabulary.
**Value:** high | **Cost:** medium

### H-C: MADA Intrinsic Dimensionality Cross-Check
**Hypothesis:** The TwoNN decreasing ID trend (layer 0→11) is replicated by the MADA estimator (which uses a different geometric assumption — local PCA tangent dimension rather than nearest-neighbor ratios), confirming that compression is a real geometric feature rather than a TwoNN artifact.
**Test:** Apply MADA (or equivalent local PCA estimator) to the same 400-gene subsets per layer, seed42. Compare ID profile shape to TwoNN results; test correlation between layer-level ID estimates from both methods.
**Expected signal if true:** MADA ID also decreases layer 0→11, Pearson r > 0.6 with TwoNN profile; z-scores replicated.
**Null/control:** Same feature-shuffle null as H01.
**Value:** medium | **Cost:** low

### H-D: Effective Rank / Spectral Compression Per Layer
**Hypothesis:** The effective rank (exp(entropy of normalized singular values)) of the gene embedding matrix decreases across transformer depth, indicating spectral compression of the representation even when CKA remains near 1.0.
**Test:** For each layer, compute full SVD of the 400 × 512 embedding matrix. Compute effective rank = exp(-Σ p_i log p_i) where p_i = σ_i² / Σσ_j². Track across layers 0–11. Null: same feature-shuffle replicates.
**Expected signal if true:** Effective rank decreases monotonically across depth; observed effective rank is lower than shuffled null at later layers.
**Null/control:** Feature-shuffle replicates of gene embedding matrix.
**Value:** medium | **Cost:** low

### H-E: Kernel CKA (RBF) Cross-Layer — Nonlinear Geometry Test
**Hypothesis:** While linear CKA ≈ 1.0 for all layer pairs, kernel CKA (RBF kernel) reveals subtle nonlinear divergence across layers, particularly between early (0–2) and late (9–11) layers.
**Test:** Compute kernel CKA with RBF kernel (bandwidth = median pairwise distance) for all 12×12 layer pairs. Compare to linear CKA matrix. Report off-diagonal kernel vs linear CKA discrepancy.
**Expected signal if true:** Kernel CKA matrix has more structure than linear (non-uniform, gradient across depth); early-vs-late layers show greater discrepancy.
**Null/control:** Feature-shuffle null; compare against linear CKA baseline.
**Value:** medium | **Cost:** low

### H-F: STRING PPI Distance Correlation
**Hypothesis:** Pairs of genes with high STRING protein-protein interaction scores have smaller embedding distances than gene pairs with low/absent STRING scores, in scGPT embedding space.
**Test:** For all gene pairs in the scGPT vocabulary with STRING score data (>400 confidence), compute Spearman correlation between STRING combined score and embedding cosine similarity per layer. Null: permute STRING scores while preserving embedding distances.
**Expected signal if true:** Spearman r > 0.1 (positive correlation between PPI strength and embedding similarity) in ≥3 layers; signal strongest at later layers.
**Null/control:** Permutation of STRING gene labels; random gene pair selection.
**Value:** high | **Cost:** medium

### H-G: Cell-Type Marker Gene Clustering
**Hypothesis:** Known cell-type marker genes (from PanglaoDB or CellMarker) cluster by cell type identity in scGPT embedding space, even when embeddings are derived from a lung single-cell dataset.
**Test:** Extract marker gene sets for 5–10 major cell types present in lung (T cell, B cell, macrophage, epithelial, endothelial, fibroblast). For each cell type pair, test whether within-type embedding distances < across-type distances. Use Kolmogorov-Smirnov test per layer.
**Expected signal if true:** Cell-type marker clusters are significantly separated (KS p < 0.01) in ≥6 layers; immune cell types separate most cleanly.
**Null/control:** Random reassignment of marker gene labels to cell types of matched size.
**Value:** high | **Cost:** medium

### H-H: H0 Betti (Connected Components) vs Layer
**Hypothesis:** The number of connected components in the Vietoris-Rips filtration (H0 Betti) decreases across transformer layers, indicating progressive integration of gene representation clusters.
**Test:** Use existing Ripser infrastructure from iter_0003 to compute H0 Betti across layers 0–11 on same 350-gene subsets, 3 seeds. Compare to feature-shuffle null. Track H0 at fixed filtration radii.
**Expected signal if true:** H0 component count decreases from layer 0 to layer 11; later layers show tighter connectivity at fixed radius.
**Null/control:** Feature-shuffle null (same as H1 persistence tests in iter_0003).
**Value:** medium | **Cost:** low

### H-I: Cross-Model CKA — scGPT vs Geneformer Matched Genes
**Hypothesis:** scGPT and Geneformer gene embeddings (for matched gene sets) show non-trivial CKA alignment, indicating shared representational geometry despite different training objectives and architectures.
**Test:** Materialize Geneformer gene-level residual embeddings for genes present in both models' vocabularies (~2000 matched genes). Compute linear CKA between scGPT layer representations and Geneformer layer representations of the same genes. Null: CKA between scGPT and randomly permuted Geneformer genes.
**Expected signal if true:** Cross-model CKA > 0.3 for ≥1 layer pair; significantly above gene-permutation null.
**Null/control:** Permute gene indices in Geneformer embedding matrix before CKA.
**Value:** high | **Cost:** high (requires materializing Geneformer residuals)

### H-J: Local Intrinsic Dimensionality by Gene Functional Class
**Hypothesis:** Genes in dense functional modules (e.g., ribosomal proteins, core splicing factors) have lower local intrinsic dimensionality in embedding space than genes in sparse or pleiotropic functions.
**Test:** Segment genes into functional classes (ribosomal, TF, signaling receptor, metabolic enzyme) using GO annotations. For each class, compute TwoNN ID on the subset across layers. Test whether class-specific ID is lower than matched random subsets.
**Expected signal if true:** Ribosomal/core complex genes show ID 2–3 points lower than random; TFs show intermediate ID; pleiotropic genes near random.
**Null/control:** Random size-matched gene subsets.
**Value:** medium | **Cost:** low

### H-K: Geodesic vs Euclidean Distance Ratio
**Hypothesis:** The ratio of geodesic distance (approximate kNN-graph shortest path) to Euclidean distance varies systematically across the gene embedding manifold, identifying non-linearly curved regions.
**Test:** Construct kNN graph (k=10) on gene embeddings per layer. For each gene pair, compute graph geodesic vs Euclidean distance ratio. Report mean, variance, and skewness of this ratio per layer. Null: same ratio distribution on feature-shuffled embeddings.
**Expected signal if true:** Real embeddings have significantly higher geodesic/Euclidean ratio (more curved, non-Euclidean structure) than shuffled; ratio increases at later layers.
**Null/control:** Feature-shuffle null.
**Value:** medium | **Cost:** low

### H-L: Persistence Entropy Per Layer
**Hypothesis:** The entropy of the H1 persistence diagram (sum of -p log p over normalized bar lengths) decreases across transformer layers, indicating progressive regularization of the loop structure.
**Test:** Compute H1 persistence entropy per layer for 3 seeds (using existing Ripser infrastructure). Track across layers 0–11. Null: shuffle replicates from iter_0003.
**Expected signal if true:** Persistence entropy decreases layer 0 → 11; later layers have more uniform, lower-entropy H1 diagrams.
**Null/control:** Feature-shuffle replicates.
**Value:** medium | **Cost:** low

### H-M: Layer-Dependent Neighborhood Overlap (Jaccard)
**Hypothesis:** For each gene, the top-10 neighbors in layer 0 overlap significantly with top-10 neighbors in layer 11, quantifying how stable the local gene neighborhoods are despite potential representation drift.
**Test:** Compute mean Jaccard similarity of k=10 nearest neighbors between layer 0 and layer k, for k=1..11, for all 400 sampled genes. Null: Jaccard for randomly reassigned gene neighborhood memberships.
**Expected signal if true:** Jaccard similarity is high (>0.7) for adjacent layers, decreasing gradually to layer 11; higher than null but not 1.0, indicating stable but not identical neighborhoods.
**Null/control:** Random neighbor assignment preserving degree distribution.
**Value:** medium | **Cost:** low

---

## Top 3 for Immediate Execution

### Slot 1 — High-Probability Discovery Candidate
**H-B: GO Term Embedding Clustering (Full Vocabulary)**

Rationale: The TRRUST test failed purely due to gene pool size (209 genes, correlated null). GO BP annotations cover thousands of genes with well-validated groupings. This is the same biological anchoring question with a properly powered design. If scGPT embeddings carry functional signal, GO clustering is the most likely place to find it given the strong topology results. Immediate, clearly interpretable.

Execution plan:
1. Load full scGPT gene vocabulary list from cycle1 embedding metadata.
2. Download/load GO BP gene sets (geneontology.org slim sets or msigdb C5 collection).
3. Filter to terms with 10–100 member genes present in scGPT vocabulary.
4. For each term × 12 layers, compute mean pairwise cosine distance.
5. Null: 50 random groups of matched size from full vocabulary.
6. Z-score per term-layer; report fraction significant, top terms, layer-profile.

### Slot 2 — High-Risk / High-Reward Candidate
**H-A: Per-Gene Layer Displacement (Residual Drift Outliers)**

Rationale: CKA ≈ 1.0 at population level is a strong result, but per-gene displacement can still vary enormously even when the mean is near zero. This is a mechanistic probe: which genes are actually processed by the transformer blocks vs passed through unchanged? If high-drift genes are enriched for known regulators, this is a direct link between geometry and mechanistic function. Risk: most genes may have near-zero drift and the signal is too subtle for GO enrichment.

Execution plan:
1. For each of 400 sampled genes, compute L2 norm of (layer_k − layer_0) for k=1..11; sum for cumulative drift score.
2. Rank genes by cumulative drift.
3. Top-10% and bottom-10% sets: test GO enrichment using Fisher exact test against background.
4. Null: permute gene-to-drift assignment, re-run GO enrichment.
5. Report: drift distribution, top 20 high-drift genes, GO enrichment p-values.

### Slot 3 — Cheap Broad-Screen Candidate
**H-D + H-C composite: Effective Rank + MADA ID Cross-Check**

Rationale: Both are SVD/PCA-based, run in minutes on existing embeddings, and directly extend the TwoNN findings. Running both together costs <1 compute-hour and either confirms (MADA replicate TwoNN trend) or refutes (MADA flat profile → TwoNN was artifact). Effective rank adds a complementary spectral view. These are fast, require no new data, and tighten the intrinsic dimensionality narrative.

Execution plan:
1. For each layer 0–11, load embedding matrix (all 209 named genes or 400 sample).
2. Compute MADA-style local PCA ID (or scikit-dimension equivalent).
3. Compute effective rank from SVD.
4. Both vs feature-shuffle null (15 replicates each).
5. Output CSV: layer, twonn_id, mada_id, effective_rank, null stats.
