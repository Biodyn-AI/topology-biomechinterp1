# Brainstormer Hypothesis Roadmap — iter_0030

---

## Retire / Deprioritize

| Hypothesis family | Status | Reason |
|---|---|---|
| Diverging-14 bifurcation geometry | **retire_now** | OOV artifact. All 14 are zero-embedding tokens; their "divergence" is trivially explained. |
| Diverging-14 cluster stability (bootstrap) | **retire_now** | OOV tokens are 100% stable by definition (they are identically zero). Meaningless. |
| Layer-of-bifurcation anatomy | **retire_now** | Measures inter/intra ratio of zero-vector set vs non-zero set. Not learned geometry. |
| Early-layer separation of 14 diverging genes | **retire_now** | OOV tokens have zero L0 norm — trivially pre-separated at L0 for the wrong reason. |
| PC1 binary biological axis | **retire_now** | Previously retired; PC1 bifurcation was OOV-driven. Confirmed closed. |
| TRRUST regulatory load of diverging genes | **retire_now** | Characterizing OOV genes' biology is a model-vocabulary finding, not geometry. Done. |
| Cross-model diverging gene test (Geneformer) | **retire_now** | Testing whether OOV genes in scGPT diverge in Geneformer is not informative about learned geometry. |
| STRING degree vs trajectory slope (209 genes) | **rescue_once** | Rerun on 195 in-vocab only. If still null, retire permanently. |

---

## New Hypothesis Portfolio

### P01: Clean Geometric Baseline — 195 In-Vocab Gene kNN Graph Topology
**Hypothesis:** The 195 in-vocabulary named genes form a connected kNN graph at deep scGPT layers (L10–L11) with topology that reflects STRING functional network structure.
**Test:** Build kNN graph (k=5,10,15) on 195 genes at each of 12 layers. Compute: (a) number of connected components, (b) spectral gap (λ2), (c) Fiedler vector as biological partition. Test if kNN graph components correspond to GO functional families or STRING network communities. Permutation null: random k-regular graphs on same node set.
**Expected signal:** Fewer components at deep layers (convergence of in-vocab genes); spectral gap increases L0→L11; Fiedler partition recovers a functional family boundary.
**Null:** Same analysis on random Gaussian embeddings.
**Value:** high | **Cost:** low

---

### P02: STRING Score → Embedding Distance (195 In-Vocab Genes, Clean)
**Hypothesis:** For the 195 in-vocabulary genes, pairwise STRING interaction score is negatively correlated with pairwise L2 distance in embedding space at deep layers, establishing structure-function correspondence.
**Test:** Extract all 195×195 pairwise distances at each of 12 layers. Extract STRING combined scores for all 195-gene pairs (0 if not connected). Spearman rho(STRING_score, −distance) at each layer. Kruskal-Wallis test across STRING score tiers (unconnected, 0.4–0.6, 0.6–0.8, ≥0.8). 500-permutation null on STRING scores.
**Expected signal:** rho increases L0→L11; deep-layer rho > 0.15, p < 0.05 by permutation. Tier ordering preserved at L10–L11.
**Null:** rho ≈ 0; tier means not ordered.
**Value:** high | **Cost:** low

---

### P03: Intrinsic Dimension Trajectory (195 In-Vocab Genes)
**Hypothesis:** The 195 in-vocabulary genes show monotonically decreasing intrinsic dimension (TWO-NN or PCA knee) from L0 to L11, consistent with progressive manifold compression.
**Test:** Compute TWO-NN intrinsic dimension estimator on 195-gene embedding matrices at all 12 layers. Also compute PCA knee (number of PCs explaining 90% variance). Report both metrics per layer. Fit linear regression to ID(layer) and report slope. Compare to random Gaussian control (same shape per layer).
**Expected signal:** ID decreases from ~10–15 at L0 to ~2–5 at L11. PCA knee decreases similarly. Slope significantly negative vs random control.
**Null:** Flat ID trajectory (no compression); ID ≈ random Gaussian.
**Value:** high | **Cost:** low

---

### P04: GO Functional Family Clustering in 195-Gene Embedding Space
**Hypothesis:** Genes sharing GO Biological Process annotation clusters occupy tighter embedding neighborhoods at deep layers than gene pairs without shared GO terms.
**Test:** For the 195 in-vocab genes, fetch GO BP annotations (mygene API). Build pairwise GO Jaccard similarity matrix (term overlap). Spearman rho(GO_jaccard, −embedding_L2) at each layer. Mann-Whitney test: within-GO-cluster pairwise distances vs cross-cluster distances at L11.
**Expected signal:** rho(GO_jaccard, −distance) > 0.1 at L10–L11; within-cluster distances < cross-cluster (MW p < 0.05).
**Null:** Permutation of GO cluster labels across 195 genes.
**Value:** high | **Cost:** medium

---

### P05: TF vs Non-TF Embedding Geometry (195 In-Vocab)
**Hypothesis:** Transcription factors among the 195 in-vocabulary genes occupy a geometrically distinct region of embedding space at deep layers (TFs cluster together independently of STRING connectivity).
**Test:** Classify 195 genes as TF/non-TF using TRRUST or AnimalTFDB. At each layer, compute: (a) centroid distance between TF cluster centroid and non-TF centroid; (b) within-TF vs within-nonTF mean pairwise distance; (c) Mann-Whitney AUROC for TF vs non-TF distance-to-TF-centroid. Permutation null: shuffle TF/non-TF labels.
**Expected signal:** AUROC > 0.65 at L10–L11; TF centroid drifts further from non-TF centroid at deep layers.
**Null:** AUROC ≈ 0.5; no centroid separation.
**Value:** medium | **Cost:** low

---

### P06: Layer-Resolved Persistent Homology on 195 In-Vocab Genes
**Hypothesis:** The 195 in-vocabulary genes exhibit growing H0 persistence (merging of components = convergence) from L0 to L11, and the H1 birth-death landscape shifts toward shorter-lived cycles at deep layers (fewer topological holes = more compact manifold).
**Test:** Compute Vietoris-Rips persistence (H0, H1) for 195-gene point clouds at all 12 layers using Ripser. Track: (a) max H0 persistence per layer, (b) number of H1 generators with persistence > 1 SD of H0 persistence. Compare bottleneck distance between consecutive layers (rate of topological change).
**Expected signal:** Max H0 persistence decreases L0→L11 (components merge); H1 count decreases (holes close). Bottleneck distance peaks at L5–L8 (fastest topological change).
**Null:** Same analysis on shuffled embeddings per layer; no systematic trend.
**Value:** high | **Cost:** medium

---

### P07: Per-Layer STRING Community Detection Alignment (195 Genes)
**Hypothesis:** Louvain communities detected on the kNN embedding graph of 195 genes align with STRING network communities (computed independently on STRING data), and this alignment increases L0→L11.
**Test:** (a) Run Louvain community detection on kNN embedding graphs (k=10) at each layer for 195 genes. (b) Run Louvain on STRING subgraph (195 genes, threshold ≥0.4). (c) Compute normalized mutual information (NMI) between embedding communities and STRING communities per layer.
**Expected signal:** NMI increases from L0 to L11; NMI(L11) > NMI(L0) with p < 0.05 by permutation.
**Null:** NMI flat or decreasing; NMI(L11) not greater than permutation null.
**Value:** high | **Cost:** medium

---

### P08: Trajectory Slope vs STRING Degree (195 In-Vocab, Clean Rerun)
**Hypothesis:** After removing the 14 OOV genes, STRING degree predicts trajectory slope (convergence rate) for the 195 in-vocabulary genes.
**Test:** Per-gene trajectory slope = Spearman rho(layer, dist_to_195-centroid) for each of 195 genes. Spearman rho(slope, STRING_degree) across 195 genes. Quartile comparison (Q75 vs Q25 degree). 500-permutation null.
**Expected signal:** rho < −0.15 (higher degree → steeper convergence); p < 0.05. Previous null result (iter_0030 H03) was confounded by OOV artifact.
**Null:** rho ≈ 0, same as iter_0030 H03.
**Value:** high | **Cost:** low

---

### P09: Geodesic Distance vs STRING Score (195 In-Vocab, Graph-Distance Metric)
**Hypothesis:** Geodesic distance on the kNN embedding graph (rather than raw L2 distance) better captures STRING interaction score correlation at deep layers.
**Test:** Build kNN graph (k=10) on 195 genes at L11. Compute all-pairs shortest path (geodesic) distances on the graph. Spearman rho(STRING_score, −geodesic_distance). Compare to rho(STRING_score, −L2_distance) from P02. Mann-Whitney test on whether geodesic outperforms L2 for predicting STRING connectivity.
**Expected signal:** rho(geodesic) > rho(L2) by ≥0.05; geodesic distance better captures STRING signal because it respects manifold structure.
**Null:** rho(geodesic) ≈ rho(L2).
**Value:** medium | **Cost:** low

---

### P10: Cross-Functional Pair Repulsion (195 Genes) — Antipodal Gene Pairs
**Hypothesis:** Gene pairs from maximally different GO functional families (e.g., DNA repair vs immune response) are pushed to antipodal positions in the embedding sphere at deep layers, while within-family pairs cluster.
**Test:** Define GO-distant pairs (Jaccard GO similarity < 0.05) and GO-close pairs (Jaccard > 0.3) among 195 genes. At each layer, compute cosine distance (not L2) for within-close vs between-distant pairs. Mann-Whitney AUROC for cosine separation. Test if cosine of cross-family pairs approaches 1 (orthogonal/antipodal) at L11.
**Expected signal:** AUROC > 0.6 for cross-family vs within-family cosine distances at L11. Mean cosine distance of far-GO pairs increases L0→L11.
**Null:** No cosine separation; AUROC ≈ 0.5.
**Value:** medium | **Cost:** low

---

### P11: Curvature Profile of 195-Gene Manifold per Layer
**Hypothesis:** The 195 in-vocab gene manifold has heterogeneous Ricci curvature at deep layers, with gene hubs (high STRING degree) at high-curvature locations.
**Test:** Compute Ollivier-Ricci curvature on kNN graph (k=10) for 195 genes at L0, L5, L10, L11. Spearman rho(STRING_degree, edge_curvature) at each layer. Report mean and variance of curvature distribution per layer.
**Expected signal:** Mean curvature increases (manifold becomes more positively curved = more sphere-like) L0→L11. High-degree nodes are at higher curvature positions.
**Null:** No correlation between degree and curvature; flat curvature trajectory.
**Value:** medium | **Cost:** medium

---

### P12: Cell Lineage / Cell Ontology Anchor — Do Lineage Markers Cluster?
**Hypothesis:** Genes annotated as cell-lineage markers (T-cell, B-cell, monocyte, NK cell from Cell Ontology or PBMC marker lists) cluster in distinct regions of the 195-gene embedding at deep layers.
**Test:** Extract cell-type marker annotations from a standard PBMC marker database (e.g., HuBMAP cell ontology or CellMarker). Classify the 195 in-vocab genes by cell-type marker status. Compute inter-lineage vs intra-lineage pairwise distances at L11. Permutation null.
**Expected signal:** Intra-lineage gene distances < inter-lineage (MW p < 0.05 at L11). Lineage marker clusters map to distinct embedding neighborhoods.
**Null:** No lineage clustering; distances random relative to marker status.
**Value:** high | **Cost:** low

---

### P13: Spectral Gap Layer Trajectory — Does 195-Gene Graph Become More Connected?
**Hypothesis:** The algebraic connectivity (Fiedler value λ2) of the kNN graph on 195 in-vocab genes increases monotonically from L0 to L11, reflecting progressive geometric consolidation.
**Test:** Build kNN graphs (k=5,10) on 195 genes at each layer. Compute normalized Laplacian spectral gap (λ2) per layer. Spearman rho(layer, λ2). Compare L0 vs L11 spectral gap by paired test. Control: random Gaussian embeddings.
**Expected signal:** rho(layer, λ2) > 0.8, p < 0.01. Spectral gap increases by ≥2× from L0 to L11.
**Null:** Flat or decreasing spectral gap.
**Value:** high | **Cost:** low

---

### P14: Embedding Norm Distribution Across Layers (195 In-Vocab)
**Hypothesis:** The L2 norm of 195 in-vocab gene embeddings shows decreasing variance (normalization effect) across layers, approaching a near-unit-sphere at L11.
**Test:** Compute L2 norm per gene per layer. Report: mean, variance, min, max of norms per layer. Test if variance decreases monotonically. Compute coefficient of variation (CV = std/mean) per layer.
**Expected signal:** CV decreases L0→L11; embeddings approach a sphere at deep layers (consistent with layer norm + attention).
**Null:** Random variation in CV with no layer trend.
**Value:** low | **Cost:** low

---

## Top 3 for Immediate Execution

### Candidate 1 — High-probability discovery: P01 + P13 (Clean Geometric Baseline + Spectral Gap)
**Rationale:** The single most urgent need is establishing the correct geometric baseline on 195 in-vocab genes. P01 (kNN graph topology + connected components + Fiedler partition) and P13 (spectral gap trajectory) are the direct clean replications of prior analyses now known to have been confounded by OOV genes. These provide the foundation for all subsequent hypotheses. They reuse existing embedding data and take <50 lines of code. High probability of producing clean positive results that anchor the paper's new results section.
**Code path:** Load 195-gene embedding matrices (filtered from existing data). Build kNN graphs at each layer. Compute spectral gap + Fiedler vector. Check if Fiedler partition correlates with STRING community or GO family.

### Candidate 2 — High-risk/high-reward: P02 + P07 (STRING Score → Embedding Distance + Community Alignment)
**Rationale:** The core structure-function hypothesis — STRING predicts embedding geometry — has never been cleanly tested on 195 in-vocab genes. P02 tests the direct continuous correlation; P07 tests it via community structure alignment (NMI of embedding vs STRING communities). If either is positive, this is the main publishable claim of the project. Risk: the signal may still be null even with correct inputs, which would be a strong negative result. Reward: a confirmed STRING → embedding correlation would be a landmark finding about scGPT's biological encoding.
**Code path:** Extract 195×195 pairwise L2 distances per layer. Pull STRING scores for 195-gene pairs. Spearman rho per layer. Louvain on both graphs. NMI per layer.

### Candidate 3 — Cheap broad screen: P03 + P08 + P14 (Intrinsic Dimension + Slope vs Degree Rerun + Norm Distribution)
**Rationale:** Three computationally cheap tests that together characterize the basic geometry of the 195-gene set and rescue the one result (STRING degree vs slope) that may have been confounded by OOV. P03 characterizes manifold compression. P08 reruns the null result from iter_0030 H03 on the correct 195-gene subset (potential positive flip). P14 is a quality-check on embedding scale that takes <10 lines. Together <100 lines of code.
**Code path:** P03: TWO-NN ID estimator per layer. P08: per-gene slope on 195-gene centroid → Spearman vs STRING degree. P14: L2 norm stats per layer.
