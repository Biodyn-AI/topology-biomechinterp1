# Brainstormer Hypothesis Roadmap — iter_0037

---

## Retire / Deprioritize

| Direction | Reason | Action |
|-----------|--------|--------|
| intrinsic_dimensionality SVD AUROC | ≥3 negative iterations, no variant rescued it | **retire_now** |
| TRRUST/GO proximity (gene-pair annotation) | Confirmed negative, mechanistic mismatch | **retire_now** |
| STRING PPI Euclidean proximity | iter_0031 H02 negative, no biological signal | **retire_now** |
| graph_topology Laplacian spectral gap | Untargeted, no follow-up signal | **retire_now** |
| module_structure GO BP enrichment | Negative, annotation-based approach doesn't capture geometry | **retire_now** |
| T-cell as primary clustering target | 4 iterations negative; keep only as null control | **deprioritize** |
| DC clustering | Negative, heterogeneity of panel likely cause; insufficient vocab | **deprioritize — rescue_once if DC-specific marker set (≥8 genes) becomes available** |

---

## New Hypothesis Portfolio

### H-A: Cross-Model B-Cell Precision@10 in Geneformer
**Hypothesis**: Geneformer (a separate transformer-based foundation model) will reproduce B-cell geometric clustering with z>3, confirming the finding is model-class general, not scGPT-specific.
**Test**: Load Geneformer token embeddings; identify B-cell marker genes in its vocabulary; compute precision@10 at equivalent transformer depths; compare z-scores to scGPT results.
**Expected signal**: z>3 in at least one Geneformer layer, with layer-specific peak matching scGPT (middle layers).
**Null/control**: T-cell markers in Geneformer should not cluster (mirrors scGPT specificity pattern).
**Value**: high | **Cost**: medium

### H-B: GC-TF vs Plasma-TF Subcluster Centroid Distance + Layer Progression
**Hypothesis**: The distance between germinal-center TF cluster centroid (BATF, SPIB, BACH2) and plasma-cell TF cluster centroid (IRF4, XBP1/PRDM1) grows monotonically from L0→L11, reflecting progressive encoding of the B→plasma differentiation axis.
**Test**: Compute per-layer centroid for GC-TFs and plasma-TFs; plot distance L0→L11; test monotonicity (Spearman r vs layer index).
**Expected signal**: Spearman r>0.7, distance at L11 > 2× distance at L0.
**Null/control**: Shuffle gene identities; compare layer-distance profile.
**Value**: high | **Cost**: low

### H-C: Persistent Homology (H0/H1) of B-Cell Cluster
**Hypothesis**: The B-cell marker cluster at L2 has non-trivial persistent H1 (loops) or unusually high H0 persistence relative to random gene sets, indicating a topologically structured rather than amorphous submanifold.
**Test**: Compute Vietoris-Rips PH of 5-gene B-cell panel embeddings vs 200 random 5-gene sets at L2. Compare max H0 persistence and number of H1 generators.
**Expected signal**: B-cell H0 max persistence z>2 vs random; or detectable H1 loop absent in null sets.
**Null/control**: 200 random 5-gene sets; also T-cell marker PH as negative control.
**Value**: high | **Cost**: medium

### H-D: PAX5 / EBF1 / BCL6 Vocabulary Scan + Centroid Distance
**Hypothesis**: Master B-cell TFs PAX5 and EBF1 (not in original panel) are in the full 4803-gene scGPT vocabulary and cluster at >90th pctile proximity to the B-cell centroid, extending the germinal-center cluster.
**Test**: Search full vocabulary (4803 gene indices → gene names) for PAX5, EBF1, BCL6, IKZF1, IKZF3; compute centroid distance for each found gene; rank vs all-195-gene null.
**Expected signal**: PAX5 and EBF1 ≥90th pctile (strong prior from known biology); BCL6 ≥80th pctile.
**Null/control**: All 195 named genes as null; compare to IRF4 (20th pctile) as negative reference.
**Value**: high | **Cost**: low

### H-E: STRING Full PPI Neighborhood Enrichment (Fisher's Exact)
**Hypothesis**: Top-50 B-cell geometric neighbors contain significantly more STRING-connected pairs (to anchor B-cell TFs) than random gene sets, confirmed with Fisher's exact test on full STRING human PPI.
**Test**: Download STRING human v12 protein network (filtered ≥400 score); for top-50 L2 neighbors, count edges to MS4A1/CD19/CD79A/BLK/PRDM1 anchor set; Fisher's exact vs random-50-gene baseline.
**Expected signal**: OR>2.5, Fisher p<0.05.
**Null/control**: 200 random 50-gene sets.
**Value**: medium | **Cost**: medium

### H-F: Layer-Resolved B-Cell Cluster Intrinsic Dimensionality (MLE)
**Hypothesis**: The local intrinsic dimensionality of the B-cell neighborhood (top-50 neighbors at each layer) decreases from L0→L2→L11, reflecting progressive compression of the B-cell program into a lower-dimensional manifold.
**Test**: At each layer L0,L2,L5,L8,L11, compute MLE intrinsic dimension (Levina-Bickel) for the top-50 B-cell neighbors; compare trend to all-gene background.
**Expected signal**: B-cell neighborhood ID decreases from ~8 (L0) to ~3 (L11); background stays flat.
**Null/control**: Random 50-gene sets at each layer; T-cell 50-NN as negative.
**Value**: medium | **Cost**: low

### H-G: NK Cell / Myeloid Clustering with Extended Vocabulary
**Hypothesis**: NK cells (NCAM1, NKG7, GNLY, KLRD1, PRF1) and myeloid cells (LYZ, CD14, CSF1R, FCGR3A, S100A8) are accessible with ≥4 in-vocab genes and will show no clustering (like T-cell/DC), further confirming B-cell specificity is not an artifact of panel size.
**Test**: Same precision@10 protocol; bootstrap null at L2; compare to B-cell z=7.55 baseline.
**Expected signal**: Both NK and myeloid z<2, confirming B-cell uniqueness.
**Null/control**: B-cell panel as positive control.
**Value**: medium | **Cost**: low

### H-H: Geodesic vs Euclidean Distance Divergence in B-Cell Region
**Hypothesis**: Geodesic distance (k-NN graph shortest path) between B-cell markers is substantially smaller than Euclidean distance relative to background gene pairs, indicating the B-cell cluster lies on a curved manifold (not flat).
**Test**: Build k-NN graph (k=10) on all 195 named genes at L2; compute geodesic distance between all B-cell marker pairs vs random same-size gene pairs; ratio geodesic/Euclidean.
**Expected signal**: B-cell geodesic/Euclidean ratio < 0.7 relative to background, Wilcoxon p<0.05.
**Null/control**: Random gene-pair sets matched in Euclidean distance; T-cell markers.
**Value**: medium | **Cost**: low

### H-I: B-Cell Centroid Drift Across Layers (Directional Vector)
**Hypothesis**: The B-cell centroid traces a non-random trajectory in 512-D embedding space from L0→L11; the direction of this trajectory is correlated with the B→plasma differentiation axis (i.e., centroid moves toward plasma cell region).
**Test**: Compute B-cell centroid at L0,L2,L5,L8,L11 (512-D vectors); compute displacement vectors; measure cosine similarity of centroid displacement L0→L11 with vector from B-cell centroid to plasma-cell centroid at L2.
**Expected signal**: Cosine similarity >0.3 (centroid drifts toward plasma region as layers deepen), Wilcoxon vs random directions p<0.05.
**Null/control**: 200 random 5-gene centroids tracked across layers.
**Value**: high | **Cost**: low

### H-J: Local Linear Approximation Quality (PCA Explained Variance) of B-Cell Neighborhood
**Hypothesis**: The top-50 B-cell neighbors at L2 lie in a lower-dimensional linear subspace (first 3 PCs explain >70% variance) compared to random 50-gene sets, indicating the B-cell cluster is approximately planar/linear in the ambient 512-D space.
**Test**: Run PCA on top-50 B-cell L2 neighbors; compute cumulative explained variance at k=1,2,3,5; compare to 200 random 50-gene sets.
**Expected signal**: B-cell top-3-PC explained variance >70% vs background ~50%.
**Null/control**: 200 random 50-gene sets; T-cell 50-NN set.
**Value**: medium | **Cost**: low

### H-K: Cross-Cell-Type Sub-Lineage Separation (Intra-B Differentiation vs Inter-Cell-Type Distance)
**Hypothesis**: The within-B-cell sub-lineage distance (GC-TF centroid ↔ plasma-TF centroid = ~5 units) is comparable to the between-cell-type distance (B-cell ↔ T-cell = 4.82), meaning scGPT encodes differentiation stages as geometrically distinct as separate cell types.
**Test**: Compute GC-TF centroid (BATF, SPIB, BACH2) and plasma-TF centroid (IRF4, PRDM1 if plasma-TF) at L2; compute intra-B differentiation distance; compare to B↔T, B↔DC distances from centroid matrix.
**Expected signal**: Intra-B GC↔plasma distance > 4 (i.e., comparable to inter-cell-type separation).
**Null/control**: Random gene-pair centroid distances.
**Value**: high | **Cost**: low

### H-L: Filtration Stability of B-Cell Cluster (Epsilon-Ball Precision Curve)
**Hypothesis**: The B-cell cluster is stable across filtration radii — precision remains above null at multiple epsilon thresholds (not just k=10), indicating a compact, well-separated cluster rather than a lucky k-NN artifact.
**Test**: Compute epsilon-ball precision for B-cell markers at L2: for each epsilon in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0], count fraction of epsilon-ball neighbors that are B-cell markers; compare to bootstrap null.
**Expected signal**: B-cell precision exceeds null z>2 across ≥4 epsilon values; T-cell shows no such pattern.
**Null/control**: Bootstrap null at each epsilon; T-cell markers.
**Value**: medium | **Cost**: low

---

## Top 3 for Immediate Execution

### Candidate 1 (High-Probability Discovery): H-B — GC-TF vs Plasma-TF Layer Progression
**Why**: The GC/plasma split is already demonstrated (iter_0037 H02). This is a direct extension: quantify it precisely and test whether it grows across layers. Low cost (just centroid distances at 5 layers), high payoff (direct mapping of differentiation axis). Data already in hand.

**Concrete plan**:
- Define GC-TF set: BATF, SPIB, BACH2 (at L2)
- Define plasma-TF set: IRF4, PRDM1 (both in vocab — PRDM1 is in reference panel)
- Compute per-layer centroid for each group, L0→L11
- Compute L2-norm between GC and plasma centroids at each layer
- Spearman test for monotonic increase; plot distance profile
- Also compute combined "differentiation axis" vector and project all B-cell program genes onto it

### Candidate 2 (High-Risk/High-Reward): H-I — B-Cell Centroid Directional Drift Toward Plasma Region
**Why**: If the B-cell centroid drifts in the direction of the plasma cell centroid as layers deepen, this would be a direct mechanistic claim: scGPT encodes the B→plasma differentiation trajectory as a geometric progression. Novel, high-impact, and testable with data already in memory.

**Concrete plan**:
- Compute B-cell centroid (5 reference genes) at L0, L2, L5, L8, L11 (512-D vectors)
- Compute displacement vector: L0→L11 centroid trajectory
- Compute reference direction: B-cell centroid → plasma-cell centroid (PRDM1 as proxy) at L2
- Cosine similarity between displacement trajectory and reference direction
- Bootstrap null: 500 random 5-gene centroid displacements
- Also test: does the GC-TF cluster centroid drift differently from the plasma-TF cluster centroid across layers?

### Candidate 3 (Cheap Broad Screen): H-D + H-G — Vocabulary Scan (PAX5/EBF1) + NK/Myeloid Panel
**Why**: Two sub-tests, both cheap (just centroid distance lookups and kNN precision). PAX5/EBF1 scan directly extends the B-cell geometric finding with master TFs. NK/myeloid panel definitively closes the "is B-cell really unique?" question.

**Concrete plan**:
- Scan full 4803-gene vocabulary for PAX5, EBF1, BCL6, IKZF1, IKZF3
- For found genes: compute distance to B-cell centroid at L2, rank vs 195-gene null
- Build NK panel (NKG7, GNLY, KLRD1, PRF1, NCAM1) and myeloid panel (LYZ, CD14, CSF1R, FCGR3A, S100A8)
- Find in-vocab genes for each panel; compute precision@10 at L2; bootstrap null
- Report z-scores vs B-cell z=7.55 as reference
