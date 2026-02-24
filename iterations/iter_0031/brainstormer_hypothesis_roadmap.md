# Brainstormer Hypothesis Roadmap — iter_0031

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| STRING score → L2 embedding distance | rho≈−0.015, AUROC≈0.494 across all layers, both gene sets tested; zero effect | `retire_now` |
| GO BP enrichment in SVD poles | Zero of 591 terms significant (iter_0011); already retired | `retire_now` |
| Repression anti-pole in SV2 | mean z=−1.41, 0/12 layers sig (iter_0013); already retired | `retire_now` |
| Bootstrap CI on act/rep co-pole differential | Chronically underpowered (n=64 repression pairs); CI never excludes zero | `retire_now` |
| TwoNN intrinsic dim with feature-shuffle null | Redundant with participation ratio (PR is stronger, PR rho=−1.000 vs TwoNN z~−2.29); PR supersedes this | `deprioritize` |
| Cross-layer CKA monotonicity | Already fully characterized (iter_0004: all pairwise CKA>0.99); no new signal expected | `retire_now` |

---

## New Hypothesis Portfolio

### HN01 — kNN Community → GO/Cell-type Annotation (High probability; high value)
**Hypothesis:** kNN graph communities at L11 (most modular layer) align with known biological groupings (GO BP families, cell-type marker sets, STRING hub clusters).
**Test:** Build k=10 kNN graph on 195 in-vocab genes at L11. Apply greedy modularity (networkx). Compute NMI between detected communities and: (a) GO BP Jaccard-derived family clusters; (b) cell-type marker genes (T-cell, B-cell, fibroblast, epithelial from iter_0023); (c) STRING high-degree hub genes. Sweep L8–L11 for stability.
**Expected signal:** NMI > 0 for ≥1 biological grouping; at least one community recoverable as a known cell-type or GO family.
**Null/control:** Random community assignments of same size distribution; shuffle gene labels.
**Value:** high | **Cost:** low

### HN02 — GO Jaccard Similarity → Embedding L2 Distance (Medium-high probability)
**Hypothesis:** Gene pairs with higher GO BP Jaccard functional similarity are closer in L2 embedding space at late layers.
**Test:** Compute pairwise GO BP Jaccard for 195 in-vocab genes using annotation cache (iter_0023/0024). Spearman rho(GO_Jaccard, L2_dist) per layer. Separate analysis for high-specificity GO terms only (size 3–30 genes).
**Expected signal:** Spearman rho < −0.05 with p < 0.01 at ≥1 layer; signal strongest at L11 where PR collapse concentrates variance.
**Null/control:** Permuted GO term assignments across genes; compare rho to null distribution.
**Value:** high | **Cost:** low

### HN03 — PC1 at L11 Biological Alignment (Medium probability)
**Hypothesis:** The dominant PC1 of L11 embeddings (explains 26% of variance) aligns with a known biological axis: cell-type identity, transcription factor vs. target status, or STRING hub degree.
**Test:** Project 195 in-vocab genes onto PC1 of L11 SVD. Test separation of: (a) T-cell vs B-cell markers; (b) TF genes (TRRUST regulators) vs target genes; (c) high STRING degree vs low degree. One-sided Wilcoxon rank-sum per group.
**Expected signal:** AUC > 0.65 for ≥1 partition; largest effect likely TF vs target (given SV2 regulatory encoding found earlier).
**Null/control:** PC1 from randomized embedding (feature-shuffle); same partition test on null PC1.
**Value:** high | **Cost:** low

### HN04 — TRRUST Co-pole Replication on 195 In-Vocab Genes
**Hypothesis:** The SV2 activation co-pole signal (12/12 layers, iter_0011/0012) survives OOV correction — i.e., it holds on the 195-gene set, not just the original 209-gene set.
**Test:** Rerun iter_0011 co-pole protocol (SV2 top-52/bottom-52 K poles, N=500 shuffles) on 195 in-vocab genes. TRRUST activation pairs filtered to 195 in-vocab gene pairs. Report n_pairs retained and layer-wise emp_p.
**Expected signal:** Activation co-pole significant at ≥10/12 layers even after OOV removal (signal should be robust given 14 OOV genes represent <7% of the set).
**Null/control:** Gene-label shuffle on 195 in-vocab genes.
**Value:** high | **Cost:** low (validation experiment)

### HN05 — Spectral Gap Layer Profile as Predictor of Biological Annotation Density
**Hypothesis:** Among 195 in-vocab genes, genes in tightly connected kNN components (high local spectral contribution) have denser GO annotation or higher STRING degree than genes in loose fringe components.
**Test:** At L11, compute per-gene kNN graph degree centrality and clustering coefficient. Spearman rho(degree_centrality, GO_annotation_count) and rho(degree_centrality, STRING_degree). Compare hub vs. peripheral gene sets.
**Expected signal:** rho > 0.1 for GO annotation count; well-annotated genes cluster in core components.
**Null/control:** Random degree assignment; compare against null distribution.
**Value:** medium | **Cost:** low

### HN06 — Dual Metric (PR × Spectral Gap) Joint Layer Predictor
**Hypothesis:** The product PR × spectral_gap_ratio is a better predictor of layer depth than either metric alone, forming a joint compression-fragmentation index.
**Test:** Compute PR and spectral_gap per layer (already in H01/H03 artifacts). Fit linear regression of layer depth on PR, spectral_gap, and PR×spectral_gap. Compare R² of single-feature vs joint model. Also test log(PR) + log(spectral_gap) as additive predictor.
**Expected signal:** R² > 0.99 for joint model; single metrics already achieve rho=−1.000 so joint model confirms independence/additivity.
**Null/control:** Random layer permutation; R² of null model.
**Value:** medium | **Cost:** low (uses existing artifact data)

### HN07 — Persistent Homology H0 (Connected Component) Count per Layer
**Hypothesis:** As layers deepen and the kNN graph becomes more fragmented, the number of connected components (H0 Betti number) increases, providing a topological count of distinct gene clusters.
**Test:** Build kNN graph (k=5, k=10) per layer on 195 genes. Count connected components. Also compute H0 persistence (0-dimensional persistent homology via Ripser on distance matrix) to track birth/death of components across filtration.
**Expected signal:** Component count increases from L0 to L11; H0 persistence diagram shifts to larger death values at late layers.
**Null/control:** Gaussian random embeddings same shape; compare component counts.
**Value:** medium | **Cost:** low

### HN08 — Cross-Model Replication: PR and Spectral Gap in Geneformer
**Hypothesis:** The dual compression+fragmentation (PR collapse + spectral gap decrease) is not scGPT-specific but a general property of single-cell transformer gene encoders.
**Test:** Load Geneformer residual-stream embeddings for the same or overlapping named gene set. Compute PR and kNN spectral gap per Geneformer layer. Spearman rho(layer, PR) and rho(layer, spectral_gap).
**Expected signal:** Both rho < −0.7 in Geneformer, mirroring scGPT's rho = −1.000.
**Null/control:** Random embedding same size; compare PR profile.
**Value:** high | **Cost:** medium (requires loading Geneformer embeddings)

### HN09 — Manifold Curvature Proxy: Deviation from Local Linearity Across Layers
**Hypothesis:** Residual-stream manifolds become less locally linear as layers deepen (consistent with increasing modularity and dimension collapse), detectable via local PCA reconstruction error.
**Test:** For each gene at each layer, fit a local PCA on its k=15 nearest neighbors. Record reconstruction error (fraction of variance not explained by top-2 local PCs). Mean across genes per layer. Spearman rho(layer, mean_recon_error).
**Expected signal:** Reconstruction error increases across layers (manifold becomes less flat); peak curvature at L11.
**Null/control:** Feature-shuffled embeddings have constant reconstruction error across layers.
**Value:** medium | **Cost:** medium

### HN10 — Filtration-Based Topology: Vietoris-Rips vs Alpha Complex Comparison
**Hypothesis:** The topological signature (H1 persistent homology) of the gene embedding manifold depends qualitatively on filtration type, revealing whether cycles emerge at characteristic length scales or diffusely.
**Test:** Compute H1 persistence via both Vietoris-Rips and alpha complex (Gudhi) on 195 in-vocab genes at L0, L6, L11. Compare: (a) total H1 persistence, (b) maximum bar length, (c) number of H1 generators. Report filtration-type differential.
**Expected signal:** Alpha complex reveals fewer but longer-lived H1 bars (more topologically significant loops) than Vietoris-Rips; layer-depth trend in both.
**Null/control:** Gaussian null embeddings same size; compare H1 bar counts.
**Value:** medium | **Cost:** medium

### HN11 — Embedding Trajectory Curvature: Do Genes Move in Curved Arcs Across Layers?
**Hypothesis:** Individual genes trace curved (non-linear) paths through 12-layer embedding space, and the curvature of these paths correlates with the gene's biological role (e.g., TFs vs housekeeping genes have different trajectory shapes).
**Test:** For each of 195 in-vocab genes, extract its 12-dimensional trajectory (one point per layer). Fit total path length and path curvature (turn angle sum). Spearman rho(curvature, STRING_degree), rho(curvature, is_TF). Compare TF vs non-TF path curvature distributions.
**Expected signal:** TFs have higher trajectory curvature (more layer-dependent geometry); housekeeping genes have flatter paths.
**Null/control:** Random shuffled gene assignments for TF/non-TF; curvature distribution under null.
**Value:** medium | **Cost:** medium

### HN12 — Bottleneck Distance Between Layer Persistence Diagrams (PH Stability)
**Hypothesis:** The bottleneck distance between consecutive-layer persistence diagrams (H1) decreases at early layers and increases at later layers, tracking where the most geometric transformation occurs.
**Test:** Compute H1 persistence diagrams for all 12 layers (Ripser). Compute bottleneck distance (persim or Gudhi) for all 11 consecutive pairs. Plot and report max-change layer interval.
**Expected signal:** Bottleneck distance peaks in L3–L6 window (mid-network transformation zone); stabilizes at early and late layers.
**Null/control:** Consecutive Gaussian noise diagrams (expected small bottleneck distance by stability theorem).
**Value:** medium | **Cost:** medium

### HN13 — Dimensionality by Biological Subgroup: Do TFs Occupy a Lower-Dimensional Submanifold?
**Hypothesis:** Transcription factor genes (TRRUST regulators, n~30) occupy a lower-dimensional subspace within the full gene embedding than non-TF genes, especially at late layers where PR collapses.
**Test:** At each layer, compute PR separately for: (a) TRRUST TF genes (in-vocab subset), (b) non-TF genes, (c) STRING hub genes, (d) cell-type markers. Compare PR trajectories across depth.
**Expected signal:** TF PR trajectory lower and/or faster-decreasing than non-TF; consistent with regulatory genes occupying a specialized, compact subspace.
**Null/control:** Random same-size subsets drawn from 195-gene pool; bootstrap CI on subgroup PR difference.
**Value:** medium | **Cost:** low

### HN14 — Embedding Space Density: Does Neighbor Distance Distribution Change Layer-Wise?
**Hypothesis:** Mean kNN neighbor distance decreases across layers (consistent with dimension collapse), and the distribution narrows (genes become more equidistant from their neighbors at late layers).
**Test:** For 195 in-vocab genes at each layer, compute mean and CV of k=10 nearest-neighbor L2 distances. Spearman rho(layer, mean_knn_dist) and rho(layer, CV_knn_dist).
**Expected signal:** Mean kNN distance decreases across layers; CV also decreases (more uniform neighbor distances — consistent with modular clustering where intra-cluster distances equalize).
**Null/control:** Gaussian null same shape; kNN distance profile under null.
**Value:** low-medium | **Cost:** low

---

## Top 3 for Immediate Execution

### Candidate A — High-Probability Discovery
**HN01: kNN Community → GO/Cell-type Annotation**
- Directly follows from the two positive iter_0031 findings (spectral gap, PR collapse)
- Community structure at L11 is physically real (spectral gap = 0.023, strongly modular)
- Biological annotation of these communities is the highest-value next experiment
- Low cost: uses existing embeddings + GO cache from iter_0024
- Success would establish the first biologically-labeled modular decomposition of the gene embedding at late layers

### Candidate B — High-Risk / High-Reward
**HN08: Cross-Model Replication in Geneformer**
- If PR collapse + spectral gap decrease replicate in Geneformer, this becomes a universal law for single-cell transformers
- Cross-model replication is the single largest validity boost available at this stage
- Risk: Geneformer residual stream may require non-trivial loading work; different gene vocabulary may reduce overlap
- Reward: Cross-model law transforms the finding from scGPT-specific curiosity to general transformer geometry

### Candidate C — Cheap Broad Screen
**HN04 + HN03 combined: TRRUST Co-pole OOV Replication + PC1 Alignment at L11**
- HN04 validates that the most important prior finding (signed regulatory geometry, SV2 co-pole) survives OOV correction — cheap validation that protects all iter_0011–0013 claims
- HN03 tests PC1 biological alignment — direct follow-on from H03, uses existing L11 artifact, one pass of Wilcoxon tests
- Both together form a single compact script (~50 lines each); negligible compute cost
- Either positive result provides a direct narrative link between the geometric findings (PR, spectral gap) and biological function
