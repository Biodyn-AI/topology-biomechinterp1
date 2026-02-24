# Brainstormer Hypothesis Roadmap — iter_0028

---

## Retire / Deprioritize

| Hypothesis family | Status | Reason |
|------------------|--------|--------|
| PC1 binary biological axis | **retire_now** | Tested 2× (iter_0027, iter_0028), all AUROCs near 0.5, 1D collapse defeats binary separation |
| TRRUST activation/repression polarity | **retire_now** | Null confound confirmed iter_0027 H02; no rescue path without fundamentally different null |
| Dorothea confidence-tier pairwise distance | **retire_now** | Inconclusive across 2 iterations; low distinguishing power |
| PC1 T-cell vs APC cell-type axis | **retire_now** (H02 iter_0028) | Both families collapse to identical PC1 value; structural impossibility |

---

## New Hypothesis Portfolio

### H-A: Fiedler Vector Biological Partition
**Hypothesis:** The Fiedler eigenvector (corresponding to lambda_2 of the normalized graph Laplacian) partitions genes into two biologically coherent communities that deepen in separation across scGPT layers.
**Test:** At each layer, extract Fiedler vector for the k=10 kNN gene graph. Rank 209 genes along Fiedler dimension. Test top-30 vs bottom-30 for TRRUST TF enrichment, immune gene family membership (HLA-I, AP1, etc.), and STRING module annotation. Report AUROC per family.
**Expected signal:** Top/bottom partition separates known functional modules (e.g., transcriptional regulators vs structural/surface genes); effect strengthens with depth as spectral gap narrows.
**Null:** Permuted gene labels assigned to same Fiedler vector rank positions.
**Value:** high | **Cost:** low

---

### H-B: Hub Centrality via TRRUST Regulatory Load
**Hypothesis:** Distance-to-centroid at L10 correlates with TRRUST transcriptional regulatory target count (independently of STRING PPI degree).
**Test:** For 69 TFs in the 209-gene vocab with TRRUST records, compute Spearman rho(TRRUST_n_targets, −dist_to_centroid_L10). Compare with STRING degree correlation (already rho=+0.225). Test significance with 500 permutations of TRRUST target counts.
**Expected signal:** rho > 0.15 with permutation null_p < 0.05; if both STRING and TRRUST measures independently predict centrality, the mechanism is model-agnostic network-hubness → geometric centrality.
**Null:** Permuted TRRUST target counts.
**Value:** high | **Cost:** low

---

### H-C: Spectral Gap k-Robustness
**Hypothesis:** The monotonic spectral gap decrease (rho=−0.993) is invariant to the choice of k in the kNN graph (k=5,10,15,20,30).
**Test:** Rerun H03 at five k values; compute Spearman rho(layer, spectral_gap) per k. Plot gap profiles and assess whether monotonicity and slope are k-stable.
**Expected signal:** rho < −0.9 at all k values; gap magnitude scales with k but direction is invariant.
**Null:** Same k sweep on random Gaussian embeddings of same shape.
**Value:** medium | **Cost:** low (direct code reuse of H03)

---

### H-D: Connected Component Gene Identity
**Hypothesis:** The two kNN graph connected components at L11 (identified structurally in H03) correspond to biologically distinct gene sets.
**Test:** Extract component assignments for 209 genes at L11 using the k=10 kNN graph. Run Fisher exact test for enrichment of TRRUST TFs, HLA family, AP1 family, STRING module, and cell-surface vs nuclear gene annotation in each component.
**Expected signal:** One component enriched for transcriptional regulators/nuclear genes; other for surface receptors/HLA/structural; consistent with immune gene family AUROC results.
**Null:** 1000 random 2-partitions of same component sizes.
**Value:** high | **Cost:** low

---

### H-E: Hub Gene Trajectory Clustering
**Hypothesis:** Individual gene embedding trajectories across layers (L0→L11) cluster by biological function — hub genes converge toward centroid while peripheral genes diverge or remain stable.
**Test:** Compute per-gene trajectory as a 12-vector of dist_to_centroid across layers. Apply k-means (k=3–5) to 209 gene trajectories. Test cluster-function association: AUROC per cluster for TF-hood, STRING degree quartile, HLA family, AP1 family.
**Expected signal:** At least one cluster shows strong convergent-to-centroid trajectory and is enriched for high STRING degree genes.
**Null:** Same k-means on shuffled trajectories (permute layer order per gene independently).
**Value:** high | **Cost:** medium

---

### H-F: Cross-Model Hub Centrality (Geneformer)
**Hypothesis:** Hub genes (high STRING degree) are also more central in Geneformer deep-layer embeddings, establishing the STRING-degree → centrality relationship as model-agnostic.
**Test:** Load Geneformer embeddings for matched genes. At the deepest available layer, compute Spearman rho(STRING_degree, −dist_to_centroid) with 500-permutation null.
**Expected signal:** rho > 0.15 and null_p < 0.05 in Geneformer, matching the scGPT L10 finding.
**Null:** 500 permutations of STRING degree labels.
**Value:** high | **Cost:** medium (requires Geneformer embedding access)

---

### H-G: Persistence Diagram Layer Distance (Bottleneck)
**Hypothesis:** The bottleneck distance between consecutive-layer persistence diagrams (H0+H1) peaks at the layer where the PR collapse accelerates most sharply.
**Test:** Compute H0+H1 persistence diagrams (Ripser) for the 209-gene subspace at all 12 layers. Compute pairwise bottleneck distances between consecutive layer diagrams. Plot bottleneck distance vs layer. Compare with PR profile.
**Expected signal:** Bottleneck distance peaks at L7–L9 (where PR drop from ~10 to ~2 occurs most sharply), confirming topological and dimensional transitions co-occur.
**Null:** Same calculation on random Gaussian embeddings (same shape per layer).
**Value:** high | **Cost:** medium

---

### H-H: Collapse Anisotropy Profile
**Hypothesis:** The near-1D collapse at L11 (PC1_var=76.7%) is anisotropic — the model compresses into a 1D line rather than a low-D hyperplane — and the collapse trajectory (PC1 share vs layer) follows a sigmoidal rather than linear profile.
**Test:** Compute cumulative PCA variance for 209 genes at all 12 layers (how many PCs needed for 80%, 90%, 95% variance). Fit sigmoidal and linear models to PC1_var(layer). Report AIC difference.
**Expected signal:** Sigmoidal fits better; inflection point at L7–L9 (where PR also collapses fastest). This characterizes the "compression event" geometry.
**Null:** Linear fit; same analysis on random Gaussian control.
**Value:** medium | **Cost:** low

---

### H-I: Intra-family Spectral Gap
**Hypothesis:** Gene families with high clustering AUROC (HLA-I, AP1, RUNX) show smaller internal kNN spectral gaps at deep layers than families with low AUROC (BCL2fam, TNFSF), reflecting tighter community structure.
**Test:** Build per-family kNN graphs (k=3 or k=5 due to small family size) for each of 9 immune families at each layer. Compute Fiedler value per family-layer. Compare Fiedler trajectory of high-AUROC vs low-AUROC families.
**Expected signal:** HLA-I and AP1 have smaller (more modular) within-family Fiedler values at L8–L11; BCL2fam and TNFSF have larger values (looser internal structure).
**Null:** Same per-family graphs with permuted inter-gene connections.
**Value:** medium | **Cost:** medium

---

### H-J: STRING Subgraph Density vs Geometric Cohesion
**Hypothesis:** Gene pairs with STRING interaction score ≥ 0.7 (high-confidence) are closer in embedding space than pairs with score 0.4–0.7 (medium), which are closer than unconnected pairs.
**Test:** Bin 209-gene STRING pairs by score tier (0.4–0.6, 0.6–0.8, ≥0.8, none). Compute mean pairwise L2 distance per tier at each layer. Kruskal-Wallis test across tiers; Spearman rho(score, −distance) for all connected pairs.
**Expected signal:** rho > 0.1 at deep layers; tier ordering preserved (higher score → shorter distance). This establishes a continuous STRING score → embedding distance relationship beyond binary connected/disconnected.
**Null:** 500 permutations of STRING scores across gene pairs.
**Value:** medium | **Cost:** low

---

### H-K: Gene Family Hub Centrality Stratification
**Hypothesis:** Within the top-centrality hub genes (top-30 most central at L10), the enrichment for specific immune families (HLA-I, AP1, RUNX) mirrors the clustering AUROC ranking from iter_0026.
**Test:** Rank 209 genes by centrality at L10 (dist to centroid). Take top-30 and bottom-30. Run Fisher exact test for each of the 9 immune families. Compare enrichment ranks with iter_0026 per-family AUROC ranks.
**Expected signal:** Top-centrality genes enriched for HLA-II, AP1, RUNX (the high-AUROC families); bottom-centrality genes enriched for TNFSF, BCL2fam (low-AUROC families). Spearman rho(AUROC_rank, enrichment_odds_ratio) > 0.6.
**Null:** Fisher exact p-values per family; also 500 random top-30 subsets.
**Value:** high | **Cost:** low

---

### H-L: Geometric Centrality Predicts Regulatory Breadth Beyond Degree
**Hypothesis:** After controlling for STRING degree, genes with high TRRUST out-degree (number of downstream regulatory targets) are closer to the centroid at L10, establishing a second independent network predictor of geometric centrality.
**Test:** Partial correlation of dist_to_centroid_L10 with TRRUST_out_degree, controlling for STRING_degree. Bootstrap CI (1000 resamples) on partial rho. Compare model fits (STRING-only vs STRING+TRRUST) by R².
**Expected signal:** Partial rho(TRRUST_out, −dist) > 0.1 after controlling STRING degree; R² increases with TRRUST added.
**Null:** 500 permutations of TRRUST out-degree (with STRING degree held fixed).
**Value:** high | **Cost:** low

---

## Top 3 for Immediate Execution

### Candidate 1 — High-probability discovery: H-D (Connected Component Gene Identity)
**Rationale:** The spectral gap data already flags 2 connected components at every layer. Component assignment for 209 genes is a direct O(n) computation from the existing kNN adjacency matrix. Biological enrichment of a data-driven geometric partition is maximally interpretable and likely to yield a positive result given the strong family-clustering signal from iter_0026. Low cost, high expected information.
**Code path:** Reuse H03 kNN graph construction; extract `nx.connected_components(G)` per layer; Fisher exact per immune family; 1000-permutation null on component sizes.

### Candidate 2 — High-risk/high-reward: H-E (Hub Gene Trajectory Clustering)
**Rationale:** Individual gene trajectories across 12 layers encode the full history of how the model transforms gene representations. Clustering trajectories by shape (convergent vs stable vs divergent) provides a novel organizing principle not yet tested. If trajectory clusters map to known biology, this establishes a "how the model moves genes" characterization that is unique in the field. Medium cost, potentially very high interpretability value.
**Code path:** Compute 209×12 distance-to-centroid matrix; apply k-means (k=3,4,5); AUROC per cluster for TF, STRING quartile, family membership; test vs shuffled trajectories.

### Candidate 3 — Cheap broad screen: H-C + H-J combined (k-Robustness + STRING score gradient)
**Rationale:** Two cheap tests that directly address the two open validation questions from iter_0028: (a) is H03 spectral gap k-invariant? and (b) does STRING interaction score continuously predict embedding distance beyond binary connectivity? Both are ≤50 lines of new code, reuse existing data, and have clear pass/fail criteria. Package as a single script for efficiency.
**Code path:** H-C: loop H03 over k=[5,10,15,20,30]; H-J: extract STRING score tiers from existing edge list, compute pairwise distance per tier, Kruskal-Wallis.
