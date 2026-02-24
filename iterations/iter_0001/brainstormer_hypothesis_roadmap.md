# Brainstormer Hypothesis Roadmap — iter_0001
Date: 2026-02-22

---

## Retire / Deprioritize

| Direction | Status | Reason |
|---|---|---|
| Rewiring-null survival (degree-preserving geodesic) | `retire_now` | 0/24 significant across 3 consecutive iterations (iter_0006–0008); quantile-constrained variant also failed; strongly negative deltas show no sign of rescue. |
| Distance-permutation null | `retire_now` | Over-adversarial (mean delta −850); not biologically interpretable; replaced by rewiring null which also failed. |
| Bridge-conditioned topology (H11) | `rescue_once_with_major_change` | Requires bridge-identifiable kNN schedule (e.g., k sweep 5→20 with bridge logging per run) to avoid split confound; low priority vs new families. |
| Cross-model feature-effect vector cosine (H02) | `rescue_once_with_major_change` | Weak significance despite high raw similarity. Rescue only with CKA or Procrustes alignment on full residual tensors — if those artifacts are not available, deprioritize. |

---

## New Hypothesis Portfolio

### N01 — kNN Graph Topology Surrogates (clustering coefficient + transitivity)
- **Hypothesis:** scGPT layer gene-embedding kNN graphs show significantly higher clustering coefficient and transitivity than feature-shuffle nulls, indicating non-random local neighborhood structure.
- **Test:** Execute `run_graph_topology_screen.py` (already written): PCA(20) → kNN(k=10) → CC + transitivity → 15 feature-shuffle null replicates; 3 seeds × 12 layers; lung domain.
- **Expected signal if true:** Mean CC/transitivity deltas positive and significant (empirical p < 0.05) in ≥ 8/12 layers.
- **Null/control:** Feature-shuffle (column-wise permutation) replicates.
- **Value:** `high` | **Cost:** `low` (code ready, data accessible, ~10 min runtime)

### N02 — TRRUST Co-regulatory Distance vs Embedding Geodesic Correlation
- **Hypothesis:** Pairs of genes with fewer TRRUST intermediary regulatory steps (closer in the regulatory graph) have shorter geodesic distances in scGPT residual embedding space than gene pairs with many steps.
- **Test:** Build TRRUST gene-gene regulatory adjacency; compute BFS shortest-path distance; compute embedding geodesic distance on kNN graph per layer; Spearman correlation vs permutation null (500 gene-pair permutations). Test top layer + 2 random layers, lung domain, 3 seeds.
- **Expected signal if true:** Significant negative Spearman (closer regulatory = shorter geodesic) in top layers; permutation p < 0.05.
- **Null/control:** Gene-label permutation (random relabeling of regulatory distance matrix).
- **Value:** `high` | **Cost:** `medium` (requires TRRUST parsing + graph ops)

### N03 — GO Term Embedding Cohesion Screen
- **Hypothesis:** Genes sharing a GO slim term cluster more tightly in embedding space (higher mean pairwise cosine similarity) than random same-size gene sets.
- **Test:** For each GO slim term with 10–200 member genes present in the embedding, compute mean pairwise cosine similarity of embeddings at top 3 layers; compare against 500 random same-size gene set draws (permutation null); report enriched terms at FDR < 0.1.
- **Expected signal if true:** ≥ 20% of tested GO terms significantly enriched; biological process terms enriched in top/biologically-active layers.
- **Null/control:** Random gene set size-matched permutation.
- **Value:** `high` | **Cost:** `medium` (needs GO annotation parsing; 200–400 term tests)

### N04 — CKA Cross-Model Layer Alignment
- **Hypothesis:** Corresponding layers of scGPT and Geneformer (matched by relative depth) show higher Centered Kernel Alignment than non-corresponding layer pairs.
- **Test:** Compute linear CKA between scGPT layer L_i and Geneformer layer L_j for all (i,j) pairs on matched genes; test if diagonal (same relative depth) CKA > off-diagonal under permutation null; 3 seeds.
- **Expected signal if true:** Diagonal CKA significantly higher than row/column mean off-diagonal; permutation p < 0.05 for ≥ 3 corresponding layer pairs.
- **Null/control:** Gene-index permutation within one model's embedding matrix.
- **Value:** `high` | **Cost:** `medium` (requires Geneformer residual embeddings — check artifact availability first)

### N05 — Persistence Diagram Seed Stability (Wasserstein Consistency)
- **Hypothesis:** Same-layer persistence diagrams from different seeds are more similar (lower Wasserstein distance) than cross-layer pairs, indicating topology is layer-specific rather than noise-driven.
- **Test:** Compute Wasserstein-1 distance between H1 persistence diagrams for (seed_a, seed_b, same layer) vs (seed_a, layer_i, seed_a, layer_j) pairs; test with permutation of seed-layer labels.
- **Expected signal if true:** Same-layer cross-seed Wasserstein distance significantly lower than cross-layer same-seed distance (permutation p < 0.05).
- **Null/control:** Random permutation of layer labels before distance computation.
- **Value:** `high` | **Cost:** `medium` (requires persim/gudhi Wasserstein; manageable with 12 layers × 3 seeds)

### N06 — H0 (Connected Components) Persistence by Layer
- **Hypothesis:** H0 lifetime sums (connected components) in scGPT layer embeddings differ significantly from feature-shuffle null, and H0 structure has a distinct layer depth profile compared with H1.
- **Test:** Run ripser with `maxdim=0`; compute H0 lifetime sum per seed-layer; compare against 20 feature-shuffle replicates using same protocol as H01. Add this as a parallel track to existing H1 runs.
- **Expected signal if true:** H0 positive in early/intermediate layers (clustering structure); possibly negative near top layers if embedding becomes more linear.
- **Null/control:** Feature-shuffle (same as H01 null).
- **Value:** `medium` | **Cost:** `low` (trivial extension of H01 ripser calls)

### N07 — Topological Lifetime Entropy by Layer
- **Hypothesis:** Shannon entropy of the H1 persistence lifetime distribution decreases monotonically with layer depth (fewer, longer cycles dominate in top layers), indicating progressive topological simplification.
- **Test:** Compute entropy of normalized lifetime histogram per layer per seed; fit layer-depth linear regression; test slope significance against seed-permuted layer-depth assignments.
- **Expected signal if true:** Significant negative slope (p < 0.05 permutation) in ≥ 2/3 domains.
- **Null/control:** Permutation of layer index assignments per seed.
- **Value:** `medium` | **Cost:** `low` (post-hoc on existing H01 persistence diagrams)

### N08 — Filtration Radius Sensitivity of H1 Signal
- **Hypothesis:** The H1 persistence vs feature-shuffle effect size is stable across a range of Vietoris-Rips filtration cutoffs (not an artifact of a single scale), indicating genuine multi-scale topological structure.
- **Test:** Run H1 computation with 5 filtration cutoff percentiles (20th–80th of pairwise distance); compute effect size (delta/null_std) at each cutoff; test if plateau region exists (≥ 3 consecutive cutoffs with delta > 2σ).
- **Expected signal if true:** Flat/plateau effect size curve over ≥ 3 cutoff levels in top layers.
- **Null/control:** Same feature-shuffle null at each cutoff level.
- **Value:** `medium` | **Cost:** `low` (extension of existing ripser calls with cutoff parameter)

### N09 — Ollivier-Ricci Curvature Distribution by Layer
- **Hypothesis:** kNN graph edges in top scGPT layers show more negative mean Ollivier-Ricci curvature (hyperbolic-like, tree-structured) than feature-shuffle nulls, indicating hierarchical manifold structure.
- **Test:** Compute OR curvature for each edge in kNN graph using earth-mover distance (Wasserstein-1) on 1-hop neighbor distributions; average over edges per layer; compare against 15 feature-shuffle null replicates; 3 seeds × 4 layers (top 2, bottom 2) for cost control.
- **Expected signal if true:** Mean OR curvature significantly more negative in observed vs null (permutation p < 0.05) in ≥ 2/4 tested layers.
- **Null/control:** Feature-shuffle kNN graph null.
- **Value:** `medium` | **Cost:** `medium` (OR curvature requires per-edge Wasserstein, expensive for large k; limit n_genes and k)

### N10 — STRING PPI Distance vs Embedding Geodesic Correlation
- **Hypothesis:** Gene pairs with higher STRING combined interaction score (closer in PPI space) have shorter embedding-space geodesic distances in scGPT top layers.
- **Test:** Subset to genes present in STRING and embedding; compute STRING-distance (1/score normalized) for all pairs; compute kNN geodesic distance per layer; Spearman correlation; permutation null (500 gene-label permutations); 3 layers × 3 seeds × lung domain.
- **Expected signal if true:** Significant negative Spearman in top layers (high STRING score = low geodesic distance), permutation p < 0.05.
- **Null/control:** Gene-label permutation on STRING distance matrix.
- **Value:** `high` | **Cost:** `medium` (requires STRING data parsing; correlation is fast)

### N11 — Gene Neighborhood Consistency Across Seeds (Topological Stability)
- **Hypothesis:** Genes that contribute to high-H1-persistence cycles in one seed also appear as neighbors in persistent cycles in other seeds (topological stability beyond seed-specific noise).
- **Test:** For each layer, identify genes in top-20% H1 lifetime cycles per seed; compute Jaccard overlap of these gene sets across seed pairs; compare against random same-size gene sets.
- **Expected signal if true:** Significantly higher Jaccard overlap than random for ≥ 8/12 layers; overlap concentrates in top layers.
- **Null/control:** Random gene set of same size drawn independently per seed.
- **Value:** `high` | **Cost:** `low` (post-hoc on existing persistence diagrams — requires cycle generator extraction from ripser)

### N12 — Louvain Community Modularity vs Feature-Shuffle Null
- **Hypothesis:** kNN graph community structure (Louvain modularity) in scGPT layer embeddings is significantly higher than feature-shuffle null, indicating clustered rather than random neighborhood organization.
- **Test:** Build kNN graph; run Louvain community detection; compute modularity Q; compare against 20 feature-shuffle null replicates; 3 seeds × 12 layers × lung domain.
- **Expected signal if true:** Mean modularity delta positive and significant (p < 0.05) in ≥ 8/12 layers.
- **Null/control:** Feature-shuffle column-wise permutation.
- **Value:** `medium` | **Cost:** `low` (scipy/igraph Louvain is fast; same data pipeline as N01)

### N13 — Intrinsic Dimension TWO-NN Estimator by Layer with GO Group Stratification
- **Hypothesis:** Genes from the same GO biological process slim term occupy lower intrinsic-dimensionality submanifolds than the full gene set in top scGPT layers.
- **Test:** Compute TWO-NN intrinsic dimensionality estimator for (a) full gene set and (b) GO-term stratified subsets (≥ 30 genes); test if GO-group ID < full-set ID using permutation of GO labels.
- **Expected signal if true:** GO-group ID significantly lower than matched-size random gene set in ≥ 20% of tested GO terms; effect stronger in top layers.
- **Null/control:** Random gene subset of same size (500 permutations per term).
- **Value:** `medium` | **Cost:** `medium` (TWO-NN fast; GO parsing adds cost)

### N14 — Persistence Image Vector Layer Classifier
- **Hypothesis:** Persistence images (flattened H1 persistence diagrams) for each layer form linearly separable clusters in a low-dimensional projection, enabling accurate layer-depth prediction without biological labels.
- **Test:** Compute persistence images for all seed-layer combinations; PCA to 10D; train linear SVM (leave-one-seed-out CV); report layer classification accuracy vs chance.
- **Expected signal if true:** Layer prediction accuracy > chance (1/12) by ≥ 3×; top vs bottom layer separation accuracy > 80%.
- **Null/control:** Shuffled layer labels (permutation test).
- **Value:** `medium` | **Cost:** `medium` (requires persim persistence images; manageable scope)

### N15 — Betti Number Layer Profile Comparison Across Domains
- **Hypothesis:** The layer profile of mean H1 Betti number (count of significant cycles at fixed filtration) is domain-specific and reflects biological context differences (lung vs immune vs external-lung).
- **Test:** Compute mean H1 Betti number per layer per domain at a single calibrated filtration cutoff (50th percentile pairwise distance); compare domain profiles via Kolmogorov-Smirnov test; control: domain-label permutation.
- **Expected signal if true:** KS test significant (p < 0.05) for ≥ 1 domain pair; immune/lung profiles differ in shape or peak layer.
- **Null/control:** Permutation of domain labels across seed-layer rows.
- **Value:** `medium` | **Cost:** `low` (post-hoc on existing per-layer H1 diagrams)

---

## Top 3 for Immediate Execution

### Candidate A — High-Probability Discovery: N01 (kNN Graph Topology Surrogates)
- **Rationale:** Code is already written and validated in `run_graph_topology_screen.py`. This is a materially new family (graph topology vs persistent homology). Feature-shuffle null has been consistently positive for H1, so CC/transitivity are likely to show similar structure. Gate recovery and new-family contribution in a single run.
- **Execution:** `python iterations/iter_0001/run_graph_topology_screen.py`. Write output CSV/JSON. Populate `executor_hypothesis_screen.json`.
- **Minimum pass bar:** ≥ 6/12 layers with CC empirical p < 0.05 (concordant with H1 positivity rate).

### Candidate B — High-Risk/High-Reward: N02 (TRRUST Co-regulatory Distance vs Embedding Geodesic)
- **Rationale:** Direct biological validation of the positive H1 signal. If embedding geodesic distances correlate with regulatory graph distances in TRRUST, this elevates the entire topology branch from statistically significant to biologically anchored. High upside; moderate implementation cost. Negative result is also informative (positive topology signal without biological specificity = artifact risk).
- **Execution:** Parse TRRUST human TF-gene table; build TF-gene bipartite graph; BFS shortest-path between gene pairs; correlate with embedding kNN geodesic per layer; 500 gene-pair permutations.
- **Minimum pass bar:** Spearman correlation significant at permutation p < 0.05 in ≥ 1 layer-domain combination.

### Candidate C — Cheap Broad Screen: N06 + N07 + N15 combo (H0/Entropy/Betti post-hoc)
- **Rationale:** All three can be computed from existing H01 ripser output files (iter_0003 and iter_0004 CSVs contain lifetime sums; raw diagrams can be re-extracted with a short additional ripser pass). Combined, they cost a single script (~60 lines) and produce 3 new hypothesis results. Low cost, broad coverage, no new data needed.
- **Execution:** Write and run a post-hoc analysis script over existing per-seed-layer persistence diagram data; compute H0 lifetimes, entropy of lifetime distributions, and Betti counts at fixed filtration; output summary CSV per hypothesis.
- **Minimum pass bar:** At least 1 of the 3 sub-hypotheses produces a directional pattern with p < 0.1 across layers.

---

*End of roadmap.*
