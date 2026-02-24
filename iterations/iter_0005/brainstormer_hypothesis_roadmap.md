# Brainstormer Hypothesis Roadmap — iter_0005

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| TRRUST co-target clustering (209-gene) | `retire_now` | Inconclusive x2, gene pool too small; reschedule only after full vocab recovery |
| Cross-model scGPT/Geneformer alignment | `retire_now` | No Geneformer residual embeddings available |
| Rewiring-null PH survival | `retire_now` | Uniformly negative across 4 variants (iters 5-8); no path to rescue |
| Distance-permutation null PH | `retire_now` | Over-adversarial, no biological interpretation |
| Cross-layer CKA as primary question | `retire_now` | Fully answered (CKA≈1.0); demote to supporting evidence |
| GO clustering (209-gene version) | `deprioritize` | Directionally correct; revisit only after full vocab |

---

## New Hypothesis Portfolio

### Cluster A — Intrinsic Dimensionality (confirmed direction, mechanistic depth)

**A1 — MADA intrinsic dimensionality per layer**
- Hypothesis: MADA (manifold-adaptive dimension estimation) will independently confirm the layer-wise ID reduction seen by TwoNN and ER, and produce a third consistent ordering.
- Test: Run MADA on the 209-gene embedding submatrix at each of 12 layers. Compute Spearman rank correlation with TwoNN and ER orderings.
- Expected signal: Spearman >0.7 with both existing estimators; MADA value in 5–15 range.
- Null/control: Random permutation of layer assignments breaks the ordering.
- Value: medium | Cost: low (uses existing data)

**A2 — ER trajectory curvature / breakpoint detection**
- Hypothesis: The 6× ER collapse is not linear — there is a phase transition (steepest drop) at a specific layer range, corresponding to where scGPT most aggressively compresses representations.
- Test: Compute first and second derivative of ER across layers 0–11. Identify breakpoint layer using change-point detection (e.g., Pettitt test or CUSUM). Compare to CKA matrix to see if low-CKA consecutive pairs correlate with high ER drop.
- Expected signal: A clear breakpoint at layers 4–8 (mid-depth), where the model transitions from feature mixing to feature compression.
- Null: Smooth monotone model (no breakpoint); breakpoint layer is not enriched for CKA transitions.
- Value: medium | Cost: low

**A3 — Per-gene local intrinsic dimensionality**
- Hypothesis: High-drift genes (immune/transcription) occupy local embedding neighborhoods with lower intrinsic dimensionality (more constrained local geometry) than low-drift genes.
- Test: For each of 209 genes, compute TwoNN ID in the 20-NN local neighborhood across cells (not across genes). Compare ID distributions between top-50 and bottom-50 drift genes using Mann-Whitney U.
- Expected signal: Top-drift genes have lower local ID (more concentrated neighborhood), p < 0.05.
- Null: Random gene assignment to top/bottom groups.
- Value: medium | Cost: medium

---

### Cluster B — Biological Anchoring (most actionable direction)

**B1 — Full-vocab residual drift + GO enrichment [CRITICAL]**
- Hypothesis: With the full 4803-gene scGPT vocabulary, the transcription-regulation and immune enrichment of high-drift genes will replicate at FDR-corrected significance.
- Test: Recover scGPT vocab.json from model checkpoint. Map all 4803 gene indices to gene symbols. Rerun L2(layer_11 - layer_0) per gene. Split top/bottom 500 genes by drift. Fisher's exact test for 500+ GO BP terms. Apply BH FDR correction.
- Expected signal: ≥3 GO terms FDR < 0.05; transcription regulation and immune families survive FDR.
- Null: Random gene relabeling preserves no enrichment.
- Value: high | Cost: medium (vocab recovery is the bottleneck)

**B2 — Full-vocab GO term embedding clustering [CRITICAL]**
- Hypothesis: With 4803 genes, GO term gene sets will cluster significantly in cosine space (z < -1.96 in ≥3 layers).
- Test: Same cosine-distance clustering method as H01 (iter_0005), but with full vocabulary. Use 200 sampled GO BP terms with 5–200 members in the full vocab. 30 null replicates.
- Expected signal: Mean z < -2.0, at least 6/12 layers significant.
- Null: Random gene relabeling.
- Value: high | Cost: medium (depends on B1 vocab recovery)

**B3 — STRING PPI network density of high-drift gene subsets**
- Hypothesis: Genes that undergo the most representational drift (immune/transcription) form a denser STRING protein-protein interaction subnetwork than low-drift genes, indicating that scGPT's representational dynamics track known co-functional relationships.
- Test: Download STRING v12 interactions at confidence ≥700. For top-50 and bottom-50 drift genes (209-gene set), compute subgraph edge density. Compare to 10,000 random permutations of gene assignment.
- Expected signal: Top-drift subgraph edge density > 95th percentile of random null (p < 0.05).
- Null: Random gene sets of same sizes.
- Value: medium | Cost: medium (STRING download, manageable)

**B4 — Layer-by-layer drift trajectory shape per gene class**
- Hypothesis: The majority of residual drift does not accumulate uniformly across layers — it is front-loaded (layers 0–4) or back-loaded (layers 8–11), and the pattern differs between high-drift and low-drift gene classes.
- Test: Compute L2(layer_k - layer_0) for each gene at each layer k. Compute the mean cumulative drift curve for top-50 vs bottom-50 drift genes. Permutation test for difference in curve shape (e.g., area under drift-profile difference).
- Expected signal: High-drift genes accumulate drift in a specific layer range; the two curves are statistically different.
- Null: No difference in curve shape between high/low-drift gene classes.
- Value: medium | Cost: low (all data in hand)

---

### Cluster C — Topology (new family, bridge biology + PH)

**C1 — PH Betti curves on high-drift vs low-drift gene subsets**
- Hypothesis: The topological structure (H0, H1 Betti numbers) of high-drift gene embedding subsets differs from low-drift gene subsets across layers, with high-drift genes forming more complex 1-cycles in middle layers.
- Test: Split 209 genes into top-50 / bottom-50 by drift. At each layer, compute H0/H1 Betti curves (Vietoris-Rips, 5 radius values). Compare Betti numbers between gene groups using permutation test (reshuffle gene-to-group labels).
- Expected signal: H1 Betti count higher for high-drift genes at ≥4 layers, p < 0.05.
- Null: Random gene-to-group assignment (permutation of drift labels).
- Value: high | Cost: medium (PH computation is runtime-bounded)

**C2 — Persistent homology across layers (layer filtration)**
- Hypothesis: There exist topological features (H1 loops) in the gene embedding space that are "born" at early layers and "die" at later layers, marking a topological phase transition in scGPT's representational processing.
- Test: Treat the 12 layers as a parameter (not just a set). Compute a "layer filtration": for each pair of consecutive layers, check if H1 features born in one layer survive to the next (using correspondence of gene embedding positions). Track Betti number change from layer 0→11 as a persistence profile.
- Expected signal: H1 features decrease monotonically with layer depth (topological simplification mirrors ER compression).
- Null: Random ordering of layers produces equally monotone Betti profiles.
- Value: high | Cost: medium-high

**C3 — Mapper graph topology of gene embedding per layer**
- Hypothesis: Applying the Mapper algorithm to the gene embedding (projected by UMAP) at each layer will reveal topologically distinct clusters corresponding to biological gene classes (by GO annotation).
- Test: At each of 3 representative layers (0, 5, 11), compute 2D UMAP of 209-gene embeddings, apply Mapper (cover: 10 intervals, 50% overlap). Extract Mapper graph node-level GO annotation purity. Compare purity to random assignment.
- Expected signal: Mapper nodes at layer 11 are more GO-pure than at layer 0 (biological compression increases purity).
- Null: Random GO-label assignment within Mapper nodes.
- Value: high | Cost: high (Mapper is slow and needs tuning)

---

### Cluster D — Manifold Geometry

**D1 — Ollivier-Ricci curvature on kNN gene graph**
- Hypothesis: High-drift genes (immune/transcription) lie in negatively curved regions of the kNN gene graph, while low-drift genes (metabolic/structural) are in positively curved regions, reflecting different local manifold geometry.
- Test: Build kNN graph (k=10) on gene embeddings at each layer. Compute Ollivier-Ricci curvature (ORC) for each edge using GraphRicciCurvature. Correlate per-gene mean edge ORC with residual drift rank.
- Expected signal: Spearman correlation between ORC and drift magnitude < -0.2 at ≥6 layers.
- Null: Random drift-to-gene assignment.
- Value: medium | Cost: medium (ORC is O(k^2 × n))

**D2 — Singular vector biology (dominant ER subspace)**
- Hypothesis: The dominant singular vector of the gene embedding matrix at layer 11 (which captures nearly all variance given ER≈1.28) is strongly aligned with a known biological axis (e.g., immune vs non-immune gene expression, or transcription factor status).
- Test: At each layer, compute top-5 singular vectors of the gene embedding matrix (209 × 512). For SV1 at layer 11, compute gene loadings. Test whether top-20 loaded genes are enriched in GO:0002376 (immune) or GO:0006357 (transcription) using hypergeometric test.
- Expected signal: Top-loaded genes on SV1 at layer 11 are significantly enriched in immune/transcription GO terms (p < 0.05).
- Null: Random gene loading ordering.
- Value: high | Cost: low (SVD on 209×512 matrix is instant)

**D3 — CKA vs ER-drop correlation by layer pair**
- Hypothesis: Consecutive layer pairs with the steepest ER drop also show the lowest CKA similarity, indicating that representational compression and functional change are coupled.
- Test: Use existing CKA matrix (iter_0004) and ER values (iter_0005). Compute delta_ER for consecutive pairs (layer k+1 - layer k). Compute off-diagonal CKA for those pairs. Spearman correlation of |delta_ER| vs CKA for 11 consecutive pairs.
- Expected signal: Negative Spearman correlation (bigger ER drop = lower CKA), p < 0.05.
- Null: Random permutation of layer ordering.
- Value: medium | Cost: low (existing data only)

---

## Top 3 for Immediate Execution

### 1 — High-probability discovery candidate
**B1: Full-vocab residual drift + GO enrichment**

Rationale: The H02 signal (OR=3.58, p=0.004) at 209 genes is compelling. With 4803 genes the enrichment should reach FDR significance. This would establish the first robust functional biological finding of the project. The main execution dependency is recovering `vocab.json` from scGPT model checkpoint — this is a one-time data engineering step that unlocks multiple downstream hypotheses.

Execution plan:
1. Locate scGPT model files in subproject_38 or official scGPT repo.
2. Load `vocab.json`, map gene indices 0–4803 to gene symbols.
3. Rerun L2(layer_11 − layer_0) for all 4803 genes.
4. Run Fisher's exact test across 500+ GO BP terms (top vs bottom 500 genes by drift).
5. Apply BH FDR. Report FDR < 0.05 hits and compare to iter_0005 nominal hits.

---

### 2 — High-risk / high-reward candidate
**D2: Singular vector biology of the dominant ER subspace**

Rationale: The ER≈1.28 at layer 11 means the gene embedding has effectively collapsed to a single dominant direction. If that direction encodes biologically interpretable information (immune axis, TF status), this is a paper-worthy mechanistic finding — it would directly characterize what scGPT "decides" gene representations mean. Risk: the dominant SV may not have a clean biological interpretation.

Execution plan:
1. Compute full SVD of the 209×512 gene embedding matrix at each of 12 layers.
2. Extract top-5 singular vectors and gene loadings.
3. For SV1 at layer 11, rank genes by loading magnitude.
4. Hypergeometric test: top-20 loaded genes enriched in immune/transcription GO terms?
5. Plot SV1 loading vs drift magnitude (scatter, Spearman r).

---

### 3 — Cheap broad-screen candidate
**A2 + D3: ER trajectory curvature + CKA-ER correlation**

Rationale: Both analyses use exclusively existing data (ER from iter_0005, CKA matrix from iter_0004). Zero new computation beyond arithmetic. Together they characterize the coupling between spectral compression and functional layer-transition. Can be done in a single script in <30 minutes.

Execution plan:
1. Load `h03_effective_rank_per_layer.csv` and `h03_cka_matrix.npy`.
2. Compute delta_ER per consecutive layer pair (11 values).
3. Compute consecutive CKA (diagonal+1 of CKA matrix).
4. Spearman correlation |delta_ER| vs CKA; report p-value.
5. Change-point detection on ER curve (Pettitt test).
6. Identify breakpoint layer; annotate on ER/TwoNN combined plot.
