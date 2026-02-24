# Brainstormer Hypothesis Roadmap — iter_0013

---

## Retire / Deprioritize

| Hypothesis | Reason | Status |
|---|---|---|
| Repression anti-pole (H02, iter_0013) | Definitive negative across all 12 layers; z=−1.41 mean | **RETIRED** |
| Bootstrap CI act/rep differential (H01, iter_0012) | Negative, underpowered, superseded by cleaner tests | **RETIRED** |
| SV1 as PPI encoding axis | z=0.39, 2/12 layers — noise-level | **RETIRED** |
| Any direct regulatory-sign-as-spatial-opposition hypothesis | Consistent negatives across multiple formulations | **RETIRED for this direction** |

---

## New Hypothesis Portfolio

### Group A: Characterizing SV2+SV3 Joint Geometry

**A1: SV2+SV3 2D scatter reveals discrete PPI community structure**
- Hypothesis: Genes in the same STRING community (Louvain or MCL cluster on the PPI graph) form localized clusters in the SV2×SV3 2D projection, beyond what either axis alone shows.
- Test: Run community detection on STRING PPI graph (score≥0.7). Plot 209 genes in SV2×SV3 space colored by community. Compute inter- vs intra-community SV2+SV3 Euclidean distance ratio, compare to null (random community relabeling).
- Signal if true: Intra-community distance significantly less than null (z>3), visible structure in 2D scatter.
- Null: Distance ratio = 1.0 under random relabeling.
- Value: **high** | Cost: **low**

**A2: SV3 pole genes are enriched for specific GO biological process categories distinct from SV2**
- Hypothesis: SV3 top/bottom poles capture a biologically distinct gene functional grouping (e.g., signaling vs. metabolic) relative to SV2, revealing a second interpretable axis in scGPT's representation.
- Test: Identify top-K=52 and bottom-K=52 SV3 pole genes (averaged across layers). Run GO biological process enrichment (Fisher's exact or hypergeometric, using goatools or g:Profiler API). Compare SV3 enriched terms to SV2 enriched terms.
- Signal if true: SV3 poles have GO terms distinct from SV2, suggesting orthogonal functional decomposition.
- Null: No GO enrichment, or same terms as SV2.
- Value: **high** | Cost: **low**

**A3: SV2+SV3 trajectory across layers shows progressive functional specialization**
- Hypothesis: The angle between the SV2 and SV3 axes (as vectors in embedding space) changes systematically across layers 0→11, reflecting increasing specialization of functional geometry in deeper layers.
- Test: Compute cosine angle between SV2 and SV3 (or their correlation in gene-rank space) at each of 12 layers. Test for monotonic trend (Spearman ρ across layers).
- Signal if true: Spearman ρ > 0.6 (or < −0.6) with p<0.05, indicating consistent layer-wise change.
- Null: No layer trend (ρ ≈ 0).
- Value: **medium** | Cost: **low**

### Group B: Hub Confound Resolution & Expanded PPI Tests

**B1: Degree-matched null for PPI co-pole — does signal survive?**
- Hypothesis: The SV2 co-pole enrichment for STRING PPI pairs is not explainable by degree alone; it persists when null pairs are matched to the same degree distribution.
- Test: For each STRING pair (gene_i, gene_j), construct null pairs by sampling gene_k, gene_l where degree(k)≈degree(i) and degree(l)≈degree(j) (±2 degree). Run co-pole test with this degree-matched null vs observed.
- Signal if true: z>3 even with degree-matched null (effect not explained by hubs).
- Null: z<1.96 with degree-matched null (pure hub artifact).
- Value: **high** | Cost: **medium**

**B2: Expanded STRING at score≥0.4 resolves hub/non-hub question**
- Hypothesis: At lower confidence threshold, sufficient non-hub pairs exist to test whether low-degree PPI partners also co-localize in SV2.
- Test: Download/load STRING at score≥0.4. Re-split into hub (degree>median) and non-hub. Run SV2 co-pole test on N>200 non-hub pairs. Compare to degree-matched null.
- Signal if true: Non-hub co-pole z>2 in ≥6/12 layers.
- Null: Non-hub z<1.5 (signal is hub-only).
- Value: **high** | Cost: **low** (just lower threshold, data likely cached)

### Group C: Cross-Regulatory-Network Tests

**C3: TRRUST indirect targets (TF→TF→target chains) show weaker SV2 co-pole than direct**
- Hypothesis: Direct TF-target pairs are more co-localized in SV2 than indirect (2-hop) TF-target pairs, showing that scGPT encodes regulatory proximity not just co-expression.
- Test: Extract all 2-hop chains in TRRUST (TF1→TF2→target). Compute SV2 co-pole for direct pairs (N=113) vs 2-hop pairs. Compare using z-scores against matched null.
- Signal if true: Direct pairs z > 2-hop pairs z, consistent across layers.
- Null: No difference in z between direct and 2-hop.
- Value: **medium** | Cost: **low**

**C4: GO co-annotation predicts SV2 co-pole better than STRING PPI**
- Hypothesis: Genes sharing GO biological process terms (even without direct PPI) are co-localized in SV2, suggesting scGPT learned functional neighborhoods not just physical interaction.
- Test: For each gene pair in 209-gene set, define GO_BP co-annotation as sharing ≥2 GO BP terms. Compare SV2 co-pole rate for GO-co-annotated pairs vs non-annotated vs STRING PPI pairs.
- Signal if true: GO-co-annotated pairs show z>3, comparable or exceeding STRING PPI.
- Null: GO co-annotation provides no SV2 co-pole enrichment.
- Value: **high** | Cost: **medium**

### Group D: Topological / Manifold Structure

**D1: Persistent homology (H0) of SV2+SV3 2D embedding — does PPI community structure leave a barcode signature?**
- Hypothesis: The SV2×SV3 2D embedding of 209 genes has topological features (H0 persistent components) that correspond to STRING communities, detectable via persistent homology.
- Test: Compute distance matrix in SV2+SV3 space. Run VR persistent homology (H0). Compare barcodes to a community-shuffled null. Test if H0 persistence landscape differs.
- Signal if true: Significantly longer H0 barcodes for community-labeled clusters than random.
- Null: PH barcodes match random-point Rips filtration.
- Value: **high** | Cost: **medium**

**D2: Local intrinsic dimension of PPI neighborhoods is lower than non-PPI neighborhoods**
- Hypothesis: Genes that form tight PPI cliques occupy lower-dimensional local manifolds in the 512-dim embedding space compared to less-connected gene sets.
- Test: For each STRING community (≥5 members), compute local intrinsic dimension using TwoNN or PCA-variance-explained method within that gene cluster. Compare to same-size random gene sets.
- Signal if true: Mean local ID for PPI communities significantly less than random sets (Wilcoxon p<0.05).
- Null: No difference in local ID.
- Value: **high** | Cost: **medium**

**D3: Geodesic distance in kNN graph of embedding vs Euclidean — which better predicts PPI?**
- Hypothesis: Geodesic distances (on the kNN graph of gene embeddings) are more predictive of PPI links than straight Euclidean distances, suggesting curved/manifold geometry matters.
- Test: Build kNN graph (k=10) on 209-gene embeddings. Compute shortest-path geodesic distance. Compare AUC for STRING PPI prediction: geodesic distance vs Euclidean distance vs SV2 projection distance.
- Signal if true: Geodesic AUC > Euclidean AUC > SV2 AUC (or geodesic > Euclidean with p<0.05).
- Null: No difference between distance metrics.
- Value: **medium** | Cost: **medium**

### Group E: Cross-Model / Generalization

**E1: Random embedding (untrained scGPT) shows no SV2 PPI co-pole signal — confirming training-dependence**
- Hypothesis: The SV2 PPI signal is specific to trained scGPT weights; randomly initialized weights of the same architecture produce z≈0.
- Test: Load scGPT architecture, randomize weights. Compute 209-gene embeddings. Run SV2 co-pole test for STRING PPI. Compare to trained model z=10.04.
- Signal if true: Random model z<2 (training-dependent signal confirmed).
- Null: Random model also shows z>3 (signal is architecture/data artifact, not learned).
- Value: **high** | Cost: **medium** (requires model init code)

**E2: scGPT SV2 co-pole vs. standard scRNA-seq PCA SV2 — is scGPT special?**
- Hypothesis: The SV2 PPI co-pole enrichment is stronger in scGPT embeddings than in standard PCA of raw scRNA-seq counts for the same genes/cells, confirming scGPT adds biological structure.
- Test: From the same cell data, compute PCA and take SV2. Run identical co-pole test for STRING PPI. Compare z-scores scGPT-SV2 vs PCA-SV2.
- Signal if true: scGPT z=10 >> PCA z (e.g., PCA z<4), confirming language model advantage.
- Null: PCA-SV2 gives similar or better z (scGPT adds nothing over linear methods).
- Value: **high** | Cost: **medium**

### Group F: Mechanistic Alignment

**F1: Layer-depth curve of SV2 z-score follows attention depth gradient**
- Hypothesis: The SV2 PPI z-score increases monotonically with layer depth (0→11), reflecting progressive encoding of biological structure in deeper transformer layers.
- Test: Extract per-layer z-scores from h01_sv_axis_specificity.json (already computed). Compute Spearman ρ between layer index and SV2 z-score.
- Signal if true: Spearman ρ > 0.7, p<0.05 (later layers encode PPI structure more strongly).
- Null: No monotonic trend; z-scores plateau or are flat.
- Value: **medium** | Cost: **low** (data already exists in h01 artifact)

---

## Top 3 for Immediate Execution

### #1 — High-Probability Discovery Candidate
**A2: GO biological process enrichment of SV3 poles**

Rationale: SV3 carries robust signal (z=7.18, 12/12 layers) but is biologically uncharacterized. GO enrichment is cheap (one API call or local annotation), high interpretability payoff, directly publishable, and has high prior probability of finding something (SV3 must encode *something* real given the z-score). This is the next natural characterization step.

### #2 — High-Risk / High-Reward Candidate
**D2: Local intrinsic dimension of PPI communities vs random sets**

Rationale: If PPI communities occupy genuinely lower-dimensional local manifolds, this would be a novel geometric finding about how biological networks are encoded. No one has tested this for gene language models. Risk: small gene set (209 genes) may not have enough members per community to estimate local ID reliably. Reward: manifold geometry finding would be a strong contribution.

### #3 — Cheap Broad-Screen Candidate
**F1: Layer-depth Spearman trend on already-computed SV2 z-scores**

Rationale: The per-layer SV2 z-scores are already computed in h01_sv_axis_specificity.json. This test costs ~5 lines of code and <1 minute. If layer depth correlates with z-score, it reveals that scGPT progressively builds PPI structure across layers — a mechanistic insight requiring zero new data.
