# Hypothesis Roadmap: Post-iter_0050

## Retire / Deprioritize

### RETIRE: "SV2-4 encodes regulatory proximity specifically"
- **Reason**: Falsified by iter_0050 H01. Co-expression pairs are always closer in SV2-4 than TRRUST pairs. The specificity claim does not hold.
- **What survives**: SV2-4 encodes co-expression structure enriched for TF modules. Keep as revised claim.
- **Status**: `retire_now` (the original framing); replace with co-expression framing.

### DEPRIORITIZE: Further TRRUST proximity testing in raw SV2-4
- **Reason**: The effect is now explained (co-expression confound). Additional proximity tests in SV2-4 without controlling for co-expression are uninformative.
- **Status**: `rescue_once_with_major_change` — only revisit if residual-regression approach is applied.

### DEPRIORITIZE: Generic PH on full 512-d embedding
- **Reason**: Multiple early iterations tested full-space PH and found no significant signal over background. The signal lives in spectral subspaces.
- **Status**: `retire_now` for full-space variant. Keep only for SV2-4 subspace variant.

---

## New Hypothesis Portfolio

### H-A: SV1 Rotation Driver Gene Identity + Enrichment
**Hypothesis**: The semantic inversion of SV1 (anti-norm at L0, pro-norm at L11) is driven by a minority of genes with extreme sign-flip in SV1 loading, enriched for housekeeping genes and/or TFs.
**Test**: Load h01_sv1_vectors.npy [12, 2039]. Compute per-gene |SV1_L11 - SV1_L0| (signed delta). Rank top-50 positive (L11-dominant) and top-50 negative (L0-dominant) movers. Map to gene names. Test enrichment vs: (a) Eisenberg housekeeping gene list, (b) TRRUST TF list, (c) immune cell markers, (d) ribosomal genes.
**Expected signal**: L0-dominant movers = low-expression genes (mitochondrial, ribosomal) that get pushed into SV1 by early layers; L11-dominant movers = high-expression cell-type markers that anchor SV1 after rotation.
**Null/control**: Permuted gene labels; top-50 random genes vs. same lists.
**Value**: high | **Cost**: low (all data exists as artifacts)

### H-B: Residual TF Enrichment After Co-expression Correction in SV2-4
**Hypothesis**: After regressing out embedding cosine similarity (co-expression proxy) from SV2-4 Euclidean distance, a residual TF-target proximity signal survives — indicating that beyond co-expression, regulatory pairs are more proximate.
**Test**: At each layer, fit a linear model: SV2-4_distance ~ embedding_cosine_similarity. Compute residuals. Run Mann-Whitney on residuals: TRRUST-positive vs. TRRUST-negative. Report rbc per layer.
**Expected signal**: If TF-target pairs cluster due to regulatory topology beyond co-expression, residual rbc > 0 in mid/late layers. If the entire SV2-4 signal is co-expression mediated, residuals are flat.
**Null**: rbc = 0 at all layers in residuals.
**Value**: high | **Cost**: low (regression is trivial, data exists)

### H-C: SV5-7 Orthogonal TF Enrichment (Independence Test)
**Hypothesis**: SV5-7 TF enrichment (rbc=0.152, iter_0050 H03) is partially orthogonal to SV2-4 — projecting gene loadings onto SV5-7 AFTER removing the SV2-4 component yields residual TF enrichment.
**Test**: For each gene at each layer, compute loading on SV5-7 space then project out the SV2-4 component (Gram-Schmidt or regression). Run Mann-Whitney on residual SV5-7 loading magnitudes: TF genes vs. non-TF. Also run for target genes. Report rbc.
**Expected signal**: If SV5-7 is independent, TF enrichment survives residualization (rbc > 0.05 in ≥5 layers).
**Null**: rbc collapses to ≈0 after SV2-4 removal.
**Value**: high | **Cost**: low (extends iter_0050 H03 code)

### H-D: SV5-7 GO Biological Process Enrichment Screen
**Hypothesis**: The top genes by SV5-7 projection magnitude are enriched for distinct GO terms from SV2-4 top genes, indicating SV5-7 captures a different biological dimension.
**Test**: At peak layer (L8), rank genes by projection magnitude onto SV5-7 vs. SV2-4. Take top-100 from each. Run hypergeometric GO-BP enrichment (gseapy) for each set. Compare enriched terms.
**Expected signal**: SV2-4 top genes = co-expression hubs (TFs, signaling genes); SV5-7 top genes = different functional category (e.g., cell cycle, metabolism, RNA processing).
**Null**: Same GO terms enriched in both sets.
**Value**: medium | **Cost**: low

### H-E: TF Hub Degree Explains SV2-4 Loading
**Hypothesis**: In the co-expression network derived from gene embeddings (edge = top-k cosine similarity), TF genes have higher degree than non-TFs, and hub degree predicts SV2-4 projection magnitude.
**Test**: Build gene co-expression graph: connect pairs with cosine similarity ≥ threshold (e.g., ≥0.8) from full-512-d embeddings at L8. Compute degree for all 2039 genes. Spearman correlation: degree vs. SV2-4 projection magnitude. Also: Mann-Whitney degree for TF vs. non-TF genes.
**Expected signal**: Positive correlation (hub genes → high SV2-4 loading); TF degree > non-TF degree. If so, the SV2-4 TF enrichment is fully explained by TFs being hubs.
**Null**: No correlation between degree and SV2-4 loading.
**Value**: high | **Cost**: low

### H-F: 0-Dim Persistent Homology on SV2-4 Distance Matrix (TRRUST Merge Test)
**Hypothesis**: In the Vietoris-Rips filtration on SV2-4 distance matrix at peak layer L8, TRRUST TF-target pairs merge into the same connected component at lower filtration radius than random gene pairs.
**Test**: Build n=2039 gene × gene SV2-4 Euclidean distance matrix at L8. Run 0-dim PH (Gudhi or scikit-tda). For each TRRUST pair: record the birth radius of their connected component merge. Compare distribution to shuffled pairs. Mann-Whitney.
**Expected signal**: TRRUST pairs merge at lower radius (topological proximity). rbc > 0.1.
**Null**: Merge radii identical to random pairs.
**Value**: high | **Cost**: medium (PH computation on 2039×2039 is tractable)

### H-G: Per-Gene SV1 Loading Trajectory Autocorrelation
**Hypothesis**: Genes with stable SV1 loading across all 12 layers (high autocorrelation of loading trajectory) are enriched for housekeeping/constitutively expressed genes; genes with low autocorrelation are enriched for context-dependent/regulated genes.
**Test**: For each gene, compute the 12-point SV1 loading trajectory. Compute lag-1 autocorrelation. Rank by autocorrelation. Mann-Whitney: housekeeping gene set (Eisenberg) vs. background. Also test: TFs vs. background.
**Expected signal**: Housekeeping genes show high autocorrelation (stable SV1 identity); TFs show lower autocorrelation (context-dependent).
**Null**: No difference between gene classes.
**Value**: medium | **Cost**: low

### H-H: SV2-4 Cluster Structure and Regulatory Module Coherence
**Hypothesis**: Hierarchical clustering of genes in SV2-4 space at L8 produces clusters with GO term coherence that exceed full-512-d clustering.
**Test**: At L8, cluster n=2039 genes using agglomerative clustering (Ward linkage) in SV2-4 subspace (3 dims) vs. full 512-d space. For each cluster, compute mean pairwise GO-BP semantic similarity (using GOSemSim or simonsays). Compare coherence scores between subspace and full-space clustering.
**Expected signal**: SV2-4 clustering shows higher within-cluster GO coherence, validating that the subspace captures biologically organized structure.
**Null**: Full-512-d clustering outperforms or equals SV2-4.
**Value**: medium | **Cost**: medium (GOSemSim computation)

### H-I: STRING Co-expression Network Signal in SV2-4 (Replication + Specificity)
**Hypothesis**: STRING high-confidence edges (score ≥700) show SV2-4 proximity similar to or stronger than TRRUST — quantifying the co-expression explanation.
**Test**: Download STRING human network. Map to gene symbols. Identify STRING edges among the 2039-gene set. Run SV2-4 proximity test (Mann-Whitney, same as H01 but using the actual STRING download, not cosine proxy). Compare STRING rbc vs. TRRUST rbc per layer.
**Expected signal**: STRING rbc ≥ TRRUST rbc at most layers, consistent with co-expression explaining the TRRUST signal.
**Null**: TRRUST rbc > STRING rbc (regulatory-specific signal).
**Value**: medium | **Cost**: medium (STRING download required)

### H-J: Cross-Layer Geodesic Length as Gene-Level Predictor
**Hypothesis**: The total geodesic length of a gene's embedding trajectory (sum of L2 distances between consecutive layer embeddings) predicts expression level and TF status — TFs have longer trajectories, reflecting more processing/transformation.
**Test**: For each of 2039 genes, compute path length = sum_{l=0}^{10} ||emb[l+1] - emb[l]||. Correlate with mean expression. Mann-Whitney: TF genes vs. non-TF.
**Expected signal**: High-variance/high-expression genes show longer paths; if TFs are transformed more during forward pass, TF path length > non-TF.
**Null**: No correlation between path length and biological category.
**Value**: medium | **Cost**: low (simple computation from existing embeddings)

### H-K: Spectral Decay Rate Quantification
**Hypothesis**: The gradient of TF enrichment across SV groups (SV2-4 > SV5-7 > SV8-10) follows a geometric decay in rbc that is predictable, indicating a single underlying variance dimension for regulatory information.
**Test**: Fit linear model: mean_rbc ~ log(SV_group_rank). Test R² and slope. Also: compute rbc for SV11-13, SV14-16 at L8 to extend the decay curve.
**Expected signal**: Geometric (log-linear) decay in rbc with SV group. If so, the regulatory information is spread continuously across SV dimensions in a predictable way.
**Null**: Non-monotonic or step-function decay.
**Value**: medium | **Cost**: low (extends H03)

### H-L: SV1 Early vs. Late Gene Sets — Cell-Type Specificity
**Hypothesis**: Genes with high SV1 loading at L0 (anti-norm population) are enriched for cell-type-specific markers; genes with high SV1 loading at L11 (pro-norm population) are enriched for broadly expressed housekeeping genes.
**Test**: Define L0-dominant genes: top-200 by |SV1 loading| at L0. Define L11-dominant genes: top-200 by |SV1 loading| at L11. Overlap these sets with known cell-type marker gene lists (e.g., PanglaoDB or CellMarker for immune cell types) and housekeeping gene lists.
**Expected signal**: L0-dominant genes = cell-type-specific markers; L11-dominant = broadly expressed. This would explain the anti-correlation at L0 (cell-type markers have variable/low mean expression) and pro-correlation at L11 (housekeeping genes are high-norm).
**Null**: No differential enrichment.
**Value**: high | **Cost**: low (gene name mapping needed once)

### H-M: TRRUST Module Cohesion in SV2-4 (Within-TF-Module Distance)
**Hypothesis**: Genes regulated by the same TF (a "TF module") form tighter clusters in SV2-4 space than genes regulated by different TFs — even controlling for total pairwise count.
**Test**: For each TF with ≥5 targets in the 2039-gene set, compute mean pairwise SV2-4 distance among its target set at L8. Compare to mean pairwise distance among random gene sets of same size. Permutation test.
**Expected signal**: TF modules show lower within-module SV2-4 distance than random. This would confirm that co-expression structure is organized by TF regulon, not just generic similarity.
**Null**: No difference vs. random.
**Value**: high | **Cost**: low

---

## Top 3 for Immediate Execution

### #1 — High-Probability Discovery: H-A (SV1 Rotation Driver Gene Identity)
**Rationale**: All data exists (h01_sv1_vectors.npy + gene names from embedding indices). This is the most direct follow-up to the iter_0050 H02 positive. The result will either confirm housekeeping/cell-type identity of the two SV1 populations (strong narrative anchor) or reveal unexpected biology (high value either way). Zero external data dependency.
**Cost**: Low. Expected runtime: <5 min.
**Deliverable**: Named gene lists + enrichment table for paper's biological validation section.

### #2 — High-Risk/High-Reward: H-F (0-Dim PH on SV2-4 at L8)
**Rationale**: The co-expression confound showed SV2-4 proximity is not TF-specific by distance. But persistence (topology) asks a different question: do TF-target pairs belong to the same connected component at lower filtration threshold? If yes, this is a genuinely new result — TF pairs are topologically co-clustered even when individual distances are ambiguous. If no, we can fully retire the SV2-4 regulatory topology claim.
**Cost**: Medium. VR on 2039×2039 matrix is tractable (<30 min with Gudhi).
**Deliverable**: Merge-radius distribution plot; TF pair vs. random comparison.

### #3 — Cheap Broad Screen: H-B (Residual TF Enrichment After Co-expression Correction)
**Rationale**: The simplest possible test to determine whether any TF-specific signal survives the co-expression confound. Linear regression residuals take <1 min. If residual rbc > 0.05 in multiple layers, the TF signal is partially orthogonal to co-expression and the paper's regulatory claim can be partially rescued with stronger framing. If rbc ≈ 0, we retire the regulatory claim cleanly.
**Cost**: Low. Expected runtime: <5 min.
**Deliverable**: Table of residual rbc per layer; binary yes/no for residual regulatory signal.
