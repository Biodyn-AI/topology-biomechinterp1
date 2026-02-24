# Brainstormer Hypothesis Roadmap — iter_0024

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| TF activation hub centrality (H02 iter_0024) | Clear null: 0/12 sig, no effect. Degree-centrality encoding disproven. | retire_now |
| GO BP enrichment in SV2 poles | Retired iter_0011 H03: 0/591 terms sig. | retired_prev |
| TRRUST repression anti-pole | Retired iter_0013 H02: confirmed null. | retired_prev |
| Repression proximity (standalone) | 0/12 sig, 0.459 AUROC (iter_0023). Accept as null. | deprioritize |
| Pure SV pole tests (no new axis) | SV2/3/4 exhaustively characterized iter_0013-0017. | deprioritize |

---

## New Hypothesis Portfolio

### H-A: Multi-Predictor Joint Model (STRING + Dorothea + GO CC + GO BP)
**Hypothesis**: STRING score, Dorothea confidence, GO CC Jaccard, and GO BP Jaccard each contribute independently to predicting scGPT L2 embedding distance.
**Test**: OLS/ridge regression: L2_distance ~ STRING_score + Dorothea_conf_binary + GO_CC_jaccard + GO_BP_jaccard for all 21,736 named-gene pairs. Compute partial R², VIF, and semi-partial Spearman for each predictor. Run at each of 12 layers.
**Expected signal if true**: Each predictor has significant partial coefficient; total R² substantially exceeds any single predictor; VIF < 5 (low multicollinearity).
**Null/control**: Permuted-label regression at layer 8 (500 shuffles); also test random Gaussian embeddings.
**Value**: high | **Cost**: medium

### H-B: GO CC Subcellular Compartment Decomposition
**Hypothesis**: Specific GO CC compartments (nucleus, cytoplasm, plasma membrane, mitochondria, ER) drive the CC proximity signal at different rates across transformer layers.
**Test**: Assign each of 209 named genes to its most specific GO CC compartment. For each compartment with ≥10 genes, compute within-compartment AUROC (within-compartment pairs vs across-compartment). Compare AUROC magnitudes and layer curves across compartments.
**Expected signal if true**: Mitochondria and nucleus compartments have highest AUROC; layer curves differ (mito: early, nucleus: late, following the CC-dominant early vs BP-dominant late pattern from iter_0024).
**Null/control**: Random same-size gene partitions (500 replicates); chromosomally-matched gene sets.
**Value**: high | **Cost**: medium

### H-C: Single-TF Target Cluster Analysis (JUN/TP53 archetypes)
**Hypothesis**: The activation targets of a high-degree TF (JUN: 17 targets, TP53-related: high-degree) form a compact geometric cluster beyond what permutation expects.
**Test**: Extract JUN activation targets from TRRUST (N≈17). Compute mean pairwise L2 distance. Compare to 1000 random same-size gene sets from the 209-gene pool. Compute cluster z-score and AUROC at each of 12 layers.
**Expected signal if true**: JUN targets z < -2 (more compact than random) at mid/late layers; effect stronger than repressor targets.
**Null/control**: Same-size random gene sets; JUN repression target set as matched control.
**Value**: medium | **Cost**: low

### H-D: Chromosomal Proximity Negative Control
**Hypothesis**: Genes on the same chromosome do NOT cluster more in scGPT embeddings than genes on different chromosomes (functional encoding, not genomic artifact).
**Test**: Map 209 named genes to chromosomes via mygene API. Compute AUROC: same-chromosome pairs vs different-chromosome pairs at each of 12 layers. Also test same chromosome arm (p/q) as stricter grouping.
**Expected signal if true**: AUROC ≈ 0.5 (null), confirming the functional proximity signals are not confounded by chromosomal co-localization.
**Null/control**: The expected null IS the hypothesis; p-value > 0.05 confirms specificity.
**Value**: medium | **Cost**: low

### H-E: Geneformer Cross-Model GO CC Alignment
**Hypothesis**: GO CC proximity signal in scGPT SV geometry is recapitulated in Geneformer static token embeddings, confirming cross-model encoding of subcellular localization.
**Test**: Use the 207 Geneformer gene token embeddings from iter_0019. Compute pairwise cosine similarity. Compute Spearman(GO_CC_Jaccard, cosine_sim). Compare to scGPT L8 Spearman(GO_CC_Jaccard, -L2_distance).
**Expected signal if true**: Geneformer also shows positive Spearman ≥ 0.07 for GO CC (lower than scGPT's 0.106 is acceptable, given static vs contextual); cross-model consistency confirmed.
**Null/control**: Random Gaussian 1152-dim embeddings for same 207 genes; permuted gene labels.
**Value**: high | **Cost**: low

### H-F: Topological Persistence of GO CC vs STRING Clusters (Comparative PH)
**Hypothesis**: GO CC-defined gene clusters have higher persistent homology H0 lifetime than STRING-defined clusters, reflecting more compact subcellular groupings.
**Test**: Build gene sets by GO CC compartment (nucleus, cytoplasm, mito, membrane, ≥8 genes each). Run Ripser H0+H1 on each set's embedding submatrix at layers 4, 8, 11. Compare H0 lifetimes to STRING quintile subsets of same sizes.
**Expected signal if true**: GO CC groups show H0 persistence ≥ STRING Q4/Q5 groups, indicating subcellular localization produces tighter manifold topology than PPI confidence alone.
**Null/control**: Randomly sampled same-size gene sets (500 replicates per layer-compartment).
**Value**: high | **Cost**: high

### H-G: Cell-Type Marker GO CC Enrichment (Integration Hypothesis)
**Hypothesis**: The strong cell-type marker clustering (AUROC=0.851) is partially explained by GO CC co-compartment membership (immune markers share surface/membrane GO CC terms; fibroblast markers share ECM/cytoplasm terms).
**Test**: Compute GO CC Jaccard for the 17-gene expanded cell-type marker panel (iter_0023). Partial Spearman(GO_CC_Jaccard, -L2_distance) after controlling for cell-type same/different membership. Test if cell-type AUROC drops after partialing out GO CC.
**Expected signal if true**: GO CC explains ~30-50% of the cell-type clustering signal; partial AUROC drops from 0.851 toward 0.7.
**Null/control**: Permuted GO CC assignments.
**Value**: medium | **Cost**: low

### H-H: Layer-Resolved Encoding Timeline (CC early, BP late, STRING mid)
**Hypothesis**: Different biological anchors peak at different transformer layers, revealing a layer-resolved biological encoding timeline: GO CC peaks early (L4-5), STRING PPI mid (L6-8), GO BP late (L7-9), cell-type identity monotonically deepens.
**Test**: Collect Spearman(anchor_jaccard, -L2_dist) for GO CC, GO MF, GO BP, STRING confidence per layer for all 12 layers. Fit peak layer per anchor. Test if the ordering CC < STRING < BP is significant (permutation of layer assignments).
**Expected signal if true**: Peak layers significantly differ across anchors (permutation test). Timeline: CC peaks ~L5, STRING ~L7-8, BP ~L7-8 (already suggested by iter_0024 data).
**Null/control**: Permuted layer ordering (500 replicates).
**Value**: high | **Cost**: low (uses existing computed data)

### H-I: Dorothea Directionality (Activation-Enriched Confidence Tiers)
**Hypothesis**: Within Dorothea A/B pairs, those with known activation direction are geometrically closer than those with repression direction (extending the TRRUST activation asymmetry to an independent database).
**Test**: If Dorothea confidence tiers encode direction (some versions do), split A/B pairs by direction. Compute AUROC(activation_A/B vs background) vs AUROC(repression_A/B vs background). If no direction info, use regulatory mode from TRRUST for TFs that appear in both databases.
**Expected signal if true**: Dorothea activation pairs show higher AUROC than repression pairs, replicating TRRUST asymmetry (iter_0023 H03: act=0.640, rep=0.459) in an independent regulatory database.
**Null/control**: Permuted direction labels within A/B tier.
**Value**: high | **Cost**: low

### H-J: STRING Score Threshold Sensitivity (0.4 vs 0.5 vs 0.7 vs 0.9)
**Hypothesis**: The STRING AUROC increases monotonically as the score threshold increases (0.4→0.9), confirming confidence-graded encoding extends to the full precision-recall curve.
**Test**: Compute binary AUROC at thresholds 0.4 (N=3092), 0.5, 0.7 (N=1022), 0.9 at each of 12 layers. Test if AUROC increases monotonically with threshold (Spearman across 4 thresholds).
**Expected signal if true**: Spearman ≥ 0.9 across layers; AUROC at 0.9 threshold > 0.65.
**Null/control**: Permuted STRING scores.
**Value**: medium | **Cost**: low

### H-K: Intrinsic Dimension per Biological Anchor Subset
**Hypothesis**: Gene subsets defined by tight biological anchors (same GO CC compartment, same STRING Q5 pairs, same cell-type) occupy lower-dimensional manifolds than random gene subsets.
**Test**: Apply TwoNN estimator to per-compartment gene subsets (≥15 genes) at each layer. Compare ID to random same-size subsets (100 replicates). Also compare ID of STRING Q5 gene set vs Q1 gene set.
**Expected signal if true**: GO CC compartment subsets and STRING Q5-connected genes have significantly lower TwoNN ID than random; consistent with tight functional clustering.
**Null/control**: Random same-size gene subsets.
**Value**: medium | **Cost**: medium

### H-L: Multi-Predictor AUROC vs Dorothea Tier in Joint STRING-Dorothea Model
**Hypothesis**: After controlling for STRING PPI score, Dorothea A/B confidence adds independent predictive power (partial AUROC gain > 0.01) because regulatory network confidence encodes non-PPI co-program structure.
**Test**: For the subset of gene pairs with both STRING and Dorothea annotations, build a logistic regression predictor: STRING_score + Dorothea_tier → L2_distance class (close/far). Compute partial AUROC for each predictor. Test at layer 8.
**Expected signal if true**: Partial AUROC for Dorothea > 0.01 after STRING controlled; VIF < 3 (low collinearity since Dorothea and STRING are different evidence types).
**Null/control**: Permuted Dorothea tier assignments; STRING-only baseline.
**Value**: high | **Cost**: medium

### H-M: GO CC + Cell-Type Composite Predictor for Embedding Distance
**Hypothesis**: A composite predictor combining GO CC Jaccard and cell-type same/different identity outperforms either alone (AUROC > 0.85 for the 17-gene marker panel).
**Test**: For the 17 cell-type markers, compute: (a) GO CC Jaccard only, (b) cell-type identity only, (c) GO CC Jaccard + cell-type identity logistic regression. Measure AUROC against L2 distance class at each layer.
**Expected signal if true**: Composite AUROC > max(individual AUROCs) by ≥ 0.02.
**Null/control**: Random GO CC + random cell-type assignments.
**Value**: medium | **Cost**: low

---

## Top 3 for Immediate Execution

### #1 — High-Probability Discovery: Multi-Predictor Joint Model (H-A)
This is the most pressing consolidation experiment in the project. Five independent positive signals (STRING, TRRUST, GO BP, GO CC, Dorothea) have never been tested together. The partial R² and VIF analysis will either confirm that these represent genuinely orthogonal biological dimensions of scGPT geometry — the core publishable claim — or reveal that some are confounded. This is the keystone experiment for the paper's main thesis.

**Implementation**: OLS regression at each of 12 layers, partial Spearman per predictor, VIF calculation. Requires joining all pair-level annotations already computed.

### #2 — High-Risk/High-Reward: Layer-Resolved Encoding Timeline (H-H)
The iter_0024 data already hints at CC peaking earlier than BP. If this resolves into a clean biological encoding timeline (GO CC early → STRING mid → GO BP late → cell-type deepening), this is a novel mechanistic finding about how transformer layers progressively encode biological abstraction levels. Uses existing data: no new API calls needed, just cross-iteration data assembly. The risk is that the layer differences are not significant.

**Implementation**: Compile Spearman-per-layer curves for all anchors from existing JSON artifacts. Fit peak layers. Permutation test on layer ordering. Cost: low.

### #3 — Cheap Broad-Screen: Chromosomal Proximity Negative Control + Dorothea Directionality (H-D + H-I combined)
H-D confirms that all positive signals are functionally specific (not genomic artifact) — a necessary specificity control for the paper. H-I extends the activation asymmetry discovery to Dorothea. Both are cheap (mygene chromosome lookup + Dorothea direction info if available). Running both in the same iteration maximizes coverage.

**Implementation**: mygene API for chromosome annotations, same-chr vs diff-chr AUROC at 12 layers. Dorothea direction field check. Combined script, ~50 lines.
