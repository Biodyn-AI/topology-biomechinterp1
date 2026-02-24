# Brainstormer Hypothesis Roadmap — iter_0029

---

## Retire / Deprioritize

| Hypothesis family | Status | Reason |
|------------------|--------|--------|
| PC1 binary biological axis | **retire_now** | Carried from iter_0028 retirement. No rescue path. |
| TRRUST activation/repression polarity | **retire_now** | Null confound confirmed iter_0027. Retired. |
| Dorothea confidence-tier pairwise distance | **retire_now** | Inconclusive, low power. Retired. |
| PC1 T-cell vs APC axis | **retire_now** | Structural impossibility confirmed iter_0028. |
| Fiedler vector biological partition (H-A) | **deprioritize** | Superseded by H01 connected component finding. The 2-component split is already the biological partition. Fiedler analysis would be derivative. |
| Spectral gap k-robustness (H-C) | **retire_now** | Completed and confirmed by H03 iter_0029. Done. |

---

## New Hypothesis Portfolio

### H-N01: Layer of Bifurcation — When Do the 14 Genes Break Away?
**Hypothesis:** The diverging-vs-converging trajectory separation has a sharply defined onset layer; before this layer all 209 genes move together, after it the 14 outlier genes diverge.
**Test:** Compute per-gene distance-to-centroid slope in a sliding 3-layer window (L0–L2, L1–L3, ..., L9–L11). For each gene, identify the window where slope sign flips. Compare flip-layer distributions: diverging-14 vs converging-195. Wilcoxon test.
**Expected signal:** Diverging genes show slope sign flip at L6–L8 (mid-network); converging genes show consistent negative slope throughout. Tight peak in flip-layer distribution for the 14.
**Null:** Same analysis on random Gaussian embeddings (same shape).
**Value:** high | **Cost:** low

---

### H-N02: GO Enrichment of the 14 Diverging Genes
**Hypothesis:** The 14 diverging genes (FOS, JUNB, HLA-A, HLA-DPB1, KLF6, LDHA, LGALS1, NCAM1, NCOA3, NR4A3, PAX5, PTGS2, TBXAS1, TNF) are significantly enriched for GO Biological Process terms related to stress response, immune activation, and context-dependent gene regulation.
**Test:** Run Fisher exact test (or use goatools) for GO BP term membership of the 14 genes vs the background of all 209 genes. Report Bonferroni-corrected p-values. Separately test GO CC (cellular compartment) and MF (molecular function).
**Expected signal:** Enrichment in "response to cytokine," "inflammatory response," "response to stress," or "cell activation" (GO BP). Depletion of "metabolic process" or "cell cycle." No specific cell-compartment bias expected a priori.
**Null:** Same Fisher test on 1000 random 14-gene draws from 209 genes.
**Value:** high | **Cost:** low

---

### H-N03: Bootstrap Stability of the 14-Gene Cluster
**Hypothesis:** The 14-gene diverging cluster identity is stable under bootstrap resampling of the layer dimension.
**Test:** Bootstrap the 12-layer distance-to-centroid matrix (resample layers with replacement, 1000 iterations). For each bootstrap replicate, run k-means (k=2) on the 209-gene trajectory matrix. Compute Jaccard stability coefficient of the diverging cluster across bootstraps.
**Expected signal:** Jaccard > 0.8 for the core diverging set (FOS, JUNB, TNF at minimum appear in >90% of replicates).
**Null:** Same bootstrap on shuffled trajectory matrix (permute layer indices per gene).
**Value:** high | **Cost:** low

---

### H-N04: Cross-Model Diverging Gene Test (Geneformer)
**Hypothesis:** The 14 genes that diverge in scGPT embeddings also show anomalously large distance-to-centroid at deep Geneformer layers, establishing model-agnostic outlier geometry.
**Test:** Load Geneformer embeddings for the 14 diverging + 195 converging genes (or the closest available matched set). At the deepest Geneformer layer, rank genes by distance-to-centroid. Compute AUROC for predicting diverging-14 membership from centroid distance rank. Permutation null (1000 shuffles of diverging label).
**Expected signal:** AUROC > 0.7; the 14 genes are preferentially in the top quartile of centroid-distance in Geneformer as well.
**Null:** AUROC ≈ 0.5 (random); permutation p > 0.05.
**Value:** high | **Cost:** medium (requires Geneformer embedding access)

---

### H-N05: Intrinsic Dimension of L11 Converging Manifold
**Hypothesis:** The 195-gene converging cluster at L11 lies on a low-dimensional (d ≤ 3) manifold, while the full 209-gene set at L0 has higher intrinsic dimension.
**Test:** Compute intrinsic dimension (TWO-NN estimator or PCA knee) for (a) all 209 genes at each layer, (b) the 195 converging genes only at each layer. Report how intrinsic dimension evolves from L0 to L11 for each subset.
**Expected signal:** Intrinsic dimension of all-209 drops monotonically; 195-gene subset at L11 approaches d=1–2 (consistent with near-1D collapse); the 14 outliers may maintain higher local dimension.
**Null:** Same ID estimator on random Gaussian embeddings.
**Value:** high | **Cost:** low

---

### H-N06: Cell-Type Expression Variance Predicts Diverging Identity
**Hypothesis:** The 14 diverging genes have higher variance of expression across cell types in typical scRNA-seq data compared to the 195 converging genes.
**Test:** Use a reference scRNA-seq dataset (PBMC 68k or similar public dataset). For each of the 209 genes, compute coefficient of variation (CV) of mean expression across cell types. Mann-Whitney test: CV(diverging-14) vs CV(converging-195).
**Expected signal:** CV of the 14 diverging genes is significantly higher (MW p < 0.05), consistent with scGPT having learned context-dependent representations for these genes that it cannot stably localize in embedding space.
**Null:** MW p > 0.05; or CV difference is explained by expression level alone (control for mean expression rank).
**Value:** high | **Cost:** medium

---

### H-N07: STRING Score Gradient — Interaction Score vs Embedding Distance
**Hypothesis:** Gene pairs with higher STRING interaction score are closer in embedding space at deep layers, establishing a continuous (not binary) STRING-score → embedding-distance relationship.
**Test:** Bin 209-gene STRING pairs by score tier (0.4–0.6, 0.6–0.8, ≥0.8, unconnected). Compute mean pairwise L2 distance per tier at each layer. Kruskal-Wallis across tiers at L10/L11. Spearman rho(score, −distance) for all connected pairs.
**Expected signal:** rho > 0.1 at deep layers; tier ordering preserved (higher score → shorter distance). Trend strengthens from L0 to L11.
**Null:** 500 permutations of STRING scores across gene pairs.
**Value:** medium | **Cost:** low

---

### H-N08: Trajectory Slope Predicts STRING Degree (Continuous)
**Hypothesis:** The continuous trajectory slope (dist-to-centroid change per layer) correlates with STRING degree across all 209 genes, not just at the binary 14-vs-195 level.
**Test:** For each gene, compute trajectory slope (linear regression slope of dist_to_centroid over L0–L11). Spearman rho(slope, STRING_degree) across all 209 genes. 500-permutation null on STRING degree labels.
**Expected signal:** rho < −0.2 (higher STRING degree → steeper convergence / more negative slope); p < 0.05 by permutation.
**Null:** rho ≈ 0; permutation p > 0.05.
**Value:** high | **Cost:** low

---

### H-N09: Persistence of the Outlier Boundary — Rips Persistence at Each Layer
**Hypothesis:** The topological boundary between the 14 diverging genes and 195 converging genes creates a detectable H1 persistent cycle that grows in persistence from L0 to L11.
**Test:** Compute Rips persistence (H0+H1) on 209 genes at all 12 layers using Ripser. Track the longest-lived H1 generator across layers. Test if its persistence (birth-death interval) monotonically increases L0→L11. Compare bottleneck distance between consecutive layer diagrams to detect structural transition.
**Expected signal:** Longest H1 cycle persistence grows from L0 to L11, consistent with the 14 outlier genes pulling further from the main component. Bottleneck distance peaks at L6–L8 (where bifurcation starts).
**Null:** Same computation on random Gaussian embeddings of same shape per layer.
**Value:** high | **Cost:** medium

---

### H-N10: TRRUST Regulatory Load of Diverging vs Converging Genes
**Hypothesis:** The 14 diverging genes include a disproportionate share of TFs with high TRRUST out-degree (number of downstream targets), independently of STRING degree.
**Test:** For TF-hood: Fisher exact test (TF vs non-TF) for diverging-14 vs converging-195. For TFs with TRRUST records: Mann-Whitney test of TRRUST out-degree between diverging and converging TFs.
**Expected signal:** Diverging-14 enriched for TFs (NR4A3, PAX5, NCOA3, KLF6 are all TFs — 4/14 = 29% vs typical TF rate in vocab); TFs in diverging-14 have higher TRRUST out-degree.
**Null:** Fisher exact null; 500 permutations of TF/non-TF labels across the 14/195 split.
**Value:** medium | **Cost:** low

---

### H-N11: Early-Layer Separation — Is the Diverging Set Pre-Separated at L0?
**Hypothesis:** The 14 diverging genes are already geometrically distinct from the 195 converging genes at L0 (before any attention layers process the input), suggesting the divergence is partly driven by input token encoding rather than attention-layer processing.
**Test:** At L0, compute centroid of all 209 genes and separately of the 195 converging genes. Compute distance-to-centroid for both the 14 and 195 sets at L0. Mann-Whitney test for L0 distance difference. Also compute kNN graph (k=10) at L0 and test whether the 14 genes form a separate component already.
**Expected signal:** If L0 separation exists (MW p < 0.05), the 14 genes start farther and stay farther — mechanism is partly input-driven. If L0 separation is absent (MW p > 0.1), the divergence is a learned/architectural effect.
**Null:** MW p > 0.1 at L0 (no early separation) vs growing separation at L5–L11.
**Value:** high | **Cost:** low

---

### H-N12: Anisotropy of Manifold Collapse — Sigmoidal vs Linear PCA Profile
**Hypothesis:** The fraction of variance explained by PC1 (already 76.7% at L11) follows a sigmoidal rather than linear growth profile from L0 to L11, with an inflection point at L6–L8 where the dimensional collapse accelerates.
**Test:** Compute PCA on 209 genes at all 12 layers. Record PC1 variance share per layer. Fit sigmoidal (logistic) and linear models to PC1_share(layer). Report AIC difference and inflection point.
**Expected signal:** Sigmoidal fits better (ΔAIC > 2); inflection at L6–L8, coinciding with the bifurcation onset.
**Null:** Linear fit; same analysis on random Gaussian control (no structure expected).
**Value:** medium | **Cost:** low

---

### H-N13: Intra-Family Spectral Gap Trajectory for High-AUROC vs Low-AUROC Families
**Hypothesis:** Gene families with high clustering AUROC (HLA-I, AP1) have tighter within-family kNN graphs (smaller Fiedler value) at deep layers compared to low-AUROC families (BCL2fam, TNFSF).
**Test:** Build per-family kNN graphs (k=3 or 5) for each of the 9 immune families at each layer. Compute Fiedler value per family per layer. Compare trajectories between high-AUROC (HLA-I, AP1, RUNX) and low-AUROC (BCL2fam, TNFSF) families.
**Expected signal:** High-AUROC families show decreasing Fiedler value (tighter structure) at L8–L11; low-AUROC families show flat or increasing Fiedler values.
**Null:** Same per-family graphs with permuted gene embeddings.
**Value:** medium | **Cost:** medium

---

## Top 3 for Immediate Execution

### Candidate 1 — High-probability discovery: H-N02 + H-N10 (GO Enrichment + TRRUST Regulatory Load of 14 Diverging Genes)
**Rationale:** The 14-gene set is now fixed and biologically coherent. GO enrichment is the most direct, interpretable, and publishable characterization of what makes these genes special. TRRUST TF analysis is <20 additional lines. Together they provide the full biological annotation of the key finding. High expected information yield, very low cost. Package as a single script.
**Code path:** goatools or scipy Fisher exact for GO BP/CC/MF (14 genes vs 209 background); TRRUST out-degree Mann-Whitney for the TF subset; bootstrap-null for GO enrichment.

### Candidate 2 — High-risk/high-reward: H-N11 + H-N01 (Early-Layer Separation + Layer of Bifurcation)
**Rationale:** These two tests together characterize the mechanism of divergence. H-N11 tests whether the separation is input-driven (L0 already separates the 14) or a learned effect (separation emerges mid-network). H-N01 pins down the layer where the slope-sign flip occurs. Together they define the temporal anatomy of the bifurcation — a result with direct implications for how scGPT encodes context-dependent genes. Medium mechanistic risk (the result could be ambiguous), high reward if clean.
**Code path:** H-N11: MW test on L0 centroid distances; kNN graph at L0. H-N01: sliding 3-layer window slope computation, Wilcoxon on flip-layer distribution.

### Candidate 3 — Cheap broad screen: H-N08 + H-N07 + H-N12 (Trajectory Slope vs STRING Degree + STRING Score Gradient + Sigmoidal PCA)
**Rationale:** Three cheap tests (<20 lines each, all reuse existing data structures) that expand the dataset without requiring new data. H-N08 tests whether the 14-vs-195 binary actually reflects a continuous STRING degree → trajectory slope correlation across all 209 genes. H-N07 (carried from iter_0028 H-J, still pending) tests the STRING score gradient. H-N12 characterizes PCA anisotropy. All three together take less time than a single medium-cost experiment and cover three different aspects of the data.
**Code path:** H-N08: linear regression per gene → Spearman rho with STRING degree; H-N07: STRING score tier bins → KW test per layer; H-N12: PCA per layer → logistic fit on PC1_share.
