# Brainstormer Hypothesis Roadmap — iter_0034 → iter_0035+

---

## Retire / Deprioritize

| Direction | Reason | Decision |
|-----------|--------|----------|
| STRING PPI distance gradient (195-gene set) | AUROC=0.494 after OOV correction (iter_0031 H02) | **retire_now** |
| GO BP enrichment in SV2 poles | 0/591 terms significant (iter_0011 H03) | **retire_now** |
| Repression anti-pole (cross-pole TF geometry) | z=-1.41, 0/12 layers (iter_0013 H02) | **retire_now** |
| Monotonic B-cell community emergence (Spearman rho framing) | rho=0.081, p=0.804 (iter_0034 H01) | **retire_now** — keep early-layer finding |
| PC2/PC3 T-cell/myeloid axes at current marker count | All p>0.14 with n≤14 markers (iter_0034 H02) | **rescue_once_with_major_change** — expand to ≥15 markers per lineage, test PC2–PC6 |
| Dorothea hub-uncorrected proximity (population null) | Collapses to AUROC<0.50 with population null (iter_0027 H03) | **retire_now** — OOV-corrected early-layer signal (iter_0032) is the valid claim |
| Chromosomal proximity | AUROC=0.515, non-significant (iter_0025 H03) | **retire_now** (confirmed null, used as control only) |

---

## New Hypothesis Portfolio

### H-A: TRRUST Regulatory Proximity on 195 In-Vocab Genes
**Hypothesis:** TRRUST activation TF-target pairs (n≈100 among 195 in-vocab genes) remain geometrically closer than random pairs, replicating the pre-OOV-correction result (iter_0023 AUROC=0.640) on the clean gene set.
**Test design:** Extract TRRUST activation/repression pairs restricted to 195 in-vocab gene names. Per-layer L2 distance AUROC vs 5000 random pairs (from same gene pool). Activation vs repression comparison.
**Expected signal if true:** Activation AUROC 0.55–0.65 at early-mid layers, decaying to ~0.50 at L11 (consistent with Dorothea decay pattern in iter_0032). Repression AUROC ~0.50 throughout.
**Null:** Random pairs from same 195-gene pool; repression pairs serve as internal control.
**Value:** high | **Cost:** low

---

### H-B: GO Jaccard Similarity Gradient on 195 In-Vocab Genes
**Hypothesis:** Gene pairs sharing more GO terms (higher Jaccard) are geometrically closer in the 195-gene embedding space, replicating iter_0024 rho=0.106 on the clean gene set.
**Test design:** Compute GO CC Jaccard for all ~19000 pairs among 195 genes. Spearman rho(GO Jaccard, L2 distance) per layer. Also test GO BP Jaccard.
**Expected signal if true:** Spearman rho -0.05 to -0.15 for GO CC (negative = more similar terms → shorter distance). Peak at mid layers (L4–L7).
**Null:** Shuffled GO labels; compare to TRRUST/STRING curves.
**Value:** high | **Cost:** low

---

### H-C: Cross-Model B-cell PC1 Validation in Geneformer
**Hypothesis:** Geneformer gene embeddings (residual stream or final layer) also place B-cell marker genes at a distinct pole of the leading PCA axis, confirming that B-cell geometry is a property of the scRNA-seq pretraining objective rather than scGPT-specific architecture.
**Test design:** Load Geneformer residual embeddings for the 195 in-vocab genes (intersection with Geneformer vocabulary). Compute SVD. Test B-cell markers (n≤9) for PC1 pole enrichment using Mann-Whitney AUROC. Compare PC1 explained variance and AUROC to scGPT L11 values.
**Expected signal if true:** B-cell AUROC > 0.20 at Geneformer final layer; PC1 explains >15% variance; B-cell mean PC1 significantly negative.
**Null:** Random 9-gene sets (bootstrap as in iter_0034 H03); shuffled Geneformer embeddings (Gaussian null at same norm).
**Value:** high | **Cost:** medium (requires Geneformer embeddings availability check)

---

### H-D: L2 Community Biological Characterization
**Hypothesis:** The kNN community at L2 that shows OR=16 for B-cell markers (iter_0034 H01) has the same B-cell/non-B-cell split as the L11 community — i.e., B-cell geometry is topologically equivalent across L2 and L11, differing only in community count (3 vs 2).
**Test design:** Run greedy modularity at L2 on 195-gene k=10 kNN. Identify PC1-negative community. Fisher OR for 9 B-cell markers. Compute ARI with L11 2-community partition. Test T-cell markers in PC1-positive community.
**Expected signal if true:** B-cell OR≥10 in one L2 community; ARI with L11 ≥0.50; the "extra" community at L2 splits a portion of the non-B-cell cluster.
**Null:** Null modularity partition (iter_0032 style, 100 random partitions, z-score).
**Value:** medium | **Cost:** low

---

### H-E: Expanded Cell-Type Marker Test on PC2–PC6
**Hypothesis:** With larger marker sets (T-cell n≥20, myeloid n≥12 from TRRUST/curated scRNA-seq), a T-cell or myeloid axis emerges in PC2–PC6 of the L11 embedding.
**Test design:** Curate expanded T-cell (LCK, CD3D, CD3E, CD3G, CD8A, CD8B, GZMB, PRF1, IFNG, IL7R, TCF7, etc.) and myeloid (LYZ, CSF1R, CD68, MRC1, ITGAM, S100A8, S100A9, FCGR3A, etc.) marker lists, filtering to 195-gene vocabulary. Test Mann-Whitney enrichment for each marker set at each PC pole from PC2 to PC6 at L11.
**Expected signal if true:** At least one PC (likely PC2) separates T-cell from myeloid at L11; AUROC≥0.65 with expanded markers.
**Null:** Bootstrap null (1000 random size-matched gene sets); same as iter_0034 H03.
**Value:** medium | **Cost:** low

---

### H-F: Persistent Homology H0/H1 on 195 In-Vocab Genes Across Layers
**Hypothesis:** The number of H0 components (connected components in Rips complex) at an intermediate filtration radius decreases monotonically from L0 to L11, consistent with convergence, while H1 loop count reflects intermediate community structure at mid-layers.
**Test design:** Run Ripser on 195-gene L2 distance matrices per layer (subsample to avoid memory issues). Track H0 count at a fixed radius (50th percentile of pairwise distances) and maximum H1 lifetime. Spearman rho(layer, H0/H1 metrics).
**Expected signal if true:** H0 count decreases monotonically (rho≈-1.0); H1 lifetime peaks at mid layers (L4–L8) then falls.
**Null:** Gaussian null Ripser (same shape tensor, random embeddings); compare H0/H1 against null.
**Value:** medium | **Cost:** medium

---

### H-G: B-cell Gene Conditional Nearest-Neighbor Profile
**Hypothesis:** The set of genes with smallest mean L2 distance to B-cell markers (conditional nearest neighbors) is biologically coherent — enriched for B-cell-biology GO terms, BCR signaling, or B-cell TF targets (PAX5, IRF4 targets) — providing a "B-cell functional neighborhood" in the model's geometry.
**Test design:** At L11, compute mean L2 distance from each non-B-cell gene to the 9 B-cell marker genes. Rank all 186 non-B-cell genes by proximity. Test top-20 and top-40 nearest neighbors for GO BP enrichment and TRRUST TF-target enrichment (especially PAX5/IRF4 targets).
**Expected signal if true:** Top-20 nearest neighbors enriched for B-cell-relevant functions (BCR signaling, B-lymphocyte differentiation, etc.); PAX5/IRF4 target genes over-represented.
**Null:** Nearest-neighbor rank profile for 9 randomly selected gene sets (bootstrap); compare enrichment.
**Value:** high | **Cost:** low

---

### H-H: B-cell Geometry Robustness — Gene Context Shuffle Null
**Hypothesis:** The B-cell PC1 signal in scGPT is NOT reproduced when gene names are randomly shuffled within the input (breaking co-expression structure while preserving marginal distributions), confirming that the geometry arises from learned co-occurrence patterns, not tokenization artifacts.
**Test design:** If scGPT context-level embeddings are available: run scGPT on the same cells but with gene indices shuffled within each cell (permuting which genes map to which expression values). Recompute 195-gene embeddings and test B-cell PC1 AUROC. If not available (static gene embeddings): test whether gene name permutation (reassigning embedding vectors to random gene names) destroys the B-cell signal — if it does, the signal is name-specific, not biology-specific.
**Expected signal if true (model is learning biology):** After gene name permutation, B-cell AUROC drops to ~0.5. After within-cell expression shuffle (if possible), geometry loses B-cell axis.
**Null:** Original ordering vs shuffled; empirical delta-AUROC.
**Value:** high | **Cost:** medium

---

### H-I: B-cell vs Non-B-cell Inter-Community Distance Trajectory
**Hypothesis:** The mean inter-community L2 distance (B-cell community centroid vs non-B-cell community centroid) increases monotonically from L0 to L11, providing a continuous geometric measure of lineage separation that complements the discrete community/OR analysis.
**Test design:** At each layer: assign each of the 195 genes to B-cell (the 9 markers) vs non-B-cell group. Compute centroid distance (B-cell centroid vs non-B-cell centroid). Spearman rho(layer, centroid distance). Compare to within-group distances (cohesion).
**Expected signal if true:** Centroid distance increases monotonically (rho≈1.0); ratio (inter/intra distance) increases monotonically — B-cell genes become increasingly separated from background.
**Null:** 1000 random 9-gene set centroid distances per layer (bootstrap null band).
**Value:** medium | **Cost:** low

---

### H-J: SV2 Co-pole Test for 195 In-Vocab Genes (TRRUST Activation Re-check)
**Hypothesis:** TRRUST activation TF-target pairs remain co-localized in SV2 poles across all 12 layers when restricted to 195 in-vocab genes, replicating the iter_0011/0012 finding on the clean gene set.
**Test design:** Extract SV2 (second left singular vector) at each layer for 195-gene matrix. Assign top/bottom K=49 genes to poles (195 × 0.25 ≈ 49). Test TRRUST activation pairs (195-gene vocab intersection) for co-pole Fisher enrichment vs 500 gene-label shuffles.
**Expected signal if true:** Activation co-pole rate significantly above null at 8–12 layers; consistent with pre-OOV result (12/12 layers in iter_0011).
**Null:** Gene-label shuffle null; repression pairs as internal control.
**Value:** high | **Cost:** low

---

### H-K: Layer 2 Mechanistic Probe — What Changes at L2?
**Hypothesis:** The transition from L1 to L2 is the critical layer for B-cell geometric separation (OR jumps from 1.25/4.31 at L1 to OR=16 at L2). This reflects a specific architectural operation at transformer layer 2, detectable as a change in the distribution of pairwise distances or in SV spectrum.
**Test design:** Compare: (1) pairwise distance distribution at L1 vs L2 (KS test; which gene pairs change most?); (2) SV spectrum (singular value profile) at L1 vs L2; (3) genes with largest L1→L2 displacement — are they enriched for B-cell markers or B-cell-proximal genes?
**Expected signal if true:** B-cell markers show anomalously large displacement from L1 to L2; the B-cell cluster compresses (within-group distance drops) while inter-group distance increases.
**Null:** Layer-to-layer displacement null: compare L1→L2 displacement to mean displacement across all adjacent layer pairs.
**Value:** medium | **Cost:** low

---

### H-L: Intrinsic Dimension Estimates by Community
**Hypothesis:** The B-cell community (9 genes) has lower local intrinsic dimensionality than the non-B-cell community (186 genes) at L11, consistent with B-cell genes occupying a tightly constrained sub-manifold within the larger embedding space.
**Test design:** Two-nearest-neighbor (TwoNN) estimator or participation ratio applied separately to the 9 B-cell markers vs 40 randomly-selected non-B-cell genes (repeated 100 times). Compare mean effective dimensionality. Also test at L2 (where OR is also high) to verify.
**Expected signal if true:** B-cell effective dimensionality < non-B-cell (e.g., PR_Bcell < 3 vs PR_nonBcell > 10 at L11).
**Null:** Bootstrap null: same estimator on 100 random 9-gene sets.
**Value:** medium | **Cost:** low

---

### H-M: TRRUST Hub-Corrected Proximity (Within-TF Null)
**Hypothesis:** TRRUST activation pairs retain geometric proximity when tested with a within-TF null (comparing each TF's targets to each other vs. comparing the TF to non-targets), eliminating hub-gene centrality confound identified in iter_0027.
**Test design:** For each TF with ≥3 named targets in 195-gene vocab: compute mean L2 distance among targets (within-TF). Compare to mean distance from same TF to n random non-targets (matched set size). Per-TF AUROC, aggregate across TFs.
**Expected signal if true:** Majority of TFs show AUROC > 0.55 (targets closer to each other than non-targets), even after hub correction.
**Null:** Random target assignment (shuffle which genes are "targets" for each TF); same within-TF structure.
**Value:** high | **Cost:** medium (requires per-TF analysis)

---

## Top 3 for Immediate Execution

### #1 — High-Probability Discovery: TRRUST + GO Jaccard on 195 In-Vocab Genes (H-A + H-B combined)
Run both TRRUST regulatory proximity (H-A) and GO CC/BP Jaccard gradient (H-B) together in one script, since they share the same per-layer L2 distance matrix infrastructure and both are long-overdue OOV-corrected retests. Together they close the two largest open validation gaps in the cumulative findings list.
- **Expected outcome:** TRRUST activation AUROC 0.55–0.65 at early layers (decay by L11); GO CC Jaccard rho -0.06 to -0.12.
- **Risk:** Low. Prior positive results at n=209 had strong effect sizes; after OOV correction the signal should diminish somewhat but not vanish.

### #2 — High-Risk/High-Reward: Cross-Model B-cell PC1 Validation in Geneformer (H-C)
If Geneformer residual embeddings are available, test whether B-cell markers are displaced to a PC1 pole in Geneformer in the same way as scGPT. This would elevate the B-cell geometry claim from a scGPT-specific curiosity to a property of the scRNA-seq pretraining objective class.
- **Expected outcome:** Geneformer B-cell AUROC 0.15–0.25 (weaker than scGPT if architecture differs; stronger if biology-driven).
- **Risk:** High. Requires confirmed Geneformer embedding artifacts at the residual-stream level. If artifacts not available, fall back to final-layer embeddings.

### #3 — Cheap Broad Screen: B-cell Conditional Nearest-Neighbor Profile + L2→L11 Centroid Trajectory (H-G + H-I)
Two cheap experiments that characterize the "B-cell neighborhood" (what genes live near B-cell markers) and the quantitative separation trajectory (centroid distance per layer). Together they build the biological narrative around the already-confirmed B-cell geometry finding.
- **Expected outcome:** Top-20 B-cell nearest neighbors enriched for B-cell GO terms; centroid distance rho≈1.0 with monotonic increase.
- **Risk:** Low. Both experiments run on existing distance matrices with no new data loading.
