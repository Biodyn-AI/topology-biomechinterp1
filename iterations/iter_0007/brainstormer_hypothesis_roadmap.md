# Brainstormer Hypothesis Roadmap — iter_0007 → iter_0008

---

## Retire / Deprioritize

| Direction | Reason | Disposition |
|-----------|--------|-------------|
| Drift TF/RNAPolII GO enrichment | Fails gene-label shuffle null (empirical p=0.124). Effect too weak relative to random gene assignment. | **retire_now** — drop as primary evidence. Keep in paper as "nominally positive but inconclusive" footnote. |
| Feature-column shuffle as null for L2 statistics | Degenerate null — permutation-invariant. Already documented as a methodological finding. | **retire_now** — replaced by gene-label shuffle in all future tests. |
| Rewiring-based PH null (from prior iterations) | Confirmed degenerate in multiple prior iterations. | **retire_now** — fully closed. |
| Whole-vocabulary drift analysis (4803 genes) | Without gene-label null, all prior drift GO enrichment results are suspect. Computationally expensive and not confirmed. | **deprioritize** — only use 209 named genes with label-shuffle null going forward unless specifically revisiting. |

---

## New Hypothesis Portfolio

### H-A: Per-Layer Gene-Label Shuffle Null for SV1 Axis (Layer 0–11)
**Hypothesis:** The extracellular space enrichment of SV1 top-quartile named genes survives gene-label shuffle null at each of the 12 layers, and specificity is strongest at later layers where SV1/SV2 ratio is highest.
**Test:** Run gene-label shuffle (N=200 reps) at each of 12 layers. Compute empirical p-value for best GO term per layer. Compare specificity (empirical p) to SV1/SV2 ratio trajectory.
**Expected signal if true:** Empirical p decreases monotonically from L0→L11, correlating with rising SV1/SV2 ratio. Strong specificity (p<0.05) at layers 9-11, weak at L0-2.
**Null:** Shuffled gene-label assignment; p ≥ 0.05 at all layers = no specificity anywhere.
**Value:** HIGH | **Cost:** MEDIUM (12 × 200 reps SVD + Fisher, ~10 min)

### H-B: SV2/SV3 Biological Axes at Layer 11
**Hypothesis:** SV2 and SV3 of the layer-11 embedding encode distinct biological axes beyond the SV1 extracellular/cytosol contrast.
**Test:** Project 209 named genes onto SV2 and SV3 at layer 11. GO BP enrichment (top-25% vs bottom-25%) with gene-label shuffle null (N=500). Check for non-redundant GO terms relative to SV1.
**Expected signal if true:** SV2 encodes a second axis (e.g., immune vs. metabolic, or proliferative vs. quiescent); SV3 may encode a third. At least one passes gene-label shuffle null.
**Null:** No significant GO enrichment in SV2 or SV3, or terms overlap SV1 completely.
**Value:** HIGH | **Cost:** LOW (already have infrastructure from H01)

### H-C: SV1 Bottom-Pole Cytosol Label-Shuffle Null Validation
**Hypothesis:** The SV1 bottom-quartile cytosol enrichment (GO:0005829, p=0.010, OR=2.96) at layer 11 is specific to the actual gene-to-embedding assignment.
**Test:** Gene-label shuffle null (N=500) targeting the SV1 bottom-pole enrichment specifically. Empirical p < 0.05 = confirmed.
**Expected signal if true:** Empirical p < 0.05; the extracellular/cytosol axis is confirmed as a true bipolar biological gradient.
**Null:** Bottom-pole cytosol enrichment does not survive null (empirical p ≥ 0.05) = SV1 encodes only a "secreted protein" signal, not a true axis.
**Value:** HIGH | **Cost:** LOW (< 5 min, re-use existing code with minor modification)

### H-D: Transient Layer Deviation Specificity (L3 Mitochondrion)
**Hypothesis:** The layer-3 mitochondrial enrichment (GO:0005739, p=0.000119, OR=inf) is a statistically specific transient encoding, not a random artifact of the small named-gene set.
**Test:** Gene-label shuffle null (N=500) at layer 3, specifically for mitochondrial GO term. Also check spatial clustering: do the top-quartile L3 SV1 genes form a tight spatial cluster, or are they dispersed within the L3 embedding?
**Expected signal if true:** Empirical p < 0.05 for L3 mitochondrial enrichment; top-quartile genes cluster more tightly in L3 than permuted genes.
**Null:** Mitochondrial enrichment at L3 is chance (empirical p ≥ 0.05); the transient is noise.
**Value:** HIGH | **Cost:** LOW

### H-E: Layers 7–8 ER Lumen Transient Specificity
**Hypothesis:** The ER lumen enrichment at layers 7-8 (GO:0005788, OR=inf) represents a real transient encoding of secretory pathway genes in an intermediate layer before final compression.
**Test:** Gene-label shuffle null (N=500) at layers 7 and 8 for ER lumen term. Check overlap between L7-8 ER-enriched genes and known ER-resident or co-translationally sorted proteins.
**Expected signal if true:** Empirical p < 0.05 at L7 or L8; ER-enriched gene set overlaps signal peptide database (e.g., SignalP predictions from UniProt).
**Null:** ER lumen enrichment does not survive null.
**Value:** HIGH | **Cost:** LOW (re-use existing code + optional SignalP annotation lookup)

### H-F: STRING PPI Network Distance vs. SV1 Axis Proximity
**Hypothesis:** Gene pairs closer on the SV1 axis at layer 11 are more connected in the STRING protein-protein interaction network.
**Test:** For all 209×208/2 = 21,736 named-gene pairs: compute |SV1_score_i - SV1_score_j|. Retrieve STRING co-expression scores for these pairs. Spearman correlation between SV1 distance and STRING score (inverted, since closer = more connected). Gene-label shuffle null (N=500) to test specificity.
**Expected signal if true:** Spearman r < -0.1 (negative: smaller SV1 distance = higher STRING score), p < 0.05; survives gene-label null.
**Null:** No significant correlation (|r| < 0.05 or p > 0.05).
**Value:** HIGH | **Cost:** MEDIUM (STRING API or local data retrieval needed)

### H-G: Signal Peptide Presence vs. SV1 Score
**Hypothesis:** Genes with predicted signal peptides (secreted proteins) have significantly higher SV1 scores at layer 11 than genes without signal peptides.
**Test:** For 209 named genes, annotate with UniProt "Signal peptide" presence/absence. Mann-Whitney U test on SV1 scores: signal-peptide-positive vs. negative. Gene-label shuffle null.
**Expected signal if true:** Signal peptide (+) genes have higher mean SV1 score; U-test p < 0.01 survives gene-label null.
**Null:** No SV1 score difference between SP(+) and SP(-) genes.
**Value:** HIGH | **Cost:** LOW (UniProt annotation lookup for 209 genes is fast)

### H-H: Effective Rank Per Biological Class (Extracellular vs. Cytosol)
**Hypothesis:** The subspace spanned by extracellular space genes in layer 11 has lower effective rank than the subspace of cytosolic genes, indicating tighter geometric clustering.
**Test:** Split 209 named genes into extracellular (GO:0005615 annotated) vs. cytosolic (GO:0005829 annotated). Compute SVD of each submatrix. Compare effective rank, cumulative variance at SV1. Label-shuffle null for the ratio.
**Expected signal if true:** Extracellular gene subspace has markedly lower effective rank (tighter cluster) than cytosolic subspace.
**Null:** Similar effective ranks for both subsets.
**Value:** MEDIUM | **Cost:** LOW

### H-I: SV2 Layer-Wise Trajectory (Is There a Second Persistent Axis?)
**Hypothesis:** SV2 also encodes a biologically specific axis across multiple layers, distinct from SV1, reflecting a secondary organizational principle in scGPT representations.
**Test:** Run full layer-wise analysis (layers 0-11) for SV2: GO enrichment at each layer, SV2/SV3 ratio trajectory. Gene-label shuffle null (N=200) at selected layers (L0, L5, L11).
**Expected signal if true:** SV2 encodes a consistent second biological axis (e.g., metabolic vs. signaling) in ≥6 layers; GO enrichment survives null.
**Null:** SV2 top-hit GO term changes randomly across layers with no persistence.
**Value:** HIGH | **Cost:** MEDIUM

### H-J: Persistent Homology on 209 Named Genes at Layer 11 (Fresh Attempt)
**Hypothesis:** The 209 named gene embeddings at layer 11 form a topologically non-trivial point cloud (non-zero H1 Betti number, i.e., loops) that is not present in random gene-label shuffles.
**Test:** Compute Vietoris-Rips persistent homology (H0, H1) on the 209-gene submatrix at layer 11 (Euclidean distance). Compare observed persistence diagram to 200 gene-label shuffles. Metric: total H1 persistence mass.
**Expected signal if true:** Non-trivial H1 features with longer lifetime than 95% of shuffles.
**Null:** H1 persistence not distinguishable from shuffled gene assignments.
**Value:** HIGH | **Cost:** MEDIUM (209-gene PH is computationally cheap; prior attempts used full vocab which was expensive)

### H-K: Intrinsic Dimensionality (TWO-NN) of Named Gene Subspace vs. Full Vocabulary
**Hypothesis:** The 209 named genes occupy a lower-intrinsic-dimensional manifold than the full gene vocabulary at layer 11, reflecting their specialized functional coverage.
**Test:** Estimate TWO-NN (two-nearest-neighbor) intrinsic dimension for (a) 209 named gene embeddings and (b) 500 randomly sampled full-vocabulary genes at layer 11. Bootstrap confidence intervals. Gene-label shuffle null for named gene ID (tests whether annotation depth, not identity, drives result).
**Expected signal if true:** Named gene ID < full-vocab intrinsic dimension; confidence intervals non-overlapping.
**Null:** Similar intrinsic dimension.
**Value:** MEDIUM | **Cost:** LOW

### H-L: Layer-Wise Correlation Between SV1/SV2 Ratio and GO Enrichment p-Value
**Hypothesis:** The SV1/SV2 ratio is a predictor of the biological specificity (GO enrichment p-value) of the SV1 axis across layers — higher dominance predicts stronger enrichment.
**Test:** Linear/Spearman correlation across 12 layers: SV1/SV2 ratio vs. -log10(best GO p-value). If extended with per-layer nulls (H-A), also correlate SV1/SV2 ratio vs. empirical p-value.
**Expected signal if true:** Spearman r < -0.5 (higher ratio = stronger enrichment), p < 0.05.
**Null:** No significant correlation.
**Value:** MEDIUM | **Cost:** LOW (data already exists from H02)

### H-M: Annotation Completeness Confounder Test
**Hypothesis:** The SV1 extracellular enrichment is partly driven by annotation completeness — secreted proteins may have more GO BP annotations per gene, biasing the Fisher exact test.
**Test:** For 209 named genes, compute annotation count per gene. Test whether SV1 top-quartile genes have significantly more GO BP annotations than bottom-quartile. If confounded, re-run GO enrichment stratified by annotation depth, or use annotation-count-matched controls.
**Expected signal if true (confounder present):** Top-quartile genes have significantly more annotations; enrichment is inflated.
**Expected signal if confound is absent:** No annotation count difference; current enrichment results are unconfounded.
**Value:** MEDIUM | **Cost:** LOW (annotation counting is fast)

### H-N: Cross-Layer SV1 Score Consistency (Individual Gene Tracking)
**Hypothesis:** Individual genes maintain their relative SV1 rank across layers — genes at the top pole at layer 11 were already at the top at layer 0.
**Test:** For each of 209 named genes, extract SV1 projection at all 12 layers. Compute Spearman rank correlation between layer-0 and layer-11 SV1 ranks. Also compute rank correlation between all consecutive-layer pairs.
**Expected signal if true:** High rank stability (Spearman r > 0.7 between L0 and L11); the biological axis is preserved through transformation.
**Null:** Low rank stability (r < 0.3) = axis identity changes through layers.
**Value:** HIGH | **Cost:** LOW (pure computation on existing data)

### H-O: SV1 Score vs. Protein Hydrophobicity / MW
**Hypothesis:** SV1 score correlates with biophysical properties associated with secretion: high hydrophobicity (GRAVY score) and molecular weight, consistent with signal-peptide-bearing extracellular proteins.
**Test:** For 209 named genes, retrieve UniProt MW and compute GRAVY index from canonical sequence. Spearman correlation with layer-11 SV1 projection. Gene-label shuffle null.
**Expected signal if true:** Positive Spearman r with hydrophobicity and/or MW (secreted proteins tend to be larger and more hydrophobic); p < 0.05 survives null.
**Null:** No significant correlation.
**Value:** MEDIUM | **Cost:** LOW (UniProt data retrieval)

### H-P: CKA Cross-Model Alignment vs. SVD Trajectory Correlation
**Hypothesis:** The peak CKA cross-model alignment between scGPT layers (from iter_0004) coincides with the layer where SV1/SV2 ratio crosses 5x, suggesting that model convergence and axis emergence are the same phenomenon.
**Test:** Re-load iter_0004 CKA alignment scores across layers. Correlate with SV1/SV2 ratio trajectory from H02. Compute Pearson r. Also check if both curves show inflection near layer 2.
**Expected signal if true:** r > 0.7 between CKA alignment and SV1/SV2 ratio trajectories; both peak around layer 9-11.
**Null:** Trajectories uncorrelated (r < 0.3).
**Value:** HIGH | **Cost:** LOW (data already exists from iter_0004)

---

## Top 3 for Immediate Execution

### Execution Candidate 1: HIGH-PROBABILITY DISCOVERY
**H-C (SV1 Bottom-Pole Label-Shuffle Null) + H-G (Signal Peptide Annotation) combined as one screen**

Rationale: Both are <1 hour of work using existing code. H-C validates the cytosol pole to make the extracellular/cytosol axis fully confirmed. H-G directly explains *why* the SV1 axis encodes extracellular biology — signal peptide presence is the mechanistic link. Together they convert the current "suggestive" finding into a confirmed, mechanistically grounded claim.

Test design:
1. H-C: gene-label shuffle (N=500) on SV1 bottom-quartile vs. bottom-of-complement; compute empirical p for cytosol term
2. H-G: annotate 209 genes with UniProt signal peptide status; Mann-Whitney on SV1 scores; gene-label shuffle null

Expected result: H-C empirical p < 0.05 (the axis is bipolar); H-G Mann-Whitney p < 0.01 (SP+ genes cluster at SV1 top). This would be the clearest mechanistic result in the project so far.

---

### Execution Candidate 2: HIGH-RISK / HIGH-REWARD
**H-D + H-E: Transient Layer Deviation Specificity (L3 Mitochondrion + L7-8 ER Lumen)**

Rationale: If the layer-3 mitochondrion and layer-7/8 ER lumen transients survive label-shuffle nulls, this is a major finding — it would mean scGPT encodes distinct organelle-compartment biology at specific processing layers, with the secretory pathway as the final attractor. This would be a new interpretability result for transformer biology. If they fail, the transients are just noise and we close them cleanly.

Test design:
1. Gene-label shuffle (N=500) at L3: test GO:0005739 (mitochondrion) empirical p
2. Gene-label shuffle (N=500) at L7 and L8: test GO:0005788 (ER lumen) empirical p
3. For any confirmed transient: check if top-quartile genes overlap with known organelle residents (MitoCarta for mitochondria, HPA for ER)
4. Spatial test: compare within-group pairwise distance for top-quartile transient genes vs. gene-label shuffle

Expected result: At least one transient (L3 or L7/L8) survives null, providing evidence for layer-specific subcellular compartment encoding.

---

### Execution Candidate 3: CHEAP BROAD-SCREEN
**H-B (SV2/SV3 axes) + H-N (cross-layer SV1 rank stability) + H-L (ratio vs enrichment correlation)**

Rationale: All three use existing data and code. H-B extends the SVD analysis to SV2/SV3 — if a second biological axis exists, it's a major new finding. H-N tests rank stability, which either validates the layer-trajectory story or reveals it's more fluid than reported. H-L tests whether the SV1/SV2 ratio mechanistically explains the enrichment trajectory — a clean quantitative claim.

Test design:
1. H-B: GO enrichment + gene-label shuffle (N=500) for SV2 and SV3 at layer 11
2. H-N: Compute Spearman(L0 SV1 rank, L11 SV1 rank) for 209 genes; all consecutive-layer pairs
3. H-L: Spearman(SV1/SV2 ratio, -log10(GO p)) across 12 layers using existing H02 data
Total cost: LOW (< 30 min, mostly pure computation on already-computed embeddings)
