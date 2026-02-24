# Brainstormer Hypothesis Roadmap — iter_0017

---

## Retire / Deprioritize

| Hypothesis direction | Reason | Decision |
|---|---|---|
| TRRUST signed regulation on SV3/SV4 | 0/12 and 2/12 layers sig; SV2-specific property confirmed | retire_now |
| GO BP enrichment in SV2 poles | 0/591 terms at any layer (iter_0011) | retire_now |
| Repression anti-pole geometry | 0/12 layers sig (iter_0013) | retire_now |
| Bootstrap CI act/rep differential | Underpowered with n=64 repression; per-axis emp_p supersedes | retire_now |
| SV5 PPI co-pole | 6/12 layers non-monotonically; weak signal, low marginal value | deprioritize |
| Additional SV4 GO profiling | SV4 biological identity established (iter_0016 H03 + iter_0017 H02) | retire_now |
| Additional quintile gradient axes (SV5+) | Three-axis story complete; diminishing returns | deprioritize |

---

## New Hypothesis Portfolio

### H-F: Geneformer SV2/SV3 PPI Co-pole
**Hypothesis:** Geneformer residual-stream gene embeddings also show STRING PPI co-pole enrichment on their leading SVD axes, confirming that confidence-graded PPI geometry is a general property of trained biological LMs, not scGPT-specific.
**Test:** Extract Geneformer layer embeddings for the 209 named genes. Compute SVD. Run co-pole enrichment test (STRING score ≥ 0.4, 3092 pairs, K=52, N=300 shuffles) on SV1–SV4 across all Geneformer layers.
**Expected signal:** ≥1 Geneformer SVD axis shows 12/12 layers sig, z>3.
**Null:** Feature-shuffle. Also check that Geneformer SV1 (random) shows no signal.
**Value:** high | **Cost:** medium
**Priority:** high-probability discovery

---

### H-G: Geometric Distance as Continuous STRING Score Predictor
**Hypothesis:** Euclidean distance in the SV2–SV3–SV4 joint embedding space is a significant predictor of STRING interaction score as a continuous variable (R²>0, Spearman rho significant after FDR correction).
**Test:** For all 3092 pairs, compute mean per-layer Euclidean distance in (SV2, SV3, SV4) coordinates. Regress against STRING score. Report R², Spearman rho, p-value. Compare to SV1 baseline.
**Expected signal:** Spearman rho ≥ 0.10 (negative: closer = higher score); SV1 baseline near zero.
**Null:** SV1 distance regression, feature-shuffle distance.
**Value:** high | **Cost:** low

---

### H-H: Out-of-Sample PPI Prediction (Precision@k)
**Hypothesis:** For each gene in the 209-gene set, its top-k geometric neighbors (by SV2–SV4 Euclidean distance) are enriched for known STRING interaction partners at above-chance rates, suggesting the embedding generalizes beyond the pairs explicitly tested.
**Test:** For each gene g, rank all other genes by distance in SV2–SV4 space. Compute precision@10 and precision@20 using STRING score ≥ 0.7 as ground truth. Compare to precision@k from random ranking.
**Expected signal:** Mean precision@10 > random (p<0.05 by permutation).
**Null:** Random gene ranking, SV1 distance ranking.
**Value:** high | **Cost:** low

---

### H-I: Random Initialization Control
**Hypothesis:** A randomly initialized (untrained) scGPT model does NOT show STRING PPI co-pole enrichment, establishing that the geometric structure requires training on biological data.
**Test:** Load scGPT architecture with random weights. Extract gene embeddings for the same 209 genes. Run SV2 co-pole test (STRING ≥ 0.7, N=300 shuffles) across all 12 layers. Compare z-scores to trained model.
**Expected signal:** Random model: z ≈ 0, 0/12 layers sig. Trained model: z ≈ 10, 12/12 layers.
**Null:** This IS the null control.
**Value:** high | **Cost:** low

---

### H-J: Cross-Seed Stability of Multi-Axis Structure
**Hypothesis:** The three-axis (SV2/SV3/SV4) PPI geometry is stable across data sampling seeds: genes that are axis-dominant in seed42 are also axis-dominant in seed43 and seed44.
**Test:** For seeds 42, 43, 44: compute axis-dominant assignments (iter_0017 method). Measure Cohen's kappa for axis assignment across seed pairs. Compute STRING co-pole z-scores per seed.
**Expected signal:** Kappa > 0.5 across seed pairs; z-score profiles highly correlated.
**Null:** Random axis reassignment.
**Value:** medium | **Cost:** low

---

### H-K: scGPT Immune Dataset Replication
**Hypothesis:** The SV2/SV3 PPI co-pole signal replicates in scGPT immune (not lung) embeddings, confirming tissue-independence.
**Test:** Extract scGPT immune residual embeddings for the 209 named genes. Run SVD. Test STRING co-pole (score ≥ 0.4, K=52, N=300 shuffles) on SV1–SV4.
**Expected signal:** SV2 and/or SV3 12/12 layers sig in immune embeddings.
**Null:** Feature-shuffle, SV1.
**Value:** medium | **Cost:** low

---

### H-L: Attention Head Alignment with SV2 Geometry
**Hypothesis:** Specific scGPT attention heads preferentially attend to STRING interaction partners (measured by attention weight between co-pole gene pairs vs. random pairs), and the heads with strongest PPI alignment map to layers where SV2 z-scores are highest (L0–L1, L11).
**Test:** Extract attention weight matrices for all 12 layers × n_heads. For each head, compute mean attention weight for STRING pairs (score ≥ 0.7) vs. random pairs. Rank heads by PPI specificity. Correlate with layer-wise SV2 z-scores.
**Expected signal:** ≥1 head per layer shows significant PPI attention enrichment; strongest heads in layers 0–1 and 11.
**Null:** Attention weights for randomly permuted gene identity labels.
**Value:** high | **Cost:** medium

---

### H-M: STRING Complex Co-clustering (Beyond Pairwise)
**Hypothesis:** Known protein complexes (STRING clusters or CORUM data) are more geometrically compact in SV2–SV4 joint space than equivalently-sized random gene sets, extending the pairwise PPI result to higher-order interaction structure.
**Test:** Load CORUM complexes or STRING MCL clusters with ≥4 named genes. For each complex, compute mean pairwise Euclidean distance in SV2–SV4 space. Compare to 1000 size-matched random gene sets from the 209-gene pool.
**Expected signal:** Known complexes are significantly more compact (lower mean distance) than random; z < -2 for ≥50% of complexes.
**Null:** Random same-size gene groups from 209-gene pool.
**Value:** high | **Cost:** medium

---

### H-N: SV2 Projection as TF Regulatory Activity Score
**Hypothesis:** SV2 projection of a TF gene predicts the fraction of its known TRRUST targets that are activation (vs repression) targets: TFs with high positive SV2 projection should have more activating targets, providing a mechanistic read-out.
**Test:** For each TF with ≥3 TRRUST targets, compute SV2 projection at layer 11. Correlate with that TF's activation fraction (n_act / (n_act + n_rep)). Report Spearman rho.
**Expected signal:** rho > 0.3, p<0.05.
**Null:** Shuffle TF-to-SV2 projection assignment.
**Value:** medium | **Cost:** low

---

### H-O: PH (Persistent Homology) on PPI-Stratified Gene Sets
**Hypothesis:** High-confidence STRING interaction genes (top-20% score) form topologically more complex neighborhoods (higher H1 persistence) than low-confidence STRING genes, revealing that PPI structure is encoded in loop topology not just polarity.
**Test:** Partition the 209 genes into high-PPI (genes with ≥1 partner at score ≥ 0.8) and low-PPI sets. Compute H1 persistence (ripser, K=50 per set, all 12 layers). Compare mean top persistence against feature-shuffle null.
**Expected signal:** High-PPI set has significantly higher H1 persistence than low-PPI set.
**Null:** Feature-shuffle, comparison to low-PPI set.
**Value:** medium | **Cost:** medium

---

### H-P: Geneformer SV2 Signed Regulatory Geometry
**Hypothesis:** If Geneformer also shows PPI co-pole geometry (H-F positive), its leading SVD axis should also encode TRRUST activation preferentially over repression, confirming regulatory sign geometry as a general BioLM property.
**Test:** (Contingent on H-F positive.) On Geneformer's top PPI axis: run TRRUST activation (n=116) and repression (n=64) co-pole test across all Geneformer layers. Report n_layers_sig for each.
**Expected signal:** Activation ≥ 8/n_layers sig; repression < 3/n_layers sig.
**Null:** Label shuffle.
**Value:** high | **Cost:** low (contingent on H-F)

---

### H-Q: Layer-Resolved Confidence Gradient Emergence
**Hypothesis:** The Q5 spike in the quintile gradient (the dominant feature of rho ≈ 0.90 on all three axes) is absent in early layers and emerges progressively, with a specific layer of onset that differs across axes, suggesting layer-specific biological organization.
**Test:** Run the quintile gradient test per-layer (not averaged across layers) for SV2, SV3, SV4. Report rho and Q5 z-score at each of 12 layers. Find layer of emergence (first layer with Q5 z > 3).
**Expected signal:** Layer-of-emergence differs across axes (e.g., SV2 early, SV3 later); matches known SV2 variance peak at layer 8.
**Null:** No layer-specific pattern (rho stable across layers).
**Value:** medium | **Cost:** low

---

### H-R: GO Term Prediction from Axis Projection (Decoding)
**Hypothesis:** A simple linear classifier (logistic regression) trained on SV2–SV4 projections can predict gene membership in specific GO categories (e.g., innate immune, membrane, apoptosis) above chance, quantifying how much biological information is recoverable from the geometric coordinates.
**Test:** For each of the 3 axis-dominant GO categories (iter_0017 H-B top terms), train leave-one-out logistic regression on (SV2, SV3, SV4) projections at layer 11. Report AUROC. Compare to SV1 baseline and shuffled-label baseline.
**Expected signal:** AUROC ≥ 0.70 for ≥2 GO categories; SV1 baseline ≈ 0.50.
**Null:** SV1 projection, label shuffle.
**Value:** medium | **Cost:** low

---

## Top 3 for Immediate Execution

### 1. High-Probability Discovery Candidate: H-F (Geneformer Cross-Model Validation)
This is the single most important missing experiment. Every iteration since iter_0010 has listed it as "highest priority." The paper's core claim—that trained biological LMs encode PPI geometry—cannot be credibly published with only one model. The experiment is structurally identical to iter_0012 H03 and should execute cleanly if Geneformer embeddings are accessible. Risk: Geneformer residual embeddings may need special extraction. Mitigation: test embedding accessibility first.

### 2. High-Risk/High-Reward Candidate: H-H (Out-of-Sample PPI Prediction)
Precision@k is a publishable quantitative benchmark: "can the model's geometry predict which proteins physically interact?" This reframes the finding from "PPI geometry is present" to "PPI geometry is predictive." If precision@10 is significantly above chance, this is the strongest single result in the paper. Risk: effect size may be modest given the 209-gene pool's density. Mitigation: use STRING score ≥ 0.5 threshold for broader coverage.

### 3. Cheap Broad-Screen Candidate: H-G + H-I bundle (Continuous Prediction + Random Init Control)
Both are one-script additions. H-G (Euclidean distance → STRING score regression) directly quantifies the gradient story as a continuous prediction. H-I (random init control) is the cleanest possible demonstration that training matters. Together these two cheap experiments would substantially strengthen the paper in one executor run.

---

## Portfolio Coverage Summary

| Family | Hypotheses |
|---|---|
| Cross-model alignment | H-F, H-P |
| Quantitative/predictive geometry | H-G, H-H, H-R |
| Controls/confounds | H-I, H-J |
| Cross-dataset replication | H-K |
| Mechanistic (attention) | H-L |
| Higher-order topology | H-M, H-O |
| Biological anchoring | H-N |
| Layer-resolved dynamics | H-Q |
