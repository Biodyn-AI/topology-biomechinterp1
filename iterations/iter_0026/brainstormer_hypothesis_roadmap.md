# Brainstormer Hypothesis Roadmap — iter_0026

**Date:** 2026-02-22

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| Chromosomal proximity confound | **RETIRED** (iter_0025) | AUROC 0.515, ruled out. |
| GO BP vs GO CC redundancy | **RETIRED** (iter_0025) | Question answered — partially independent. |
| Further univariate STRING/Dorothea sweeps | **DEPRIORITIZE** | 12+ iterations of characterization; diminishing returns on new univariate results. |
| TRRUST solo activation AUROC | **DEPRIORITIZE** | Characterized at L5=0.659. No further solo TRRUST unless testing new sub-hypothesis. |
| Dorothea MOR direction test via Dorothea file | **DEPRIORITIZE** | Blocked — MOR column absent. Pivot to TRRUST A/R split. |
| BCL2fam / TNFSF family analysis | **LOW PRIORITY** | Both negative in H01 (AUROC 0.48–0.54). Likely reflects diverse regulatory context in training data. Not worth further investment without a clear biological rescue hypothesis. |
| KLF family analysis (without outlier removal) | **RESCUE ONCE** | AUROC=0.31 (anti-clustered). Remove KLF14 from KLF group and retest. If KLF2/4/6 cluster well without KLF14, this confirms KLF14 is a tissue-context outlier — one useful data point for the paper. |

---

## New Hypothesis Portfolio

### H-N1: TRRUST Activation vs. Repression Geometric Polarity
**Hypothesis:** TF→target pairs annotated as activating (A) in TRRUST are geometrically closer in scGPT embeddings than repressing (R) pairs, encoding regulatory direction geometrically.
**Test:** Load TRRUST file, filter to pairs in 209-gene vocab. Split by TF_type=A vs TF_type=R. Compute pairwise L2 distances for each class. AUROC(−d_activator, −d_repressor) at each of 12 layers. Also test AUROC(each class vs null random pairs).
**Expected signal:** AUROC(activator < repressor) > 0.55 at L5–L8; activators closer than repressors, reflecting co-expression vs anti-correlation in training data.
**Null/control:** Shuffle A/R labels within TRRUST set.
**Value:** HIGH — directly resolves the directional polarity question blocked in H02. Uses available data.
**Cost:** LOW — TRRUST already available, same pipeline as Dorothea.

---

### H-N2: PC1 Biological Identity at Layer 11 (Dominant Axis Decoder)
**Hypothesis:** The dominant principal component (PC1) at layer 11 — which accounts for ~50%+ of all variance after dimensional collapse — encodes a biologically interpretable axis (e.g., immune vs. non-immune, TF vs. non-TF, hub vs. peripheral).
**Test:** At L11, compute PCA on 208 gene embeddings. Extract PC1 loadings per gene. Correlate loadings with binary gene attributes: is_TF (known TF from TRRUST/Dorothea), is_HLA, is_hub (STRING degree top quartile), GO_CC_nucleus, GO_CC_membrane. Compute Spearman ρ between PC1 loading and each attribute. Test significance (permutation).
**Expected signal:** At least one biological attribute shows |ρ| > 0.25 with PC1 loadings (p < 0.01 by permutation). Prediction: is_TF or is_hub has the strongest correlation.
**Null/control:** Permute gene labels within each attribute; also test PC2, PC3 to check specificity to PC1.
**Value:** HIGH — this is the mechanistic keystone for the dimensional collapse story. If PC1 encodes TF identity or network centrality, it directly connects H03 to the regulatory encoding results.
**Cost:** LOW — pure computation on existing embeddings + existing gene metadata. No new data needed.

---

### H-N3: Per-Family Layer Curves (HLA-I vs. AP1 Peak Layer Test)
**Hypothesis:** Structural gene families (HLA-I, sequence-similar paralogs) reach peak family AUROC at earlier layers than regulatory families (AP1/RUNX, co-regulated TF modules), reflecting the geometric hierarchy where structural identity precedes functional identity.
**Test:** For each of the 9 immune gene families, compute AUROC at all 12 layers (already have global layer curves but not per-family). Find the peak layer for each family. Test: structural families (HLA-I, HLA-II) peak at earlier layer than regulatory families (AP1, RUNX, IL2_path)?
**Expected signal:** HLA-I peaks at L5 or earlier; AP1/RUNX peak at L7–L8. Layer difference ≥ 2.
**Null/control:** Permuted family labels; bootstrap AUROC confidence intervals.
**Value:** HIGH — if confirmed, this is a clean mechanistic mapping from family type to layer. Directly extends H01 into the layer-hierarchy story.
**Cost:** LOW — data already collected (h01_immune_family_auroc.json has layer-level global but not per-family layer curves; need small re-run).

---

### H-N4: PR Collapse Rate per Biological Gene Subset
**Hypothesis:** Biologically coherent gene subsets (HLA-I, AP1, GO CC nucleus-annotated genes) undergo faster dimensional collapse across layers than random same-size gene subsets, because biological constraint compresses the representational space more aggressively.
**Test:** For each of 5 gene subsets (HLA-I, AP1+RUNX, GO_CC_nucleus, STRING_top_degree_quartile, random), compute PR at each of 12 layers using subset-specific PCA. Compute collapse rate = (PR_layer0 − PR_layer11) / PR_layer0. Compare subset collapse rates to random baseline.
**Expected signal:** HLA-I and AP1+RUNX collapse faster (collapse rate > 0.85) than random subsets (collapse rate ~0.7–0.75).
**Null/control:** Same-size random gene subsets (10 bootstrap samples).
**Value:** MEDIUM-HIGH — connects dimensional collapse (H03) to biological structure (H01), creates a unified narrative.
**Cost:** LOW — re-use existing embeddings; compute subset PCA at 12 layers.

---

### H-N5: STRING Hub vs. Peripheral Embedding Density
**Hypothesis:** High-degree STRING network hubs (top 25% by degree among 209 genes) occupy denser local neighborhoods in embedding space than low-degree peripheral genes.
**Test:** Compute STRING degree for 209 named genes from existing cache. Split into quartiles. For each gene, compute mean L2 distance to k=10 nearest neighbors at each layer. Mann-Whitney test hub (Q4) vs peripheral (Q1) local density per layer.
**Expected signal:** Hub genes have lower mean k-NN distance (denser) at L7–L11; Mann-Whitney z > 2.
**Null/control:** Random hub assignment (shuffle degree labels).
**Value:** MEDIUM — links network topology to embedding geometry; cheap characterization.
**Cost:** LOW — reuses STRING cache + existing embeddings.

---

### H-N6: STRING Score Threshold Sensitivity Sweep
**Hypothesis:** The STRING Spearman ρ signal degrades monotonically as the interaction score threshold increases from 0.4 to 0.9, indicating low-confidence interactions carry genuine geometric signal.
**Test:** Recompute STRING Spearman ρ at L8 for 5 thresholds: 0.4, 0.5, 0.6, 0.7, 0.9. Plot ρ vs threshold. Spearman rank correlation of ρ-vs-threshold as test statistic.
**Expected signal:** Monotonic decrease (Spearman rho negative, p < 0.05). ρ drops from ~0.15 at threshold 0.4 to ~0.05 at threshold 0.9.
**Null/control:** Same analysis with shuffled STRING scores.
**Value:** MEDIUM — validates signal quality and methods robustness; useful for methods section.
**Cost:** LOW — reuses STRING cache; simple threshold filter on existing pipeline.

---

### H-N7: Oncogene / Tumor Suppressor Polarity
**Hypothesis:** Oncogenes and tumor suppressors from COSMIC/CGC among the 209 named genes are linearly separable in scGPT embedding space at late layers.
**Test:** Pull COSMIC Gene Census list; label 209 named genes as oncogene / tumor suppressor / neither. Require ≥10 genes per class. Compute LDA classification AUROC at each layer (leave-one-out). Test peak AUROC layer.
**Expected signal:** AUROC > 0.65 at L8–L11; better than early layers.
**Null/control:** Random class assignment of same size.
**Value:** MEDIUM — clinically compelling narrative; depends on having enough labeled genes.
**Cost:** LOW-MEDIUM — COSMIC data available; need label construction + LDA.

---

### H-N8: Persistent Homology of GO CC Subcellular Compartment Clusters
**Hypothesis:** GO CC-defined subcellular compartments (nucleus, cytoplasm, membrane) form topologically distinct 0-dimensional barcodes in scGPT embedding space, with nuclear genes showing longer-lived clusters than cytoplasmic.
**Test:** For each major GO CC class among 209 named genes (require ≥15 genes), extract embedding submanifold per layer. Compute 0D PH (VR filtration). Measure mean and max persistence per class per layer. Compare across classes.
**Expected signal:** Nuclear genes show mean persistence 2× longer than cytoplasm genes at L5 (peak GO CC layer).
**Null/control:** Random same-size gene subsets.
**Value:** HIGH — topological validation of the pairwise distance signal; topologically-grounded finding if confirmed.
**Cost:** MEDIUM — requires ripser; subsets are small (≤50 genes), so PH is fast.

---

### H-N9: Geneformer Cross-Model Joint Predictor Replication
**Hypothesis:** The multi-predictor joint model (STRING + Dorothea + GO CC) replicates in Geneformer embeddings with comparable total partial R².
**Test:** Load Geneformer gene embeddings for the same 209 named genes (check subproject_38 outputs). Compute identical partial Spearman R² joint model across Geneformer layers.
**Expected signal:** Total partial R² > 0.01; similar predictor ranking (STRING dominant or tied with GO CC).
**Null/control:** Shuffled gene identities.
**Value:** HIGH — confirms generalizability; the most impactful validation experiment available.
**Cost:** MEDIUM — requires Geneformer embedding access; check if available in subproject_38.

---

### H-N10: Disease Gene Manifold Clustering (DisGeNET)
**Hypothesis:** Genes associated with the same disease in DisGeNET cluster in scGPT embedding space more than gene pairs from different diseases.
**Test:** Fetch DisGeNET gene-disease associations for 209 named genes via API or local file. Classify pairs as same-disease vs. different-disease (requiring ≥5 genes per disease). AUROC at each layer.
**Expected signal:** AUROC > 0.58 at mid-late layers.
**Null/control:** Random disease assignment.
**Value:** HIGH — clinical relevance angle; strong narrative fit.
**Cost:** MEDIUM — requires DisGeNET data download + new pairwise test.

---

### H-N11: KLF Outlier Rescue — Remove KLF14 and Retest
**Hypothesis:** KLF14 is a tissue-context outlier (liver/adipose) causing the KLF family AUROC to collapse; removing it reveals proper clustering of KLF2/4/6 (ubiquitous TFs with similar regulatory targets).
**Test:** Remove KLF14 from the KLF family group. Recompute AUROC for {KLF2, KLF4, KLF6} as a 3-gene cluster. Check if AUROC rises above 0.65 at L7–L8.
**Expected signal:** AUROC(KLF without KLF14) > 0.65. This confirms KLF14 is the outlier.
**Null/control:** Compare to original KLF AUROC=0.31.
**Value:** MEDIUM — confirms interpretation of H01 result; single computation.
**Cost:** LOW — trivial modification to existing script.

---

### H-N12: Elbow Layer Alignment — PR Collapse Rate vs. GO/Dorothea Signal Peaks
**Hypothesis:** The layer at which the PR collapse rate is steepest (the "geometric elbow") corresponds to the layer where Dorothea/GO BP signals peak, linking geometric compression to information encoding milestones.
**Test:** Compute per-layer PR derivatives (ΔPR_l = PR_l − PR_{l+1}). Find the elbow (layer of maximum ΔPR). Compare to peak Dorothea AUROC layer and peak GO BP layer from existing data. Test if they co-occur within ±1 layer.
**Expected signal:** Peak ΔPR at L7–L8, matching Dorothea peak (L7) and GO BP peak (L6–L7).
**Null/control:** Random permutation of layer ordering.
**Value:** MEDIUM-HIGH — provides mechanistic narrative linking compression to encoding; uses entirely existing data. Zero new computation required except derivative calculation.
**Cost:** LOW — entirely computed from existing h03 and prior AUROC data.

---

### H-N13: Expanded TF Family Screen (ETS, STAT, NFkB, Nuclear Receptor, Homeodomain)
**Hypothesis:** The per-family clustering AUROC signal extends beyond immune gene families to broader TF superfamilies present in the 209-gene vocab.
**Test:** Define TF superfamily groups from the 209-gene vocab: ETS (ETS1, ETS2, GABPA, FLI1), STAT (STAT1/2/3/5/6), NFkB (RELA, RELB, NFKB1/2), nuclear receptors (NR3C1, NR4A1/2, RARA), homeodomain (HOXA/B members if present). Compute AUROC within-family vs cross-family at each layer.
**Expected signal:** STAT and NFkB show AUROC > 0.75 at L7–L8 (tight regulatory co-regulation). Nuclear receptors may be lower (diverse ligand contexts).
**Null/control:** Random family assignment of same family sizes.
**Value:** MEDIUM-HIGH — extends H01 to non-immune context; tests whether gene family encoding is a general principle.
**Cost:** LOW — same pipeline as H01 immune family AUROC; only need to curate family lists.

---

### H-N14: Co-Expression Mediation Control (STRING Signal Independent of Expression)
**Hypothesis:** The STRING PPI geometric signal is not primarily mediated by gene co-expression; adding co-expression as a predictor does not substantially reduce STRING partial R².
**Test:** Obtain a co-expression reference matrix for the 209 genes (from GEO or existing single-cell data). Add co-expression as a 5th predictor in the partial Spearman R² joint model. Test if STRING partial R² drops by > 50% when co-expression is included.
**Expected signal:** STRING partial R² stays ≥ 0.005 after co-expression inclusion (not mediated by expression).
**Null/control:** Replace co-expression with shuffled co-expression.
**Value:** HIGH — critical confound control; rules out "scGPT learned expression correlation" alternative.
**Cost:** HIGH — requires external co-expression reference dataset; significant data fetch.

---

## Top 3 for Immediate Execution

### Candidate 1 — High-Probability Discovery: H-N2 (PC1 Biological Identity at L11)
**Rationale:** The dimensional collapse (H03) is already a major finding. Identifying what PC1 encodes at L11 transforms it from a geometric observation into a mechanistic one. If PC1 correlates with TF status, network degree, or HLA identity, this is a direct mechanistic result with no parallel in prior scGPT literature. Fully computable from existing embeddings + gene metadata already assembled across 26 iterations. Near-zero marginal cost.
**Predicted outcome:** PC1 loadings at L11 correlate with is_TF (ρ > 0.3, p < 0.01) and/or STRING hub degree (ρ > 0.25). Possibly: HLA genes cluster at one extreme of PC1.
**Fallback:** If no single attribute dominates PC1, report the null finding clearly (PC1 is a mixed axis) and pivot to PC1 vs PC2 decomposition.

---

### Candidate 2 — High-Risk / High-Reward: H-N1 (TRRUST Activator vs. Repressor Polarity)
**Rationale:** Directional regulatory encoding — whether the model distinguishes gene activation from repression geometrically — is the single most novel biological finding we haven't confirmed yet. TRRUST has the A/R field. Iter_0026 H02 was blocked because Dorothea lacks MOR; TRRUST is the correct data source. If AUROC(activator < repressor distance) > 0.55, this is a first-in-kind demonstration that a genomic language model encodes regulatory sign. High reward. Low cost.
**Predicted outcome:** AUROC 0.55–0.62 at L5–L8 for activator pairs being closer than repressor pairs.
**Fallback:** If TRRUST has too few pairs (< 50 in each class), combine with Dorothea high-confidence A/B pairs from any available full Dorothea file with MOR.

---

### Candidate 3 — Cheap Broad Screen: H-N6 + H-N12 combined (STRING Threshold Sensitivity + PR Elbow Alignment)
**Rationale:** H-N6 (STRING threshold sweep) characterizes signal robustness at essentially zero marginal cost — filters existing STRING cache at 5 thresholds, runs in < 5 minutes. H-N12 (PR elbow vs. biological peak alignment) requires no new computation at all — purely derives from existing h03 and prior AUROC data. Together they provide two methods-validation results and one mechanistic narrative linkage in a single cheap pass.
**Predicted outcome:** H-N6: monotonic ρ decrease with threshold (rho < -0.9); H-N12: peak ΔPR at L7–L8 matches Dorothea peak at L7.
**Fallback:** Either analysis is self-contained and can be run independently.

---

## Summary Table

| ID | Hypothesis | Value | Cost | Priority |
|----|-----------|-------|------|----------|
| H-N1 | TRRUST activator/repressor polarity | HIGH | LOW | **Execute now (C2)** |
| H-N2 | PC1 biological identity at L11 | HIGH | LOW | **Execute now (C1)** |
| H-N3 | Per-family layer curves (structural vs regulatory) | HIGH | LOW | Tier 1 |
| H-N4 | PR collapse rate per bio subset | MEDIUM-HIGH | LOW | Tier 1 |
| H-N5 | Hub vs. peripheral embedding density | MEDIUM | LOW | Tier 2 |
| H-N6 | STRING threshold sensitivity | MEDIUM | LOW | **Execute now (C3)** |
| H-N7 | Oncogene/tumor suppressor polarity | MEDIUM | LOW-MED | Tier 2 |
| H-N8 | PH of GO CC subcellular clusters | HIGH | MEDIUM | Tier 2 |
| H-N9 | Geneformer cross-model replication | HIGH | MEDIUM | Tier 1 |
| H-N10 | Disease gene clustering (DisGeNET) | HIGH | MEDIUM | Tier 2 |
| H-N11 | KLF outlier rescue | MEDIUM | LOW | Tier 3 (trivial add-on) |
| H-N12 | PR elbow alignment with bio peaks | MEDIUM-HIGH | LOW | **Execute now (C3)** |
| H-N13 | Expanded TF family screen | MEDIUM-HIGH | LOW | Tier 1 |
| H-N14 | Co-expression mediation control | HIGH | HIGH | Tier 3 |
