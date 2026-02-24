# Brainstormer Hypothesis Roadmap — iter_0028

**Date:** 2026-02-22

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| Dorothea regulatory proximity, same-gene-pool null | **RETIRE NOW** | H03 confirmed null-sensitivity; AUROC=0.658 in iter_0026 was hub-centrality artifact. |
| TRRUST activation/repression population-level polarity | **RETIRE NOW** | H02: AUROC=0.514 across all layers and metrics. Clean null. |
| PC1 identity via binary TF/HLA/AP1 features | **RETIRE NOW** | Three binary features tested, all AUROC 0.43–0.49. Binary approach exhausted. |
| STRING Spearman ρ without degree-matched null | **DEPRIORITIZE** | Likely same centrality confound as Dorothea. Must use matched null before reporting. |
| BCL2fam / TNFSF family analysis | **RETIRED** (iter_0026) | Negative; no rescue. |
| ETS1 within-TF polarity | **LOW PRIORITY** | Only 2 repression targets; insufficient statistical power. Revisit only if STRING MOR data available. |
| KLF family outlier rescue (KLF14) | **LOW PRIORITY** | Tier-3 single data point; execute only as a cheap add-on. |

---

## New Hypothesis Portfolio

### N1: Hub Centrality as Primary Signal — Direct Characterization
**Hypothesis:** High-STRING-degree hub genes occupy significantly more central positions in scGPT embedding space than low-degree genes, and this centrality increases across layers (hubs become more central in deeper layers).
**Test:** For each of 209 named genes, compute embedding centrality = mean L2 distance to all other 208 genes at each layer. Correlate with STRING degree (Spearman ρ). Test per-layer. Also compute Mann-Whitney for degree Q4 vs Q1 at each layer.
**Expected signal:** Spearman ρ(degree, −centrality_distance) > 0.25 at L7–L11 (p < 0.01). Centrality-degree correlation increases monotonically with layer.
**Null/control:** Shuffle STRING degree assignments across genes.
**Value:** HIGH | **Cost:** LOW
**Justification:** This converts the H03 confound into a first-class result. "Hub genes become geometrically central in deep layers" is a concrete, novel finding about how scGPT encodes connectivity — publishable in its own right.

---

### N2: Matched-Degree Null for Regulatory Proximity
**Hypothesis:** Dorothea A/B regulatory pairs are geometrically closer than degree-matched random pairs (controlling for hub centrality), indicating specificity of regulatory encoding beyond mere hub clustering.
**Test:** For each Dorothea A/B pair (TF_i → target_j), construct 10 matched null pairs by sampling: a random gene with STRING degree within ±2 of TF_i, and a random gene with STRING degree within ±2 of target_j. Compute AUROC(d_regulatory < d_matched_null) at each layer.
**Expected signal if true:** AUROC > 0.55 at L5–L8 (genuine regulatory specificity beyond hub effect).
**Expected signal if false (hub confound only):** AUROC ≈ 0.5 (null is flat once degree is matched).
**Null/control:** Shuffle TF→target assignments while preserving degree distribution.
**Value:** HIGH | **Cost:** MEDIUM (requires degree-matching logic + more null pairs)
**Justification:** This is the definitive experiment to determine whether prior positive results have any real signal. If AUROC collapses to ~0.5, the regulatory proximity story is dead. If it stays >0.55, we have a genuine specificity finding.

---

### N3: PC1 Cell-Type Marker Axis at L11
**Hypothesis:** The dominant PC1 axis at L11 encodes a T-cell activation vs. antigen-presentation cell-type polarity, reflecting the dominant cell-type composition in PBMC training data.
**Test:** Obtain canonical PBMC cell-type marker lists (e.g., T-cell effectors: CD8A, GZMB, PRF1, CD3D, LCK, CCL5; antigen-presenting: HLA-A, HLA-DPB1, CIITA, CD74, HLA-DRA). Score each of 209 named genes by "T-cell score" (overlap or enrichment with T-cell marker list) vs "APC score." Compute Spearman ρ(PC1 loading, T-cell_score − APC_score) at L11. Permutation test.
**Expected signal:** ρ > 0.3, p < 0.01 by permutation; JUN/LCK/CCL5 at one end, HLA genes at other.
**Null/control:** Random gene-score assignments.
**Value:** HIGH | **Cost:** LOW (marker lists are curated text, no external API needed)
**Justification:** Direct mechanistic interpretation of the dominant geometric axis. JUN-top/FOS-bottom already suggests this; LCK (T-cell kinase) and HLA genes confirm the polarity.

---

### N4: Gene Family Clustering Confound Audit
**Hypothesis:** Prior immune family AUROC results (HLA-I AUROC ~0.75, AP1/RUNX ~0.70) are partially confounded by differential hub status within family vs cross-family comparisons.
**Test:** For each immune family tested, compute mean STRING degree within-family vs. mean STRING degree in null pairs. If within-family mean degree > null mean degree, hub confound is present. Compute "hub-corrected AUROC" by restricting cross-family null to genes with matching degree distribution.
**Expected signal:** HLA-I AUROC is robust (HLA genes are not hubs); AP1/JUN/FOS family AUROC is partially confounded.
**Null/control:** Original null from iter_0026 family analysis.
**Value:** HIGH | **Cost:** LOW (uses existing results; adds degree check)
**Justification:** Critical for knowing which prior positive results to keep in the paper. HLA-I result probably survives; AP1 result needs scrutiny.

---

### N5: Layer Eigenspectrum Shape — Beyond PR
**Hypothesis:** Across layers, the PCA eigenspectrum transitions from approximately power-law (early layers, many components contribute) to rank-1 or rank-2 dominated (late layers after collapse), and this transition is sharp (concentrated at a specific layer), not gradual.
**Test:** At each of 12 layers, compute the full PCA eigenspectrum for all 209 genes. Normalize. Fit to power-law and exponential distributions. Compute "spectral entropy" (−Σ p_i log p_i where p_i = var_i / total_var). Track spectral entropy across layers. Find the sharpest drop in entropy.
**Expected signal:** Spectral entropy drops sharply at L8–L10; not linear. Sharp transition layer ≠ early/late but concentrated in a few layers.
**Null/control:** Shuffled embeddings (expected flat entropy profile).
**Value:** MEDIUM-HIGH | **Cost:** LOW (pure eigenspectrum computation on existing embeddings)

---

### N6: PBMC Cell-Type Specificity Score vs Embedding Centrality
**Hypothesis:** Genes with high cell-type specificity (expressed in 1–2 cell types only) occupy peripheral positions in embedding space at late layers, while broadly expressed "housekeeping" genes are central.
**Test:** Use existing PBMC single-cell expression data (from scGPT training data or GEO) to compute a "specificity score" for each of 209 named genes (e.g., Gini coefficient of expression across cell types, or τ-score). Correlate with L11 embedding centrality (mean distance to all others). Spearman ρ per layer.
**Expected signal:** Spearman ρ(specificity, centrality_distance) > 0.25 at L11 (more specific → more peripheral).
**Null/control:** Shuffle specificity scores.
**Value:** HIGH | **Cost:** MEDIUM (needs PBMC expression reference; check if scGPT training data accessible in subproject_38)

---

### N7: Persistent Homology of Hub vs. Peripheral Gene Subsets
**Hypothesis:** Hub genes (top 25% STRING degree) form topologically denser 0D clusters (longer-lived connected components) in VR filtration than peripheral genes, and this density increases in deeper layers.
**Test:** Split 209 genes into hub (Q4 degree) and peripheral (Q1 degree) subsets (~52 genes each). At each layer, compute 0D PH (Ripser) on each subset. Compare mean 0D persistence and number of components at death threshold. Test: hub subset shows more cohesive topology (fewer, longer-lived components) at L8–L11.
**Expected signal:** Mean 0D persistence for hub subset > peripheral subset at L10–L11 (p < 0.05 by bootstrap).
**Null/control:** Random same-size subsets.
**Value:** MEDIUM-HIGH | **Cost:** MEDIUM (Ripser on 52-gene subsets is fast; needs PH setup)

---

### N8: STRING Score Threshold Sweep with Matched-Degree Null
**Hypothesis:** The STRING geometric signal (Spearman ρ between STRING score and pairwise embedding distance) is robust to degree-matching but degrades with increasing score threshold.
**Test:** At L8, for 5 STRING score thresholds (0.4, 0.5, 0.6, 0.7, 0.9), compute Spearman ρ between binary STRING interaction and L2 distance, using degree-matched null. Compare to ρ without degree matching.
**Expected signal:** Degree-matched ρ is lower than unmatched ρ at all thresholds (hub correction matters), but monotonic decrease with threshold is preserved.
**Null/control:** Unmatched ρ (same as prior iterations); shuffled STRING scores.
**Value:** MEDIUM | **Cost:** LOW

---

### N9: Geneformer Cross-Model Hub Centrality Replication
**Hypothesis:** The hub centrality signal (high-degree genes are more central in embedding space) replicates in Geneformer embeddings, indicating this is a property of genomic language model training rather than scGPT-specific architecture.
**Test:** Load Geneformer embeddings for 209 named genes (check subproject_38 outputs). Compute mean pairwise L2 distance per gene at each Geneformer layer. Correlate with STRING degree. Compare Spearman ρ profile across Geneformer layers to scGPT layers.
**Expected signal:** Geneformer shows Spearman ρ(degree, −centrality) > 0.2 at mid-late layers.
**Null/control:** Shuffled Geneformer gene identities.
**Value:** HIGH | **Cost:** MEDIUM (depends on Geneformer availability in subproject_38)

---

### N10: Oncogene / Tumor Suppressor Polarity (COSMIC)
**Hypothesis:** Oncogenes and tumor suppressors annotated in COSMIC among the 209 genes are linearly separable in scGPT embedding space at late layers.
**Test:** Label 209 genes with COSMIC Gene Census (oncogene / tumor suppressor / other). Require ≥8 per class. LDA AUROC leave-one-out at each layer. Also test if tumor suppressors are more peripheral (higher mean distance to all others) than oncogenes.
**Expected signal:** LDA AUROC > 0.65 at L8–L11. Tumor suppressors more peripheral than oncogenes.
**Null/control:** Random class assignments of same size.
**Value:** MEDIUM-HIGH | **Cost:** LOW-MEDIUM

---

### N11: Layer-by-Layer Effective Connectivity (kNN Graph Clustering Coefficient)
**Hypothesis:** The embedding k-NN graph (k=10) becomes more clustered (higher global clustering coefficient) in deeper layers, reflecting progressive formation of tight gene neighborhoods as representation consolidates.
**Test:** At each of 12 layers, build k-NN graph (k=10) on 209 genes using L2 distance. Compute global clustering coefficient. Test monotonic increase.
**Expected signal:** Clustering coefficient increases from ~0.15 at L0 to ~0.40 at L11 (monotonic or near-monotonic).
**Null/control:** Random permutation of embedding positions; expected flat clustering profile.
**Value:** MEDIUM | **Cost:** LOW

---

### N12: AP1 JUN/FOS Pair Geometric Trajectory — Layer-by-Layer Distance
**Hypothesis:** JUN and FOS, which are top vs. bottom PC1 genes at L11 (opposing extremes of the dominant axis), start geometrically CLOSE in early layers (co-regulation signal) and diverge in late layers as the PC1 axis crystallizes.
**Test:** Compute L2(JUN, FOS) at each of 12 layers. Test for monotonic increase from L0 to L11. Compare to mean L2 between random gene pairs (matched by early-layer distance).
**Expected signal:** L2(JUN, FOS) increases from ~L5 onward, crossing above mean random distance by L9–L11.
**Null/control:** 20 random gene pairs with similar L0 distance to JUN-FOS.
**Value:** MEDIUM-HIGH | **Cost:** LOW (single pair trajectory, near-zero cost)
**Bonus:** Can extend to all AP1 family pairs and to HLA-A/HLA-DPB1 pair for comparison.

---

### N13: Expanded Non-Immune TF Superfamily Clustering (ETS, STAT, NFkB, Nuclear Receptor)
**Hypothesis:** TF superfamilies (ETS, STAT, NFkB, nuclear receptors) present in the 209-gene vocab show within-family AUROC > 0.65 at L7–L8, confirming that TF family clustering is a general property beyond immune genes.
**Test:** Curate TF superfamily groups from 209-gene vocab: ETS (ETS1, ETS2, GABPA, FLI1), STAT (STAT1/2/3/5/6), NFkB (RELA, RELB, NFKB1/2), NR (NR3C1, NR4A1/2/3, RARA). Compute within-family AUROC at each layer vs cross-family same-size null.
**Expected signal:** STAT and NFkB AUROC > 0.70 at L7–L8. Nuclear receptors lower (~0.60) due to diverse ligand contexts.
**Null/control:** Random family assignment of same family sizes.
**Value:** MEDIUM-HIGH | **Cost:** LOW

---

### N14: PR Collapse Rate per Biologically-Defined Gene Subset
**Hypothesis:** Biologically coherent gene subsets (HLA-I family, STAT family, STRING hub Q4) show faster PR collapse across layers than size-matched random subsets.
**Test:** For each of 4 subsets (HLA-I, STAT, hub Q4, random), compute PR at each layer (subset-specific PCA). Compute collapse rate = (PR_L0 − PR_L11) / PR_L0. Mann-Whitney test subset collapse rate vs random bootstrap.
**Expected signal:** HLA-I and STAT collapse faster (rate > 0.90) than random (rate ~0.75–0.80).
**Null/control:** 10 random same-size gene subsets per biological subset.
**Value:** MEDIUM | **Cost:** LOW

---

## Top 3 for Immediate Execution

### Candidate 1 — High-Probability Discovery: N1 (Hub Centrality as Primary Signal)
**Rationale:** The H03 finding revealed that hub centrality is driving apparent regulatory proximity signals. Rather than treating this as an embarrassment, characterize it head-on: hub genes become progressively more central in deeper layers. This is directly testable with existing embeddings and STRING data (no new data), runs in minutes, and is interpretable — high-degree network genes occupy geometric centers in scGPT representation, likely because they co-occur across many cell types in training data. This is a novel, positive, well-controlled finding that reframes the iter_0022–0026 results constructively.
**Test:** Spearman ρ(STRING degree, −mean_L2_to_all_others) per layer. Per-layer correlation profile is the main deliverable.
**Fallback:** If STRING degree data unavailable, use TRRUST target count (already computed as "hub degree" in H01).

---

### Candidate 2 — High-Risk / High-Reward: N3 (PC1 Cell-Type Marker Axis at L11)
**Rationale:** H01 showed PC1 at L11 is not separable by binary features. But the gene pattern (JUN, LCK, BIRC3 at top; FOS, HLA-A, HLA-DPB1 at bottom) is too biologically coherent to be noise. A T-cell activation vs APC cell-type axis is the most parsimonious biological interpretation. If a simple curated PBMC marker gene score correlates with PC1 loadings (ρ > 0.3), this transforms the collapse finding into a mechanistic narrative: the model compresses immune gene representations onto a dominant cell-type axis. This is a high-reward narrative pivot. Risk: marker list curation requires some care; only 209 genes, so marker overlap may be thin.
**Fallback:** If PBMC marker overlap < 15 genes, use GenAge / COSMIC binary annotation as a broader continuous proxy, or try gene expression variance from a PBMC reference.

---

### Candidate 3 — Cheap Broad Screen: N12 + N5 combined (JUN/FOS Divergence + Layer Spectral Entropy)
**Rationale:** N12 (JUN-FOS distance trajectory across 12 layers) requires computing L2 for a single gene pair per layer — 5 minutes of code. N5 (spectral entropy of eigenspectrum per layer) requires computing PCA eigenspectra across layers — also existing computation from prior PCA analyses, just not extracted as entropy. Together they provide (a) a concrete single-pair narrative anchor for the PC1 divergence story, and (b) a clean geometric characterization of when the dimensional collapse is sharpest. Both are pure computation on existing embeddings, no new data.

---

## Summary Table

| ID | Hypothesis | Value | Cost | Priority |
|----|-----------|-------|------|----------|
| N1 | Hub centrality as primary signal | HIGH | LOW | **Execute now (C1)** |
| N2 | Matched-degree null for regulatory proximity | HIGH | MEDIUM | Tier 1 |
| N3 | PC1 cell-type marker axis at L11 | HIGH | LOW | **Execute now (C2)** |
| N4 | Gene family clustering confound audit | HIGH | LOW | Tier 1 |
| N5 | Layer eigenspectrum shape / spectral entropy | MEDIUM-HIGH | LOW | **Execute now (C3 part)** |
| N6 | PBMC cell-type specificity vs embedding centrality | HIGH | MEDIUM | Tier 1 |
| N7 | PH of hub vs peripheral gene subsets | MEDIUM-HIGH | MEDIUM | Tier 2 |
| N8 | STRING threshold sweep, matched-degree null | MEDIUM | LOW | Tier 2 |
| N9 | Geneformer cross-model hub centrality replication | HIGH | MEDIUM | Tier 1 |
| N10 | Oncogene/tumor suppressor polarity (COSMIC) | MEDIUM-HIGH | LOW-MED | Tier 2 |
| N11 | kNN graph clustering coefficient per layer | MEDIUM | LOW | Tier 3 |
| N12 | JUN/FOS pair divergence trajectory | MEDIUM-HIGH | LOW | **Execute now (C3 part)** |
| N13 | Expanded TF superfamily clustering | MEDIUM-HIGH | LOW | Tier 1 |
| N14 | PR collapse rate per bio subset | MEDIUM | LOW | Tier 2 |
