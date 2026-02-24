# Brainstormer Hypothesis Roadmap — iter_0025

**Date:** 2026-02-22

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| Chromosomal proximity confound | **RESOLVED — retire** | AUROC 0.515, 10× below functional signals, fades at L11. Definitively ruled out. No further work. |
| GO BP vs GO CC redundancy question | **RESOLVED — retire** | Both are significant and partially independent (partial R²: CC=0.0052, BP=0.0001 after joint residualization). Question answered. |
| Further univariate STRING/Dorothea sweeps | **DEPRIORITIZE** | Univariate signal well-characterized across 12+ iterations. Marginal returns now low. |
| TRRUST solo analysis | **DEPRIORITIZE** | TRRUST activation AUROC characterized (L5=0.659). Further solo TRRUST work has diminishing returns unless testing new TRRUST-specific sub-hypotheses. |
| Topological screen without biological anchor | **DEPRIORITIZE** | Pure topology without biological grounding produced weak/noisy signal in early iterations. Only revisit if anchored to specific biology. |

---

## New Hypothesis Portfolio

### H-A: Geneformer Cross-Model Joint Predictor Replication
**Hypothesis:** The multi-predictor joint model (STRING + Dorothea + GO CC) replicates in Geneformer embeddings with comparable or higher total partial R².
**Test:** Load Geneformer gene embeddings for the same 209 named genes. Compute identical partial Spearman R² joint model across Geneformer layers.
**Expected signal:** Total partial R² > 0.01 (comparable to scGPT's 0.017); similar predictor ranking.
**Null/control:** Shuffled gene identities → R² collapses to 0.
**Value:** HIGH. Confirms generalizability; this is the single most impactful validation experiment available.
**Cost:** MEDIUM (requires Geneformer embedding infrastructure; may already exist in subproject_38).

---

### H-B: Protein Complex Co-Membership (CORUM)
**Hypothesis:** Genes co-belonging to the same protein complex (CORUM database) are geometrically closer in scGPT embeddings than gene pairs not sharing a complex.
**Test:** Download CORUM human complex annotations. For 209 named genes, classify pairs as same-complex vs. different-complex. Mann-Whitney AUROC at each layer.
**Expected signal:** AUROC > 0.60 at mid layers (complexes are tightly coordinated; should be stronger signal than generic PPI STRING).
**Null/control:** Randomly re-assign complex memberships → AUROC → 0.5.
**Value:** HIGH. CORUM is a cleaner functional signal than generic STRING (direct physical co-complex vs. indirect functional association).
**Cost:** LOW (CORUM download + same infrastructure as STRING AUROC test).

---

### H-C: Attention Head Attribution — Which Heads Drive GO BP Late Encoding
**Hypothesis:** Specific attention heads in layers L6–L7 drive the late GO BP signal; ablating or masking these heads reduces GO BP Spearman ρ while leaving GO CC (L5) intact.
**Test:** Extract per-head attention weights from scGPT at L5–L7. For each head, compute correlation between attention weight to a gene and GO BP Jaccard similarity of that gene to its neighbors. Rank heads by GO BP attribution score. Optionally: zero out top-k heads and re-measure GO BP Spearman ρ.
**Expected signal:** 1–3 attention heads with GO BP attribution > 2× mean; ablation reduces GO BP ρ by ≥ 30%.
**Null/control:** Same analysis for shuffled GO BP labels.
**Value:** HIGH (mechanistic; publication-grade finding if confirmed).
**Cost:** HIGH (requires attention weight extraction from scGPT forward pass; non-trivial code).

---

### H-D: Regulatory Direction Encoding (Dorothea Activators vs. Repressors)
**Hypothesis:** Dorothea-annotated activator TF→target pairs are geometrically closer in scGPT embeddings than repressor pairs, encoding regulatory polarity.
**Test:** From Dorothea, separate pairs by direction (activation vs. repression). Compute L2 distances for each group. Mann-Whitney test: activator pairs closer than repressor pairs?
**Expected signal:** AUROC (activator < repressor distance) > 0.55 at mid-late layers.
**Null/control:** Shuffle activation/repression labels within the Dorothea set.
**Value:** HIGH (if positive, shows directional regulatory information is geometrically encoded — not just interaction strength).
**Cost:** LOW (Dorothea has mode_of_regulation field; reuses existing pipeline).

---

### H-E: Disease Gene Manifold Clustering (OMIM/DisGeNET)
**Hypothesis:** Genes associated with the same disease (from OMIM or DisGeNET) cluster in scGPT embedding space more than gene pairs from different diseases.
**Test:** Fetch disease-gene associations for the 209 named genes. Classify pairs as same-disease vs. different-disease. AUROC at each layer.
**Expected signal:** AUROC > 0.58 in mid-late layers for same-disease gene pairs.
**Null/control:** Random disease assignment.
**Value:** HIGH (clinical relevance angle; strong for paper narrative).
**Cost:** MEDIUM (API fetch + new pairwise classification).

---

### H-F: Persistent Homology of GO CC Subcellular Compartment Clusters
**Hypothesis:** GO CC-defined subcellular compartments (nucleus, cytoplasm, membrane, mitochondria) form topologically persistent 0-dimensional barcodes in the scGPT manifold, with persistence lengths proportional to compartment specificity.
**Test:** For each major GO CC class among 209 named genes, extract the gene embedding submanifold. Compute 0D persistent homology (VR filtration). Measure mean persistence of connected components. Compare across compartment types and across layers.
**Expected signal:** Nucleus/membrane genes show distinct, persistent clusters (long barcodes) vs. cytoplasm (more dispersed). Peak persistence at L5 (matches GO CC peak).
**Null/control:** Randomly labeled gene subsets of same size.
**Value:** HIGH (topological validation of the pairwise distance signal; publishable if PH barcodes differentiate compartments).
**Cost:** MEDIUM (requires ripser/gudhi; subsets small enough for PH computation).

---

### H-G: Intrinsic Dimension per Biological Anchor Subset
**Hypothesis:** Genes from the same GO CC compartment or STRING PPI cluster occupy a lower-dimensional submanifold than randomly selected gene sets, reflecting biological constraint on representational degrees of freedom.
**Test:** For each biological subset (same compartment, same complex, top-K STRING neighbors), estimate intrinsic dimension via TwoNN or MLE. Compare to random same-size subsets. Test if biologically defined subsets have significantly lower intrinsic dimension.
**Expected signal:** GO CC same-compartment subsets: intrinsic dimension 20–40% lower than random subsets of same size.
**Null/control:** Random gene subsets of matching size.
**Value:** MEDIUM-HIGH (new geometric characterization; fills topological gap).
**Cost:** MEDIUM (TwoNN is fast; need clean subset extraction).

---

### H-H: Layer-Resolved Timeline for Geneformer (Cross-Model Timeline Comparison)
**Hypothesis:** Geneformer encodes the same biological hierarchy (compartment before process) but with different peak layers, reflecting architectural differences from scGPT.
**Test:** Replicate H02 (layer timeline) on Geneformer: compile per-layer STRING Spearman, GO CC Spearman, GO BP Spearman across Geneformer layers. Compare peak layer ordering to scGPT.
**Expected signal:** Same qualitative ordering (structural < functional) but shifted peak layers; peak span ≥ 3 layers.
**Null/control:** Random label shuffling per anchor.
**Value:** HIGH (cross-model mechanistic comparison; direct follow-up to H02).
**Cost:** LOW-MEDIUM (reuses H02 pipeline; only needs Geneformer embeddings).

---

### H-I: STRING Score Threshold Sensitivity — Signal at Low-Confidence Interactions
**Hypothesis:** The STRING Spearman ρ signal degrades monotonically as score threshold increases from 0.4 to 0.9, indicating that low-confidence interactions contain genuine geometric signal (not noise).
**Test:** Recompute STRING Spearman ρ at L8 for 5 thresholds: 0.4, 0.5, 0.6, 0.7, 0.9. Plot ρ vs threshold. Test for monotonic decrease using Spearman rank correlation of thresholds vs ρ.
**Expected signal:** ρ decreases from ~0.152 (threshold 0.4) to ~0.05 at threshold 0.9. Monotonic decrease (Spearman ρ of ρ-vs-threshold = negative, p < 0.05).
**Null/control:** Same analysis with shuffled STRING scores.
**Value:** MEDIUM (characterizes signal quality; useful for methods section; may reveal threshold-dependent biology).
**Cost:** LOW (reuses STRING cache; simple threshold filter).

---

### H-J: Network Hub vs. Peripheral Gene Embedding Density
**Hypothesis:** High-degree STRING network hubs (top 10% by degree) occupy a denser region of embedding space (more nearest neighbors within fixed radius) than low-degree peripheral genes.
**Test:** Compute STRING degree for 209 named genes using the existing cache. Split into quartiles. For each gene, compute mean L2 distance to its k nearest neighbors in embedding space (local density). Test for hub-vs-peripheral density difference via Mann-Whitney at each layer.
**Expected signal:** Hub genes have lower mean k-NN distance (denser local neighborhood) at mid-late layers.
**Null/control:** Random hub assignment.
**Value:** MEDIUM (links network topology to embedding geometry; new angle).
**Cost:** LOW (reuses STRING cache + existing embedding).

---

### H-K: Pathway-Specific Layer Timeline (KEGG/Reactome Sub-Analysis)
**Hypothesis:** Specific biological pathways (e.g., cell cycle, DNA repair, immune response) have sharper, more concentrated layer-resolved peaks than the aggregate GO BP signal, because they represent more coherent functional modules.
**Test:** Fetch KEGG/Reactome pathway memberships for 209 named genes. For the 5 largest pathways represented (≥10 genes), compute within-pathway mean pairwise L2 distance at each layer. Find peak layer per pathway. Compare peak sharpness (peak/mean ratio) to aggregate GO BP.
**Expected signal:** Cell-cycle pathway peaks sharply at specific layer; peak/mean ratio > 1.5 vs. GO BP ratio ~1.1.
**Null/control:** Random gene set of same size as each pathway.
**Value:** MEDIUM-HIGH (refinement of H02 finding; pathway-specific mechanistic insight).
**Cost:** MEDIUM (API fetch for KEGG/Reactome; loop over pathways).

---

### H-L: Co-Expression Mediation Test (Is STRING Signal Mediated by Co-Expression?)
**Hypothesis:** The STRING PPI geometric signal is NOT primarily mediated by gene co-expression; PPI-interacting genes are geometrically similar even after controlling for expression correlation.
**Test:** Obtain co-expression correlation matrix for the 209 named genes from a reference scRNA-seq dataset. Add co-expression as a 5th predictor in the partial Spearman model. Test whether STRING partial R² drops substantially (>50%) when co-expression is added.
**Expected signal:** STRING partial R² stays ≥ 0.006 after co-expression inclusion (i.e., not primarily mediated by expression).
**Null/control:** Replace co-expression with shuffled co-expression.
**Value:** HIGH (critical confound control; rules out the "scGPT just learned expression correlation" alternative).
**Cost:** HIGH (requires appropriate co-expression reference dataset).

---

### H-M: Oncogene / Tumor Suppressor Polarity in Embedding Space
**Hypothesis:** Oncogenes and tumor suppressors occupy distinct, separable regions of scGPT embedding space at late layers, reflecting their opposing functional roles.
**Test:** Curate lists of known oncogenes and tumor suppressors among the 209 named genes (from CGC/COSMIC). Project embeddings to 2D (UMAP/PCA). Compute AUROC for linear separation (oncogene vs. tumor suppressor classification using LDA on embeddings). Test layer by layer.
**Expected signal:** AUROC > 0.65 at L8–L11 for oncogene/tumor-suppressor separation.
**Null/control:** Random oncogene/tumor suppressor label shuffling.
**Value:** MEDIUM (clinically compelling narrative; depends on having enough labeled genes).
**Cost:** LOW-MEDIUM (CGC list available; requires ≥20 labeled genes per class in the 209 set).

---

## Top 3 for Immediate Execution

### Candidate 1 — High-Probability Discovery: H-B (Protein Complex Co-Membership, CORUM)
- **Rationale:** Protein complexes are the strongest possible functional signal (direct physical interaction, curated database). Infrastructure is already in place. Expected AUROC > 0.60, stronger than generic STRING. Low cost, high probability of strong positive result.
- **Predicted outcome:** AUROC 0.62–0.72 at L5–L8; cleaner, sharper signal than STRING PPI.
- **Fallback:** If CORUM genes are underrepresented in the 209 set, pivot to STRING subnetworks at threshold 0.9.

### Candidate 2 — High-Risk / High-Reward: H-C (Attention Head Attribution)
- **Rationale:** Mechanistic attribution of which attention heads drive the GO BP late-layer signal would transform the paper from observational to mechanistic. High implementation cost but the payoff (specific circuit identification) is publication-defining.
- **Predicted outcome:** 1–2 heads in L6–L7 with significantly elevated GO BP attribution scores.
- **Fallback:** If full ablation is infeasible, settle for correlation between per-head attention patterns and GO BP Jaccard similarity (no ablation required — pure correlation analysis is cheap).

### Candidate 3 — Cheap Broad Screen: H-D + H-I combined (Regulatory Direction + STRING Threshold Sensitivity)
- **Rationale:** Both reuse existing infrastructure with minimal new code. H-D (Dorothea activator/repressor split) adds directionality signal at near-zero cost. H-I (STRING threshold sweep) characterizes robustness of the core result. Together they fill important methods validation + mechanistic refinement gaps in one pass.
- **Predicted outcome:** H-D: AUROC ~0.55–0.60 for activator < repressor distance; H-I: monotonic ρ decrease with threshold (ρ_Spearman of ρ-vs-threshold < -0.9).
- **Fallback:** If Dorothea direction field is missing for most pairs, do H-I alone.

---

## Summary Table

| ID | Hypothesis | Value | Cost | Priority |
|----|-----------|-------|------|----------|
| H-A | Geneformer joint model replication | HIGH | MEDIUM | Tier 1 |
| H-B | CORUM protein complex clustering | HIGH | LOW | **Execute now** |
| H-C | Attention head attribution (GO BP) | HIGH | HIGH | **Execute now** |
| H-D | Dorothea activator/repressor polarity | HIGH | LOW | **Execute now** |
| H-E | Disease gene clustering (OMIM) | HIGH | MEDIUM | Tier 2 |
| H-F | PH of GO CC subcellular clusters | HIGH | MEDIUM | Tier 2 |
| H-G | Intrinsic dimension per bio subset | MEDIUM-HIGH | MEDIUM | Tier 2 |
| H-H | Geneformer layer timeline | HIGH | LOW-MED | Tier 1 |
| H-I | STRING threshold sensitivity | MEDIUM | LOW | **Execute now** |
| H-J | Hub vs. peripheral embedding density | MEDIUM | LOW | Tier 2 |
| H-K | Pathway-specific layer timeline | MEDIUM-HIGH | MEDIUM | Tier 2 |
| H-L | Co-expression mediation test | HIGH | HIGH | Tier 3 |
| H-M | Oncogene/tumor suppressor polarity | MEDIUM | LOW-MED | Tier 2 |
