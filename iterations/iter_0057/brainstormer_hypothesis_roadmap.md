# Brainstormer Hypothesis Roadmap — iter_0057

**Date**: 2026-02-23

---

## Retire / Deprioritize

| ID | Description | Status | Reason |
|----|-------------|--------|--------|
| H02-iter57 | SV basis permutation hypothesis | `retire_now` | PA=74.8° rules it out; no rescue value |
| Subspace rotation angle as AUROC predictor | (iter_0056 H03) | `retire_now` | Definitively negative (rho=-0.27); superseded by H03-iter57 which gives the right predictor (energy fraction, not rotation) |
| STRING PPI → embedding proximity | (iter_0031) | `retire_now` | Previously retired; do not revisit |
| GO co-regulation → embedding distance | (iter_0032) | `retire_now` | Previously retired |
| TDA loops on circuit genes | (iter_0052) | `retire_now` | Previously retired |
| SV8-14 secondary signal | (rescue candidate iter_0056) | `deprioritize` | Cross-seed validation of 6D is the priority; SV8-14 has been consistently weak |

**No rescue candidates this round** — the SV8-14 rescue has been repeatedly weak and the main signal is already validated.

---

## New Hypothesis Portfolio

### A. Boundary Interpretation (Direct Follow-Ups from H01 Validation)

**A1. TF/Target Boundary Consistency Across Seeds**
- *Hypothesis*: The same TF genes rank highest in classifier confidence across all 3 seeds at L3, indicating the 6D boundary geometry converges to a seed-stable solution.
- *Test*: At L3, for each seed, fit LR on 6D SV2-7. Extract per-gene classifier probability scores. Compute Spearman rank correlation of scores across seeds (main vs s43, main vs s44, s43 vs s44). Permutation null: shuffle gene-to-seed assignment.
- *Expected signal*: Spearman rho > 0.7 for all three pairs; top-20 TFs by confidence overlap ≥ 50% across seeds.
- *Null*: rho < 0.3; top-20 sets share < 5 genes.
- *Value*: high | *Cost*: low

**A2. TF Family Enrichment at 6D Boundary**
- *Hypothesis*: TF genes that score highest in the 6D classifier (top-20 by LR probability) are enriched for specific TF families (e.g., zinc-finger, bHLH, bZIP) based on JASPAR/TFdb annotation, compared to bottom-ranked TFs.
- *Test*: At L3, extract per-TF classifier scores across 3 seeds. Rank TFs. Cross-reference top-20 vs bottom-20 against JASPAR TF family annotation. Fisher's exact test per TF family (≥ 3 members required). FDR correction.
- *Expected signal*: At least one TF family enriched at p_FDR < 0.05 in the high-confidence boundary region.
- *Null*: No TF family enriched above chance in top-20 boundary genes.
- *Value*: high | *Cost*: low

**A3. Non-Linear Boundary Structure (Kernel SVM vs LR)**
- *Hypothesis*: The TF/target separation in 6D SV2-7 is non-linear; RBF-kernel SVM outperforms LR by ≥ 0.03 AUROC at L2-L3.
- *Test*: At L0, L2, L3, L8: compare 5-fold stratified LR vs RBF-SVM (GridSearchCV C and gamma). 100-perm null. Report AUROC delta.
- *Expected signal*: SVM exceeds LR by ≥ 0.03 at peak layers; delta collapses at L8 (where signal is weak), indicating the structure becomes trivially low-dimensional there.
- *Null*: SVM AUROC ≤ LR + 0.02 at all layers.
- *Value*: medium | *Cost*: low

### B. Effective Rank / Intrinsic Dimension Hypotheses

**B1. Effective Rank as Layer Discriminability Predictor (Cross-Seed)**
- *Hypothesis*: The effective rank (exp(entropy of normalized singular value distribution) of the full embedding at each layer correlates positively with joint AUROC, confirming that "richer geometry = better discrimination."
- *Test*: For all 3 seeds at all 12 layers, compute effective rank of the [n_nz × 512] embedding matrix. Spearman(effective_rank, joint_AUROC) pooled across seeds. Compare to rho from H03-iter57 (SV energy fraction, rho=-0.93).
- *Expected signal*: rho > 0.85 (positive, since effective rank is inverse of concentration). Stronger predictor than raw SV fraction.
- *Null*: rho < 0.5 or sign inconsistent across seeds.
- *Value*: high | *Cost*: low

**B2. TwoNN Intrinsic Dimension: TF vs Target Gene Manifolds**
- *Hypothesis*: TF genes occupy a lower intrinsic-dimensional manifold in full 512D embedding space than target-only genes, reflecting concentrated functional role.
- *Test*: At each of 12 layers, estimate TwoNN intrinsic dimension separately for TF genes (n~73) and target-only genes (n~222) in 512D. Track ID ratio (TF/target) across depth. Bootstrap CI (100 resamples of each subpopulation).
- *Expected signal*: TF ID < target ID at L0-L5; ratio may flip at L9 (integration point).
- *Null*: TF ID ≈ target ID; no systematic manifold dimension difference.
- *Value*: medium | *Cost*: low

**B3. Participation Ratio of SV2-7 Subspace Across Seeds**
- *Hypothesis*: The participation ratio (inverse Simpson index of SV energy distribution within SV2-7) is conserved across seeds even as individual SV directions drift, explaining why 6D span is stable while axes are not.
- *Test*: For each seed and layer, compute PR = (sum S[1:7])² / sum(S[1:7]²). Track PR across layers and seeds. Test: PR coefficient of variation across seeds < CV across layers.
- *Expected signal*: Within-layer seed SD of PR < 0.05; PR is a seed-stable descriptor of subspace structure.
- *Null*: PR varies as much across seeds as across layers (no stability advantage).
- *Value*: medium | *Cost*: low

### C. Topological / Geometric Structure

**C1. Persistent Homology on TF vs Target Gene Sets at L3**
- *Hypothesis*: TF genes form a topologically distinct cluster in 6D SV2-7 (lower H0 fragmentation, distinct H1 structure) vs target-only genes at L3 (peak AUROC).
- *Test*: At L3 for main seed: VR persistence on 73 TF genes and 222 target genes in 6D. Compare: H0 component count at epsilon=median_nn_distance; H1 bar lifetime distribution. 200-bootstrap null (random label assignment preserving sizes).
- *Expected signal*: TFs have fewer H0 components and shorter H1 bars (more cohesive, more planar) than targets; bootstrap p < 0.01.
- *Null*: No PH profile difference beyond label-shuffle null.
- *Value*: high | *Cost*: medium

**C2. TRRUST Graph Laplacian Spectral Alignment with SV5-7 Gene Loadings**
- *Hypothesis*: The top eigenvectors of the TRRUST regulatory graph Laplacian are geometrically encoded in the SV5-7 right singular vectors (gene loadings), indicating the embedding directly reflects network topology.
- *Test*: Build TRRUST interaction graph over 295 in-vocab named genes. Compute graph Laplacian; extract top-5 eigenvectors. Compute max cosine similarity between each eigenvector and SV5-7 right singular vectors at L0, L2, L3. 1000-perm null (shuffle gene identities on eigenvectors).
- *Expected signal*: At least one Laplacian eigenvector achieves cosine ≥ 0.15 with an SV5-7 direction; p_perm < 0.05.
- *Null*: Max cosine ≤ random; p_perm > 0.1.
- *Value*: high | *Cost*: medium

**C3. Persistence Landscape Stability of TRRUST Genes Across Seeds**
- *Hypothesis*: The persistence landscape of TRRUST circuit genes in 6D SV2-7 at L0 is stable across seeds (low inter-seed Wasserstein distance) relative to a random gene-set bootstrap.
- *Test*: Compute persistence landscapes (H0, H1) for the 295 TRRUST genes in 6D SV2-7 at L0 for main, s43, s44. Pairwise Wasserstein-2 distance. Compare to 100-bootstrap random 295-gene landscapes (same procedure). Test: circuit gene inter-seed W2 < 5th percentile of random-set W2.
- *Expected signal*: Circuit genes show significantly smaller landscape variance across seeds vs random gene sets.
- *Null*: No stability advantage; circuit gene W2 within random-set distribution.
- *Value*: medium | *Cost*: medium

### D. Generalization and Cross-Model

**D1. OOD Generalization: Train on Immune, Classify on Lung**
- *Hypothesis*: A classifier trained on immune SV2-7 at L3 generalizes to lung embeddings for TF/target prediction, indicating the 6D regulatory geometry is cell-type-agnostic.
- *Test*: Train LR on cycle4_immune SV2-7 at L3 (all 3 seeds, majority-vote ensemble). Apply to lung embeddings (same scGPT, matched gene set). Measure AUROC on lung set vs lung-trained LR upper bound and chance.
- *Expected signal*: Cross-tissue AUROC > 0.60; gap to lung-trained < 0.10.
- *Null*: Cross-tissue AUROC ≤ 0.55.
- *Value*: high | *Cost*: medium

**D2. Variance Spread — AUROC Relationship in Lung/Other Cell Types**
- *Hypothesis*: The rho=-0.93 (H03-iter57) between SV energy concentration and AUROC is reproduced in lung/other cell-type embeddings, confirming it is a general architectural property.
- *Test*: Repeat H03 protocol on lung scGPT embeddings: compute sv_joint_frac and joint_AUROC at each of 12 layers. Spearman correlation. Compare sign and magnitude to immune rho=-0.93.
- *Expected signal*: rho < -0.7 in lung set; confirms cross-tissue generality of the geometry-AUROC relationship.
- *Null*: rho > -0.3 or sign reverses.
- *Value*: high | *Cost*: low

### E. Mechanistic / Attention-Level

**E1. L8 Attention Sink Hypothesis**
- *Hypothesis*: The AUROC minimum at L8 is caused by a dominant attention sink head that transiently flattens regulatory geometry.
- *Test*: Extract attention weights at L8 vs L3 and L5. Per-head attention entropy (lower = more concentrated). Identify minimum-entropy heads at L8. Compare entropy distributions across the three layers. Correlation: minimum head entropy at each layer vs joint AUROC.
- *Expected signal*: L8 has at least one near-degenerate attention head (entropy < 1 bit) not present at L3.
- *Null*: Attention entropy profiles similar at L3 and L8; no layer-AUROC correlation.
- *Value*: high | *Cost*: high

**E2. Net Displacement Conservation in 6D Space**
- *Hypothesis*: The 6D mean TF→target displacement vector's L2 norm is approximately conserved across layers (CV < 20%) even as the displacement migrates between SV2-4 and SV5-7.
- *Test*: At each layer, compute mean TRRUST-pair displacement in full 6D SV2-7. Compute L2 norm of that displacement. Track norm across 12 layers; compute CV. Also decompose: how much variance is in the SV2-4 vs SV5-7 subspaces.
- *Expected signal*: 6D displacement norm CV < 20%; per-subspace contributions trade off (one rises as other falls).
- *Null*: 6D norm varies > 40% (no conservation).
- *Value*: medium | *Cost*: low

**E3. GO Pathway Enrichment of Top vs Bottom SV2-7 Scorers**
- *Hypothesis*: TF genes with highest classifier probability in 6D SV2-7 at L3 are enriched for GO:BP terms related to transcriptional regulation (RNA Pol II, chromatin remodeling) vs bottom-ranked TFs.
- *Test*: Rank all TF genes by mean cross-seed LR probability. Run gProfiler or gseapy GO:BP enrichment on top-20 vs bottom-20. FDR < 0.05 threshold.
- *Expected signal*: Top-20 enriched for "regulation of transcription by RNA polymerase II" or chromatin-related GO terms.
- *Null*: No GO term enriched at FDR < 0.05 in either set.
- *Value*: medium | *Cost*: low

---

## Top 3 for Immediate Execution

### Priority 1: High-Probability Discovery Candidate
**A1 — TF/Target Boundary Consistency Across Seeds + A2 — TF Family Enrichment (combined script)**

Rationale: H01 validated *that* the 6D boundary exists and is seed-stable; the next urgent question is *which genes* define it and whether they have biological meaning. A1 (per-gene probability correlation across seeds) + A2 (TF family enrichment) can be computed in a single script reusing the already-validated 6D embeddings. Expected result: either a clean list of "anchor TFs" with enriched family identity, or evidence the boundary is distributed (both publishable). Cost: very low; impact: directly advances the biological interpretation section of the paper.

### Priority 2: High-Risk / High-Reward Candidate
**C2 — TRRUST Graph Laplacian Spectral Alignment with SV5-7**

Rationale: If the Laplacian eigenvectors of the TRRUST regulatory network are directly encoded in the SV5-7 gene loadings, this is a mechanistic bridge between graph topology and manifold geometry — the strongest possible interpretability claim for this project. Never tested in any prior iteration. Requires building the Laplacian (1–2 hours setup) but uses existing SVD outputs. If positive, this is the headline result of the paper.

### Priority 3: Cheap Broad-Screen Candidate
**B1 — Effective Rank as Layer Discriminability Predictor + E2 — Net Displacement Conservation (combined)**

Rationale: Both tests use existing per-layer SVD outputs with minimal additional code. B1 tests whether effective rank (the correct functional inverse of H03-iter57's energy fraction) is a universal layer descriptor, cross-seed. E2 tests whether the 6D displacement vector is conserved — mechanistically motivated by the crossover finding. Together they form a single analysis pass (~30 min) that either confirms or refutes the "richer geometry = better discrimination" principle and the conservation hypothesis. Both results are reportable regardless of sign.
