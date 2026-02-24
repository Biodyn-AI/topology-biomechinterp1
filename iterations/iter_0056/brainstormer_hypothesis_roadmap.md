# Brainstormer Hypothesis Roadmap — iter_0056

**Date**: 2026-02-23

---

## Retire / Deprioritize

| ID | Description | Reason |
|----|-------------|--------|
| H03-iter56 | Subspace rotation angle as AUROC predictor | Definitively negative (rho=-0.27, p=0.42); retire permanently |
| STRING PPI → proximity | (retired iter_0031) | — |
| Annotation/GO/co-reg → embedding distance | (retired iter_0032) | — |
| TDA loops on circuit genes | (retired iter_0052) | — |
| BCL6 metabolic divergence | (retired iter_0041) | — |
| SV2-4 as regulatory axis | (revised to co-expression iter_0051) | — |

**Rescue-once candidates** (allowed one more test with substantial change):
- SV8-14 secondary signal: previously seen weak signal at L8 only; rescue by testing SV8-14 in joint classifier (SV2-14) vs SV2-7, measuring whether higher-order SVs add discriminative power.

---

## New Hypothesis Portfolio

### A. Direct Follow-Ups (High Priority)

**A1. Joint 6D Cross-Seed Validation**
- *Hypothesis*: The joint SV2-7 AUROC profile (mean=0.744) replicates in seeds 43 and 44 with consistent layer-wise pattern.
- *Test*: Run the same iter_0056 H01 protocol on seed43 and seed44 embeddings. Extract SV2-7, concatenate, 5-fold LR, 100-perm null. Compare mean AUROC and layer-profile shape across 3 seeds.
- *Expected signal*: Seeds match main ± 0.03; all three show AUROC ≥ 0.70 at most layers.
- *Null*: AUROC ≤ 0.55 or high seed-to-seed variance (SD > 0.05 mean cross-layer).
- *Value*: high | *Cost*: low

**A2. Seed43 SV Basis Alignment (Crossover Failure Investigation)**
- *Hypothesis*: Seed43's SV5-7 right singular vectors (gene loading space) are aligned with main's SV2-4 (i.e., the labeling of spectral subspaces is permuted), which would explain why sv57_mag never exceeds sv24_mag.
- *Test*: At each layer, compute principal angles between seed43's SV2-4 loading matrix and main's SV5-7 loading matrix. If principal angles are small (< 30°), the subspaces are aligned differently (permuted). Also compare seed43 SV2-4 AUROC and SV5-7 AUROC trajectories to main.
- *Expected signal*: Mean principal angle between seed43 SV2-4 and main SV5-7 < 30° at early layers, confirming basis permutation.
- *Null*: Principal angles uniform ~55-60° (random alignment in 512D space).
- *Value*: high | *Cost*: low

**A3. Subspace Mutual Information Across Layers**
- *Hypothesis*: SV2-4 and SV5-7 representations of the same genes decouple over depth — low MI at early layers (independent encoding) and higher MI at late layers (convergence at integration point L9).
- *Test*: At each layer, estimate MI between 3D SV2-4 projection and 3D SV5-7 projection of nonzero genes using kNN-based estimator (e.g., MINE or KSG). Track MI across 12 layers.
- *Expected signal*: MI low L0-L3, rises L4-L9, providing a continuous measure of the "integration" that Claim 49 posits at L9.
- *Null*: MI flat or random across layers.
- *Value*: medium | *Cost*: medium

### B. Topological / Geometric Depth

**B1. Persistent Homology on SV2-7 TF vs Target Gene Sets**
- *Hypothesis*: TF genes form a topologically distinct cluster in 6D SV2-7 space (lower H0 fragmentation, distinct H1 structure) compared to target-only genes, detectable via Vietoris-Rips PH.
- *Test*: At L3 (peak AUROC 0.789): project 73 TF genes and 222 target-only genes into 6D SV2-7. Run separate PH (maxdim=1) on each set. Compare: (a) number of H0 components at given epsilon, (b) H1 bar count and mean lifetime. 100-bootstrap null (random label assignment of same N).
- *Expected signal*: TFs have fewer H0 components (tighter cluster) and shorter/fewer H1 bars (more planar arrangement); targets more fragmented.
- *Null*: No PH-profile difference between TF and target gene sets.
- *Value*: high | *Cost*: medium

**B2. Effective Rank / Intrinsic Dimension of TF vs Target Gene Subpopulations**
- *Hypothesis*: TF genes occupy a lower-dimensional manifold in residual space than target-only genes, reflecting the concentrated functional role of TFs.
- *Test*: TwoNN intrinsic dimension separately for TF genes (n=73) and target-only genes (n=222) in full 512D embedding at each of 12 layers. Track ID ratio (TF/target) across depth.
- *Expected signal*: TF ID < target ID across most layers; ratio compresses or reverses at integration point.
- *Null*: TF ID ≈ target ID at each layer (no manifold distinction).
- *Value*: medium | *Cost*: low

**B3. TRRUST Graph Laplacian Spectral Alignment with SV5-7**
- *Hypothesis*: The top eigenvectors of the TRRUST regulatory graph Laplacian are geometrically encoded in the SV5-7 embedding subspace; cosine alignment between graph eigenvectors and SV5-7 gene loadings should be significantly above null.
- *Test*: Build the TRRUST interaction graph over the 295 in-vocab named genes. Compute graph Laplacian, extract top-5 eigenvectors (over gene dimensions). Compute cosine similarity between each eigenvector and SV5-7 right singular vectors at L0 (where regulatory signal peaks). Permutation null: shuffle gene identities 1000x.
- *Expected signal*: One or more Laplacian eigenvectors achieves cosine ≥ 0.15 with SV5-7 directions, significantly above permutation null.
- *Null*: Max cosine ≤ random; permutation p > 0.1.
- *Value*: high | *Cost*: medium

**B4. Filtration Stability of Regulatory Geometry (Persistence Landscape)**
- *Hypothesis*: The persistence landscape of TRRUST circuit genes in SV5-7 is stable across seeds (low Wasserstein distance between seed landscapes), while non-circuit gene landscapes are more variable.
- *Test*: Compute persistence landscapes for TRRUST 295 genes in SV5-7 at L0 for main, seed43, seed44. Compute pairwise Wasserstein distances (H1) between seeds. Compare to non-circuit gene control (random 295 genes, 100 bootstrap).
- *Expected signal*: Circuit gene landscape inter-seed Wasserstein distance significantly smaller than non-circuit gene landscape variance.
- *Null*: No stability advantage for circuit gene landscapes.
- *Value*: medium | *Cost*: medium

### C. Cross-Architecture / Generalization

**C1. SV5-7 Regulatory Geometry in Geneformer (Cross-Model)**
- *Hypothesis*: The early-layer SV5-7 regulatory proximity signal (Claim 42) is not scGPT-specific but appears in Geneformer residual embeddings for the same immune genes at the corresponding early layers.
- *Test*: Extract Geneformer residual embeddings for matched TRRUST immune genes. Run the same SVD procedure, extract SV5-7 at early layers, compute rbc (TRRUST pairs vs null). Compare to scGPT result at matched depth fractions (e.g., L0 = first layer regardless of absolute depth).
- *Expected signal*: Geneformer SV5-7 rbc > 0.10 at early layers, confirming architectural generality.
- *Null*: Geneformer SV5-7 rbc ≤ 0.05 or not significant.
- *Value*: high | *Cost*: high

**C2. Non-Linear Classification on 6D SV2-7 (Kernel SVM)**
- *Hypothesis*: The TF/target boundary in 6D SV2-7 space is non-linear; kernel SVM outperforms logistic regression by ≥ 0.03 AUROC at peak layers.
- *Test*: At L3 and L0 (peak layers), compare LR vs RBF-kernel SVM (5-fold stratified CV, 100-perm null) on the 6D joint representation.
- *Expected signal*: SVM AUROC > 0.81 at L3, suggesting non-linear separability.
- *Null*: SVM AUROC ≤ LR AUROC + 0.02.
- *Value*: medium | *Cost*: low

**C3. OOD Generalization: Train on Immune, Classify on Lung**
- *Hypothesis*: A classifier trained on immune scGPT SV2-7 generalizes to lung embeddings for TF/target prediction, indicating the 6D regulatory geometry is cell-type-agnostic.
- *Test*: Train LR on cycle4_immune SV2-7 at L3. Apply to lung embeddings (same scGPT, matched gene set, TF/target labels from TRRUST). Measure AUROC on lung set. Compare to lung-trained LR (upper bound) and random (lower bound).
- *Expected signal*: Cross-tissue AUROC > 0.60, significantly above chance.
- *Null*: Cross-tissue AUROC ≤ 0.55 (geometry is tissue-specific).
- *Value*: high | *Cost*: medium

### D. Mechanistic / Attention-Level

**D1. L8 Bottleneck Mechanism — Attention Sink Hypothesis**
- *Hypothesis*: The AUROC minimum at L8 (SV5-7: 0.543, SV2-4: 0.732) coincides with a dominant attention sink (one attention head captures > 50% of total attention weight), which flattens regulatory geometry transiently.
- *Test*: Extract attention weights at L8 vs L3 for the immune embedding run. Compute per-head attention entropy; identify minimum-entropy heads. Compute correlation between attention sink strength and AUROC drop. Compare to random layer (L5).
- *Expected signal*: L8 has at least one near-degenerate attention head (entropy < 1 bit) absent at L3.
- *Null*: Attention entropy profiles similar at L3 and L8.
- *Value*: high | *Cost*: high

**D2. SV5-7 Singular Value Magnitude as Layer-Wise Information Proxy**
- *Hypothesis*: The sum of singular values S5+S6+S7 (variance explained by the regulatory subspace) tracks the SV5-7 AUROC trajectory across layers, providing a data-free proxy for regulatory encoding strength.
- *Test*: Compute sum(S[4:7]) at each of 12 layers. Spearman rho with SV5-7 AUROC trajectory. Compare to SV2-4 (sum S[1:4] vs SV2-4 AUROC).
- *Expected signal*: Spearman rho > 0.6 for SV5-7 singular value sum vs AUROC; provides interpretable variance-based predictor.
- *Null*: rho < 0.3 (singular value magnitude does not track AUROC).
- *Value*: medium | *Cost*: low

**D3. Net Displacement Conservation in SV2-7 Space**
- *Hypothesis*: The net directional signal (TF→target mean displacement) in the full 6D SV2-7 space is conserved across layers even as its distribution across SV2-4 vs SV5-7 shifts (the crossover).
- *Test*: Compute the full 6D mean displacement vector for valid TRRUST pairs at each layer. Compute its L2 norm. Compare to the sum of sv24_mag + sv57_mag. Plot conservation ratio across 12 layers.
- *Expected signal*: 6D norm is approximately conserved (CV < 20% across layers), suggesting the regulatory direction is encoded in the joint space even as it migrates between subspaces.
- *Null*: 6D norm varies substantially (CV > 40%), no conservation.
- *Value*: medium | *Cost*: low

---

## Top 3 for Immediate Execution

### Priority 1: High-Probability Discovery Candidate
**A1 — Joint 6D Cross-Seed Validation**

Rationale: Claim 50 is the strongest new claim in the paper. Cross-seed validation is the direct path to hardening it from "promising" to "validated." All data and code are available; cost is minimal. Expected result: 3-seed AUROC agreement, advancing Claim 50 to fully validated status.

### Priority 2: High-Risk / High-Reward Candidate
**B3 — TRRUST Graph Laplacian Spectral Alignment with SV5-7**

Rationale: If the regulatory graph's spectral structure is literally encoded in the embedding subspace, this would be the deepest mechanistic finding in the project — a direct correspondence between graph topology and manifold geometry. Never tested; requires building the Laplacian and running alignment. If positive, this is a paper-defining result.

### Priority 3: Cheap Broad-Screen Candidate
**A2 + D2 combined** — Seed43 SV Basis Alignment Investigation + SV Singular Value as AUROC Proxy

Rationale: Both tests reuse existing data and the crossover CSV. A2 resolves the Claim 49 conditionality question. D2 provides a data-free proxy for AUROC that, if confirmed, is a useful analytical tool for future iterations. Combined, they are one short script and one CSV analysis — cost < 1 hour compute.
