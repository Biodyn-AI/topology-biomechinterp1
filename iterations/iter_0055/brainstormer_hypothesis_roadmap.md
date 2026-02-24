# Brainstormer Hypothesis Roadmap — iter_0055

**Date**: 2026-02-23
**Research phase**: Consolidation + mechanistic dissection of dual-subspace regulatory architecture

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| H1 Betti loops on circuit genes in SV2-4 | **RETIRED** | Negative result iter_0052, no loops found |
| OOV bifurcation narrative (14-gene cluster) | **RETIRED** | Artifact, iter_0030 |
| TCR circuit attractor | **RETIRED** | Insufficient vocabulary (2-3 genes), underpowered |
| Lineage centroid cosine orthogonality | **RETIRED** | Indistinguishable from null due to SV1 dominance, iter_0046 |
| BCL6 metabolic divergence as targeted attractor | **RETIRED** | Co-movement confound, iter_0044 |
| Dorothea TF-target proximity (hub-corrected) | **DEPRIORITIZED** | Requires matched within-TF null that has never been run; low rescue odds given co-expression confound history |
| Cross-lineage TwoNN compression law | **DEPRIORITIZED** | Methodological incompatibility between estimators shown iter_0044; revisit only with global-manifold estimator |
| SV2-4 as direct regulatory proximity encoder | **RETIRED** | Explained by co-expression, confirmed iter_0051 |

---

## New Hypothesis Portfolio

### 1. Joint SV2-7 (6D) Layer-Stable Classifier
**Hypothesis**: Combining SV2-4 and SV5-7 into a 6D joint feature vector yields a logistic regression AUROC that is higher and more uniform across all 12 layers than either subspace alone.
**Test**: At each of 12 layers, train 5-fold stratified logistic regression on 6D joint projection (SV2-7 coordinates of 2039 nonzero genes). Compare layer profile to H02 results (SV5-7 and SV2-4 separately). 100-perm null.
**Expected signal if true**: Joint AUROC ≥0.74 across most layers with <0.05 variance; "fills in" the L8 SV5-7 minimum.
**Null/control**: Permuted gene labels (same null as H01/H02 iter_0055).
**Value**: high | **Cost**: low

### 2. L9 Crossover — Directionality Cross-Seed Replication + Shuffle Null
**Hypothesis**: The L9 SV5-7 > SV2-4 directionality crossover replicates on seed43 and seed44, and both displacement magnitudes exceed label-shuffle null at every layer.
**Test**: Run the H03 iter_0055 displacement protocol (mean TF→target displacement magnitude per layer) on seed43 and seed44 embeddings. Add 500-label-shuffle null for magnitude. Report crossover layer per seed.
**Expected signal if true**: Crossover at L8-L10 for all seeds; magnitude >null at all layers.
**Null/control**: Label-shuffle null (targets assigned to random non-TF genes).
**Value**: high | **Cost**: low

### 3. GO/Biological Enrichment of L9 Top-Displacement Genes
**Hypothesis**: Genes that move most in SV5-7 between L8 and L9 (i.e., drive the directionality amplitude peak) are enriched for specific biological processes, potentially signaling pathway components.
**Test**: For each gene, compute |SV5-7 embedding vector at L9| - |SV5-7 embedding vector at L8|. Rank all 2039 genes by this delta. Run mygene GO BP enrichment on top-50 and bottom-50. Test whether TF fraction differs between top/bottom.
**Expected signal if true**: Top-moving genes enriched for transcriptional regulation, signal transduction, or immune activation terms; TF fraction >background.
**Null/control**: Random same-size gene sets.
**Value**: high | **Cost**: low

### 4. SV5-7 Directionality Vector Biological Identity
**Hypothesis**: The mean displacement vector (mean of target-TF in SV5-7) at L0 points in a biologically meaningful direction — i.e., genes co-aligned with this vector share a functional annotation (e.g., target-like properties).
**Test**: Compute mean displacement unit vector at L0 (from H03). For each gene, compute its scalar projection onto this vector. Rank all 2039 genes. Test GO enrichment of top/bottom 50. Compare TF fraction in high vs low projection groups.
**Expected signal if true**: High-projection genes enriched for regulatory targets; low-projection for TFs; specific GO terms associated with either pole.
**Null/control**: Random 50-gene sets; gene-label shuffled direction vector.
**Value**: high | **Cost**: low

### 5. L8 SV5-7 AUROC Minimum — Gene-Level Anatomy
**Hypothesis**: The SV5-7 AUROC minimum at L8 (0.543, near-random) reflects a specific set of TF genes that have migrated out of the SV5-7 regulatory subspace at L8, not a global blurring.
**Test**: At each of L0, L3, L8, L11: extract SV5-7 coordinates of TF genes (n=73) and target-only genes (n=222). Compare distribution overlap metrics (KL divergence, Mahalanobis distance). Identify which TF genes are misclassified at L8 but classified correctly at L0/L9.
**Expected signal if true**: A stable subset of TFs (~20-30) shifts to target-like SV5-7 coordinates at L8, driving the AUROC drop; these TFs may have specific functional properties (e.g., mid-layer active TFs, broadly expressed TFs).
**Null/control**: Same set of misclassified vs correctly classified genes at L0 (baseline).
**Value**: high | **Cost**: medium

### 6. SV5-7 Early-Layer Encoding — Within-TF-Family Null
**Hypothesis**: The SV5-7 TF/target distinction at L0 is not driven purely by TF-family co-expression clusters but survives a within-TF-family null where negative pairs are drawn from the same TF's target set.
**Test**: For each TRRUST TF with ≥5 named targets, construct within-TF positive (TF→target) and within-TF negative (target-target from same TF) pairs. Compare SV5-7 distances. 100-perm null.
**Expected signal if true**: Even within same TF regulon, TF is geometrically separated from its targets in SV5-7, ruling out co-expression family confound.
**Null/control**: Within-TF label-shuffled pairs.
**Value**: high | **Cost**: medium

### 7. Mutual Information / Cross-Subspace Information Routing
**Hypothesis**: The SV5-7 coordinate of a gene at L0-L2 is predictive of its SV2-4 coordinate at L5-L8 (information routing hypothesis: early SV5-7 encoding rotates into SV2-4 encoding at mid-depth).
**Test**: For n=2039 genes, compute Spearman rank correlation between SV5-7 coordinates at L0 and SV2-4 coordinates at L6-L8. Test whether TRRUST gene pairs maintain correlated cross-subspace coordinates (i.e., pairs close in SV5-7 at L0 are close in SV2-4 at L6).
**Expected signal if true**: Positive Spearman correlation for same-pair distances across subspaces; AUROC > 0.5 for SV5-7-L0 distance predicting SV2-4-L6 distance.
**Null/control**: Randomized gene-index pairing between L0-SV5-7 and L6-SV2-4 coordinates.
**Value**: high | **Cost**: medium

### 8. Persistent Homology on SV5-7 Circuit Genes at L0 vs L9
**Hypothesis**: TRRUST circuit genes in SV5-7 space have distinct persistent homology (H0/H1) topology at L0 (regulatory encoding peak) vs L9 (directionality peak), compared to L8 (AUROC minimum).
**Test**: Extract SV5-7 coordinates of 295 TRRUST circuit genes at L0, L8, L9. Run Ripser (maxdim=1). Compare H0 and H1 diagrams: number of components, persistence lifetimes, Wasserstein distance between diagrams.
**Expected signal if true**: L0 and L9 show more compact/structured topology (fewer long-lived H0 bars = fewer clusters) than L8 where circuit gene geometry is maximally diffuse.
**Null/control**: Size-matched non-circuit gene sets at each layer.
**Value**: medium | **Cost**: medium

### 9. SV8-SV14 Regulatory Proximity — Layer Profile
**Hypothesis**: The weak secondary regulatory signal in SV8-14 (iter_0053) follows a distinct layer profile from both SV5-7 and SV2-4, providing a third independent encoding regime.
**Test**: At each of 12 layers, compute co-expression-residualized rbc for SV8-10 and SV11-14 (TRRUST positive vs negative). Compare layer profile to SV5-7 (early-dominant) and SV2-4 (L8-dominant).
**Expected signal if true**: SV8-14 signal peaks at L6-L7 (intermediate), distinct from both primary subspaces.
**Null/control**: SV21-30 (iter_0053 showed near-zero residualized rbc).
**Value**: medium | **Cost**: low

### 10. Cross-Model SV5-7 Test in Geneformer
**Hypothesis**: Geneformer's layer representations contain an analogous early-layer subspace where TF vs target-only AUROC is elevated, despite the Geneformer architecture being different from scGPT.
**Test**: Load Geneformer immune gene embeddings (if available as contextual layer representations). SVD of centered gene embeddings. Test SV5-7 TF/target AUROC at each layer.
**Expected signal if true**: Geneformer shows elevated AUROC in some SV range at early layers, suggesting cross-model architectural convergence.
**Null/control**: Permuted gene labels; SV1 as distal control.
**Value**: high (if positive) | **Cost**: medium (data access dependent)

### 11. SV5-7 Subspace Rotation Angle Trajectory
**Hypothesis**: The SV5-7 subspace (as a 3D linear subspace) rotates continuously across layers, and the rotation rate peaks at L8 (where AUROC is minimum), explaining the AUROC trajectory as a subspace misalignment effect.
**Test**: Compute the principal angle between the SV5-7 subspace (spanned by left singular vectors u4,u5,u6) at consecutive layers. Report principal angles L0→L1, L1→L2, ..., L10→L11. Test correlation with AUROC trajectory.
**Expected signal if true**: Principal angle largest at L7→L8 or L8→L9 transition; anticorrelation with AUROC trajectory (more rotation = lower AUROC).
**Null/control**: SV2-4 subspace rotation angles as comparison; random Gaussian matrices as absolute null.
**Value**: high | **Cost**: low

### 12. TF-Family Geometric Diversity in SV5-7 at L0
**Hypothesis**: Different TF families (AP1, ETS, IRF, RUNX, KLF, STAT, SP1) occupy geometrically distinct regions in SV5-7 space at L0, and TF targets cluster near their cognate TF family.
**Test**: At L0, extract SV5-7 coordinates for all 73 TFs and 222 target-only genes. Use TRRUST to label TFs by family. AUROC for TF-family separation (LR multi-class). Compute TF→nearest-target-in-SV5-7 precision@k.
**Expected signal if true**: TF families form distinct clusters in SV5-7; targets cluster near cognate TF family at L0.
**Null/control**: Random TF-family label shuffle; compare to SV1 as distal control.
**Value**: high | **Cost**: medium

---

## Top 3 for Immediate Execution

### [1] HIGH-PROBABILITY DISCOVERY CANDIDATE
**Hypothesis 1: Joint SV2-7 (6D) Layer-Stable Classifier**

Rationale: Direct synthesis of the last two iterations. H02 showed complementary AUROC trajectories — joining both subspaces should yield a layer-flat, high-AUROC classifier that constitutes the strongest single claim in the paper. Computationally trivial (add 6D feature vector). Very likely positive given the complementarity already demonstrated.

Primary metric: AUROC at each of 12 layers for joint 6D LR classifier (SV2-7 combined). Compare to SV5-7 alone and SV2-4 alone.

### [2] HIGH-RISK / HIGH-REWARD CANDIDATE
**Hypothesis 11: SV5-7 Subspace Rotation Angle Trajectory**

Rationale: Explains the mechanistic origin of the AUROC minimum at L8 and the directionality crossover at L9. If subspace rotation rate peaks at L7-L8 and anticorrelates with AUROC, it establishes a rotation-AUROC correspondence that is mechanistically novel and publishable. High-risk because the correlation may not hold, but computationally cheap.

Primary metric: Spearman correlation between inter-layer SV5-7 principal angles and SV5-7 AUROC trajectory.

### [3] CHEAP BROAD-SCREEN CANDIDATE
**Hypothesis 2: L9 Crossover Cross-Seed Replication + Shuffle Null**

Rationale: The crossover at L9 is the most striking result of iter_0055 but has only been seen on the main seed. Replication on seed43/seed44 is required before any publication claim. Also adds a label-shuffle null for displacement magnitude. Cost is low (same code as iter_0055, applied to existing seed embeddings). Gates whether Claim 49 can be promoted to a primary claim.

Primary metric: Crossover layer on seed43 and seed44; displacement magnitude vs shuffle null p-value at all 12 layers.
