# Brainstormer Hypothesis Roadmap — iter_0035 → iter_0036+

---

## Retire / Deprioritize

| Direction | Reason | Decision |
|-----------|--------|----------|
| STRING PPI → embedding proximity | AUROC=0.494 (iter_0031, OOV-corrected) | **retire_now** |
| TRRUST/GO co-annotation → proximity | AUROC=0.506/0.473 (iter_0035 H01) — 2nd failure | **retire_now** |
| GO enrichment at SV2 poles | 0/591 terms significant (iter_0011 H03) | **retire_now** |
| Repression anti-pole geometry | z=-1.41, 0/12 layers (iter_0013 H02) | **retire_now** |
| Monotonic B-cell community emergence (Spearman rho) | rho=0.081, null (iter_0034 H01) | **retire_now** — subsumed by kNN precision |
| Chromosomal proximity | AUROC=0.515, null (iter_0025 H03) | **retire_now** — use as negative control only |
| Hub-uncorrected Dorothea proximity | Collapses with population null (iter_0027 H03) | **retire_now** |
| PC2/PC3 T-cell/myeloid axes (n<14 markers) | Inconclusive (iter_0034 H02) | **rescue_once_with_major_change** — ≥20 markers, test PC2–PC6 |

---

## New Hypothesis Portfolio

### H-A: Expanded B-cell kNN precision@10 (n=15-20 markers)
**Hypothesis:** Expanding the B-cell marker set from 5 to 15-20 in-vocab genes (adding CD79B, IGHM, IGHD, BANK1, BLNK, BACH2, EBF1, PAX5 etc.) strengthens precision@10 signal and enables per-gene stability analysis.
**Test:** Run same kNN precision@10 protocol (k=10, 500-bootstrap null) with expanded marker set at L2, L5, L8, L11. Compare per-layer z-scores to iter_0035 baseline.
**Expected signal:** precision@10 increases with more markers (true positive rate), z-scores ≥6 at L2.
**Null:** Bootstrap null on 500 random size-matched gene sets.
**Value:** high | **Cost:** low

---

### H-B: Geneformer Cross-Model B-cell kNN Validation
**Hypothesis:** Geneformer gene embeddings also show elevated B-cell kNN precision@10, confirming that B-cell geometric coherence is a property of the scRNA-seq pretraining task class, not scGPT architecture-specific.
**Test:** Load Geneformer residual embeddings (or final layer if residual not available) for in-vocab genes. Intersect with B-cell marker set. Run kNN precision@10 at available layers. Compare z-scores to scGPT iter_0035 H02.
**Expected signal:** Geneformer precision@10 z ≥ 2.5 at ≥1 layer.
**Null:** Same bootstrap null; scGPT result as positive comparator.
**Value:** high | **Cost:** medium (artifact availability check needed)

---

### H-C: B-cell Neighborhood Functional Characterization
**Hypothesis:** The non-B-cell genes with smallest mean L2 distance to B-cell markers at L11 are functionally enriched for B-cell biology (BCR signaling, PAX5/IRF4 targets, B-lymphocyte differentiation GO terms) — i.e., the model's B-cell geometric cluster extends functionally into known B-cell biology.
**Test:** At L11, rank all 190 non-B-cell genes by mean L2 to B-cell marker centroids. Test top-20 nearest for GO BP enrichment (Fisher exact, Bonferroni) and TRRUST TF-target enrichment (PAX5, IRF4, BCL6). Compare to random 20-gene sets (bootstrap).
**Expected signal:** ≥2 significant GO terms (BCR signaling, B-cell activation); PAX5 or IRF4 targets over-represented.
**Null:** 1000 random size-20 sets from 190 non-B-cell genes; compare enrichment.
**Value:** high | **Cost:** low

---

### H-D: Layer 2 Mechanistic Probe — B-cell Displacement Anomaly
**Hypothesis:** B-cell markers undergo anomalously large displacement from L1 to L2 relative to all other adjacent-layer transitions, identifying L2 as the architecturally critical layer for B-cell geometric separation (peak kNN z=5.20 at L2 vs z≈2.5 at later layers).
**Test:** Compute per-gene L2 norm of displacement vector for each adjacent layer pair (L0→L1, L1→L2, ..., L10→L11). Test whether B-cell markers' L1→L2 displacement is anomalously large: (1) z-score vs all 195 genes' L1→L2 displacement; (2) L1→L2 displacement vs mean displacement across all 11 transitions (for B-cell genes specifically).
**Expected signal:** B-cell marker mean L1→L2 displacement is ≥2 SD above all genes' L1→L2 displacement; B-cell L1→L2 displacement is the maximum per-transition displacement for these genes.
**Null:** All genes' displacement per transition; mean displacement across non-B-cell genes.
**Value:** medium | **Cost:** low

---

### H-E: T-cell/Myeloid Axes with Expanded Markers (PC2–PC6, ≥20 markers)
**Hypothesis:** With T-cell (≥20 in-vocab markers: LCK, CD3D, CD3E, CD3G, CD8A, CD8B, GZMB, PRF1, IL7R, TCF7, etc.) and myeloid (≥15 markers: LYZ, CSF1R, CD68, ITGAM, S100A8, S100A9, FCGR3A, etc.) expanded marker lists, PC2 or PC3 at L11 separates T-cell from myeloid markers (AUROC ≥ 0.65 for at least one axis).
**Test:** Curate expanded marker lists filtered to 195-gene vocab. Test Mann-Whitney enrichment for each cell type at each pole of PC2–PC6 at L11. Bootstrap null (1000 random size-matched sets).
**Expected signal:** At least one PC clearly separates T-cell from myeloid; AUROC > 0.65 for ≥1 marker set.
**Null:** Bootstrap null; B-cell axis as positive control comparator.
**Value:** medium | **Cost:** low

---

### H-F: B-cell Sub-Manifold Intrinsic Dimensionality
**Hypothesis:** The B-cell marker set (5 in-vocab, or 15-20 expanded) occupies a lower-dimensional sub-manifold than a size-matched random gene set, quantifiable by participation ratio (PR) or TwoNN — consistent with the B-cell genes lying near a constrained line/plane in the 512-D embedding space.
**Test:** Compute PR (||v||_1^2 / ||v||_2^2 applied to the covariance eigenspectrum) for the B-cell marker embedding matrix at L11. Compare to 1000 random size-matched gene sets. Also test TwoNN estimator. Test at both L2 and L11.
**Expected signal:** B-cell PR < mean random PR by ≥ 1 SD (e.g., PR_Bcell ≈ 1-2 vs PR_random ≈ 4-8).
**Null:** Bootstrap null on 1000 random gene sets.
**Value:** medium | **Cost:** low

---

### H-G: kNN Precision Across All Marker Sets (Immune Cell Type Panel)
**Hypothesis:** kNN precision@10 enrichment is specific to B-cell markers (not a general cell-type property) at L2. Testing a panel of cell types (T-cell, myeloid, NK, plasma cell markers) reveals that B-cell geometric coherence is the dominant signal, with other cell types showing weaker or absent enrichment.
**Test:** For each cell type (B-cell, T-cell CD4, T-cell CD8, myeloid, NK, plasma cell): identify in-vocab markers (≥4 per type). Run precision@10 at L2 and L11 with same bootstrap null. Compare z-scores across cell types.
**Expected signal:** B-cell z-score highest (z≥5 at L2); T-cell and myeloid z-scores lower (1-3) or null; NK/plasma unclear.
**Null:** Same bootstrap null (500 random size-matched sets per cell type).
**Value:** high | **Cost:** low (panel screen on existing distance matrices)

---

### H-H: SV2 TRRUST Activation Co-Pole on 195 In-Vocab Genes
**Hypothesis:** TRRUST activation TF-target pairs remain co-localized in SV2 poles when restricted to the 195 in-vocab gene set, replicating iter_0011/0012 and completing the OOV-corrected validation.
**Test:** Extract SV2 per layer for 195-gene embedding matrix. Assign top/bottom 25% to poles. Test TRRUST activation pairs (195-gene vocab intersection) for co-pole Fisher enrichment vs 500 gene-label shuffles. Per-layer enrichment profile.
**Expected signal:** Fisher OR > 1.5, FDR < 0.05 at ≥6/12 layers; early-layer peak consistent with iter_0032 Dorothea result.
**Null:** Gene-label shuffle; repression pairs as internal control.
**Value:** high | **Cost:** low

---

### H-I: Persistent Homology on 195 In-Vocab Genes — H0/H1 Layer Trajectory
**Hypothesis:** H0 component count at fixed filtration radius decreases monotonically across layers (convergence), and H1 loop count peaks at mid-layers, on the 195 in-vocab gene distance matrix (extending the iter_0003 scGPT-lung result to the immune vocabulary).
**Test:** Run Ripser on 195-gene L2 distance matrices at L0, L2, L5, L8, L11 (subsampled). Track H0 count at 50th-percentile filtration radius and maximum H1 lifetime. Spearman rho(layer, metric). Feature-shuffle null (10 replicates).
**Expected signal:** H0 rho < -0.90; H1 lifetime peaks at L5–L8.
**Null:** Gaussian null (random 195-vector embeddings at same norm); feature-shuffle null.
**Value:** medium | **Cost:** medium (Ripser computation on full 195 × 512)

---

### H-J: Negative Control Validation — Gene Name Permutation Destroys B-cell Signal
**Hypothesis:** The B-cell kNN precision@10 signal disappears when gene names are randomly reassigned to embedding vectors (permuting which biology maps to which position in the manifold), confirming that the signal is learned from co-expression data, not a tokenization artifact.
**Test:** Randomly permute the 195 gene-name-to-embedding-vector assignments (1000 permutations). For each permutation, compute kNN precision@10 for "B-cell markers" (now pointing to random embeddings). Compare observed precision@10 to permutation distribution.
**Expected signal:** Observed precision@10 at L2 (0.14) is ≫ all 1000 permutation values (confirming signal is name/biology specific, not geometry artifact).
**Null:** The permutation distribution itself.
**Value:** high | **Cost:** low (no new data — permutation on existing embeddings)

---

### H-K: Dorothea Activation Decay Profile — OOV-Corrected Full Layer Curve
**Hypothesis:** OOV-corrected Dorothea activation pair proximity (AUROC ~0.60 at early layers, confirmed iter_0032) decays monotonically to AUROC ~0.50 by L11, with the decay rate distinguishable from random fluctuation, characterizing how regulatory geometry is progressively "washed out" by later transformer layers.
**Test:** Run per-layer AUROC for Dorothea activation pairs (OOV-corrected, 195-gene vocab) at all 12 layers. Fit linear decay model. Compare decay slope to bootstrap null (random pair AUROC per layer).
**Expected signal:** Spearman rho(layer, AUROC) < -0.80; AUROC > 0.55 at L0–L4, ~0.50 at L9–L11.
**Null:** 1000 random pair sets of same size, same per-layer AUROC distribution.
**Value:** medium | **Cost:** low

---

### H-L: Cross-Layer CKA for 195 In-Vocab Genes (Replication at Immune Gene Set)
**Hypothesis:** The near-identity CKA finding from iter_0004 (cross-layer CKA ≈ 1.0 for scGPT lung embeddings) replicates for the 195 in-vocab immune gene set, confirming residual stream stability as a model-level property rather than a gene-set-level artifact.
**Test:** Compute all 12×12 linear CKA pairs for 195-gene embedding matrix per layer. Compare adjacent CKA to non-adjacent and to feature-shuffle null. Also compare CKA curve shape (layer similarity vs distance) against the earlier lung-derived curve.
**Expected signal:** Adjacent-layer CKA ≥ 0.99 for all pairs; same structure as iter_0004 result.
**Null:** Feature-shuffle null (same as iter_0004 protocol).
**Value:** low | **Cost:** low (confirmation/replication)

---

## Top 3 for Immediate Execution

### #1 — High-Probability Discovery: Multi-Cell-Type kNN Precision Panel + Expanded B-cell Markers (H-G + H-A combined)
Run kNN precision@10 for a panel of 5-6 cell types (B-cell expanded to n=15-20 in-vocab, T-cell CD4/CD8, myeloid, NK, plasma) in one script using the existing 195-gene distance matrices. This simultaneously: (a) strengthens the B-cell finding with more markers, (b) establishes specificity by testing whether other cell types show similar or weaker enrichment, (c) provides a clean comparative figure for the paper.
- **Expected outcome:** B-cell z ≥ 6 at L2 with expanded markers; T-cell/myeloid z < 2 or null (OR shows specificity).
- **Risk:** Low. All runs on existing embeddings. Marker curation is the main effort.

### #2 — High-Risk/High-Reward: Geneformer Cross-Model B-cell Validation (H-B)
If Geneformer embedding artifacts are accessible, test B-cell kNN precision@10 on Geneformer. This is the highest-impact single experiment available: a positive result would allow the paper to claim that B-cell geometry reflects the scRNA-seq pretraining objective, not scGPT-specific architecture.
- **Expected outcome:** Geneformer precision@10 z ≥ 2.0 at ≥1 layer.
- **Risk:** High. Requires artifact availability check first; if unavailable, pivot to H-J (negative control permutation).

### #3 — Cheap Broad Screen: Negative Control Permutation + B-cell Neighborhood Characterization (H-J + H-C combined)
Both experiments run on existing embeddings with no new data loading. H-J closes the critical "is this a tokenization artifact?" question. H-C characterizes what the B-cell cluster actually represents biologically. Together they complete the evidence triangle: (signal is real) + (signal is biology-specific) + (signal has biological meaning).
- **Expected outcome:** Permutation confirms B-cell z=5.20 is in the top 0.1% of name-permuted values; top-20 B-cell neighbors enriched for BCR signaling or PAX5 targets.
- **Risk:** Low. Pure computation on existing artifacts.
