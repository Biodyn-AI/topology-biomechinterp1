# Brainstormer Hypothesis Roadmap — iter_0059
Date: 2026-02-23

---

## Retire / Deprioritize

| Direction | Reason | Action |
|---|---|---|
| eff_rank as independent AUROC predictor | Fully confounded by layer depth (partial_r=-0.045). Claim 53 retracted. | `retire_now` |
| Hub TF degree → boundary margin | r=-0.012, p=0.90. Falsified. | `retire_now` |
| Pairwise TF-target L2 distance as edge classifier (current form) | AUROC=0.565, far below threshold. Weak signal only. | `rescue_once_with_major_change` — must add asymmetric features (signed displacement, TF family, node embeddings as covariates) |
| TRRUST Laplacian spectral alignment | Near-orthogonal principal angles (~85°). No meaningful imprint. | `retire_now` |
| SV subspace rotation angle vs AUROC | Non-significant (p=0.417). No predictive structure. | `retire_now` |
| H1 Betti loops on circuit genes | Negative (circuit fewer/shorter loops than non-circuit). | `retire_now` |
| kNN transitivity | Negative (lower in real than shuffled). | `retire_now` (already retired) |
| Housekeeping gene SV1 enrichment | Inconclusive due to overlap n=7. Power too low to rescue. | `retire_now` |

---

## New Hypothesis Portfolio

### H-A: Within-Layer Sparsity as AUROC Driver
**Hypothesis**: The layer-depth AUROC decline is driven by progressive sparsification — the fraction of nonzero gene entries decreases with depth, reducing the effective information content available to the classifier.
**Test**: For each layer, compute (a) fraction of genes with nonzero embeddings, (b) mean L2-norm of nonzero genes, (c) per-gene norm variance. Regress these against AUROC across 12 layers × 3 seeds. Partial Spearman controlling for layer index.
**Expected signal**: If true, sparsity fraction will be a significant partial predictor (partial_r > 0.3) of AUROC after controlling for layer.
**Null**: Fraction nonzero is constant or not correlated with AUROC within layers.
**Value**: high | **Cost**: low

---

### H-B: TF Family-Stratified Layer AUROC Trajectory
**Hypothesis**: Different TF families peak their geometric separability at different transformer layers — STAT and Forkhead families peak early (L2-L3) matching master regulator identity, while bZIP/C2H2-ZF families with negative margins peak mid-depth.
**Test**: At each layer, for each TF family (STAT, Forkhead, ETS, bZIP, C2H2-ZF, Rel, RUNX), compute mean margin in the 6D SV2-7 space. Plot family-specific trajectories across 12 layers. One-sample t-test: is family mean margin > 0 at each layer?
**Expected signal**: STAT/Forkhead: early positive margin; bZIP/C2H2-ZF: non-positive margin throughout or flipping by depth.
**Null**: All TF families show identical layer profiles.
**Value**: high | **Cost**: low

---

### H-C: GO Biological Process Cluster Compactness in 6D
**Hypothesis**: Genes annotated to the same GO Biological Process term are significantly closer in the 6D SV2-7 subspace than random gene sets of matched size, and this compactness is layer-dependent.
**Test**: Select GO-BP terms with 10-50 member genes in the active 2039-gene set. For each term, compute mean pairwise L2 distance in 6D at each layer. Compare to 1000 bootstrap resamples of matched-size random gene sets. One-tailed z-test per term. Compute fraction of significant terms per layer.
**Expected signal**: Immune-relevant GO terms (e.g., "cytokine signaling," "T cell activation") are significantly compact at early layers where TF geometry is sharpest.
**Null**: GO-term compactness not different from random gene sets.
**Value**: high | **Cost**: medium

---

### H-D: Seed43 Crossover Anomaly — Geometry Investigation
**Hypothesis**: Seed43's failure to show the L9 SV2-4/SV5-7 crossover (iter_0056) results from a quantitatively different spectral energy distribution at late layers, not random noise. Seed43 is a mechanistically informative outlier.
**Test**: For seed43 only, at L9: compute SV energy fractions, subspace rotation angles from prior layers, TF/target mean distance in SV2-4 and SV5-7. Compare to main/seed44. Also test: does seed43 maintain the same AUROC trajectory despite lacking crossover?
**Expected signal**: Seed43 has anomalously high or low SV5-7 energy at L9, explaining the suppressed crossover while preserving AUROC.
**Null**: No systematic difference; pure sampling variance.
**Value**: medium | **Cost**: low

---

### H-E: Signed Displacement as Edge Predictor (Rescue of H02)
**Hypothesis**: Using the *signed* displacement vector (target − TF position in 6D), rather than scalar L2 distance, captures directional geometry and produces AUROC > 0.62 for TF-target edge prediction.
**Test**: For each positive and negative TRRUST edge, compute the 6D signed displacement (target_proj − TF_proj). Project onto the mean TF→target direction (from iter_0053/0055). Use this scalar projection as the edge score. AUROC vs null (1000 permutations). Also test: cosine similarity between displacement and mean direction.
**Expected signal**: AUROC > 0.62, since directional coding was established (Claims 44, 49).
**Null**: AUROC ≤ 0.565 (same as scalar distance).
**Value**: high | **Cost**: low

---

### H-F: STRING PPI Proximity in 6D
**Hypothesis**: Gene pairs with high STRING protein-protein interaction confidence scores are significantly closer in 6D SV2-7 space than pairs with low STRING confidence, independent of TRRUST regulatory proximity.
**Test**: For genes in the 2039 active set, pull STRING v12 scores (combined_score ≥ 400 vs < 150). Compute L2 distances in 6D at L2 and L8. Mann-Whitney + co-expression residualization (same approach as SV5-7 regulatory signal). Compare to TRRUST signal magnitude.
**Expected signal**: STRING signal in 6D, possibly orthogonal or additive to TRRUST signal. This would indicate the 6D subspace encodes physical interaction networks not just regulatory ones.
**Null**: STRING proximity not separable from random after co-expression control.
**Value**: high | **Cost**: medium

---

### H-G: Regulatory Motif Triangles (Feed-Forward Loops)
**Hypothesis**: Feed-forward loops (FFL: A→B, A→C, B→C) in the TRRUST network form triangular configurations in 6D space, where the intermediate TF (B) is geometrically between the master TF (A) and target (C).
**Test**: Extract all 3-gene FFLs from TRRUST (A→B directed, A→C directed, B→C directed). For each FFL, test the betweenness condition: does B project between A and C on the A→C direction vector? One-sample test against 1000 permuted FFL triples.
**Expected signal**: B position is between A and C significantly more than permuted triples (>50% of cases).
**Null**: B is random relative to A-C vector.
**Value**: high | **Cost**: medium

---

### H-H: Within-Layer Seed Variance vs eff_rank
**Hypothesis**: Although eff_rank is collinear with layer (defeating between-layer analysis), *within* a single layer, seed-to-seed AUROC variance does correlate with seed-level eff_rank differences.
**Test**: Within each layer separately, compute Pearson(seed_eff_rank, seed_AUROC) across 3 seeds. This removes the layer confound entirely. Report mean r across 12 layers and sign consistency.
**Expected signal**: Consistent positive r within layers (seeds with higher eff_rank → higher AUROC), providing partial rescue of the eff_rank→AUROC mechanistic link.
**Null**: r near zero within layers; effect was purely between-layer.
**Value**: medium | **Cost**: low

---

### H-I: Persistent Homology on 6D Active Circuit Genes (Not SV2-4 Only)
**Hypothesis**: The retraction of H1 loops in iter_0052 used SV2-4 (3D) only. In the full 6D SV2-7 subspace, circuit genes may exhibit topological loop structure (H1) that was geometrically compressed in 3D.
**Test**: For TF ∪ target genes (295 genes) at L2 and L3 (peak layers), compute Vietoris-Rips PH (maxdim=1) in 6D SV2-7 coordinates. Compare H1 bar count and lifetime to matched non-circuit genes and to 500 bootstrap resamples. Use Ripser.
**Expected signal**: More and/or longer H1 bars in circuit vs non-circuit in 6D (not seen in 3D).
**Null**: H1 structure not different from 3D result; dimension increase doesn't help.
**Value**: medium | **Cost**: medium

---

### H-J: TF Boundary Anchor GO Annotation
**Hypothesis**: The stable high-margin TFs (STAT4, BACH2, NFATC2, RUNX1, ZEB1) share biological function (immune regulation) that biologically explains their geometric extremality, while low-margin TFs are metabolically or structurally annotated.
**Test**: For the 116 TF×layer records in h01_crossseed_stability.csv, annotate each TF with top-3 GO-BP terms (manually from QuickGO/UniProt). Compute mean margin for immune-annotated vs metabolic/structural-annotated TFs at L2/L3. Fisher exact test: high-margin (top quartile) enriched for immune GO terms?
**Expected signal**: OR > 2, p < 0.05 for immune annotation in high-margin TFs.
**Null**: GO annotation independent of margin.
**Value**: medium | **Cost**: low

---

### H-K: Intrinsic Dimension Layer Profile via Two-NN on 6D Subspace
**Hypothesis**: TwoNN intrinsic dimension in the 6D SV2-7 subspace (not full 512D) follows a different layer profile than the full-embedding TwoNN, and the 6D ID trajectory predicts AUROC trajectory better than full-embedding ID.
**Test**: For each of 12 layers, compute TwoNN ID on the 295 circuit genes' 6D SV2-7 projections (and on the 2039 nonzero genes). Compare trajectory to 512D ID from iter_0047. Spearman(6D_ID, AUROC) and partial Spearman(6D_ID, AUROC | layer).
**Expected signal**: 6D ID has a non-monotone layer profile, potentially predicting the SV5-7 vs SV2-4 AUROC complement structure.
**Null**: 6D ID tracks layer depth same as full-embedding ID.
**Value**: medium | **Cost**: low

---

### H-L: Cross-Model Alignment: SV2-7 Geometry vs Geneformer
**Hypothesis**: The biological structure encoded in scGPT's 6D SV2-7 subspace (TF/target separation, directional asymmetry) partially transfers to Geneformer's representation space, indicating model-agnostic encoding of regulatory identity.
**Test**: Using matched genes (TRRUST TFs/targets present in both scGPT and Geneformer embeddings), compute CKA between scGPT SV2-7 projection and Geneformer's top-6D PCA subspace of the same genes. Test: CKA > null from random matched gene sets. Also test AUROC of Geneformer 6D for TF vs target classification.
**Expected signal**: CKA > 0.3 and Geneformer AUROC > 0.62.
**Null**: CKA ≈ 0 (orthogonal representations).
**Value**: high | **Cost**: high

---

## Top 3 for Immediate Execution

### Candidate 1 — High-Probability Discovery
**H-E: Signed Displacement as Edge Predictor**
- Rationale: Directional asymmetry (Claims 44, 49) is one of the most robust findings. The scalar distance approach (H02, AUROC=0.565) is the wrong projection. Projecting onto the established mean TF→target direction is a principled and cheap upgrade that should recover substantially more signal.
- Cost: low. Requires only re-analysis of existing h02 data with directional projection.
- Risk: low. Established direction vector exists from iter_0053/0055 data.

### Candidate 2 — High-Risk / High-Reward
**H-G: Regulatory Motif Triangles (FFL Geometry)**
- Rationale: If feed-forward loops form ordered geometric configurations, this would be the first mechanistic motif-level finding — connecting graph structure to geometric structure in a non-trivial way. No prior test of this.
- Cost: medium. Need to extract FFLs from TRRUST (small number likely), test betweenness condition.
- Risk: moderate. FFL count may be low (power issue); betweenness may require careful definition.

### Candidate 3 — Cheap Broad Screen
**H-B: TF Family-Stratified Layer AUROC Trajectory**
- Rationale: The TF family gradient is established (Claim 52: Forkhead/STAT high, bZIP/C2H2-ZF low). The natural extension is to trace each family's margin trajectory across all 12 layers. This is zero-cost (existing data), reveals layer-specific family encoding, and may identify which transformer depth processes which regulatory circuit type.
- Cost: low. Pure re-analysis of existing h01_crossseed_stability.csv + h01_tf_family_enrichment.csv data with layer resolution.
