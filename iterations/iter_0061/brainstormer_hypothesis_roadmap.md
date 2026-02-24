# Brainstormer Hypothesis Roadmap — iter_0061

**Date**: 2026-02-23

---

## Retire / Deprioritize

| Hypothesis Family | Reason | Decision |
|-------------------|--------|----------|
| FFL geometric ordering (betweenness/interpolation t) | Permutation null ≈ 0.50 by symmetry; two negative permutation tests; N=264 still not significant | **retire_now** |
| L8 CKA transition boundary | Consecutive CKA 0.977-0.981 uniformly; no drop at L7-L8 | **retire_now** |
| Dual-role gene identity ambiguity | All 38 consistently target-like; HIF1A flip is n=1 outlier | **retire_now** |
| Effective rank as causal AUROC predictor | Partial Spearman r=-0.045 after controlling for layer depth | **retire_now** |
| Directional projection vs scalar distance | Identical AUROC (0.563 vs 0.565); no gain from direction | **retire_now** |
| Parametric geometric ordering tests without permutation | Systematic null bias; all such tests must be re-done with permutation | **methodology ban** |

---

## New Hypothesis Portfolio

### TOPOLOGY

**T1: Deep-Layer (L10-L11) Persistent Homology Divergence**
- Hypothesis: L10-L11 have higher H1 Betti (topological holes) than L0-L4, consistent with cross-seed CKA divergence indicating less-constrained geometry.
- Test: Ripser on PCA-20 projections of 400 sampled genes at each layer for cycle4_immune (main, seed43, seed44). Compare H0/H1 Betti counts and persistence sums. Permutation null: gene-position shuffle. Cross-seed variance of Betti numbers.
- Expected signal: H1 Betti increases L4→L11; cross-seed Betti variance higher at L10-L11.
- Null: No difference between early and late layers in Betti counts.
- Value: **high** | Cost: **medium**

**T2: Topological Stability Profile — Which Layer Is Most Robust?**
- Hypothesis: The layer with maximum cross-seed topological consistency (minimum Wasserstein distance between persistence diagrams of different seeds) corresponds to the layer of peak biological signal.
- Test: At each layer, compute H1 persistence diagram for each of 3 seeds. Compute pairwise Wasserstein distance. Compare to per-layer AUROC (from T3/H3 below).
- Expected signal: Topological stability peak co-occurs with AUROC peak (L2-L3 or L6-L9).
- Null: No correlation between Wasserstein stability and AUROC layer profile.
- Value: **high** | Cost: **medium**

**T3: cycle4_immune Filtration Variant — VR vs Alpha Complex**
- Hypothesis: Alpha complex filtration (density-aware) reveals different topological features than VR filtration for gene embeddings, particularly in distinguishing TF vs target submanifolds.
- Test: Use GUDHI to compute alpha complex persistence for TF genes vs target genes separately at L2-L3. Compare Betti profiles.
- Expected signal: TF submanifold has different H1 count than target submanifold.
- Null: Permuted TF/target labels give same Betti difference.
- Value: **medium** | Cost: **medium**

---

### MANIFOLD GEOMETRY

**M1: cycle4_immune TF-Target AUROC at Full Scale (735 pairs)**
- Hypothesis: With 735 TRRUST pairs (vs 288 in cycle1), TF-target SV5-7 AUROC exceeds 0.70 at peak layer, establishing the signal as robust and biologically interpretable.
- Test: Run SV5-7 LR classification + pairwise L2 distance AUROC for cycle4_immune [12, 4941, 512] with 735 positive + 2205 negative pairs. All 12 layers. Permutation null (2000 shuffles).
- Expected signal: AUROC ≥ 0.72 at L2-L3; stable across seeds.
- Null: Permuted gene-to-role labels; AUROC ≈ 0.50.
- Value: **high** | Cost: **low**

**M2: Layer-Depth AUROC Mechanism — Sparsity vs Biological Separation**
- Hypothesis: The per-layer AUROC rise (L0→L9 in cycle1) is driven by sparsification (fewer active genes), not learned biological structure. Partial correlation of AUROC on sparsity (n_nonzero genes) after controlling for layer should reveal whether sparsity is an independent driver.
- Test: For each of 12 layers × 3 seeds: (a) n_nonzero genes per layer, (b) SV5-7 AUROC, (c) fraction variance in top SV. Partial Spearman(AUROC, sparsity | layer) and Partial Spearman(AUROC, layer | sparsity).
- Expected signal: If biological: partial r(AUROC, layer | sparsity) > 0.5. If confound: partial r ≈ 0.
- Null: Theoretical — if AUROC tracks sparsity, it's an artifact.
- Value: **high** | Cost: **low** (pure reanalysis of existing data)

**M3: Intrinsic Dimensionality at L10-L11 vs L0-L4 (TwoNN, cycle4_immune)**
- Hypothesis: L10-L11 embeddings have higher intrinsic dimensionality than L0-L4 for cycle4_immune, opposite to the cycle1/lung trend (iter_0004, which showed slight decrease). This would confirm that deep immune layers are more "spread out."
- Test: TwoNN estimator (Facco 2017) at each layer, 400 genes, 15 null replicates per layer.
- Expected signal: TwoNN ID at L10-L11 > TwoNN ID at L0-L4 for cycle4.
- Null: Feature shuffle; expected TwoNN ≈ 20 (noise).
- Value: **medium** | Cost: **low**

**M4: Gene Embedding Trajectory Curvature — Peak Curvature Layer**
- Hypothesis: Individual gene embeddings follow non-straight trajectories through layer space, with curvature distribution non-uniformly peaking at L7-L9 (the transition zone suggested by HIF1A flip and cross-seed CKA divergence onset).
- Test: For 500 randomly sampled genes, compute 12-layer trajectory in SV5-7. Discrete curvature = angle between consecutive displacement vectors (L_n-1→L_n vs L_n→L_n+1). Identify each gene's peak-curvature layer. Test if distribution is non-uniform (chi-squared test).
- Expected signal: Modal peak-curvature at L7-L9; distribution significantly non-uniform.
- Null: Uniform distribution across layers.
- Value: **medium** | Cost: **low**

---

### CROSS-MODEL / CROSS-CYCLE

**C1: Cross-Cycle Classifier Transfer (cycle4→cycle1)**
- Hypothesis: A TF-target classifier trained on cycle4_immune SV5-7 (larger, more constrained) transfers to cycle1_main SV5-7 after Procrustes alignment with AUROC > 0.65, demonstrating universal geometric structure.
- Test: Train LR on cycle4 L3 SV5-7 (using cycle4 TRRUST genes). Procrustes-align cycle1 L3 SV5-7 to cycle4. Apply classifier to cycle1. Compare transfer AUROC vs within-cycle1 AUROC. Permutation null: random rotation of cycle4 classifier.
- Expected signal: Transfer AUROC ≥ 0.65 (vs within-cycle1 ≈ 0.744).
- Null: Random rotation control shows transfer AUROC ≈ 0.50.
- Value: **high** | Cost: **medium**

**C2: Seed-Consensus TF-Target Subspace (Intersection Geometry)**
- Hypothesis: The TF-target separation direction common to all 3 seeds has higher AUROC than any single seed, because seed-invariant directions capture biologically determined structure.
- Test: At L2-L3, compute LR weight vector (normalized) per seed. Average (consensus). Test AUROC of consensus classifier vs individual seeds. Permutation: consensus of shuffled labels.
- Expected signal: Consensus AUROC ≥ 0.76 (≥ best single seed 0.762).
- Null: Shuffled-label consensus AUROC ≈ 0.50.
- Value: **high** | Cost: **low**

---

### BIOLOGICAL ANCHORING

**B1: HIF1A Regulatory Neighborhood Co-Transition at L8**
- Hypothesis: HIF1A's TRRUST targets (genes it regulates) also shift toward TF-like margins at L8-L11, parallel to HIF1A itself, indicating the embedding captures the regulatory program propagation.
- Test: Identify HIF1A's TRRUST targets. Track LR margin trajectory for each target across 12 layers. One-sample t-test: mean(margin_L8:L11) - mean(margin_L0:L7) > 0. Compare to random same-size gene sets.
- Expected signal: HIF1A target margins shift positive (Δ > 0.3) at L8.
- Null: Random same-size gene set shows no L8 shift.
- Value: **medium-high** | Cost: **low** (reanalysis of existing embeddings)

**B2: GO Biological Process Module Compactness**
- Hypothesis: Genes co-annotated to immune-relevant GO BP terms (immune response, cytokine production, lymphocyte activation) cluster more tightly than random at some layer, with peak at L2-L6.
- Test: For each of 10 selected immune GO terms with ≥15 genes in the 4941-gene cycle4 pool: compute mean pairwise L2 distance in SV5-7. Compare to 100 random same-size sets. Z-score. Test across all 12 layers.
- Expected signal: At least 3/10 GO terms show z < -2 at some layer.
- Null: Random gene sets; expected z ≈ 0.
- Value: **high** | Cost: **medium** (need GO annotations for gene set)

**B3: STRING PPI Partners as Geometric Neighbors**
- Hypothesis: Protein-protein interaction partners (STRING high-confidence ≥700) are closer in SV5-7 embedding space than random gene pairs, at some layer.
- Test: Intersect STRING human PPIs (≥700 confidence) with cycle4 4941 genes. Compute pairwise distance for PPI edges vs random pairs. AUROC and permutation null.
- Expected signal: AUROC > 0.58 at some layer.
- Null: Permuted gene-pair labels.
- Value: **high** | Cost: **medium** (STRING data download needed)

**B4: Temporal Activation Order Geometry (Early vs Late Response TFs)**
- Hypothesis: Known early-response TFs (FOS, JUN, EGR1 — activated within 30min) vs late-response factors (PRDM1, BCL6, IRF4 — activated at 24-72h) are geometrically ordered along the first PC of the SV5-7 subspace, with early-response at one end.
- Test: For cycle4_immune at each layer, project known early (n≈5) and late (n≈5) immune response TFs to first PC of SV5-7. Test if group PC1 means are separated (t-test). Permutation: shuffle early/late labels.
- Expected signal: mean(PC1_early) ≠ mean(PC1_late), t-test p < 0.05 at some layer.
- Null: Permuted early/late labels show no separation.
- Value: **medium-high** | Cost: **low**

---

### ALGORITHMIC / MECHANISTIC

**A1: Stable TF Anchor Configuration Consistency Across Layers**
- Hypothesis: The 6 stable-low-CV TF anchors from iter_0059 (STAT4, BACH2, ZEB1, NFATC2, RUNX1, JARID2) maintain a consistent pairwise geometric configuration across L0-L9 but diverge at L10-L11, mirroring the cross-seed CKA break.
- Test: At each layer, extract SV5-7 coordinates of 6 anchor TFs. Compute Procrustes distance from L2 configuration. Also compute vs random 6-gene Procrustes as null.
- Expected signal: Anchor Procrustes distance < 0.1 at L0-L9, then rises at L10-L11. Random-set Procrustes higher throughout.
- Null: Random 6-gene sets show similar or higher Procrustes distances.
- Value: **medium** | Cost: **low**

**A2: Effective Gradient Direction in Subspace (Layer-to-Layer Displacement Alignment)**
- Hypothesis: The dominant direction of gene displacement (layer L→L+1) is consistent across genes — i.e., there is a "global flow" in SV5-7 that all genes follow, and this flow direction is aligned with the TF-target separation axis.
- Test: For each of L0→L11 transitions, compute displacement vectors per gene in SV5-7. Find dominant direction (first PC of displacement matrix). Test cosine alignment with LR weight vector (TF-target axis). Also test if TF vs target genes have differently signed projections on this flow direction.
- Expected signal: Cosine(flow, LR_axis) > 0.5 at some layer; TF and target projections on flow differ (t-test p < 0.05).
- Null: Random displacement directions; cosine ≈ 0.
- Value: **medium** | Cost: **low**

---

## Top 3 for Immediate Execution

### Candidate 1: HIGH-PROBABILITY DISCOVERY
**M1 — cycle4_immune TF-Target AUROC at Full Scale (735 pairs)**

Rationale: All prior AUROC work (0.565-0.744) used the 288-pair cycle1 dataset. cycle4 has 735 TRRUST pairs for the same immune context, giving 2.5× the statistical power. This is the most direct scaling test. Data is already loaded in prior scripts. If AUROC reaches ≥0.72 at peak layer, this is a clear publishable result establishing the signal's robustness. Cost is minimal — the exact same code from cycle1 runs on cycle4. Add permutation null and cross-seed replication.

Concrete plan:
1. Load cycle4_immune_main/layer_gene_embeddings.npy [12, 4941, 512]
2. Load cycle4_immune/cycle1_edge_dataset.tsv (or equivalent cycle4 edge file) — confirm 735 positive pairs and their gene indices
3. SVD of centered nonzero embeddings at each layer → project to SV5-7
4. LR classifier (TF vs target) with LOO or 5-fold CV → AUROC
5. Pairwise L2 distance AUROC for all positive+negative pairs
6. Permutation null (2000 gene-label shuffles) at peak layer
7. Repeat for seed43 and seed44

---

### Candidate 2: HIGH-RISK / HIGH-REWARD
**C1 — Cross-Cycle Classifier Transfer (cycle4→cycle1)**

Rationale: If a TF-target classifier trained on cycle4 generalizes to cycle1 via Procrustes alignment, this is the strongest possible evidence for geometry being architecturally determined rather than dataset-specific. This would be the most novel finding of the project so far — universal TF-target geometric separation. The risk: Procrustes alignment may fail (if subspaces are incompatible), giving null. If it works, it's publication-grade.

Concrete plan:
1. Train LR on cycle4 L3 SV5-7 (peak AUROC layer, N=larger dataset)
2. Get cycle1 L3 SV5-7 embeddings for shared genes
3. Procrustes (scipy.spatial.procrustes or sklearn orthogonal_procrustes) align cycle1 to cycle4 using shared gene subset
4. Apply cycle4 LR to Procrustes-aligned cycle1 embeddings
5. Compute transfer AUROC vs within-cycle1 baseline AUROC
6. Permutation: random orthogonal rotation of cycle4 LR weight as null

---

### Candidate 3: CHEAP BROAD-SCREEN
**M2 — Layer-Depth AUROC Mechanism (Sparsity vs Biological Structure)**

Rationale: The single most important unresolved confound. Until we know whether the layer-depth AUROC trend reflects learned biological separation or is driven by gene sparsification, all AUROC claims are suspect. This test uses only existing data (iter_0059 artifacts + any layer's embeddings), takes <30 min to code, and gives a definitive diagnostic. It's a prerequisite for interpreting all AUROC results.

Concrete plan:
1. For cycle4_immune (main, seed43, seed44), at each of 12 layers:
   - Count n_nonzero genes (genes with non-zero embedding norms)
   - Compute SV5-7 AUROC (LR classifier, CV)
   - Compute variance explained by top SV (eigenvalue ratio λ1/sum)
2. Partial Spearman(AUROC, n_nonzero | layer_index) — both ranked
3. Partial Spearman(AUROC, layer_index | n_nonzero) — resolves which drives which
4. Report: if partial r(AUROC, layer | sparsity) > 0.5 → biological. If < 0.2 → sparsity confound.
