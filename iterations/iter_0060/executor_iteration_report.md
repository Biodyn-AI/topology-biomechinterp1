# Executor Iteration Report — iter_0060
Date: 2026-02-23

## Summary
Three hypotheses executed as directed by brainstormer: (H01) signed displacement projection AUROC, (H02) TF family-stratified margin trajectory across layers, (H03) feed-forward loop triangle geometry in SV5-7 space. All ran to completion with machine-readable artifacts.

## Command Trace

```bash
# Main experiment script
conda run -n subproject40-topology python /tmp/iter60_main.py

# Data used:
# Embeddings: .../subproject_38.../implementation/outputs/cycle1_main/layer_gene_embedding_sum.npy [12, 4803, 512]
# Also cycle1_seed43 and cycle1_seed44 (loaded but H01/H02 use main only)
# Edge dataset: .../implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv
# TF family map: .../iter_0059/h01_crossseed_stability.csv

# Output artifacts in iter_0060/:
# h01_signed_proj_auroc.csv
# h02_family_margin_trajectory.csv
# h03_ffl_geometry.csv
# h03_ffl_permutation.csv
```

---

## H01: Signed Displacement Projection AUROC (Directional Geometry)

**Novelty**: new_method in manifold_distance family. Brainstormer prediction: directional projection should exceed AUROC=0.62 (scalar distance at 0.565).

**Method**: At each of 12 layers, compute SVD of centered nonzero embeddings [~1972, 512]. Project all genes to SV5-7 (3D). For all 288 positive TRRUST edges, compute displacement vectors (target_proj - source_proj). Compute mean direction (unit vector). For all 961 test pairs, project displacement onto this mean direction (signed scalar). AUROC = signed_projection predicts positive edge. Permutation test (2000 shuffles) at peak layer.

**Results**:

| Layer | AUROC | mean_dir_mag |
|-------|-------|-------------|
| L0 | 0.563 | 45.4 |
| L1 | 0.554 | 35.7 |
| L2 | 0.553 | 24.6 |
| L3 | 0.524 | 21.9 |
| L4 | 0.535 | 27.3 |
| L5 | 0.516 | 27.9 |
| L6 | 0.504 | 50.2 |
| L7 | 0.489 | 46.5 |
| L8 | 0.502 | 57.9 |
| L9 | 0.500 | 58.2 |
| L10 | 0.516 | 48.5 |
| L11 | 0.521 | 41.2 |

Max AUROC = **0.563 at L0** (perm_p=0.001, null=0.499±0.021)
Mean AUROC across layers: 0.523

**Interpretation**: Directional projection onto mean TF→target displacement gives max AUROC=0.563, essentially identical to scalar distance AUROC=0.565 from iter_0059/H02. The brainstormer prediction (>0.62) is NOT met. The signal is real (perm_p=0.001) but weak. Directional projection does not improve over scalar distance — the information in SV5-7 pairwise geometry for edge prediction is truly limited at ~0.563.

**Decision**: negative. Directional projection is equivalent to scalar distance, not better. The 0.62 threshold is not met.

---

## H02: TF Family-Stratified Margin Trajectory Across Layers

**Novelty**: new_family approach using biological stratification (cheap re-analysis).

**Method**: At each of 12 layers, SVD of centered nonzero embeddings → project to SV5-7. Fit LR (TF vs target genes). Compute signed margin per TF gene. Stratify by TF family (bZIP, ETS, C2H2-ZF, bHLH, GATA, RUNX, Rel, Other). Report mean margin per family per layer. Also record overall AUROC.

**TF families**: 9 families across 58 TF genes. bHLH n=1 (HIF1A), bZIP n=4 (BATF/FOS/JUN/JUNB), GATA n=1, RUNX n=2, ETS n=2, C2H2-ZF n=2, Rel n=1.

**Results — Overall AUROC by layer**:

| Layer | AUROC (SV5-7, TRRUST genes only) |
|-------|----------------------------------|
| L0 | 0.570 |
| L2 | 0.494 |
| L3 | 0.589 |
| L6 | 0.640 |
| L9 | 0.725 |
| L11 | **0.739** |

Note: These AUROCs differ from prior results because here only TRRUST-subset genes are used in the classifier (not all nonzero genes), and no CV is applied.

**Key findings — Family margin trajectory**:

1. **bHLH (HIF1A, n=1)**: Dramatic sign-reversal. Consistently negative margins at L0-L7 (target-like), then flips strongly positive at L8-L11 (margin: +1.68 → +2.76 → +2.53). HIF1A acquires a TF-like representation only in the deepest transformer layers.

2. **bZIP (BATF, FOS, JUN, JUNB, n=4)**: Most consistently target-like across all layers (mean margin: -0.77 to -0.98). These immune/stress response TFs are always on the "target" side of the SV5-7 boundary.

3. **ETS (n=2)**: Least negative family at early layers (mean margin -0.35), with mean_margin=+0.36 at L11.

4. **Other family (n=45)**: Large spread, but mean stays negative throughout (around -0.37 to -0.73).

Peak positive family at each layer:
- L0-L5: bHLH (HIF1A) or C2H2-ZF (least negative)
- L6-L11: bHLH strongly flips positive while all others remain negative or mildly negative

**Caveat**: bHLH result is single-gene (HIF1A). Biological interpretation: HIF1A may transition from a regulated gene (context-dependent hypoxia response) to a master regulator representation in the deeper transformer.

**Decision**: neutral (bHLH effect is dramatic but n=1; bZIP pattern robust across 4 genes; layer-wise AUROC trend is informative but contradicts prior reports using full gene set).

---

## H03: Feed-Forward Loop Triangle Geometry in SV5-7

**Novelty**: new_family — first test of FFL motif geometry in this dataset. Novel, no prior art.

**Method**: Build directed regulatory graph from 288 positive TRRUST edges. Identify FFLs (A→B, A→C, B→C where A,B=TFs, C=target). For each FFL triplet (A,B,C) and each layer: project to SV5-7. Compute t = (B-A)·(C-A) / |C-A|² (interpolation parameter). If B is geometrically "between" A and C on the A→C axis, t ∈ (0,1). One-sample t-test vs t=0. Permutation test for frac_between at L2 (2000 shuffles of gene positions).

**FFLs found**: 22 triplets. Examples: (RUNX1→JUN→IL2), (STAT4→TBX21→IFNG), (ETS1→FLI1→TGFBR2).

**Results**:

| Layer | n_ffls | t_mean | frac_between | p_vs0 |
|-------|--------|--------|-------------|-------|
| L0 | 20 | 0.565 | 0.550 | 0.0016 |
| L1 | 20 | 0.586 | 0.600 | 0.0007 |
| **L2** | **20** | **0.598** | **0.550** | **0.0001** |
| L3 | 20 | 0.770 | 0.650 | 0.0011 |
| L4 | 20 | 0.622 | 0.400 | 0.0002 |
| L5 | 20 | 0.600 | 0.450 | 0.0016 |
| L6 | 20 | 0.752 | 0.500 | 0.0018 |
| L7 | 20 | 0.459 | 0.350 | 0.0172 |
| L8 | 20 | 0.226 | 0.350 | 0.6496 |
| L9 | 20 | 0.197 | 0.250 | 0.4459 |
| L10 | 20 | 0.116 | 0.350 | 0.5799 |
| L11 | 20 | 0.106 | 0.500 | 0.5939 |

**Permutation test at L2** (shuffle gene positions):
Real frac_between=0.550, perm_null=0.461±0.124, **p=0.3015** (not significant)

**Interpretation**:
- **t_mean signal**: Highly significant at L0-L6 (p<0.002). On average, the intermediate TF B is displaced ~0.6 of the way from A toward C along the A→C axis. This is significantly > 0 (B is in the positive direction from A toward C).
- **Signal collapse at L8-L11**: t_mean drops to ~0.1 and becomes non-significant, suggesting the geometric ordering breaks down in deeper layers.
- **Permutation test**: frac_between (binary betweenness) is NOT significant after permutation (p=0.30), meaning the fraction geometrically "between" could arise by chance.
- **Conflict**: t_mean test (mean continuous) is significant but frac_between permutation test is not. The t_mean test is higher-powered for detecting mean direction, while frac_between is a weaker binary test.
- **Power limitation**: N=22 FFLs is small. The signal is consistent across L0-L6 but underpowered for strong conclusions.

**Decision**: inconclusive. t_mean signal real at L0-L6 (p<0.0001 at L2) but betweenness permutation not significant; N=22 FFLs insufficient for strong claims. Pattern (geometric ordering at early-to-mid layers, collapse at deep layers) is biologically interesting and novel.

---

## Summary Table

| Hypothesis | Family | Decision | Direction | Key Metric |
|------------|--------|----------|-----------|------------|
| H01: Signed displacement projection | manifold_distance | negative | positive | max AUROC=0.563 (perm_p=0.001) |
| H02: TF family margin trajectory | module_structure | neutral | mixed | bHLH +2.76 at L10 (n=1); bZIP -0.94 (n=4) |
| H03: FFL triangle geometry | module_structure | inconclusive | mixed | t_mean=0.60 at L2 (p=0.0001); perm frac_between p=0.30 |
