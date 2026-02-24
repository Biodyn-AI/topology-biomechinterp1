# Executor Iteration Report — iter_0059
Date: 2026-02-23

## Summary
Three hypotheses executed: (H01) cross-seed TF boundary anchor stability, (H02) pairwise TF→target 6D distance as edge predictor, (H03) partial Spearman of effective rank controlling for layer depth. All experiments ran to completion with machine-readable artifacts.

## Command Trace

```bash
# Main experiment script
conda run -n subproject40-topology python /tmp/iter59_fixed.py

# Data paths used:
# BASE = .../subproject_38.../implementation/outputs
# Embeddings: cycle1_main/layer_gene_embedding_sum.npy  [12, 4803, 512]
# Edge dataset: cycle1_main/cycle1_edge_dataset.tsv
# Prior data: iterations/iter_0058/h03_effrank_auroc.csv
```

## Data Schema
- Embedding arrays: `[12_layers, 4803_genes, 512_dim]`
- TRRUST gene2emb_idx: 209 unique genes via `source_idx`/`target_idx` columns
- TF genes: 64 unique TFs, target genes: 165 unique targets
- Positive edges: 288, Negative edges: 864

---

## H01: Cross-Seed TF Boundary Anchor Stability + Out-Degree

**Method**: At peak layers L2 and L3 (from iter_0058), run SVD of centered nonzero embeddings [1972,512] per seed (main/seed43/seed44). Project to SV5-7. Fit logistic regression (TF vs target binary, StandardScaler, C=1.0). Compute per-TF signed margin. Measure cross-seed variance (CV = std/|mean|) and Spearman correlation with TRRUST out-degree.

**Results**:

| Layer | Seed | n_TF | n_target | AUROC |
|-------|------|------|----------|-------|
| L2 | main | 60 | 127 | 0.744 |
| L2 | seed43 | 59 | 131 | 0.729 |
| L2 | seed44 | 59 | 130 | 0.739 |
| L3 | main | 60 | 127 | 0.762 |
| L3 | seed43 | 59 | 131 | 0.741 |
| L3 | seed44 | 59 | 130 | 0.762 |

**Stability** (116 gene×layer records):

Top stable positive anchors (low CV across seeds):
- JARID2 L2: mean_margin=0.383, CV=0.010 (family=Other, out-degree=1)
- STAT4 L3: mean_margin=1.254, CV=0.067 (family=Other, out-degree=5)
- NFATC2 L3: mean_margin=0.860, CV=0.098 (family=Other, out-degree=4)
- BACH2 L2: mean_margin=1.103, CV=0.108 (family=Other, out-degree=2)
- RUNX1 L2: mean_margin=0.315, CV=0.152 (RUNX, out-degree=9)
- ZEB1 L2: mean_margin=0.810, CV=0.158 (C2H2-ZF, out-degree=2)

**Out-degree correlation**: Spearman(mean_margin, TRRUST out-degree) = **r=-0.012, p=0.90** → NO correlation. High-margin TFs are not the hub regulators (most are low-degree).

**Decision**: neutral/inconclusive. Boundary geometry is reproducible but high-margin TFs are low-degree specialized factors, not hub TFs. The mechanistic "boundary anchored by hub TFs" hypothesis is **falsified**.

---

## H02: Pairwise TF→Target 6D Distance as Edge Probability Predictor

**Method**: For each of 12 layers, compute SVD of nonzero embeddings, project to SV5-7. For all positive (288) and negative (864) TRRUST edges, compute L2 distance between TF and target projections. AUROC = edge predictor (small distance → positive edge). Mann-Whitney test. Permutation null (1000 shuffles of last layer's distances).

**Results**:

| Layer | AUROC | pos_mean | neg_mean | mw_p |
|-------|-------|----------|----------|------|
| L0 | 0.498 | 3.760 | 3.779 | 0.543 |
| L2 | 0.531 | 2.922 | 3.074 | 0.074 |
| L6 | 0.549 | 2.363 | 2.600 | 0.011 |
| L8 | 0.558 | 2.288 | 2.541 | 0.003 |
| **L9** | **0.565** | **2.073** | **2.338** | **0.001** |
| L10 | 0.554 | 1.747 | 1.918 | 0.006 |
| L11 | 0.552 | 1.362 | 1.495 | 0.007 |

Mean AUROC across layers: 0.536
Permutation null: 0.500 ± 0.021
Permutation p-value at max AUROC (0.565): **p=0.002**

**Interpretation**: Signal is statistically real but small (AUROC 0.565 vs 0.500 null). Does NOT reach the 0.62 threshold for a "qualitative leap" in regulatory inference. Positive TRRUST edges are slightly but consistently closer in SV5-7 than negative edges, particularly in later layers (L6-L11). The signal is consistent but too weak for practical use.

**Decision**: negative relative to H02's original goal (0.62 threshold), weak positive signal exists.

---

## H03: Partial Spearman (Effective Rank vs AUROC | Layer Depth) + Lipschitz

**Method**: Use iter_0058's h03_effrank_auroc.csv (36 seed×layer observations). Partial Spearman: residualize ranked eff_rank and ranked AUROC on ranked layer index, then Pearson on residuals. Also compute mean gene displacement between consecutive layers in SV5-7 (Lipschitz regularity).

**Results**:

| Metric | Value |
|--------|-------|
| Raw Spearman(eff_rank, AUROC) | r=0.855, p=3.3e-11 |
| **Partial Spearman(eff_rank, AUROC \| layer)** | **r=-0.045, p=0.795** |
| Spearman(layer, eff_rank) | r=-0.997 |
| Spearman(layer, AUROC) | r=-0.859 |
| Lipschitz Spearman(layer, displacement) | r=-0.573, p=0.066 |

**Critical finding**: The r=0.855 correlation between effective rank and AUROC is **entirely a layer-depth confound**. After partialing out layer, the correlation collapses to r=-0.045. Both eff_rank and AUROC decrease monotonically with layer depth (r=-0.997 and -0.859 respectively), creating spurious correlation.

**Lipschitz**: Mean inter-layer displacement in SV5-7 shows a non-monotone pattern (peaks at L5: 3.88, trough at L10: 0.83). Weak negative trend with layer (r=-0.573, p=0.066) — later layers make smaller geometric steps, consistent with convergence.

**Decision**: negative (for eff_rank as causal predictor). The rho=0.855 claim is a confound artifact. Effective rank is NOT an independent predictor of TF-target discriminability.

---

## Summary Table

| Hypothesis | Family | Decision | Direction | Key Metric |
|------------|--------|----------|-----------|-----------|
| H01: Cross-seed TF anchor stability | module_structure | neutral | mixed | CV<0.15 stable; out-degree r=-0.01 |
| H02: Pairwise 6D dist AUROC | manifold_distance | negative | positive-weak | AUROC=0.565 < threshold 0.62 |
| H03: Partial Spearman eff_rank | intrinsic_dimensionality | negative | negative | partial_r=-0.045 vs raw_r=0.855 |

## Key Insight for Next Iteration
The partial Spearman result (H03) is the most important finding: the eff_rank→AUROC correlation is spurious. This means the layer-depth trend in AUROC is the primary signal, not the dimensionality per se. The question becomes: what drives the layer-depth AUROC trend? Candidates: (1) feature sparsification (fewer nonzero genes), (2) feature compression (lower eff_rank as consequence not cause), (3) actual separation of TF vs target clusters in later layers.

Artifacts:
- `h01_crossseed_stability.csv` (116 rows)
- `h02_pairwise_dist_auroc.csv` (12 rows)
- `h03_partial_spearman.csv` (1 row summary)
- `h03_lipschitz.csv` (11 rows)
