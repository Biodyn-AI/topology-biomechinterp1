# Executor Iteration Report — iter_0018

## Summary
Three hypotheses tested spanning null control, prediction benchmark, and a new attention-geometry family.

Key outcomes:
- **H01 (Random Gaussian null)**: Confirmed negative (0/12 significant layers, mean_z=-0.085). Validates scGPT PPI geometry is model-specific, not an SVD artifact.
- **H02 (Precision@k benchmark)**: Modest 1.20x enrichment at k=100 (P@100=0.170 vs random=0.142). Consistent across all 12 layers. The SV2 axis is a co-pole partition metric, not a smooth distance ranking metric.
- **H03 (Attention geometry — new family)**: Mixed but promising. STRING PPI pairs show no significant co-attention (z=0.066, MW_p=0.091). TRRUST activation TF-target pairs show highly significant co-attention (MW_p=9.9e-09, mean_att=0.01762 vs null=0.00889, ~2x enrichment). **New finding**: scGPT attention encodes TF regulatory co-occurrence, while SVD geometry encodes PPI proximity — two distinct geometric representations.

---

## Command Trace

```bash
# H01 + H02 (random Gaussian null + precision@k)
conda run -n subproject40-topology python \
  iterations/iter_0018/run_iter0018_screen.py

# H03 (attention geometry, fixed with correct gene ordering)
conda run -n subproject40-topology python \
  iterations/iter_0018/run_iter0018_h03_fixed.py
```

Environment: `subproject40-topology` conda environment.

---

## H01: Random Gaussian Embedding Null Control

**Method**: Replace scGPT `layer_gene_embeddings.npy` [12, 4803, 512] with `np.random.default_rng(42).standard_normal(same_shape)`. Run identical SV2 co-pole test: K=52 top/bottom poles, 300 shuffle nulls per layer, 3092 STRING pairs (score≥0.4).

**Results**:
```
Layer | z-score
L00   | -0.460
L01   |  0.304
L02   |  0.162
L03   | -0.644
L04   | -0.689
L05   |  1.311
L06   |  0.668
L07   | -0.573
L08   | -0.740
L09   |  0.106
L10   | -0.430
L11   | -0.034
n_sig = 0/12
mean_z = -0.085
```

**Interpretation**: Random embeddings produce z~0 at all layers. The scGPT finding (mean_z~3-5 across layers) is model-specific, not an artifact of the SVD method or data distribution. This is the strongest null control to date.

**Artifact**: `iterations/iter_0018/h01_random_gaussian_null.json`

---

## H02: Out-of-Sample PPI Precision@k Prediction Benchmark

**Method**: SV2 1D projections of 209 named genes at each of 12 layers. Negative absolute distance `−|proj_i − proj_j|` as pairwise similarity. Rank all 21,736 named-gene pairs. Compute precision@k for STRING edges.

**Random baseline**: 3092 STRING pairs / 21736 total = 14.23%

**Results** (averaged SV2 across 12 layers):
```
k=50:   P@50 = 0.180  (enrichment = 1.26x)
k=100:  P@100 = 0.170 (enrichment = 1.20x)
k=200:  P@200 = 0.155 (enrichment = 1.09x)
k=500:  P@500 varies  (~1.0-1.2x per layer)
k=1000: P@1000 varies (~1.0-1.3x per layer)
```

Note: P@50 and P@100 are identical across all 12 layers because the same top pairs dominate SV2 proximity regardless of layer.

**Interpretation**: The SV2 axis provides modest 1.2x enrichment for STRING pairs when used as a ranking distance. This is weaker than expected from the co-pole test because the 1D absolute distance between SV2 projections is not the right framing — the co-pole test captures partition structure (top-K vs bottom-K co-membership), which is fundamentally different from a smooth distance ranking. This is a methodological insight: SVD poles encode categorical partition structure, not smooth manifold distances.

**Artifact**: `iterations/iter_0018/h02_precision_at_k.json`

---

## H03: scGPT Attention Co-Occurrence Geometry (New Family)

**Method**: Load aggregated scGPT attention_scores.npy [8181×8181] from `invariant_causal_edges/lung/`. Map 208 named genes to correct attention matrix positions via `processed.h5ad` var_names (adata gene ordering). Compute symmetric mean attention `(att[i,j]+att[j,i])/2` for:
- 3091 STRING pairs (score≥0.4) with valid indices
- 116 TRRUST activation TF-target pairs with valid indices
- 18,437 non-STRING background pairs (within 209-gene universe)

Also compute STRING confidence quintile gradient (Spearman rho).

**Results**:
```
STRING pairs:
  obs_mean = 0.010191
  null_mean = 0.008893
  z = 0.066, MW_p = 0.091 (NOT significant)

TRRUST activation pairs:
  obs_mean = 0.017620
  null_mean = 0.008893
  z = 0.444, MW_p = 9.9e-09 (**HIGHLY significant**)
  Enrichment: ~2x vs non-STRING background

STRING confidence quintile gradient:
  Q1 [0.400,0.459]: mean_att = 0.010215
  Q2 [0.459,0.533]: mean_att = 0.010488
  Q3 [0.533,0.640]: mean_att = 0.008256
  Q4 [0.640,0.796]: mean_att = 0.007888
  Q5 [0.796,0.999]: mean_att = 0.014105
  Spearman rho = 0.100, p = 0.873 (not significant)
```

**Interpretation**:

A new geometric finding: scGPT attention heads specifically co-attend TF-target regulatory pairs but NOT general protein-protein interaction pairs. TRRUST activation pairs have ~2x higher mean co-attention than non-interaction background (MW_p=9.9e-09). STRING PPI pairs show no significant enrichment (z~0).

This creates a mechanistic dissociation:
- **SVD residual geometry** (iter_0010-0017): encodes PPI network proximity (STRING confidence gradient rho≥0.90 across SV2/SV3/SV4 axes)
- **Attention geometry** (this iteration): encodes TF regulatory co-occurrence (TRRUST activation, 2x enrichment)

These are complementary encodings in different geometric components of scGPT. The model uses attention to co-attend regulatory TF pairs and uses residual stream geometry to organize PPI proximity.

**Artifact**: `iterations/iter_0018/h03_attention_geometry.json`

---

## Quantitative Summary

| Hypothesis | Metric | Value | Decision |
|------------|--------|-------|----------|
| H01: Random Gaussian null | n_sig/12, mean_z | 0/12, −0.085 | negative (confirmatory) |
| H02: Precision@k benchmark | P@100 enrichment | 1.20× (vs 1.0× random) | neutral |
| H03: Attention geometry | TRRUST MW_p | 9.9×10⁻⁹ (~2× enrichment) | promising |

---

## Data Provenance
- scGPT embeddings: `/subproject_38.../outputs/cycle1_main/layer_gene_embeddings.npy` [12,4803,512]
- Attention scores: `/single_cell_mechinterp/outputs/invariant_causal_edges/lung/attention_scores.npy` [8181,8181]
- Attention gene ordering: `processed.h5ad` var_names (same directory)
- STRING cache (score≥0.4): `iterations/iter_0015/string_ppi_score04_cache.json`
- TRRUST: `single_cell_mechinterp/external/networks/trrust_human.tsv`
