# Executor Iteration Report — iter_0061

**Date**: 2026-02-23
**Status**: Complete — 3 hypotheses tested

---

## Summary

1. **H01** (Layer-to-layer CKA trajectory): **Negative** — No sharp drop at L7→L8 (consecutive CKA 0.977–0.981). Cross-seed CKA declines monotonically from L0 to L11, steepest at L10–L11 (0.87→0.78).
2. **H02** (Dual-role gene margin trajectory): **Neutral** — All 38 dual-role (TF∩target) genes are consistently target-like (negative margins) across all 12 layers. No sign reversal; HIF1A flip (iter_0060) is not generalizable.
3. **H03** (FFL geometry expanded + permutation, SV2-4/SV5-7/SV2-7): **Negative (→retire)** — Permutation test shows random triplets have t_mean ~0.50, real FFL triplets 0.34–0.59. No significant geometric ordering above null.

---

## Command Trace

```bash
# Main experiment
conda run -n subproject40-topology python /tmp/iter61_full.py

# Supplementary permutation (SV5-7 peak layers)
conda run -n subproject40-topology python /tmp/iter61_perm_sv57.py
```

**Data paths used:**
- `cycle4_immune_main/layer_gene_embeddings.npy` [12, 4941, 512]
- `cycle4_immune_seed43/layer_gene_embeddings.npy`
- `cycle4_immune_seed44/layer_gene_embeddings.npy`
- `cycle4_immune_main/cycle1_edge_dataset.tsv` (735 positive TRRUST pairs)

**Environment:** `subproject40-topology` (conda)

---

## H01: Layer-to-layer CKA Trajectory

### Method
- Linear CKA on 1496 shared nonzero genes (union of mask across seeds at all layers)
- Three types: consecutive-layer CKA (LN→LN+1), cross-seed CKA (per layer), each-layer vs L0

### Results

**Consecutive-layer CKA (mean across 3 seeds):**
| Layer pair | CKA |
|-----------|-----|
| L0→L1     | 0.983 |
| L3→L4     | 0.973 |
| L6→L7     | 0.981 |
| L7→L8     | 0.977 |
| L8→L9     | 0.979 |
| L10→L11   | 0.972 |

No sharp drop at L7→L8. Representations change smoothly throughout the network.

**Cross-seed CKA (mean across 3 pairs):**
| Layer | Cross-seed CKA |
|-------|---------------|
| L0    | 0.979 |
| L4    | 0.953 |
| L8    | 0.932 |
| L10   | 0.870 |
| L11   | 0.779 |

Monotonic decline; steepest at L10–L11. Seeds diverge most in deep layers, suggesting deep-layer representations are more sensitive to random seed.

### Interpretation
The L8 boundary hypothesis (from FFL collapse and HIF1A flip in iter_0060) is **not confirmed** by CKA. Consecutive layer similarity is uniformly high. However, the dramatic cross-seed CKA drop at L10–L11 suggests these layers are computationally less constrained — potentially more "creative" or task-specific.

**Decision: negative** (for L8 boundary). Side finding: deep layers (L10-L11) warrant investigation.

---

## H02: Dual-Role Gene Margin Trajectory

### Method
- 38 dual-role genes (appear as both TF and target in TRRUST): BCL11A, GATA1, BCL6, NFKB1, STAT3, CEBPB, PRDM1, PPARG, FOS, JUN, RELA, SPI1, FOXO3, EGR1, etc.
- At each layer: SVD of nonzero embeddings → SV5-7. Fit LR (30 pure TFs vs 197 pure targets). Report LR decision margin for each dual gene.
- Pre-L8 vs post-L8 t-test for sign change.

### Results

**Mean margin by layer (main seed, negative = target-like):**
| Layer | Mean margin |
|-------|------------|
| L0    | -1.481 |
| L2    | -1.599 |
| L4    | -1.715 |
| L6    | -1.878 |
| L8    | -1.914 |
| L10   | -1.760 |
| L11   | -1.726 |

All dual-role genes are classified as target-like across all layers. The margin becomes slightly more negative at L6–L8 (deeper target-classification) but never reverses.

**Pre-L8 mean = -1.685, Post-L8 mean = -1.797** → t=1.733, p=0.084 (not significant at α=0.05).

### Interpretation
The HIF1A sign reversal at L8 (iter_0060) was an outlier — not representative of the 38-gene dual population. The general rule is: genes that appear as TRRUST targets (even if they also are TFs) are embedded in target-like space. This is a robust negative for "dual genes are intermediate."

**Decision: neutral** (confirms target-classification dominance for dual-role genes; closes the dual-role branch).

---

## H03: FFL Geometry Multi-Subspace with Expanded N and Permutation

### Method
- Enumerate all FFLs (A→B, A→C, B→C; A,B=TFs, C=target) from 735 positive TRRUST edges → **226–264 valid FFLs** (vs N=22 in iter_0060)
- For each layer (0–11) and subspace (SV2-4, SV5-7, SV2-7): project, compute interpolation t
- Permutation test (2000 shuffles of gene positions within nonzero set) at best subspace

### Results

**Parametric t-test (misleading without permutation):**
- SV2-4 at L0: t_mean=0.601, p<1e-23 (but null is also ~0.50)
- SV5-7 at L7: t_mean=0.603, p<1e-42
- SV2-7: highly significant at all layers

**Permutation test (correct null):**
| Subspace | Layer | real_t_mean | null_mean | null_std | perm_p |
|----------|-------|------------|-----------|----------|--------|
| SV2-7    | L2    | 0.345       | 0.503     | 0.129    | 0.887  |
| SV5-7    | L5    | 0.572       | 0.508     | 0.183    | 0.367  |
| SV5-7    | L6    | 0.583       | 0.502     | 0.181    | 0.325  |
| SV5-7    | L7    | 0.594       | 0.499     | 0.180    | 0.297  |
| SV5-7    | L8    | 0.555       | 0.496     | 0.175    | 0.379  |

**All perm_p > 0.29.** Not significant. The null distribution has t_mean ~0.50 because a random gene B lies geometrically between random A and C ~50% of the time (by symmetry in embedding space). FFL-specific ordering does not exceed this baseline.

### Key Insight
The parametric t-test (t_mean > 0) was systematically misleading because the correct null is t_mean ≈ 0.50, not 0. FFL triplets do not show statistically special geometric ordering compared to random triplets from the same embedding space. **Retiring this branch.**

**Decision: negative (→retire)**

---

## Artifacts Generated

| File | Content |
|------|---------|
| `h01_cka_trajectory.csv` | CKA values for all layer pairs, seeds, types |
| `h02_dual_role_margins.csv` | Per-gene per-layer margins for 38 dual-role genes |
| `h03_ffl_sv_comparison.csv` | t_mean/p_val/frac_between for all layers × subspaces |
| `h03_perm_tmeans.npy` | Permutation null distribution (SV2-7, L2) |
| `h03_perm_summary.csv` | Permutation summary (SV2-7, L2) |
| `h03_perm_sv57_peak_layers.csv` | Permutation summary (SV5-7, L5–L8) |
| `executor_hypothesis_screen.json` | Machine-readable hypothesis outcomes |

---

## Retirements This Iteration

- **FFL geometric ordering** (H03 this iteration, H03 iter_0060): Retire. Permutation null ~0.50 exceeds real in some cases; SV5-7 peak real=0.59 not significant (perm_p=0.30). Two permutation tests both negative.
- **Dual-role gene L8 flip**: Not generalizable; retire as standalone hypothesis.
- **L8 transition as CKA boundary**: Falsified.

## Emerging Directions

1. **Deep layer divergence (L10-L11)**: Cross-seed CKA drops sharply here (0.87→0.78). Test: are deep-layer subspaces less reproducible because they encode more sample-specific features? Test with persistent homology or intrinsic dimensionality at L10-L11 vs L2-L3.
2. **Cross-seed alignment stability by layer**: Which layer's geometry is most consistent across seeds? (L0 highest; L11 lowest). This may map to regulatory vs contextual representations.
3. **Pure null baseline for geometric tests**: The permutation results reveal that naive geometric tests (t_mean > 0) are biased. All future geometric tests should use permutation-corrected baselines.
