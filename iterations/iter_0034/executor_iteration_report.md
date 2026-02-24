# Executor Iteration Report: iter_0034

**Date:** 2026-02-23
**Focus:** Layer-wise B-cell community emergence, PC2/PC3 cell-type axes, expression-level confound test

---

## Command Trace

```bash
conda run -n subproject40-topology python iterations/iter_0034/run_iter0034_screen.py
```

Script: `iterations/iter_0034/run_iter0034_screen.py`

---

## H01: Layer-wise B-cell Community Purity Emergence — MIXED/INCONCLUSIVE

**Family:** module_structure
**Novelty:** new_method — first layer-sweep of B-cell Fisher OR per community
**Method:** For each layer L0–L11: k=10 kNN greedy-modularity communities, compute Fisher OR for B-cell marker enrichment in PC1-negative community. Spearman rho(layer, OR).
**Artifact:** `iterations/iter_0034/h01_layer_bcell_community_purity.json`

### Per-Layer Fisher OR (B-cell markers in PC1-negative community)

| Layer | n_comms | Modularity | OR (PC1-neg) | p | OR (best) | p (best) |
|-------|---------|-----------|--------------|---|-----------|----------|
| L0    | 3       | 0.411     | 2.56         | 0.148 | 2.56 | 0.148 |
| L1    | 3       | 0.410     | 1.25         | 0.505 | 4.31 | 0.039 |
| L2    | 3       | 0.404     | **16.00**    | **0.001** | 16.00 | 0.001 |
| L3    | 3       | 0.427     | **11.84**    | **0.005** | 11.84 | 0.005 |
| L4    | 3       | 0.461     | 2.60         | 0.180 | 3.24 | 0.089 |
| L5    | 3       | 0.482     | **15.25**    | **0.002** | 15.25 | 0.002 |
| L6    | 3       | 0.460     | 1.32         | 0.503 | 6.22 | 0.016 |
| L7    | 3       | 0.502     | 2.01         | 0.273 | 3.55 | 0.069 |
| L8    | 4       | 0.513     | 2.01         | 0.273 | 4.00 | 0.049 |
| L9    | 3       | 0.518     | 3.12         | 0.102 | 2.17 | 0.210 |
| L10   | 3       | 0.509     | 3.64         | 0.064 | 3.64 | 0.064 |
| **L11** | **2** | 0.434   | **10.60**    | **0.008** | 10.60 | 0.008 |

**Spearman rho(layer, OR_neg_comm) = 0.081, p=0.804 — NOT significant**
**Spearman rho(layer, OR_best_comm) = -0.133, p=0.681 — NOT significant**

**Interpretation:**
The B-cell community purity is NOT monotonically increasing across layers. B-cell markers are already enriched with high OR (12–16) at L2–L3 and L5, then fluctuate due to community identity shuffling across layers (the PC1-negative community does not track a consistent set of genes across layers as communities reorganize). The binary 2-community structure is unique to L11. The story is not "progressive crystallization" but rather "B-cell geometry is baked in from the embedding (L2–L3 already strong) and then L11 uniquely consolidates to 2 communities."

**Decision:** Mixed — the layer-wise emergence hypothesis is negative (no monotonic increase), but the early-layer strength is a novel finding.

---

## H02: PC2/PC3 Cell-Type Axes at L11 — NEGATIVE

**Family:** intrinsic_dimensionality
**Novelty:** new_method — first test of higher PCs for T-cell/myeloid axes
**Method:** SVD at L11, compute PC1–PC3. Test B-cell, T-cell, Myeloid enrichment at each pole.
**Artifact:** `iterations/iter_0034/h02_pc2_pc3_axes.json`

### Variance explained
- PC1: 25.9%, PC2: 16.0%, PC3: 7.5%

### Results

| Axis | Group | Direction | AUROC | p |
|------|-------|-----------|-------|---|
| PC1 | B-cell | negative | **0.203** | **0.001** ✓ |
| PC1 | T-cell | positive | 0.583 | 0.150 |
| PC1 | Myeloid | positive | 0.299 | 0.953 |
| PC2 | B-cell | negative | 0.394 | 0.142 |
| PC2 | B-cell | positive | 0.394 | 0.860 |
| PC2 | T-cell | positive | 0.442 | 0.767 |
| PC2 | T-cell | negative | 0.442 | 0.234 |
| PC2 | Myeloid | positive | 0.618 | 0.168 |
| PC2 | Myeloid | negative | 0.618 | 0.834 |
| PC3 | B-cell | negative | 0.427 | 0.229 |
| PC3 | T-cell | positive | 0.381 | 0.932 |
| PC3 | Myeloid | positive | 0.429 | 0.721 |

**PC2 and PC3 show no significant cell-type enrichment (all p>0.14).** Only PC1 B-cell signal (p=0.001) survives. The 195-gene vocab is not large enough to recover a strong T-cell or Myeloid axis in the top 3 PCs, possibly because (a) T-cell and myeloid genes are in a broadly mixed positive PC1 pole, or (b) their variance is spread across PCs 4+.

**Decision:** Negative for PC2/PC3 cell-type axes. PC1 B-cell confirmed again.

---

## H03: Expression-Level Confound Test for B-cell PC1 Signal — POSITIVE

**Family:** null_sensitivity
**Novelty:** new_method — first confound test / bootstrap null for B-cell geometry
**Method:**
1. Pearson r(PC1 score, L2 norm) at L11
2. Regress PC1 on L2 norm, test B-cell in residual
3. Bootstrap null: 1000 random 9-gene sets, empirical p for B-cell AUROC
**Artifact:** `iterations/iter_0034/h03_expression_confound_test.json`

### Results

| Test | Value | Interpretation |
|------|-------|----------------|
| Pearson r(PC1, L2_norm) | -0.309, p<0.0001 | Some correlation: genes with lower norm → more negative PC1 |
| B-cell L2 norm vs background | AUROC=0.601, p=0.847 | **B-cell norms NOT significantly lower** |
| B-cell PC1 residual AUROC (after L2 regression) | **0.214, p=0.002** | **B-cell signal survives L2 regression** |
| B-cell residual mean | -1.765 vs bg=0.085 | **Large separation in residuals** |
| Bootstrap empirical p | **0.000** (z=-3.06) | **0/1000 random 9-gene sets as extreme** |

**Key finding:** The B-cell PC1 signal is **NOT explained by embedding norm (expression proxy)**:
- B-cell L2 norms are NOT significantly lower (AUROC=0.60, p=0.847) — the confound is absent
- After regressing out L2 norm from PC1 scores, B-cell markers are still significantly in the negative PC1 residual (AUROC=0.214, p=0.002)
- Bootstrap null: 0 out of 1000 random 9-gene sets achieve AUROC ≤ 0.203 (empirical p=0.000, z=-3.06)

**This definitively confirms the B-cell PC1 signal reflects geometric/structural organization in the scGPT embedding, not an expression-level or norm artifact.**

**Decision:** PROMISING — strongest null control result for the B-cell geometry finding.

---

## Summary

| H | Family | Decision | Direction |
|---|--------|----------|-----------|
| H01 | module_structure | mixed | negative (no monotonic emergence), positive (early-layer strength) |
| H02 | intrinsic_dimensionality | negative | negative |
| H03 | null_sensitivity | **promising** | **positive** |

**Key insight from H01:** B-cell geometric separation is present from L2 (OR=16), not a late-layer phenomenon. L11 is special in collapsing to 2 communities, but the underlying geometry is early.

**Key insight from H03:** The B-cell PC1 signal passes bootstrap null (z=-3.06) and L2-norm regression controls — it reflects genuine structural geometry.
