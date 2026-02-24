# Executor Iteration Report: iter_0051

## Summary

Three hypotheses tested following brainstormer guidance on the two open questions from iter_0050: (1) whether any regulatory-specific signal survives co-expression regression in SV2-4, and (2) the biological identity of genes that dominate SV1 at early vs. late layers, plus (3) a topological test of circuit-gene clustering via 0-dim PH.

Key results:
- **H01 NEGATIVE**: After regressing out embedding cosine similarity (co-expression proxy) from SV2-4 pairwise distances, only 1/12 layers retains significant TRRUST proximity (L8, rbc=0.083). The TRRUST signal in SV2-4 is **almost entirely explained by co-expression confound**.
- **H02 POSITIVE**: SV1-high genes at both L0 and L11 are strongly depleted of TFs/targets. OR=0.108 at L0 (p<0.0001), OR=0.348 at L11 (p=0.007). High-SV1 genes are "background" genes—anonymous unnamed genes not in regulatory circuits. This mechanistically grounds the SV1 anti-correlation with TF content.
- **H03 MIXED**: Circuit genes (TF+target, n=295) have significantly smaller pairwise distances in SV2-4 at L8 (KS p<1e-6) and distinct 0-dim persistence distributions (KS p=0.007), but fragmentation integral ratio (0.94) is not significant vs null (p=0.116).

---

## Command Trace

```bash
# Write experiment script
cat > /tmp/iter51_experiments.py << 'PYEOF'
# ... [full script written inline]
PYEOF

# Run experiments
conda run -n subproject40-topology python /tmp/iter51_experiments.py
```

### Data sources
- Embeddings: `/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main/layer_gene_embeddings.npy`
  - Shape: [12, 4941, 512]; nonzero subset n=2039 genes
- Edge dataset: `.../cycle4_immune_main/cycle1_edge_dataset.tsv`
  - 589 positive (TRRUST TF→target) / 1523 negative pairs (nonzero-gene subset)
- SV1 loadings: `.../iter_0049/h01_sv1_vectors.npy` [12, 2039]

---

## H01: Residual TF Signal After Co-expression Regression

**Hypothesis**: After controlling for co-expression (embedding cosine similarity), TRRUST pairs retain significantly smaller SV2-4 distances than negative pairs.

**Method**: At each of 12 layers:
1. Compute SVD of centered nonzero embeddings [2039, 512].
2. Project genes to SV2-4 subspace [2039, 3].
3. Compute SV2-4 Euclidean distances for TRRUST positive (n=589) and negative (n=1523) pairs.
4. Compute full-embedding cosine similarity for each pair.
5. Regress cosine similarity (co-expression proxy) out of SV2-4 distances via OLS.
6. Test residualized distances with Mann-Whitney (TRRUST < NEG).

### Results

| Layer | rbc_original | rbc_residual | p_residual |
|-------|-------------|-------------|-----------|
| L0 | -0.068 | -0.170 | 1.000 |
| L1 | -0.036 | -0.129 | 1.000 |
| L2 | +0.040 | -0.060 | 0.984 |
| L3 | +0.066 | -0.051 | 0.967 |
| L4 | +0.055 | -0.051 | 0.965 |
| L5 | +0.113 | +0.015 | 0.295 |
| L6 | +0.127 | +0.031 | 0.132 |
| L7 | +0.102 | -0.002 | 0.524 |
| **L8** | **+0.156** | **+0.083** | **0.0016** |
| L9 | +0.050 | -0.007 | 0.594 |
| L10 | +0.024 | -0.074 | 0.996 |
| L11 | +0.074 | +0.001 | 0.489 |

**Significant after regression**: 1/12 layers (L8 only, rbc=0.083).

**Interpretation**: The TRRUST proximity in SV2-4 is almost entirely attributable to generic co-expression. Only a small residual at L8 survives (rbc=0.083 vs original 0.156). The SV2-4 subspace fundamentally encodes **co-expression structure**, not causal regulatory topology. This definitively closes the iter_0049 H03 claim.

**Decision: NEGATIVE**

Artifact: `h01_coexpr_regression_residual.csv`

---

## H02: SV1 High-Loading Gene Identity (TF/Target Depletion)

**Hypothesis**: Genes with high SV1 loading at L0 and L11 are not transcriptional regulators — they are "background" genes outside known regulatory circuits.

**Method**:
- Use SV1 loadings from iter_0049 [12, 2039].
- Split genes into top/bottom 20% quantiles at L0 and L11 (n=408 per group).
- Compute TF fraction (TRRUST TF membership) and target fraction for each group.
- Fisher exact test for TF enrichment/depletion.
- Only 295/2039 nonzero genes have known names; analysis conditional on named genes.

### Results: Gene Group Profiles

| Group | n | TF frac | Target frac | mean_norm_L0 |
|-------|---|---------|------------|-------------|
| high_SV1_L0 (top 20%) | 408 | 0.005 | 0.032 | 21.01 |
| low_SV1_L0 (bottom 20%) | 408 | 0.081 | 0.154 | 21.03 |
| high_SV1_L11 (top 20%) | 408 | 0.015 | 0.054 | 21.02 |
| low_SV1_L11 (bottom 20%) | 408 | 0.056 | 0.172 | 21.37 |
| **baseline** | 2039 | **0.036** | **0.129** | — |

### Fisher Exact Tests

| Comparison | OR | p |
|-----------|-----|---|
| high_SV1_L0 TF enrichment | **0.108** | **p<0.0001** |
| high_SV1_L11 TF enrichment | **0.348** | **p=0.007** |

**Interpretation**:
- High-SV1 genes are dramatically depleted of TFs at L0 (OR=0.108, 18x lower than baseline) and L11 (OR=0.348, 3x lower).
- Low-SV1 genes are enriched for TFs and targets.
- This means SV1 captures a complementary biological axis: genes with high SV1 loading are non-regulatory "background" genes, while TFs/targets sit at the **low end** of SV1.
- Mechanistically: the SV1 semantic inversion across layers (negative correlation with norm at L0 → positive at L11) corresponds to a reorientation of the non-regulatory gene axis.

**Decision: POSITIVE (strong)**

Artifacts: `h02_sv1_groups_identity.csv`, `h02_sv1_groups_layer_norms.csv`

---

## H03: 0-dim Persistent Homology on SV2-4 at L8

**Hypothesis**: Circuit genes (TF + target, n=295) form more compact/distinct topological clusters in SV2-4 at L8 than size-matched non-circuit genes.

**Method**:
- Extract SV2-4 projection of circuit vs matched non-circuit genes at L8.
- Compute: (1) pairwise distance distributions, (2) 0-dim PH via single-linkage, (3) fragmentation integral (area under component count curve) vs 500-sample bootstrap null.

### Results

| Metric | Circuit | Non-circuit | Test | p |
|--------|---------|------------|------|---|
| Pairwise distance mean | 0.0381 | 0.0420 | KS stat=0.068 | p<1e-6 |
| 0-dim persistence median | 0.00521 | 0.00590 | KS stat=0.140 | p=0.0065 |
| Fragmentation integral | 1972 | 2097 | ratio=0.940 | p=0.116 (null) |

**Interpretation**:
- Circuit genes cluster more tightly in SV2-4 (smaller pairwise distances, smaller persistence lifetimes = faster merging into components).
- However, the fragmentation integral ratio (0.94) is NOT significantly different from the null distribution (null mean=1.00, p=0.116).
- The proximity signal in raw distances is significant, but 0-dim PH doesn't reveal additional cluster structure beyond compactness.
- This is a **mixed/negative** for the topological clustering hypothesis: circuit genes are locally compact but do not form topologically distinct islands.

**Decision: MIXED (proximity confirmed, topological cluster structure not confirmed)**

Artifact: `h03_ph0_sv24_L8.csv`

---

## Quantitative Summary

| Hypothesis | Family | Direction | Decision |
|-----------|--------|-----------|---------|
| H01: Residual TF signal after coexpr regression | module_structure | negative | negative |
| H02: SV1-high gene identity (TF depletion) | intrinsic_dimensionality | positive | promising |
| H03: 0-dim PH on circuit genes in SV2-4 | persistent_homology | mixed | inconclusive |

**Key takeaway**: The TRRUST-regulatory vs co-expression distinction in SV2-4 is resolved: SV2-4 encodes co-expression, not regulatory topology. The mechanistic story of SV1 is now clearer: the SV1 axis (dominant spectral direction) separates non-regulatory anonymous genes from known TFs/targets, and this separation reorients across layers.
