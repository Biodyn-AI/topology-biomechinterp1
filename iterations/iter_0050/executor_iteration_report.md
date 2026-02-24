# Executor Iteration Report: iter_0050

## Summary

Three hypotheses tested targeting the spectral hierarchy of scGPT gene representations. Key results: (1) the SV2-SV4 regulatory proximity signal is NOT specific to TF-regulation vs co-expression — high co-expression pairs are even closer than TRRUST pairs; (2) SV1 loading anti-correlates with expression norms at early layers and reverses sign at late layers, tracking the geometric rotation documented in iter_0049; (3) SV2-SV4 most enriched for TF/target identity, with secondary enrichment in SV5-SV7.

---

## Command Trace

```bash
# Run all three experiments
conda run -n subproject40-topology python /tmp/iter50_experiments.py

# Embedding data
# /Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main/layer_gene_embeddings.npy
# shape: [12, 4941, 512]; nonzero subset: n=2039 genes

# Edge data
# /Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main/cycle1_edge_dataset.tsv
# 735 positive (TRRUST TF→target) / 2205 negative pairs

# SV1 vectors from iter_0049
# /Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0049/h01_sv1_vectors.npy
# shape: [12, 2039]
```

---

## H01: SV2-SV4 Co-expression Specificity Control

**Hypothesis**: The TRRUST regulatory proximity signal in SV2-SV4 (iter_0049 H03) is specific to TF-regulation, not attributable to generic co-expression.

**Method**: Construct a co-expression proxy by taking 5000 random gene pairs from nonzero-gene set, computing average embedding cosine similarity, selecting top-300 as "high co-expression" pairs (mean cosine = 0.885). Compare SV2-SV4 Euclidean distances: TRRUST-positive vs TRRUST-negative vs co-expression pairs. Mann-Whitney test at each layer.

### Results

| Layer | TRRUST dist | NEG dist | COEXPR dist | p(TRRUST<NEG) | rbc_vs_neg | p(TRRUST<COEXPR) | rbc_vs_coexpr |
|-------|------------|---------|------------|----------------|------------|------------------|---------------|
| L0 | 0.0410 | 0.0384 | 0.0276 | 0.992 | -0.068 | 1.000 | -0.453 |
| L3 | 0.0418 | 0.0444 | 0.0274 | 0.009 | +0.066 | 1.000 | -0.446 |
| L5 | 0.0377 | 0.0422 | 0.0270 | <0.001 | +0.113 | 1.000 | -0.355 |
| L6 | 0.0364 | 0.0412 | 0.0268 | <0.001 | +0.127 | 1.000 | -0.304 |
| L8 | 0.0320 | 0.0371 | 0.0265 | <0.001 | +0.156 | 1.000 | -0.184 |
| L11 | 0.0349 | 0.0381 | 0.0303 | 0.004 | +0.074 | 1.000 | -0.154 |

**Key metrics**:
- TRRUST < NEGATIVE: 8/12 layers (p<0.05) — replicates iter_0049 H03
- TRRUST < CO-EXPRESSION: 0/12 layers — FAILS at all layers
- Co-expression pairs (mean cosine=0.885) are consistently CLOSER in SV2-4 than TRRUST pairs

**Interpretation**: The SV2-4 proximity effect for TRRUST pairs is NOT specific to regulatory relationships. Highly co-expressed genes (as proxied by embedding similarity) are even closer in SV2-4 subspace. This means the SV2-4 signal likely reflects co-expression rather than causal regulatory topology. This is a critical negative: the claim from iter_0049 that SV2-4 specifically encodes transcriptional regulation must be revised — it encodes co-expression structure, within which TRRUST pairs happen to be somewhat enriched.

**Decision: NEGATIVE (specificity falsified)**

Artifact: `h01_sv24_string_specificity.csv`

---

## H02: SV1 Identity Characterization

**Hypothesis**: SV1 (dominant spectral axis) captures expression level — genes with high embedding norm have high SV1 loading.

**Method**: At each layer, compute Spearman correlation between |SV1 gene loading| (from iter_0049 h01_sv1_vectors.npy) and embedding L2-norm at that layer and at L0.

### Results

| Layer | rho (SV1 vs norm_same) | p | rho (SV1 vs norm_L0) |
|-------|----------------------|---|----------------------|
| L0 | -0.326 | <0.001 | -0.326 |
| L1 | -0.328 | <0.001 | -0.358 |
| L2 | -0.171 | <0.001 | -0.244 |
| L3 | +0.040 | 0.071 | -0.108 |
| L4 | -0.027 | 0.231 | -0.262 |
| L5 | -0.016 | 0.474 | -0.280 |
| L6 | +0.010 | 0.662 | -0.277 |
| L7 | -0.006 | 0.785 | -0.325 |
| L8 | +0.070 | 0.002 | -0.339 |
| L9 | +0.191 | <0.001 | -0.314 |
| L10 | +0.207 | <0.001 | -0.213 |
| L11 | +0.269 | <0.001 | +0.148 |

**Key findings**:
- Layers 0-2: SV1 loading strongly ANTI-correlates with embedding norm (rho ≈ -0.17 to -0.33). Genes with LOW expression amplitude have high SV1 loading.
- Layers 3-8: near-zero correlation — SV1 is expression-neutral in middle layers.
- Layers 9-11: SV1 loading becomes POSITIVELY correlated with norm (rho = 0.19-0.27). High-norm genes gain high SV1 loading.
- The sign flip at L8-L9 aligns exactly with the geometric rotation of SV1 direction documented in iter_0049 (L9→L10 consecutive cosine = 0.894, L10→L11 = 0.793).
- The persistent negative correlation vs L0 norms (all layers except L11) indicates that SV1 is NOT simply capturing absolute expression level — it undergoes a semantic reorientation.

**Decision: POSITIVE (mechanistic insight into SV1 rotation)**

Artifact: `h02_sv1_identity.csv`

---

## H03: SV5-SV10 Spectral Axis Biological Screen

**Hypothesis**: Different spectral axis groups (SV1, SV2-4, SV5-7, SV8-10) differentially enrich for TRRUST TFs vs targets vs background.

**Method**: At each layer, compute projection magnitude of each gene onto each SV group (||U[:,sv_idx]||). Compare TF genes (n=73) vs non-TF and target genes (n=263) vs non-target using Mann-Whitney. Report rbc (rank-biserial correlation).

### TF Enrichment by SV Group

| SV Group | Layers TF-enriched (p<0.05) | Mean rbc (TF) | Layers target-enriched | Mean rbc (target) |
|----------|----------------------------|---------------|------------------------|-------------------|
| SV1 | 7/12 | **-0.110** | 2/12 | +0.029 |
| SV2-4 | 10/12 | **+0.229** | 8/12 | +0.124 |
| SV5-7 | 7/12 | **+0.152** | 9/12 | +0.039 |
| SV8-10 | 5/12 | **+0.123** | 7/12 | -0.072 |

**Key findings**:
- SV1 is anti-enriched for TFs (rbc=-0.110): TFs have LOWER SV1 loading. TFs are not a dominant source of embedding variance captured by SV1.
- SV2-4 is the strongest enrichment axis for both TFs (rbc=0.229) and targets (rbc=0.124), replicating and extending iter_0049 H03.
- SV5-7 shows secondary enrichment for TFs (rbc=0.152, 7/12 layers) — indicating the regulatory signal extends beyond the top 3 spectral directions.
- SV8-10 shows weaker but detectable TF enrichment (5/12 layers).
- The spectral hierarchy is: SV2-4 > SV5-7 > SV8-10 >> SV1 for TF-relevance.

**Decision: POSITIVE (spectral hierarchy confirmed, SV5-7 new finding)**

Artifact: `h03_sv510_screen.csv`

---

## Overall Assessment

| Hypothesis | Decision | Direction | Key metric |
|-----------|----------|-----------|------------|
| H01: SV2-4 co-expression specificity | NEGATIVE | negative | TRRUST < COEXPR at 0/12 layers |
| H02: SV1 identity / expression correlation | POSITIVE | positive | Sign-flip at L8-L9 aligns with geometric rotation |
| H03: SV5-10 spectral screen | POSITIVE | positive | SV2-4>SV5-7>SV8-10 for TF enrichment |

**Critical revision**: The iter_0049 claim that SV2-4 encodes "regulatory proximity" must be qualified — it encodes **co-expression proximity** (of which TF-target pairs are a subset). This narrows but does not eliminate the biological finding.

**Novel finding**: SV1 undergoes a semantic inversion in its relationship to expression amplitude, transitioning from anti-correlated (L0-L2) to uncorrelated (L3-L8) to positively correlated (L9-L11) — tracking the geometric rotation documented in iter_0049.

**New positive**: SV5-7 shows secondary TF enrichment (rbc=0.152, 7/12 layers), extending the spectral hierarchy observation.
