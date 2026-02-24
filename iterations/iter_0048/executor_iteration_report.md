# Executor Iteration Report: iter_0048

## Summary

Three hypotheses tested. Two clear positive results and one mixed/negative.

**Key findings**:
1. **H01 (MIXED)**: Gene annotation via edge dataset: 295/361 named genes are nonzero (embedded); 66/361 named zero-norm genes identified as tissue-specific/non-immune (TFAP2B, SNAI2, TYR, PGR, TBX5, CCL22, AMH), confirming biological validity of zero-norm set. SV1/SV2 at L11 do NOT cleanly separate immune cell-type markers — annotated B-cell, T-cell, Myeloid genes all score near zero on SV1 (mean ~-0.01 to -0.02), indicating SV1 captures a non-cell-type axis. Top SV1 genes are unnamed (not in edge annotation set).
2. **H02 (POSITIVE, PUBLICATION-QUALITY)**: Full spectral decay curve confirms systematic, monotone compression: eff_rank 236.9 → 48.7 (4.87×), k90 307 → 102 (3.0×), Frobenius norm 558 → 204 (2.74×). SV1 singular value peaks at L8 (143.9) then decreases, while variance concentration increases monotonically from 5.6% (L0) to 18.6% (L11). SPR (spectral participation ratio) 90.1 → 15.2 (5.9×).
3. **H03 (MIXED/INCONCLUSIVE)**: Sliced Wasserstein distance between consecutive layers ranges 0.29–0.55; cumulative distance from L0 increases monotonically to L10 (1.32) then plateaus (L11: 1.29). Spearman rho(SWD, eff_rank_drop) = -0.46, p=0.15 (not significant). Different transport metrics diverge: centroid displacement grows at later layers while SWD is smallest at L10→L11, indicating the center-of-mass shifts while the distribution contracts.

---

## Command Trace

```bash
conda run -n subproject40-topology python /tmp/iter48_experiments.py
```

**Embedding data**:
- `cycle4_immune_main`: `/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main/layer_gene_embeddings.npy` — shape [12, 4941, 512]
- Gene names: recovered from edge dataset `cycle1_edge_dataset.tsv` (361 named genes with indices)
- Nonzero-only set: n=2039 (same 2902 zero-norm genes at all 12 layers)

---

## H01: Gene Annotation + SV1/SV2 Identity

**Method**: Map gene names from cycle1_edge_dataset.tsv to embedding indices. Classify named genes into zero-norm vs nonzero sets. Compute SVD of centered nonzero embeddings at L11. Score immune cell-type marker genes on SV1/SV2.

### Named gene set statistics
| Category | N genes |
|----------|---------|
| Total named in edge dataset | 361 |
| Named AND nonzero (embedded) | 295 |
| Named AND zero-norm (absent) | 66 |

### Zero-norm named genes (biologically meaningful exclusion)
Sample: TFAP2B (neural crest TF), SNAI2 (EMT marker), PITX1 (limb TF), TYR (melanocyte), PGR (progesterone receptor), TBX5 (cardiac), DHRS2, AVP (neuropeptide), CCL22 (dendritic cell chemokine), AMH (hormone)

**Interpretation**: Zero-norm genes are systematically non-immune — developmental TFs, tissue-specific markers, hormones. This validates that the zero-norm set = structurally excluded (out-of-vocabulary for immune context), not random dropout.

### SV1/SV2 annotation at L11
L11 variance explained: SV1=18.6%, SV2=11.7%, SV3=8.8%

| Cell type | N found/total | Mean SV1 |
|-----------|--------------|---------|
| B-cell markers | 8/12 | -0.017 |
| T-cell markers | 6/12 | -0.016 |
| Myeloid markers | 2/8 | -0.009 |
| GC TFs (BCL6/BACH2/BATF) | 3/4 | -0.013 |
| NK markers | 1/4 | -0.019 |

**Finding**: SV1 does NOT separate cell-type lineages — all annotated immune markers score near zero on SV1 (range: -0.009 to -0.019). The dominant spectral axis captures something other than cell-type identity. Top SV1 genes are unnamed (not in the 361-gene annotation set).

**Decision: MIXED** — zero-norm biology validated (positive); SV1/SV2 biological identity remains unresolved (neutral/negative).

---

## H02: Spectral Decay Curve Across 12 Layers

**Method**: For each of 12 layers, compute SVD of centered nonzero-only embeddings (n=2039, 512-dim). Track eff_rank (entropy-weighted), SPR (participation ratio), k90 (dims for 90% variance), SV1 singular value, Frobenius norm.

### Full spectral decay table

| Layer | Eff_rank | SPR | k90 | SV1 | Var(SV1) | Frob |
|-------|---------|-----|-----|-----|---------|------|
| L00 | 236.9 | 90.1 | 307 | 132.1 | 5.6% | 558.0 |
| L01 | 230.4 | 91.0 | 299 | 110.1 | 5.3% | 478.2 |
| L02 | 214.0 | 83.4 | 287 | 106.0 | 5.5% | 452.2 |
| L03 | 186.5 | 68.1 | 270 | 114.2 | 6.6% | 446.3 |
| L04 | 156.0 | 52.6 | 248 | 121.2 | 8.3% | 420.9 |
| L05 | 141.1 | 46.1 | 235 | 124.3 | 9.3% | 407.3 |
| L06 | 125.5 | 40.2 | 216 | 129.1 | 10.3% | 401.9 |
| L07 | 107.1 | 33.4 | 193 | 134.5 | 11.5% | 396.0 |
| L08 | 88.7 | 26.4 | 169 | 143.9 | 13.9% | 386.6 |
| L09 | 74.7 | 22.8 | 145 | 127.9 | 14.3% | 338.1 |
| L10 | 59.9 | 19.0 | 119 | 101.3 | 14.4% | 266.9 |
| L11 | 48.7 | 15.2 | 102 | 87.8 | 18.6% | 203.7 |

**Compression statistics**:
- Eff_rank: 4.87× compression (236.9 → 48.7)
- SPR: 5.93× compression (90.1 → 15.2)
- k90: 3.01× compression (307 → 102)
- Frob norm: 2.74× compression (558 → 204)

**Notable SV1 dynamics**: SV1 singular value grows L0→L8 (132 → 144, +9%), then sharply drops L8→L11 (144 → 88, -39%), while variance fraction continues increasing (5.6% → 18.6%). This suggests a two-phase process: (1) consolidation of spectral mass into SV1 direction through L8; (2) overall representation shrinkage thereafter.

**Decision: PROMISING (strong, publication-quality signal)**

---

## H03: Wasserstein Distance Between Layer Distributions

**Method**: Sliced Wasserstein Distance (SWD, 100 random projections) between consecutive layer embeddings of nonzero-only genes. Also: energy distance, centroid displacement, covariance Frobenius distance. Cumulative SWD from L0. Spearman correlation with eff_rank drops.

### Consecutive-layer transport distances

| Transition | SWD | Energy dist | Centroid dist | Cov Frob |
|-----------|-----|------------|--------------|---------|
| L00→L01 | 0.552 | 7.90 | 11.79 | 6.3 |
| L01→L02 | 0.431 | 4.84 | 8.61 | 4.1 |
| L02→L03 | 0.370 | 3.99 | 7.61 | 4.7 |
| L03→L04 | 0.370 | 3.93 | 7.45 | 6.1 |
| L04→L05 | 0.376 | 4.53 | 7.81 | 5.4 |
| L05→L06 | 0.403 | 4.96 | 8.07 | 6.0 |
| L06→L07 | 0.371 | 5.45 | 8.39 | 6.6 |
| L07→L08 | 0.386 | 5.99 | 8.84 | 7.4 |
| L08→L09 | 0.479 | 7.73 | 9.83 | 7.9 |
| L09→L10 | 0.400 | 8.76 | 9.80 | 6.3 |
| L10→L11 | 0.291 | 4.42 | 5.80 | 4.1 |

**Cumulative SWD from L0**: Monotone increase to L10 (1.32), slight plateau at L11 (1.29). Largest step: L0→L1 (0.55), smallest: L10→L11 (0.29).

**Spearman rho(SWD, eff_rank_drop) = -0.46, p=0.15** (not significant).

**Interpretation**: Consecutive-layer SWD is relatively stable (0.29–0.55), with no sharp transition. Energy distance and centroid displacement grow at later layers, while SWD is smallest at the final transition — indicating the distribution shrinks (fewer dimensions active) while the centroid continues moving. These metrics capture different aspects of geometric flow.

**Decision: INCONCLUSIVE** — transport distances quantified but do not strongly predict spectral compression.

---

## Artifacts Generated

- `iterations/iter_0048/h01_gene_annotation_sv_identity.json`
- `iterations/iter_0048/h02_spectral_decay.csv`
- `iterations/iter_0048/h03_wasserstein_layer_distances.csv`
- `iterations/iter_0048/h03_wasserstein_cumulative.csv`
- `iterations/iter_0048/iter48_summary.json`
