# Executor Iteration Report — iter_0031

## Summary

Three hypotheses tested on the corrected **195 in-vocabulary gene** universe (excluding 14 OOV genes discovered in iter_0030). Two strong positive findings and one definitive null.

**Key findings:**
1. **H01 (spectral gap):** Spectral gap ratio decreases monotonically from L0 → L11 (Spearman rho = −1.000, p = 0). Real embedding mean gap = 0.067 vs Gaussian null = 0.373. The kNN graph becomes increasingly modular/fragmented as depth increases.
2. **H02 (STRING → embedding distance):** STRING interaction score does NOT predict embedding proximity among 195 in-vocab genes (mean Spearman rho = −0.015, AUROC = 0.494). Clean null result on the corrected gene set.
3. **H03 (participation ratio):** Participation ratio (effective dimensionality) decreases monotonically from PR≈58 (L0) to PR≈9.5 (L11) with Spearman rho = −1.000. Variance in PC1 rises from 8% to 26%. The manifold collapses to a small number of dominant directions across layers.

---

## Command Trace

```bash
# Run main screen (H01 + H03)
conda run -n subproject40-topology python iterations/iter_0031/run_iter0031_screen.py

# H02 inline correction (STRING cache format fix)
conda run -n subproject40-topology python3 -c "... H02 STRING inline script ..."
```

All outputs written to `iterations/iter_0031/`.

---

## H01: kNN Spectral Gap — 195 In-Vocab Genes

**Family:** graph_topology
**Method:** k=10 kNN graph per layer on 195 in-vocab genes. Normalized Laplacian eigenvalues via `scipy.sparse.linalg.eigsh`. Spectral gap ratio = λ₂ / λ_max. Compared to Gaussian null (same shape, random embeddings).
**Artifact:** `h01_spectral_gap_195.json`

**Results per layer:**
| Layer | λ₂     | λ_max  | Ratio  |
|-------|--------|--------|--------|
| L00   | 0.2061 | 1.4221 | 0.1449 |
| L03   | 0.1120 | 1.4042 | 0.0798 |
| L06   | 0.0633 | 1.3982 | 0.0453 |
| L09   | 0.0477 | 1.3891 | 0.0343 |
| L11   | 0.0314 | 1.3576 | 0.0231 |

- Spearman rho(layer, ratio) = **−1.000**, p = 0.0 (perfect monotonic decrease)
- Real mean ratio = **0.067** vs Gaussian null mean = **0.373** (5.6× difference)
- Interpretation: The 195 in-vocab gene kNN graph has far more modular structure than random (spectral gap << null), and this modular fragmentation deepens with each layer.

**Decision:** positive

---

## H02: STRING Score → Embedding Distance (195 In-Vocab Genes)

**Family:** manifold_distance
**Method:** Load STRING PPI pairs (score ≥ 0.4) from iter_0015 cache. Filter to 195 in-vocab named genes. Compute L2 embedding distance per pair per layer. Spearman rho(STRING_score, L2_dist) and AUROC (high ≥ 0.7 vs low < 0.5 score pairs).
**Artifact:** `h02_string_dist_195.json`

**Results:**
- N valid pairs = 2592; high score (≥0.7) = 796; low score (<0.5) = 830
- Mean Spearman rho across layers = **−0.015** (range: −0.026 to −0.010)
- Mean AUROC (high vs low) = **0.494** (essentially 0.50)
- Layer trend Spearman rho = −0.755, p = 0.004 (slight improvement at late layers but effect size tiny)

**Interpretation:** STRING interaction score does not predict embedding proximity in the scGPT residual stream for these 195 named genes. The correlation is essentially zero at all layers. This is a definitive null for STRING → embedding geometry.

**Decision:** negative

---

## H03: Participation Ratio + L2 Norm Distribution

**Family:** intrinsic_dimensionality
**Method:** For 195 in-vocab genes per layer: mean-center, SVD, participation ratio = (Σs²)² / Σs⁴. Also compute L2 norm mean, std, CV across genes per layer. Spearman rho(layer, PR) and rho(layer, norm_CV).
**Artifact:** `h03_intdim_norm_195.json`

**Results:**
| Layer | PR    | var_PC1 | norm_mean | norm_CV |
|-------|-------|---------|-----------|---------|
| L00   | 58.09 | 0.080   | 21.09     | 0.017   |
| L03   | 44.25 | 0.096   | 21.13     | 0.008   |
| L06   | 28.86 | 0.120   | 21.73     | 0.009   |
| L09   | 14.84 | 0.175   | 23.61     | 0.008   |
| L11   | 9.53  | 0.259   | 22.10     | 0.004   |

- Spearman rho(layer, PR) = **−1.000**, p = 0.0 (perfect monotonic decrease)
- PR collapse: 58.09 → 9.53 (6.1× reduction across 12 layers)
- var_PC1 increases: 8% → 26% (more variance concentrated in PC1 at late layers)
- norm_CV is very small and slightly decreasing (genes become more norm-homogeneous)

**Interpretation:** The effective dimensionality of the in-vocab gene embedding manifold collapses monotonically across layers. Combined with H01 (growing modularity), the manifold is simultaneously becoming lower-dimensional AND more fragmented — consistent with the model organizing genes into distinct clusters in a low-dimensional subspace.

**Decision:** positive

---

## Artifact Paths

- `iterations/iter_0031/h01_spectral_gap_195.json`
- `iterations/iter_0031/h02_string_dist_195.json`
- `iterations/iter_0031/h03_intdim_norm_195.json`
- `iterations/iter_0031/run_iter0031_screen.py`

---

## OOV Correction Note

All analyses use the 195 in-vocabulary genes (excluding 14 OOV genes: FOS, HLA-A, HLA-DPB1, JUNB, KLF6, LDHA, LGALS1, NCAM1, NCOA3, NR4A3, PAX5, PTGS2, TBXAS1, TNF). Prior iterations (≤ iter_0029) analyzed 209 genes including OOV; results from those iterations involving the 14-gene diverging cluster should be interpreted with caution.
