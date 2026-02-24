# Executor Iteration Report: iter_0047

## Summary

Three hypotheses tested. Critical validity correction made to prior iter_0046 claims.

**Key findings**:
1. **H01 (POSITIVE, CORRECTS PRIOR)**: Zero-norm genes (2902/4941 = 58.7%) are the SAME set at ALL 12 layers — they are structurally absent from embeddings, not dropped progressively. TwoNN on nonzero-only (n=2039) is VALIDATED: L0=34.7 → L11=20.0 (ratio=0.58). SVD eff_rank CORRECTED: prior 14× collapse was an artifact of zero vectors dominating the global PCA. Real collapse (nonzero-only): L0=236.3 → L11=48.4 (4.9×, still real but much weaker than claimed).
2. **H02 (PARTIAL/NEGATIVE)**: SV1 explains only 18.6% of variance in nonzero-only genes at L11 (not 93.4%). Prior claim was inflated by the zero block collapsing onto a single axis. Gene names file missing, so biological annotation blocked.
3. **H03 (NEGATIVE)**: ID compression is approximately linear across layers (R²=0.986, slope=-1.30/layer). No sharp breakpoint. Largest per-layer drops at L3 (Δ=-2.43) and L9 (Δ=-1.94) but no evidence of phase transition.

**Critical correction from this iteration**: The SVD eff_rank finding from iter_0046 (14× collapse, SV1=93.4% variance) was an artifact of including 2902 zero-norm gene embeddings. These collapse onto the origin and artificially dominate the low-rank signal. The real (nonzero-only) analysis gives:
- SVD eff_rank: 4.9× collapse (L0=236 → L11=48), still significant but not dramatic
- SV1: 18.6% variance (not 93.4%), spectrum is distributed across ~5-10 directions

**Validated claim (intact)**: TwoNN intrinsic dimensionality compression L0=34.7 → L11=20.0 (ratio=0.58) on nonzero-only genes, consistent with prior cross-seed results.

---

## Command Trace

```bash
conda run -n subproject40-topology python /tmp/iter47_experiments.py
```

**Embedding data**:
- `cycle4_immune_main`: `/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_38_geometric_residual_stream_interpretability/implementation/outputs/cycle4_immune_main/layer_gene_embeddings.npy` — shape [12, 4941, 512]
- Gene names file: NOT FOUND at expected path (`gene_names.npy`)

---

## H01: Zero-Norm Validity Check

**Method**: Compute L2 norm of each gene's embedding vector at each layer. Identify zero-norm genes (norm < 1e-8). Re-run TwoNN and SVD on nonzero-only genes.

### Zero-norm statistics (constant across all layers)

| Layer | N_zero | Pct_zero | Mean_norm_nonzero |
|-------|--------|----------|-------------------|
| L00   | 2902   | 58.7%    | 21.17             |
| L01   | 2902   | 58.7%    | 20.69             |
| ...   | 2902   | 58.7%    | ~21-24            |
| L11   | 2902   | 58.7%    | 22.10             |

**Critical finding**: The same 2902 genes have zero norm at ALL 12 layers. This is not dropout or progressive sparsification — these genes were never embedded (out-of-vocabulary or structurally excluded). The valid embedding set is n=2039 throughout.

### TwoNN on nonzero-only genes (n=2039)

| Layer | ID_all_genes | ID_nonzero_only | ID_all_layers_nz |
|-------|-------------|-----------------|------------------|
| L00   | 32.57       | 33.67           | 34.68            |
| L01   | 31.66       | 32.86           | 33.27            |
| L02   | 29.60       | 33.01           | 32.18            |
| L03   | 27.68       | 30.69           | 30.34            |
| L04   | 26.56       | 28.43           | 27.91            |
| L05   | 26.71       | 27.10           | 26.96            |
| L06   | 27.80       | 26.68           | 26.44            |
| L07   | 23.97       | 25.40           | 25.27            |
| L08   | 21.15       | 24.81           | 24.69            |
| L09   | 22.80       | 23.48           | 22.98            |
| L10   | 20.33       | 21.46           | 21.03            |
| L11   | 17.93       | 19.89           | 19.98            |

**TwoNN compression validated**: Monotone decrease from ~34.7 → ~20.0 (ratio=0.58) holds for nonzero-only genes. Zero genes slightly inflated the all-gene ID at later layers (mixing zero-cluster with real data). Compression signal is genuine.

### SVD eff_rank: all vs nonzero-only

| Layer | Eff_rank (all genes) | Eff_rank (nonzero only) | Top1 (all) | Top1 (nonzero) |
|-------|---------------------|------------------------|------------|----------------|
| L0    | 23.22               | 236.26                 | 0.539      | 0.056          |
| L11   | 1.65                | 48.38                  | 0.934      | 0.187          |

**CORRECTION**: The extreme SVD collapse (14×, SV1=93.4%) was driven by zero vectors:
- With all genes: zeros collapse onto origin → single dominant direction → eff_rank artificially low
- With nonzero only: L0=236 → L11=48 (4.9× real collapse, top1 rises from 5.6% → 18.7%)
- Still a real effect (4.9× collapse is significant), but not the 14× dramatic claim

**Decision: PROMISING (validates TwoNN; corrects SVD claim)**

---

## H02: SV1 Identity

**Method**: SVD of centered nonzero-gene embeddings at L11. Measure SV1/SV2 variance fractions. Attempt biological annotation via gene names (blocked — file missing).

### Results

- SV1 explains **18.6%** of variance in nonzero-only genes at L11
- SV2 explains **11.7%** of variance
- Spectrum is distributed, not dominated by one direction

Prior claim (93.4% for SV1) was an artifact of including zero vectors which all project onto a single direction.

**Gene names file not found**: Biological annotation of SV1/SV2 blocked. Cannot determine if SV1 corresponds to expression level, cell-type axis, or other signal.

**Decision: NEGATIVE/PARTIAL — prior finding invalidated. Biological annotation needed.**

---

## H03: ID Compression Breakpoint Detection

**Method**: TwoNN at each of 12 layers on nonzero-only genes (n=2039). Piecewise linear fit to identify breakpoint layer. Finite differences for per-layer drop sizes.

### Results

| Layer | TwoNN ID (nonzero) | Δ from prev |
|-------|-------------------|-------------|
| L00   | 34.68             | —           |
| L01   | 33.27             | -1.41       |
| L02   | 32.18             | -1.09       |
| L03   | 30.34             | -1.83       |
| L04   | 27.91             | -2.43       |
| L05   | 26.96             | -0.94       |
| L06   | 26.44             | -0.53       |
| L07   | 25.27             | -1.17       |
| L08   | 24.69             | -0.57       |
| L09   | 22.98             | -1.72       |
| L10   | 21.03             | -1.94       |
| L11   | 19.98             | -1.05       |

- Global linear fit: **R²=0.986**, slope=-1.30 ID units/layer
- Best piecewise breakpoint at L5 (slopes: -1.61 before, -1.22 after) — marginal
- Largest drops at L3-L4 (early) and L9-L10 (late), but no clear phase transition

**Decision: NEGATIVE — compression is smooth and linear, no discrete phase change.**

---

## Artifact paths

- `iterations/iter_0047/h01_zero_norm_validity.json`
- `iterations/iter_0047/h02_sv1_identity.json`
- `iterations/iter_0047/h03_id_breakpoint.json`
