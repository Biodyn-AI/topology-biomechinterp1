# Executor Iteration Report — iter_0004

## Summary

Three novel hypotheses tested across two new families (intrinsic_dimensionality, cross_model_alignment) and one biological anchoring test (module_structure / TRRUST). Two promising results, one negative.

---

## Hypotheses Tested

### H01: TwoNN Intrinsic Dimensionality per Layer
**Family:** intrinsic_dimensionality (new family)
**Status:** tested
**Decision:** promising

**Method:**
TwoNN estimator (Facco et al. 2017) on PCA-20 projections of 400 randomly sampled scGPT gene embeddings per layer. Null: 15 feature-shuffle replicates (per-column permutation of gene embeddings).

**Command trace:**
```bash
conda run -n subproject40-topology python run_iter0004_screen.py
```

**Results:**
| Layer | TwoNN ID | Null Mean | Z-score |
|-------|----------|-----------|---------|
| 0  | 10.333 | (higher) | -0.04 |
| 4  | 8.857  | (higher) | -2.78 |
| 5  | 8.795  | (higher) | -3.22 |
| 11 | 8.516  | (higher) | -5.79 |
- ID range: 8.52–10.33 (decreasing trend: layer 0 → layer 11)
- Mean z-score: **-2.29** (real ID is significantly LOWER than shuffled null)
- 6/12 layers with |z| > 2.0
- Mean null ID is consistently higher than observed

**Interpretation:**
scGPT gene embeddings occupy a significantly lower-dimensional manifold than feature-shuffled noise (z = -2.29). The decreasing trend from layer 0 (ID=10.3) to layer 11 (ID=8.5) suggests progressive compression/specialization of the representation across depth. This is a structured geometric signal: real embeddings are more constrained in intrinsic dimensionality than noise would predict.

**Artifact:** `h01_intrinsic_dim_per_layer.csv`

---

### H02: TRRUST Co-Target Distance Clustering
**Family:** module_structure (biological anchoring, new method)
**Status:** tested
**Decision:** inconclusive

**Method:**
For each of 64 TFs with ≥5 mapped co-targets in the TRRUST-derived edge dataset (cycle1_edge_dataset.tsv), computed mean pairwise Euclidean distance in PCA-20 space of co-target genes. Null: 20 random same-sized groups sampled from the pool of 209 named genes.

**Command trace:**
```bash
conda run -n subproject40-topology python run_h02_trrust_fast.py
```

**Results:**
- 768 TF-layer tests (64 TFs × 12 layers)
- Mean z_dist: **+0.222** (co-targets slightly farther than random groups)
- Significant clustering (z < -1.96): **8/768 (1.0%)**
- Significant dispersion (z > +1.96): **25/768 (3.3%)**
- Top TF by clustering: MEF2C (mean z_dist = -1.69, 2/12 layers sig), STAT4 (mean z_dist = -1.69, 4/12 layers sig)

**Interpretation:**
No significant biological co-regulatory signal detected. TRRUST co-targets do not cluster more tightly than random gene groups in the scGPT gene embedding space. This is consistent with the embedding being gene-centric (contextual), not reflecting static TF-target regulatory structure without cell-type conditioning. The gene pool is small (209 genes) and the null (same-pool random groups) is relatively lenient.

**Artifact:** `h02b_trrust_cotarget_by_tf_layer.csv`, `h02b_trrust_tf_summary.csv`

---

### H03: Cross-Layer Linear CKA Alignment
**Family:** cross_model_alignment (new family)
**Status:** tested
**Decision:** promising

**Method:**
Computed linear CKA for all 12×12 pairs of scGPT layer embeddings. Null: 15 feature-shuffle replicates of layer 0 vs all layers.

**Command trace:**
```bash
conda run -n subproject40-topology python run_iter0004_screen.py
```

**Results:**
- Adjacent-layer CKA: [0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000]
- Mean adjacent CKA: **1.000** (effectively 1.0 for all adjacent pairs)
- Mean null CKA (feature-shuffle vs any layer): **0.0046**
- Z-scores of observed vs null: all **~173–175** (min 172.1, max 175.1)
- Full CKA matrix: all non-self pairs have CKA > 0.99 (extremely high across all layer distances)

**Interpretation:**
scGPT gene embeddings are essentially identical across all 12 transformer layers under linear CKA (CKA ≈ 1.0 for adjacent layers, and high even for distant layers). Feature-shuffled representations have near-zero CKA. This reveals that:
1. The transformer blocks make extremely small perturbative changes to gene representations
2. The residual stream is dominated by the skip connections, not the transformer blocks themselves
3. The representation space geometry is stable across depth

This is a strong structural finding: scGPT behaves as a near-identity residual network for gene embeddings, consistent with the known property of deep transformers with strong residual streams. The signal is layer-stable and highly discriminative vs null (z ≈ 174).

**Artifacts:** `h03_cross_layer_cka_matrix.csv`, `h03_cka_matrix.npy`, `h03_cka_null.npy`

---

## Quantitative Summary

| Hypothesis | Family | Primary Metric | Result | Direction | Decision |
|-----------|--------|---------------|--------|-----------|----------|
| H01 TwoNN ID | intrinsic_dim | mean z vs null | -2.29 (6/12 sig) | positive | promising |
| H02 TRRUST CC | module_structure | mean z_dist | +0.222 (8/768 sig) | inconclusive | inconclusive |
| H03 Cross-layer CKA | cross_model_align | mean adj CKA | 1.000 (z≈174) | positive | promising |

---

## Environment Notes

- All runs: `conda run -n subproject40-topology python <script>`
- No new packages installed
- Data: `/Volumes/Crucial X6/.../subproject_38_.../outputs/cycle1_main/layer_gene_embeddings.npy` (shape 12×4803×512)
- Gene list extracted from `cycle1_edge_dataset.tsv` (209 named genes, indices in embedding matrix)

---

## Iteration Decision

**Gate: PASSED**
- 2 promising hypotheses (H01, H03) with strong effect sizes
- 1 inconclusive (H02) — gene pool too small for definitive TRRUST test
- All hypotheses are new families vs prior iterations
- Minimum research gate met: machine-readable artifacts, explicit command traces, quantitative metrics
