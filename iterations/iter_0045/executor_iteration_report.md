# Executor Iteration Report: iter_0045

## Summary
Three hypotheses tested. Primary positive finding: monotone intrinsic dimensionality (ID) compression
across all 12 scGPT layers (L0=32.57 → L11=18.05), far below ambient 512D space. kNN lineage purity
test negative (underpowered). TCR signaling circuit directional signal (CD28: L0=595→L6=86) but
statistically underpowered (p=0.32).

## Command Trace

```bash
# Embeddings: cycle4_immune_main/layer_gene_embeddings.npy [12, 4941, 512]
# Vocabulary: cycle4_immune_main/cycle1_edge_dataset.tsv (361 genes with known indices)
conda run -n subproject40-topology python /tmp/iter45_experiments.py
# Then additional nulls for H02:
conda run -n subproject40-topology python -c "[inline null computation for Gaussian + feature shuffle]"
```

## H01: kNN Lineage Purity — k=10 community structure

**Method**: For in-vocab lineage markers (B-cell n=3, T-cell n=2, Myeloid n=1, NK n=1), compute
fraction of 10-nearest neighbors from same lineage. Null: 500 label-shuffle permutations.

**Lineage sets used**:
- B-cell: MS4A1, CD79A, BLK (n=3 in vocab from 5 candidates)
- T-cell: CD8A, PRF1 (n=2 from 7 candidates; CD3E/CD3D/TRAC/TRBC1 not in vocab)
- Myeloid: LYZ (n=1)
- NK: FCGR3A (n=1)

**Results (mean purity across all markers)**:
| Layer | Mean Purity | Null Mean | p-value |
|-------|-------------|-----------|---------|
| L0    | 0.214       | 0.214     | 0.958   |
| L2    | 0.214       | 0.214     | 0.958   |
| L5    | 0.214       | 0.214     | 0.958   |
| L8    | 0.214       | 0.214     | 0.958   |
| L11   | 0.214       | 0.214     | 0.958   |

**Artifact**: `h01_knn_lineage_purity.json`

**Interpretation**: Completely flat purity = null baseline, p=0.96. NEGATIVE.
Root cause: with Myeloid (n=1) and NK (n=1) having 0 possible same-lineage neighbors, and only 2
T-cell markers, the test is fundamentally underpowered. Cannot detect lineage community structure
with <3 genes per lineage. Method valid; sample size invalid for this vocab.

## H02: Full-Manifold TwoNN Intrinsic Dimensionality

**Method**: TwoNN estimator (Facco et al. 2017) applied to n=2000 randomly sampled gene vectors from
all 4941 genes at each of 12 layers. Nulls: (a) standard Gaussian in 512D, (b) unit-sphere Gaussian,
(c) L11 embeddings with features shuffled column-wise.

**Results**:
| Layer | ID Estimate |
|-------|-------------|
| L0    | 32.57       |
| L1    | 31.25       |
| L2    | 29.84       |
| L3    | 28.01       |
| L4    | 26.32       |
| L5    | 24.66       |
| L6    | 23.99       |
| L7    | 22.96       |
| L8    | 21.72       |
| L9    | 20.82       |
| L10   | 19.42       |
| L11   | 18.05       |

**Null comparison**:
- Gaussian 512D (raw): ID = 122.97
- Gaussian 512D (unit sphere): ID = 135.53
- L11 feature-shuffled: ID = 18.05 (same as actual)

**Key metrics**:
- Delta ID L0→L11: −14.52 (−44.6% compression)
- Compression ratio: 0.554
- Ambient vs actual: scGPT gene manifold (32D) ≈ 25% of Gaussian baseline (123D)
- Monotone decrease: every layer reduces ID (no increase anywhere)

**Artifact**: `h02_fullmanifold_twonn.json`

**Interpretation**: POSITIVE / PROMISING.
1. scGPT gene embeddings live on a ~32D manifold in 512D space (well below Gaussian null ~123D)
2. Monotone compression across all 12 layers — every transformer block reduces intrinsic dimensionality
3. L0→L11 compression = 44.6%; compression ratio 0.554
4. Feature shuffle preserves ID (expected — shuffling columns doesn't change pairwise distances),
   so this null doesn't discriminate. Better null: Gaussian or cross-layer comparison.
5. Consistent with iter_0042/0043 B-cell specific findings; this is the full-manifold version.

**Caveat**: n_valid=846/2000 at each layer (some points identical ratio). Subsampling introduces
noise but trend is clear and monotone.

## H03: TCR Signaling Circuit Convergence (clean, no CD247 confound)

**Method**: T-cell identity centroid = CD247, CD8A, RUNX3 (removing CD28/LCK from centroid).
Circuit genes = CD28, LCK (pure upstream signaling, no overlap with centroid). Measure rank of
circuit genes to centroid at each layer. Null: 200 random 2-gene sets.

**Results**:
| Layer | CD28 rank | LCK rank | Mean rank |
|-------|-----------|----------|-----------|
| L0    | 595       | 6        | 300.5     |
| L1    | 511       | 7        | 259.0     |
| L2    | 471       | 7        | 239.0     |
| L3    | 297       | 7        | 152.0     |
| L4    | 188       | 7        | 97.5      |
| L5    | 129       | 6        | 67.5      |
| L6    | 86        | 6        | 46.0      |
| L7    | 108       | 5        | 56.5      |
| L8    | 107       | 5        | 56.0      |
| L9    | 148       | 5        | 76.5      |
| L10   | 154       | 6        | 80.0      |
| L11   | 149       | 8        | 78.5      |

**Key metrics**:
- Rank change L0→L11: −222.0
- Null mean: −15.6 ± 507.8
- p-value (one-sided, circuit converges more): 0.320

**Artifact**: `h03_tcr_circuit_clean.json`

**Interpretation**: INCONCLUSIVE.
- CD28: strong convergence from rank 595 to minimum 86 at L6 (V-shape)
- LCK: always rank 5-8 (already maximally close at L0 — constitutive co-embedding)
- Null p=0.32: not significant; LCK's constant-proximity dominates mean rank
- The CD28-specific trajectory (595→86→149) is suggestive but single-gene, underpowered
- To test properly: need more signaling circuit genes in vocab (ZAP70, PTPRC, ITK not in vocab)

## Retirements / Carry-forwards
- H01 (kNN lineage purity with few markers): RETIRED for this vocab. Revive only with expanded vocab
  or different community metric (modularity on larger sets).
- H02 (full-manifold TwoNN): CARRY FORWARD. Needs better null (cross-layer bootstrap, layer-shuffled
  controls) and extension to other datasets.
- H03 (TCR circuit): CARRY FORWARD with rescue — need expanded signaling circuit vocab or different
  centroid definition.
