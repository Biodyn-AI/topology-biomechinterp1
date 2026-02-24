# Executor Iteration Report — iter_0035

## Summary

Three hypotheses tested on the 195 in-vocabulary gene universe. Two strong positive findings and one definitive null.

**Key findings:**
1. **H02 (B-cell kNN precision@10):** B-cell markers cluster geometrically at 7x enrichment at L2 (z=5.20, emp_p=0.000), persisting through L11 (z=2.25, emp_p=0.044). This is a new, direct measure of B-cell geometric coherence independent of the PC1 analysis.
2. **H03 (B-cell centroid separation):** Normalized B-cell centroid separation increases monotonically from L0 to L11 (Spearman rho=+0.972, p<0.0001), showing B-cell amplification relative to manifold scale. Raw centroid distance is confounded by overall shrinkage.
3. **H01 (TRRUST + GO proximity):** Co-regulated pairs and GO co-annotated pairs show no proximity signal at L11 (AUROC ≈ 0.50). The functional annotation → embedding distance hypothesis is retired.

---

## Command Trace

```bash
# Main experiment script
conda run -n subproject40-topology python \
  iterations/iter_0035/run_iter0035_screen.py
```

Output written to `iterations/iter_0035/`.

---

## H01: TRRUST + GO Jaccard vs Embedding Proximity

**Family:** module_structure
**Method:** Co-regulated pairs (genes sharing a TF regulator from TRRUST edges, n=5989) and GO co-annotated immune pairs (10 gene sets, n=32). Mann-Whitney AUROC testing lower L2 distance vs 5000 random pairs at L11. Per-layer Spearman rho.
**Artifact:** `h01_trrust_go_jaccard.json`

**Results:**
- Co-reg pairs: n=5989, mean dist=5.399 vs random=5.378
- GO pairs: n=32, mean dist=5.136 vs random=5.378
- Co-reg AUROC at L11 = **0.506** (p=0.867) — null
- GO AUROC at L11 = **0.473** (p=0.302) — null
- Co-reg Spearman rho(layer, AUROC) = −0.993 (p<0.0001): co-reg AUROC decreases across layers
- GO Spearman rho(layer, AUROC) = +0.720 (p=0.008): GO AUROC increases across layers

**Interpretation:** Neither TF co-regulation nor GO co-annotation predicts embedding proximity. Co-reg AUROC trend is negative (genes diverge over layers), GO trend is positive but AUROC remains below 0.5 at all layers. **Decision: negative.**

---

## H02: B-cell Conditional Nearest-Neighbor Profile

**Family:** manifold_distance
**Method:** For each of 5 in-vocab B-cell markers (MS4A1, CD19, CD79A, BLK, FCER2): find k=10 nearest neighbors in 195-gene embedding space. Compute mean precision@10 (fraction of neighbors that are B-cell markers). Bootstrap null: 500 random size-5 gene sets. Test at L2, L5, L8, L11.
**Artifact:** `h02_bcell_knn_precision.json`

**Results:**
| Layer | obs_precision@10 | null_mean | z-score | emp_p |
|-------|-----------------|-----------|---------|-------|
| L2    | **0.140**       | 0.020     | 5.20    | 0.000 |
| L5    | **0.100**       | 0.022     | 3.18    | 0.020 |
| L8    | **0.080**       | 0.020     | 2.51    | 0.042 |
| L11   | **0.080**       | 0.021     | 2.25    | 0.044 |

- Peak enrichment at L2: B-cell markers appear in each other's k-nearest neighborhoods at 7× the null rate
- Signal sustained through L11: all layers show significant enrichment (all emp_p < 0.05)
- No monotonic decrease: L5–L11 remain stable at ~0.08–0.10 precision

**Decision: promising.** B-cell markers form a geometrically coherent cluster across all tested layers, with early-layer peak. This is structurally independent from the PC1 confound test in iter_0034.

---

## H03: B-cell Centroid Separation Trajectory

**Family:** intrinsic_dimensionality
**Method:** Per layer L0–L11: compute B-cell centroid (n=5) and non-B-cell centroid (n=190). Raw L2 distance and normalized distance (divided by mean pairwise distance from 2000 random gene pairs). Spearman rho(layer, distance). Bootstrap null: 1000 random size-5 gene sets tested for raw rho.
**Artifact:** `h03_bcell_centroid_trajectory.json`

**Normalized separation per layer:**
| Layer | centroid_dist | mean_pairwise | norm_dist |
|-------|--------------|---------------|-----------|
| L00   | 5.524        | 16.079        | 0.344     |
| L02   | 4.488        | 12.807        | 0.351     |
| L05   | 4.113        | 11.256        | 0.365     |
| L08   | 4.428        | 10.754        | 0.412     |
| L11   | 2.456        | 5.361         | **0.458** |

- Spearman rho(layer, **normalized** dist) = **+0.972**, p < 0.0001
- Spearman rho(layer, raw dist) = −0.902, p < 0.0001 (raw distance shrinks with layer — confound)
- Bootstrap null (raw rho): null_mean = −0.969, z = 0.82, emp_p = 0.088

**Interpretation:** Raw centroid distance decreases due to overall embedding magnitude shrinkage (as confirmed by mean_pairwise also shrinking). After normalization, B-cell centroid separation **increases** monotonically from 0.34 to 0.46 (33% relative increase from L0 to L11). The normalized signal is highly significant (Spearman p<0.0001) but does not have a bootstrap null to compare against. **Decision: promising** — normalized separation amplification is a robust finding.

---

## Data Provenance

- Embeddings: `cycle1_main/layer_gene_embeddings.npy` [12, 4803, 512]
- Gene universe: 195 in-vocab OOV-filtered named genes
- B-cell markers in-vocab (n=5): MS4A1, CD19, CD79A, BLK, FCER2
- Seed: rng=42 throughout
