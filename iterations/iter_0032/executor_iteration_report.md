# Executor Iteration Report: iter_0032

**Date:** 2026-02-23
**Focus:** Community structure biology, regulatory distance decay (OOV-corrected), PC1 biological polarity

---

## Summary

Three hypotheses tested across module_structure, manifold_distance, and intrinsic_dimensionality families.

**Key finding:** Dorothea high-confidence regulatory proximity is an early-layer signal that decays to null at L11 (Spearman rho=-0.853, p=0.0004) on the OOV-corrected 195-gene set. This contrasts sharply with iter_0024's full-gene result (stable AUROC ~0.67 across layers), suggesting OOV genes with zero embeddings artificially inflated late-layer regulatory signal previously.

---

## H01: Community Detection at L11 → Biological Annotation

### Method
- Greedy modularity community detection (networkx) on k=10 kNN graph at L11
- 195 in-vocab named genes
- Null: 100 random partitions of same community sizes
- Biological enrichment: Fisher's exact test for curated gene families and cell-type markers

### Command Trace
```
conda run -n subproject40-topology python iterations/iter_0032/run_iter0032_screen.py
```

### Results
- **Communities found:** 2 (sizes 107, 88)
- **Modularity:** 0.4342
- **Null modularity:** -0.001 ± 0.013
- **Z-score:** 34.12 (extremely strong above null)
- **Biological enrichment (Fisher p):**
  - B-cell markers (CD19, MS4A1, CD79A): 3/3 in community 1, OR=∞, p=0.090
  - IL2_pathway (IL2, IL2RA, IL2RG): 3/3 in community 0, OR=∞, p=0.163
  - T-cell markers: 4/5 in community 0, OR=3.38, p=0.251
  - AP1, HLA-I/II, BCL2fam, TNFSF, RUNX: p > 0.5 (no enrichment)

### Interpretation
The 2-community structure at L11 is geometrically very robust (z=34) but does not map cleanly to any single curated family at p<0.05. B-cell markers show marginal trend (p=0.09) toward community 1. The bifurcation likely captures a broad functional axis not fully represented by the curated families tested.

**Decision:** PROMISING (strong geometry, biological signal unclear)

---

## H02: Dorothea Regulatory Pairs → Distance Decay Across Layers (195 in-vocab)

### Method
- Load Dorothea regulon (dorothea_human.tsv), filter to 195 in-vocab genes
- High-conf (A/B): 205 pairs; Low-conf (C/D/E): 693 pairs
- Compare L2 embedding distance to 5000 random gene pairs per layer
- AUROC (regulatory pairs closer than random) and Spearman rho(layer, AUROC)

### Command Trace
```
conda run -n subproject40-topology python3 -c "... [inline script]"
```
(Full script inline in run_iter0032_screen.py extended section)

### Results
| Layer | High-conf AUROC | Low-conf AUROC |
|-------|----------------|----------------|
| L0    | 0.564          | 0.513          |
| L1    | 0.566          | 0.504          |
| L2    | 0.571          | 0.499          |
| L3    | 0.546          | 0.488          |
| L4    | 0.547          | 0.500          |
| L5    | 0.537          | 0.497          |
| L6    | 0.513          | 0.492          |
| L7    | 0.501          | 0.499          |
| L8    | 0.500          | 0.501          |
| L9    | 0.509          | 0.511          |
| L10   | 0.510          | 0.509          |
| L11   | 0.505          | 0.508          |

- **High-conf Spearman rho(layer, AUROC):** -0.853 (p=0.0004)
- **Low-conf Spearman rho:** +0.175 (p=0.587, null)
- **Mean AUROC high:** 0.531; **Mean AUROC low:** 0.500

### Interpretation
High-confidence regulatory pairs are closer in embedding space in **early layers** (AUROC ~0.57 at L0-L2), but this signal completely decays by L8-L11 (AUROC ~0.50). Low-confidence pairs show no signal at any layer.

**This is a new and important finding** for the OOV-corrected gene set: regulatory geometry exists in early layers but is destroyed by the manifold compression/fragmentation seen in late layers (PR collapse to 9.5 from 58 dimensions, spectral gap halving — iter_0031).

Contrast with iter_0024 (full 209-gene set, including OOV): AUROC ~0.67 across all layers. The discrepancy suggests OOV genes (zero embedding at L0, still non-zero at later layers) created an artificial late-layer baseline that inflated AUROC in prior analysis.

**Decision:** PROMISING (novel mechanistic finding, corrects prior result)

---

## H03: PC1 at L11 vs Biological Partitions

### Method
- SVD of 195-gene embeddings at L0, L5, L11 (mean-centered)
- Mann-Whitney AUROC (two-sided) for:
  - TFs (n=13) vs non-TF (n=182)
  - T-cell markers (n=5 in vocab)
  - B-cell markers (n=3 in vocab)
  - Monocyte markers (n=0 in vocab)

### Command Trace
```
conda run -n subproject40-topology python iterations/iter_0032/run_iter0032_screen.py
```

### Results
| Layer | PC1 var% | TF AUROC (p) | T-cell AUROC | B-cell AUROC (p) |
|-------|----------|--------------|--------------|-----------------|
| L0    | 8.0%     | 0.640 (0.093)| 0.496        | 0.219 (0.099)   |
| L5    | 10.9%    | 0.633 (0.111)| 0.583        | 0.161 **(0.041)**|
| L11   | 25.9%    | 0.655 (0.062)| 0.588        | 0.215 (0.095)   |

### Interpretation
B-cell markers (CD19, MS4A1, CD79A) are consistently at the **negative pole** of PC1 (AUROC < 0.5, significant at L5 p=0.041). TFs trend toward the positive pole (AUROC 0.64-0.66, p≈0.06 at L11). PC1 variance increases 3.2x from L0→L11.

However, the biological separation is limited by very small marker sets in the 195-gene vocab. The pattern is suggestive of a TF-vs-B-cell axis but underpowered for definitive claim.

**Decision:** INCONCLUSIVE (correct direction but limited power; needs larger gene sets)

---

## Artifacts Generated

| File | Description |
|------|-------------|
| `h01_community_detection_l11.json` | Communities, modularity, null, enrichment results |
| `h02_dorothea_dist_195invocab.json` | AUROC per layer for high/low Dorothea confidence pairs |
| `h03_pc1_biological_partitions.json` | PC1 variance and biological group AUROCs per layer |

---

## Key Takeaways

1. **Regulatory geometry is early-layer**: Dorothea high-conf pairs are significantly closer at L0-L2 (AUROC ~0.57) but not at L8-L11 (AUROC ~0.50). This OOV-corrected result revises the iter_0024 finding.
2. **L11 graph is strongly bimodal**: 2 communities with z=34 above null, but biology not yet identified.
3. **PC1 polarity**: B-cell markers at negative pole, TFs at positive pole — suggestive but underpowered.
