# Executor Iteration Report: iter_0044

## Summary

Three hypotheses tested. All returned negative or inconclusive results, clarifying boundaries of prior positive findings and retiring several directions.

---

## Hypotheses Tested

### H01: Cross-lineage ID compression law
**Family**: intrinsic_dimensionality | **Novelty**: new_method | **Lineage**: iter_0043_H03

**Motivation**: iter_0043 found effector CD8 gene set shows ΔID = -45.8 (compression) at L7 while pan-T expands (+26.2). Hypothesis: this is a general law — effector cell types compress while reference types expand.

**Method**: TwoNN intrinsic dimensionality on gene subsets (5-point clouds) at each layer. Pre/post L7 split.

**Command**:
```bash
conda run -n subproject40-topology python run_iter0044_screen.py
```

**Results**:
| Gene set | In-vocab | Pre-L7 ID | Post-L7 ID | ΔID |
|---|---|---|---|---|
| NK effector | FCGR3A/KLRB1/KLRD1/GNLY (n=4) | 14.7 | 14.0 | -0.7 |
| Plasma cell | IGHG1/IGKC/JCHAIN/MZB1/DERL3 (n=5) | 9.6 | 12.4 | +2.7 |
| Neutrophil | S100A8/S100A9/CXCR1/FCGR3B/MPO (n=5) | 2.2 | 1.8 | -0.4 |
| Pan-myeloid | LYZ/CST3/AIF1/FCN1/CD14 (n=5) | 17.4 | 17.9 | +0.5 |
| Effector CD8 | CD8A/GZMB/PRF1/NKG7/GNLY (n=5) | 8.6 | 11.4 | +2.8 |
| General T | CD3E/CD3D/TRAC/TRBC1 (n=4) | 8.0 | 6.2 | -1.8 |

**Critical observation**: These ΔID values are ~100x smaller than iter_0043's results (which showed ΔID = -45.8 for effector CD8). The iter_0043 computation likely used TwoNN over the full ~4941-gene manifold and extracted local dimensionality estimates; the current 5-point TwoNN is measuring pairwise distance ratios among only 5 genes.

**Decision**: **NEGATIVE** — No cross-lineage compression pattern found. Effector CD8 shows +2.8 (expansion not compression) with this approach. Methodological inconsistency with iter_0043 confirms these are measuring different quantities. Direction retired.

**Artifact**: `h01_crosslineage_id.json`

---

### H02: TCR activation circuit attractor test
**Family**: manifold_distance | **Novelty**: new_method | **Lineage**: none

**Motivation**: The B-cell GC circuit (BATF/BACH2 converging toward B-cell centroid) is a positive result. Hypothesis: TCR signaling circuit genes (CD28/LAT/LCK/ZAP70/CD247) similarly converge toward T-cell centroid at deep layers.

**Method**: Precision@20 and mean rank of TCR circuit genes relative to T-cell centroid (CD3E/CD3D/TRAC/TRBC1/CD247) at each layer. Null: 200 random gene sets.

**In-vocab**: CD28, LCK, CD247 (3/5 found; LAT, ZAP70 absent)

**Results by layer**:
| Layer | Mean rank | Prec@20 |
|---|---|---|
| L0 | 92 | 0.67 |
| L5 | 44 | 0.67 |
| L6 | 37 | 0.67 |
| L11 | 41 | 0.67 |

- Rank change L0→L11: **-50.7** (ranks improve from 92 to 41)
- Precision@20: constant at **0.67** throughout all layers
- Null mean rank change: +0.6 ± 211
- **Null p-value = 0.275** (not significant)

**Confound**: CD247 is in both the circuit gene set AND the centroid gene set, meaning it always ranks 1st relative to the centroid (it IS the centroid). This inflates precision@20 artificially to 0.67 even at L0.

**Decision**: **INCONCLUSIVE** — The rank improvement (92→41) is directionally interesting but fails null significance. The overlap between circuit and centroid genes confounds interpretation. Would need clean separation of circuit vs centroid genes.

**Artifact**: `h02_tcr_circuit_attractor.json`

---

### H03: BCL6 rank B-cell vs metabolic co-movement
**Family**: manifold_distance | **Novelty**: refinement | **Lineage**: iter_0042_H02

**Motivation**: BCL6 was previously found to be metabolically isolated (near a metabolic gene cluster). Test if BCL6's distance to B-cell centroid and distance to metabolic centroid are anti-correlated across layers (divergence) or co-correlated (global movement).

**Method**: Track BCL6 rank-to-B-cell-centroid (MS4A1/CD79A/BLK/PAX5) and rank-to-metabolic-centroid (NAMPT/GLUL/PFKFB3/ACSL1/NIBAN1/FNDC3B/VMP1) at L0-L11. Spearman correlation.

**Results**:
- BCL6 B-cell ranks: [655, 773, 1312, 1078, 1165, 1097, 930, 1011, 791, 751, 770, 928]
- BCL6 metabolic ranks: [133, 161, 282, 209, 227, 200, 204, 238, 222, 220, 209, 202]
- **Spearman rho = 0.504, p = 0.095**
- Interpretation: **co-move** (both tracks rise and fall together)

**Decision**: **NEGATIVE** — BCL6 does not selectively diverge toward the metabolic cluster while leaving the B-cell cluster. Both distance tracks co-vary, indicating global positional shifts rather than a targeted metabolic attractor. The BCL6 metabolic isolation from iter_0042 may reflect stable proximity rather than directed movement. Retired.

**Artifact**: `h03_bcl6_rank_metabolic.json`

---

## Methodological Note

The discrepancy between H01 results and iter_0043 H03 results highlights an important issue: TwoNN on N=5 isolated gene embeddings measures the geometry of a 5-point simplex, not local intrinsic dimensionality within the embedding manifold. Future ID estimates should use TwoNN on the full gene manifold with local subsets as query points.

---

## Artifacts Generated
- `h01_crosslineage_id.json`
- `h02_tcr_circuit_attractor.json`
- `h03_bcl6_rank_metabolic.json`
- `executor_hypothesis_screen.json`
- `executor_iteration_report.md`
- `executor_next_steps.md`

## Command Trace
```bash
cd /Volumes/Crucial\ X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0044
conda run -n subproject40-topology python run_iter0044_screen.py
# Output: 3 JSON artifacts written, all experiments completed successfully
# Runtime: ~60s
```
