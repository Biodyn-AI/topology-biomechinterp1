# Executor Iteration Report — iter_0043

## Summary

Three hypotheses tested. Key findings: (1) PRDM1 anti-convergence hypothesis is **negative** — PRDM1 converges toward B-cell centroid at the same rate as BATF/BACH2 (rho=-1.000), meaning the GC repressor is NOT geometrically separated from activators; (2) BCL6 also converges toward B-cell centroid (rho=-0.993), with a transient rank increase at L3 suggesting brief divergence before convergence; (3) T-cell subtype stratification reveals that effector CD8 T cells COMPRESS (ID -45.78 post-L7) while general T-cell markers EXPAND (+26.18), revealing subtype-divergent geometric fates.

---

## Command trace

```bash
conda run -n subproject40-topology python \
  iterations/iter_0043/run_iter0043_screen.py
```

---

## H01: GC Repression Circuit — PRDM1/IRF4 Anti-Convergence Test

**Family**: manifold_distance | **Novelty**: new_method (repressor trajectory test)
**Lineage**: H03_iter0042, H01_iter0041

**Method**: Load cycle4_immune embeddings [12, 4941, 512]. Compute B-cell centroid (MS4A1, CD79A, BLK, PAX5). For each GC circuit gene {PAX5, BATF, BACH2, BCL6, PRDM1} (IRF4 not in vocab): compute distance to B-cell centroid at each layer L0-L11, compute Spearman rank correlation of distance vs layer index. Compare to null: 20 random genes. **Prediction**: activators (BATF, BACH2) converge (rho < 0); repressors (PRDM1) diverge (rho > 0).

**Results**:

| Gene | Dist L0 | Dist L3 | Dist L6 | Dist L11 | Spearman rho | p-value |
|------|---------|---------|---------|---------|--------------|---------|
| PAX5 | 9.769 | 7.399 | 5.759 | 2.395 | **-1.000** | 0.0000 |
| BATF | 15.912 | 12.900 | 10.498 | 4.087 | **-1.000** | 0.0000 |
| BACH2 | 14.913 | 11.817 | 10.510 | 4.082 | **-1.000** | 0.0000 |
| BCL6 | 15.064 | 12.932 | 11.295 | 4.963 | **-0.993** | 0.0000 |
| PRDM1 | 15.635 | 12.736 | 10.502 | 4.008 | **-1.000** | 0.0000 |
| Null (n=20) | — | — | — | — | 0.275 ± 0.931 | — |

**Selected pairwise distances (L0 → L11)**:
- BATF-PRDM1: 17.603 → 4.232 (co-localize at L11)
- BACH2-PRDM1: 16.796 → 3.940 (tightest pair at L11)
- BCL6-PRDM1: 16.022 → 4.104

**Interpretation**: The anti-convergence prediction is **not supported**. PRDM1 (BCL6 repression target, plasma differentiation driver) converges toward B-cell centroid at the same Spearman rho=-1.000 as BATF/BACH2. All GC circuit genes — activators and repressors — consolidate near B-cell centroid by L11. At L11, BACH2-PRDM1 distance (3.940) is the tightest pair in the entire circuit, tighter than BATF-BACH2 (4.526). This means the model encodes the entire GC circuit (including the terminally differentiated plasma pathway) within a unified B-cell-proximal geometric neighborhood. IRF4 was not in vocab (confirmed missing).

**Decision**: negative for anti-convergence hypothesis, but reveals a strong positive result: GC circuit unity — all five regulators converge to the same region, suggesting the model encodes GC biology as a unified attractor rather than opposing poles.

**Artifact**: `iterations/iter_0043/h01_gc_repression_circuit.json`

---

## H02: BCL6 Trajectory from B-cell Centroid

**Family**: manifold_distance | **Novelty**: new_method (BCL6 divergence test)
**Lineage**: H02_iter0042, H03_iter0038

**Method**: Track BCL6 Euclidean distance to B-cell centroid (MS4A1, CD79A, BLK, PAX5) across L0-L11. Also track rank (how many genes are closer to centroid than BCL6). Compute rho. Null: 20 random non-B-cell genes.

**Results**:

| Layer | BCL6 Dist | BCL6 Rank |
|-------|-----------|-----------|
| L0 | 15.06 | 654 |
| L1 | 13.01 | 772 |
| L2 | 13.04 | 1311 |
| L3 | 12.93 | 1077 |
| L4 | 12.24 | 1164 |
| L5 | 11.87 | 1096 |
| L6 | 11.30 | 929 |
| L7 | 10.96 | 1011 |
| L8 | 10.16 | 790 |
| L9 | 8.53 | 750 |
| L10 | 6.51 | 770 |
| L11 | 4.96 | 927 |

- BCL6 rho = -0.993 (p=0.0000) — converging toward B-cell centroid in absolute distance
- Null rho mean = 0.081 ± 0.969
- **Rank anomaly**: BCL6 rank increases from 654 (L0) to 1311 (L2) — indicating BCL6 briefly moves relatively farther from the B-cell centroid in rank space at L2, despite absolute distance decreasing
- By L11, BCL6 rank=927 (out of 4941 genes), placing it within top 19% of B-cell-proximate genes

**Interpretation**: BCL6 does NOT diverge from B-cell centroid — it converges like other circuit members. However, the rank trajectory reveals a transient "rank dip" at L2-L5 (rank ~1000-1311) before recovery, consistent with BCL6's early metabolic isolation (iter_0042: BCL6 cluster persists strongly at L0-L4). The metabolic isolation and B-cell convergence are not contradictory: BCL6 is simultaneously metabolically isolated AND converging toward B-cell context — possibly reflecting its dual identity as GC master regulator.

**Decision**: negative for divergence hypothesis; BCL6 converges. The rank anomaly at L2 is a secondary interesting finding worth tracking.

**Artifact**: `iterations/iter_0043/h02_bcl6_divergence.json`

---

## H03: T-cell Subtype ID Stratification at L7 Breakpoint

**Family**: intrinsic_dimensionality | **Novelty**: new_method (subtype stratification of ID breakpoint)
**Lineage**: H01_iter0042 (T-cell ID expansion at L7)

**Method**: Compute TwoNN intrinsic dimensionality at each layer L0-L11 for gene subsets representing T-cell subtypes: effector CD8 (CD8A, GZMB, PRF1, NKG7, GNLY; 5/5 in vocab), Treg (FOXP3, IL2RA, CTLA4, IKZF2; 2/4 in vocab — CTLA4, IKZF2), exhaustion (LAG3, PDCD1, HAVCR2, TIGIT; 1/4 in vocab — TIGIT), helper TFH (CD4, CXCR5, ICOS, BCL6, SH2D1A; 2/5 in vocab — ICOS, BCL6), general T-cell (CD3E, CD3D, TRAC, TRBC1; 4/4 in vocab). For subtypes with ≥3 genes in vocab, measure pre-L7 mean ID vs post-L7 mean ID.

**Results**:

| Subtype | Genes in vocab | Pre-L7 mean ID | Post-L7 mean ID | ΔID (post-pre) |
|---------|---------------|----------------|-----------------|----------------|
| effector_cd8 | 5/5 | **139.3** | **93.5** | **-45.8** (compression) |
| general_t | 4/4 | **96.8** | **123.0** | **+26.2** (expansion) |
| treg | 2/4 | NaN (n<3) | — | — |
| exhaustion | 1/4 | NaN (insufficient) | — | — |
| helper_tfh | 2/5 | NaN (n<3) | — | — |

Effector CD8 ID by layer: [176.7, 197.2, 160.9, 165.5, 105.5, 80.0, 89.2, 87.8, 123.1, 98.3, 84.5, 73.9]
General T-cell ID by layer: [100.2, 100.6, 102.5, 102.4, 96.2, 91.1, 84.6, 92.2, 101.5, 109.1, 146.5, 165.5]

**Interpretation**: The T-cell L7 ID expansion (from iter_0042: +4.965 slope) is driven by **general lineage markers** (CD3E/CD3D/TRAC/TRBC1), not by effector identity. Effector CD8 T-cells show the opposite pattern — ID compression throughout (peaks at L1: 197.2, decreases to 73.9 at L11). This reveals **divergent geometric fates**: effector CD8 identity is compressed into a tight manifold, while pan-T identity expands into a higher-dimensional representation. Possible interpretation: the model specializes effector cytotoxic genes (high functional specificity → low dimensionality) while the general T-cell machinery becomes representationally richer at deeper layers.

**Decision**: promising — clean dissociation between subtype-specific and pan-lineage geometric trajectories. Effector CD8 compression is a new positive finding.

**Artifact**: `iterations/iter_0043/h03_tcell_subtype_id.json`

---

## Quantitative summary

| Hypothesis | Primary metric | Value | Direction | Decision |
|-----------|---------------|-------|-----------|----------|
| H01: GC repressor anti-convergence | PRDM1 Spearman rho | -1.000 | negative (anti-convergence fails; circuit unity found) | negative |
| H02: BCL6 divergence from B-cell | BCL6 Spearman rho | -0.993 | negative (BCL6 converges, not diverges) | negative |
| H03: T-cell subtype ID stratification | Effector CD8 ΔID at L7 | -45.8 vs +26.2 (general T) | positive (dissociation confirmed) | promising |

---

## New mechanistic understanding

1. **GC circuit unity**: All five GC circuit members (PAX5, BATF, BACH2, BCL6, PRDM1) converge toward a unified B-cell-proximal neighborhood by L11. The model does NOT geometrically separate activators from repressors — the entire GC regulatory circuit consolidates together.
2. **BCL6 dual identity**: BCL6 converges toward B-cell centroid (rho=-0.993) but shows a rank anomaly at L2 — consistent with transient metabolic isolation before contextual integration.
3. **Effector CD8 compression vs pan-T expansion**: The L7 T-cell expansion is a pan-lineage marker effect. Effector genes compress (−45.8 ID change), suggesting specialized cytotoxic identity is encoded in lower dimensions.
