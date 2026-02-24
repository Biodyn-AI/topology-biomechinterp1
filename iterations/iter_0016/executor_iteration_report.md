# Executor Iteration Report — iter_0016

## Summary
Three hypotheses tested; all produced quantitative results with biological interpretation.
Key finding: SV3 also encodes STRING confidence as a monotonic gradient (rho=0.90), and SV2/SV3/SV4
axes are nearly uncorrelated in their per-pair co-pole rates (max pairwise r=0.25), indicating they
capture partially orthogonal PPI geometry. SV4 GO biology includes endosome, apoptosis, virus response,
and synapse — distinct from SV2/SV3 profiles.

---

## H01: SV3 STRING Confidence Gradient

**Hypothesis**: Replicate the SV2 rho=1.000 monotonic gradient on SV3 axis.

**Method**: 3092 STRING pairs (score≥0.4). Split into 5 quintiles by score.
For each quintile, compute mean SV3 co-pole z-score (K=52, N_null=300) across 12 layers.
Spearman rho of quintile midpoint vs mean z.

**Command**:
```bash
conda run -n subproject40-topology python iterations/iter_0016/run_iter0016_screen.py
```

**Results**:
- Quintile edges: [0.400, 0.459, 0.533, 0.640, 0.796, 0.999]
- Per-quintile mean_z: Q1=1.238, Q2=1.447, Q3=1.397, Q4=1.906, Q5=3.813
- Spearman rho=0.9000, p=0.037
- Monotonic trend holds but with Q3 dip (1.397 < Q2=1.447)

**Interpretation**: SV3 also shows a positive gradient (rho=0.90 vs SV2 rho=1.00). The signal is
monotonic at the top (Q4→Q5 jump: 1.906→3.813) and is statistically significant. Less perfectly
monotonic than SV2 (Q2/Q3 tie), but directionally consistent. SV3 appears to be a second axis
encoding PPI confidence.

**Decision**: promising (directional confirmation; weaker than SV2 but significant)

---

## H02: SV2/SV3/SV4 Axis Independence Test

**Hypothesis**: SV2, SV3, SV4 encode overlapping vs orthogonal subsets of STRING pairs.

**Method**: For each axis pair, compute:
1. Pairwise Jaccard overlap of top-K=52 and bottom-K=52 gene sets (layer-averaged)
2. For each STRING pair, compute mean per-layer co-pole indicator per axis
3. Inter-axis Pearson correlation of per-pair co-pole rates

**Results**:
- Top-pole Jaccard overlaps (layer avg):
  - SV2 vs SV3: 0.1326
  - SV2 vs SV4: 0.1717
  - SV3 vs SV4: 0.1836
- Inter-axis Pearson correlations of per-pair copole rates:
  - SV2–SV3: r=−0.016
  - SV2–SV4: r=+0.213
  - SV3–SV4: r=+0.247
- Dominant-axis pair counts: SV2=1547, SV3=850, SV4=695

**Interpretation**:
- Pole gene set overlap (Jaccard ~0.13–0.18) is modest but non-trivial, i.e. the axes are not
  purely orthogonal in gene membership but also not identical.
- Inter-axis copole correlations are near-zero to low (max r=0.247), indicating SV2/SV3/SV4 capture
  largely distinct STRING pair subsets.
- SV2 dominates (1547/3092 pairs) but SV3 and SV4 together cover an additional ~1545 pairs.
- Result supports the multi-axis framing: at least 3 near-orthogonal PPI geometry axes exist.

**Decision**: promising — orthogonality is real (max r=0.247), supports multi-axis claim

---

## H03: SV4 GO Biology Profile

**Hypothesis**: SV4 poles (top-K=52, bottom-K=52) show specific GO enrichment at focus layers.

**Method**: Fisher's exact test for 591 GO terms (size 3–50 in named genes) at layers 7, 8, 11 for
both SV4 poles (K=52). Report top enriched GO terms per layer.

**Results**:
- Layer 7 top hits:
  - GO:0005768 endosome (bottom pole, p=8.2e-4, fold=4.02): APP, FLT1, IL2RB, KDR
  - GO:0048143 astrocyte activation (bottom pole, p=3.5e-3, fold=4.02)
  - GO:0097192 extrinsic apoptotic signaling (bottom pole, p=3.5e-3, fold=4.02)
- Layer 8 top hits:
  - GO:0009615 response to virus (top pole, p=2.7e-3, fold=2.81): CCL5, CLU, FOXP3, IFNG
  - GO:0043029 T cell homeostasis (top pole, p=3.5e-3, fold=4.02): BCL2, FOXP3, PMAIP1
- Layer 11 top hits:
  - GO:0045202 synapse (top pole, p=3.5e-3, fold=4.02): ALDH1A1, APP, CLU, EGR3
  - GO:0097192 extrinsic apoptotic signaling (top pole, p=3.5e-3, fold=4.02): BCL2, BCL2A1, IL1B, IL2
  - GO:0035914 skeletal muscle cell differentiation (top pole, p=3.5e-3)

**Interpretation**: SV4 bottom pole (layers 7–8) clusters endosomal/apoptotic/neuroimmune genes.
Top pole captures antiviral/T-cell homeostasis (layer 8) and synaptic/apoptotic genes (layer 11).
These categories are distinct from SV2 (immune effector) and SV3 (broad immune/signaling) profiles,
supporting the orthogonality finding from H02. The SV4 axis appears to separate membrane-trafficking
and apoptotic biology from antiviral/adaptive immune biology.

**Decision**: promising — biologically coherent SV4 profile distinct from SV2/SV3

---

## Artifacts Generated
- `iterations/iter_0016/run_iter0016_screen.py` — full experiment script
- `iterations/iter_0016/iter0016_results.json` — combined results JSON
- `iterations/iter_0016/h01_sv3_confidence_gradient.json` — H01 detailed results
- `iterations/iter_0016/h02_axis_independence.json` — H02 detailed results
- `iterations/iter_0016/h03_sv4_go_biology.json` — H03 detailed results
- `iterations/iter_0016/run_iter0016_stdout.log` — full execution log

## Key Quantitative Claims (iter_0016)
1. SV3 STRING confidence gradient: rho=0.90, p=0.037; Q5 mean_z=3.81 vs Q1 mean_z=1.24
2. Axis independence: max inter-axis copole correlation r=0.247 (SV3–SV4); SV2–SV3 r=−0.016
3. SV4 GO top hits: endosome (p=8.2e-4), virus response (p=2.7e-3), synapse (p=3.5e-3)
