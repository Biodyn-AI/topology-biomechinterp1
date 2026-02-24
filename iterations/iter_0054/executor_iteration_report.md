# Executor Iteration Report — iter_0054

**Date**: 2026-02-23
**Focus**: TF/target classification from SV5-7 geometry; BFS hierarchy depth encoding; Layer-resolved directionality trajectory

---

## Summary

Three hypotheses tested, all yielding significant quantitative results:

- **H01** (LR classifier, SV5-7 L0): SV5-7 coordinates at L0 classify TF vs target-only genes with AUROC=0.694 (5-fold CV), vs null AUROC=0.494 (p=0.000; 0/100 permutations exceeded observed). L8 SV5-7 gives AUROC=0.499 (chance level), confirming early-layer specificity.
- **H02** (BFS depth vs SV5-7): BFS hierarchy depth (from master TFs, n=33) correlates weakly but significantly with SV5 (rho=0.167, p=0.0048) and SV7 (rho=0.169, p=0.0044) at L0. Effect is absent at L8, again confirming L0 specificity of the SV5-7 regulatory subspace.
- **H03** (Layer trajectory): Directionality signal (combined magnitude of mean TF→target displacement in SV5-7) increases monotonically from L0 (0.0057) through L9 (0.0200), then plateaus/decreases slightly at L10–L11. All 12 layers show at least one SV axis significantly displaced (p<0.05). Peak Cohen's d values reach 0.53 at L9 (SV5). This establishes a clear directional amplification trajectory.

---

## Command Trace

```bash
# Script: /tmp/iter54_main.py

# Run all experiments
conda run -n subproject40-topology python /tmp/iter54_main.py

# Data sources
# Embeddings: /Volumes/Crucial X6/.../cycle4_immune_main/layer_gene_embeddings.npy [12, 4941, 512]
# Edge dataset: .../cycle4_immune_main/cycle1_edge_dataset.tsv [2940 edges: 735 pos/2205 neg]
# Nonzero genes: 2039 (nonzero norm at all 12 layers)
# After nonzero filter: 589 pos / 1523 neg pairs
```

---

## H01: Logistic Regression TF vs Target-Only Classification (SV5-7)

**Hypothesis**: SV5-7 coordinates at L0 predict TF vs target-only gene identity (AUROC > 0.6).

**Setup**:
- TF genes: appear as source in ≥1 positive TRRUST edge (n=68 nonzero)
- Target-only genes: appear as target but never as source (n=215 nonzero)
- Features: 3D SV5-7 coordinates at L0
- Model: LogisticRegression (C=1.0, StandardScaler), 5-fold stratified CV
- Null: 100 label permutations, same CV

**Results**:

| Condition | AUROC (mean ± std) |
|-----------|-------------------|
| SV5-7, L0 | **0.694 ± 0.017** |
| SV2-4, L0 (control) | 0.664 ± 0.132 |
| SV5-7, L8 | 0.499 ± 0.075 |
| Null (permuted labels) | 0.494 ± 0.051 |

- p vs null: **0.000** (0/100 permutations ≥ 0.694)
- Fold scores: [0.713, 0.714, 0.688, 0.680, 0.674]

**Interpretation**: SV5-7 at L0 meaningfully encodes TF identity (AUROC 0.694 vs null 0.494). Effect is early-layer specific: L8 is at chance. SV2-4 gives similar point estimate but far higher variance, suggesting less stable encoding. **Decision: promising**

---

## H02: BFS Hierarchy Depth vs SV5-7 Dominant Axis

**Hypothesis**: Genes deeper in the TF→target cascade are displaced along the dominant SV5-7 directionality axis, proportional to BFS depth from master TFs.

**Setup**:
- Built directed graph from 589 positive TRRUST pairs (nonzero genes)
- Master TFs: in-degree=0 nodes (n=33)
- BFS from all master TFs simultaneously to assign hierarchy depth
- Spearman correlation: depth vs each SV5-7 coordinate at L0 and L8

**Results**:

| SV | Layer | rho | p-value |
|----|-------|-----|---------|
| SV5 | L0 | **0.167** | **0.0048** |
| SV6 | L0 | -0.057 | 0.337 |
| SV7 | L0 | **0.169** | **0.0044** |
| SV5 | L8 | -0.072 | 0.227 |
| SV6 | L8 | 0.079 | 0.184 |
| SV7 | L8 | -0.016 | 0.793 |

- n=283 nodes with assigned depth
- Depth distribution: {0: 33, 1: 57, 2: 157, 3: 36}

**Interpretation**: Modest but significant correlation between hierarchy depth and SV5/SV7 axes at L0 (rho~0.17, p<0.005). Absent at L8. Suggests SV5-7 at L0 partially encodes cascade position, not just TF/target identity. However, effect size is small. **Decision: neutral** (consistent with directed geometry but weak signal; not independently conclusive)

---

## H03: Layer-Resolved Directionality Trajectory (SV5-7, All 12 Layers)

**Hypothesis**: The TF→target directional asymmetry in SV5-7 space amplifies with layer depth, peaking mid-to-late network.

**Setup**:
- At each layer: SVD of centered nonzero embeddings [2039, 512], project to SV5-7
- For each of 589 positive pairs (nonzero both): compute displacement vector (target − source)
- Summary: combined magnitude of mean displacement, Cohen's d per SV axis

**Results**:

| Layer | Combined magnitude | SV5 d | SV6 d | SV7 d | Notable |
|-------|--------------------|--------|--------|--------|---------|
| L0 | 0.00568 | 0.120 | -0.139 | 0.106 | All 3 p<0.05 |
| L1 | 0.00677 | -0.085 | **0.224** | -0.087 | SV6 p<1e-7 |
| L2 | 0.00821 | 0.025 | **0.332** | 0.021 | SV6 p<1e-14 |
| L3 | 0.01110 | -0.195 | **-0.412** | -0.051 | SV6 p<1e-21 |
| L4 | 0.01231 | 0.238 | **0.391** | 0.187 | All p<1e-5 |
| L5 | 0.01313 | **-0.366** | -0.336 | 0.077 | 2/3 p<1e-14 |
| L6 | 0.01393 | 0.002 | -0.212 | **-0.399** | 2/3 p<1e-6 |
| L7 | 0.01418 | 0.100 | 0.152 | **-0.418** | SV7 p<1e-21 |
| L8 | 0.01717 | **-0.423** | -0.218 | **-0.390** | All 3 p<1e-7 |
| L9 | **0.01997** | **-0.529** | -0.309 | 0.403 | SV5 p<1e-33 |
| L10 | 0.01835 | -0.093 | **0.457** | **0.483** | SV6+7 p<1e-25 |
| L11 | 0.01708 | -0.020 | **-0.510** | -0.400 | SV6 p<1e-31 |

Key patterns:
1. **Monotonic amplification**: magnitude increases ~3.5× from L0→L9
2. **Peak at L9**: magnitude=0.0200, SV5 d=−0.529 (strongest single-axis d seen)
3. **All layers significant**: every layer has ≥1 SV axis with p<0.05
4. **Axis rotation**: dominant axis shifts (SV6 at L2-L4, SV5 at L5/L8-L9, SV7 at L6-L7/L10, SV6 at L11) — the space "rotates" while amplifying
5. **Plateau at L9-L11**: slight decrease after L9 peak

**Decision**: **promising** — this cleanly characterizes the trajectory of regulatory geometry through the network. The ~3.5× amplification from L0→L9 is a new, definitive result.

---

## Artifact Summary

| File | Description |
|------|-------------|
| `h01_lr_tf_target_sv57.json` | H01 AUROC results (all conditions) |
| `h02_bfs_depth_sv57_corr.csv` | H02 Spearman rho by SV and layer |
| `h03_layer_directionality_sv57.csv` | H03 per-layer directionality metrics |
| `iter54_summary.json` | Combined summary |
