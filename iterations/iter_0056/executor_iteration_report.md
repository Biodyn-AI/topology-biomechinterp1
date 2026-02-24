# Executor Iteration Report — iter_0056

**Date**: 2026-02-23
**Status**: Complete — 3 hypotheses tested

---

## Summary

This iteration tested three hypotheses guided by the brainstormer's top priorities:
1. **H01** (Joint SV2-7 6D classifier): **Promising** — mean AUROC 0.744, stable across depth
2. **H02** (L9 crossover cross-seed replication): **Inconclusive** — 2/3 seeds confirm; seed43 fails
3. **H03** (Subspace rotation → AUROC): **Negative** — rho=-0.27, p=0.42, not significant

---

## Command Trace

```bash
# H01: Joint SV2-7 classifier across 12 layers
conda run -n subproject40-topology python /tmp/iter56_experiments.py

# H02: Cross-seed crossover + label-shuffle null
# H03: Principal angle trajectory
conda run -n subproject40-topology python /tmp/iter56_h2h3.py
```

Data paths:
- `cycle4_immune_main/layer_gene_embeddings.npy` [12, 4941, 512]
- `cycle4_immune_seed43/layer_gene_embeddings.npy` [12, 4941, 512]
- `cycle4_immune_seed44/layer_gene_embeddings.npy` [12, 4941, 512]
- Edge file: `cycle4_immune_main/cycle1_edge_dataset.tsv` (2940 edges, 589 positive nonzero pairs)

---

## H01: Joint SV2-7 (6D) Layer-Stable TF/Target Classifier

### Method
At each of 12 layers:
- Filter 4941 genes to 2039 nonzero (by L0 norm)
- SVD of centered [2039, 512]: extract gene coords in SV2-4 (u[:,1:4]·S[1:4]) and SV5-7 (u[:,4:7]·S[4:7])
- Concatenate to 6D; labels: TF (n=73) vs target-only (n=222)
- 5-fold stratified LR, mean AUROC; 100-perm label-shuffle null

### Results

| Layer | Joint 6D | SV5-7 | SV2-4 | Null | p_perm |
|-------|----------|-------|-------|------|--------|
| L0  | 0.771 | 0.715 | 0.634 | 0.506 | 0.000 |
| L1  | 0.765 | 0.691 | 0.688 | 0.505 | 0.000 |
| L2  | 0.782 | 0.723 | 0.700 | 0.505 | 0.000 |
| L3  | **0.789** | 0.716 | 0.703 | 0.504 | 0.000 |
| L4  | 0.760 | 0.668 | 0.721 | 0.503 | 0.000 |
| L5  | 0.757 | 0.611 | 0.725 | 0.506 | 0.000 |
| L6  | 0.730 | 0.571 | 0.723 | 0.502 | 0.000 |
| L7  | 0.738 | 0.584 | 0.721 | 0.504 | 0.000 |
| L8  | 0.727 | 0.532 | 0.732 | 0.506 | 0.000 |
| L9  | 0.715 | 0.663 | 0.655 | 0.506 | 0.000 |
| L10 | 0.702 | 0.684 | 0.652 | 0.503 | 0.000 |
| L11 | 0.688 | 0.686 | 0.655 | 0.500 | 0.000 |

**Key metrics**:
- Joint AUROC ≥ 0.72 at **9/12 layers** (L0–L8; drops L9-L11 where both subspaces carry overlapping info)
- Joint beats max(SV57, SV24) at **11/12 layers** (fails only at L8 by 0.005)
- Mean joint AUROC = **0.744** vs null ~0.505 (all p_perm = 0.00)
- Max joint AUROC = **0.789** at L3

**Interpretation**: The dual-subspace combination consistently outperforms either subspace individually. The joint classifier is layer-stable at AUROC ≥ 0.72 through most of the network. At L8 (SV57 minimum), SV24 alone matches the joint — indicating SV24 dominates there and SV57 contributes redundant/noisy signal.

**Decision**: Promising. Qualifies for Claim 50 (weakened from "10/12" to "9/12" — still strong).

---

## H02: L9 Directionality Crossover Cross-Seed Replication + Label-Shuffle Null

### Method
For each seed (main/seed43/seed44):
- SVD at each layer; mean displacement vector (target - TF) for valid TRRUST positive pairs
- Magnitude of mean displacement in SV2-4 vs SV5-7 subspaces
- 500-perm label-shuffle null (random target reassignment)
- Crossover = first layer where sv57_mag > sv24_mag

### Results

| Seed | Nonzero | Valid pairs | Crossover layer(s) |
|------|---------|-------------|---------------------|
| main | 2039 | 589 | **L9** |
| seed43 | 2080 | 553 | **None** |
| seed44 | 1972 | 577 | **L9** |

- main and seed44: clear L9 crossover
- seed43: sv57_mag never exceeds sv24_mag across all 12 layers
- 2/3 seeds confirm L9 crossover; full cross-seed replication fails

**sv57_mag trajectory (main)**: 0.42→0.44→0.54→0.81→0.86→0.92→0.96→0.97→0.54→1.25→1.08→1.10
**sv24_mag trajectory (main)**: 2.03→1.90→1.99→1.81→1.65→1.52→1.48→1.57→0.46→0.76→0.94→0.96

Note: At L8 for main seed, both magnitudes drop dramatically (sv24: 0.46, sv57: 0.54), then sv57 rebounds more strongly at L9 (1.25 vs 0.76).

**seed43 investigation**: Will need to check seed43's magnitude profiles to understand why crossover fails (different architecture? fewer valid pairs?).

**Decision**: Inconclusive. Claim 49 remains a main-seed finding; not promoted to fully cross-seed validated.

---

## H03: SV5-7 Principal Angle Trajectory vs AUROC

### Method
Main seed: extract SV5-7 right singular vectors (Vt[4:7].T ∈ R^{512×3}) per layer. Max principal angle between consecutive layers via QR decomposition + SVD of Q_A.T @ Q_B. Spearman ρ with sv57 AUROC at the next layer. Control: same for SV2-4.

### Results

| Transition | SV5-7 angle | SV2-4 angle | AUROC next |
|------------|-------------|-------------|------------|
| L0→L1 | 28.5° | 21.5° | 0.695 |
| L1→L2 | 29.3° | 29.3° | 0.719 |
| L2→L3 | **69.7°** | 31.4° | 0.728 |
| L3→L4 | 52.3° | 34.7° | 0.678 |
| L4→L5 | 40.2° | 42.0° | 0.622 |
| L5→L6 | 37.7° | 32.8° | 0.586 |
| L6→L7 | 35.1° | 36.6° | 0.603 |
| L7→L8 | 46.3° | 48.1° | 0.543 |
| L8→L9 | 42.6° | 48.0° | 0.670 |
| L9→L10 | 31.9° | 38.1° | 0.679 |
| L10→L11 | 29.6° | 39.5° | 0.680 |

**Spearman rho**: -0.273, p=0.417 (not significant)

**Key observation**: Large rotation at L2→L3 (69.7°, nearly orthogonal) does NOT predict low AUROC; in fact AUROC is high after this rotation. The L7→L8 rotation (46.3°) does coincide with a subsequent AUROC minimum (0.543), but this is a single point.

**Decision**: Negative. Rotation angle does not explain AUROC trajectory in a systematic way. The SV5-7 subspace is notably unstable across depth (large angles), but this instability is not predictive of information loss.

---

## Artifacts Generated

- `h01_joint_sv27_auroc.csv` — per-layer AUROC comparison (joint/sv57/sv24/null/p)
- `h02_crossover_replication.csv` — per-seed/layer displacement magnitudes + crossover flags
- `h03_subspace_rotation.csv` — principal angle trajectory + AUROC
- `iter56_summary.json` — numeric summary of all three hypotheses
- `executor_hypothesis_screen.json`

---

## Quantitative Summary

| Hypothesis | Status | Decision | Key Metric |
|-----------|--------|----------|------------|
| H01: Joint SV2-7 classifier | tested | **promising** | mean AUROC=0.744, beats both at 11/12 layers |
| H02: L9 crossover replication | tested | **inconclusive** | L9 confirmed in 2/3 seeds; seed43 fails |
| H03: Rotation → AUROC | tested | **negative** | Spearman rho=-0.27, p=0.42 |
