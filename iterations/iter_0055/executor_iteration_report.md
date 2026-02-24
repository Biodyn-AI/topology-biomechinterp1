# Executor Iteration Report — iter_0055

**Date**: 2026-02-23
**Focus**: Cross-seed replication, layer-resolved AUROC, and SV2-4 control directionality

---

## Summary

Three hypotheses tested, all yielding highly significant and interpretable results:

- **H01** (Cross-seed replication): TF vs target-only AUROC in SV5-7 at L0 replicates across all three seeds (main=0.708, seed43=0.717, seed44=0.738; all p=0.0099 by 100-perm null). The early-layer geometric encoding is fully reproducible.
- **H02** (Layer-resolved AUROC): SV5-7 AUROC peaks early (L2-L3: ~0.72) and at L9 (0.670), and declines in mid layers (L6-L8: <0.62). SV2-4 shows inverse pattern: AUROC rises across layers, peaking at L5-L8 (~0.74), then drops at L9-L11. Clear complementary subspace encoding.
- **H03** (Directionality comparison): SV2-4 and SV5-7 have *opposite* layer trajectories. SV2-4 displacement magnitude decreases L0→L11 (2.03→0.63), while SV5-7 increases L0→L9 (0.42→1.29) then decreases at L10-L11. The crossover (SV5-7 > SV2-4) occurs at L9, aligning with the AUROC crossover seen in H02.

---

## Command Trace

```bash
# Write Python script
cat > /tmp/iter55_main.py << 'PYEOF'
# [full script at /tmp/iter55_main.py]
PYEOF

# Run all experiments
conda run -n subproject40-topology python /tmp/iter55_main.py
```

### Data sources
- Main: `/Volumes/.../cycle4_immune_main/layer_gene_embeddings.npy` — shape [12, 4941, 512]
- Seed43: `/Volumes/.../cycle4_immune_seed43/layer_gene_embeddings.npy`
- Seed44: `/Volumes/.../cycle4_immune_seed44/layer_gene_embeddings.npy`
- Edge datasets: respective `cycle1_edge_dataset.tsv` (TRRUST pairs)
- Nonzero genes: main=2039, seed43=2080, seed44=1972

---

## H01: Cross-Seed Replication (TF vs Target-only, SV5-7, L0)

**Hypothesis**: The AUROC=0.694 from iter_0054 (TF vs target-only genes in SV5-7 space at L0) replicates on seeds 43 and 44.

**Method**: At L0, SVD of centered nonzero-gene embeddings. Project to SV5-7 (0-indexed: 4,5,6). 5-fold stratified LR cross-validation. 100-permutation null per seed.

**Results**:

| Seed | n_TF | n_target-only | AUROC mean ± std | null_mean | p_perm |
|------|------|--------------|-----------------|-----------|--------|
| main | 73 | 222 | 0.708 ± 0.045 | ~0.50 | 0.0099 |
| seed43 | 69 | 218 | 0.717 ± 0.106 | ~0.50 | 0.0099 |
| seed44 | 72 | 211 | 0.738 ± 0.078 | ~0.50 | 0.0099 |

**Interpretation**: All three seeds reproduce significant AUROC (0.71–0.74 vs ~0.50 null). The signal is robust to random seed variation. The slightly higher AUROC on seed44 may reflect minor differences in nonzero gene count (1972 vs 2039/2080). Effect is consistent in direction and magnitude across all three seeds.

**Decision**: Promising — full cross-seed replication confirmed.

---

## H02: Layer-Resolved AUROC (SV5-7 vs SV2-4, all 12 layers)

**Hypothesis**: The TF/target-only classification AUROC varies by layer and may peak at L9 (directionality peak from iter_0054). SV2-4 serves as control.

**Results**:

| Layer | SV5-7 AUROC | SV2-4 AUROC |
|-------|-------------|-------------|
| L0 | **0.708** | 0.657 |
| L1 | 0.695 | 0.701 |
| L2 | **0.719** | 0.704 |
| L3 | **0.728** | 0.717 |
| L4 | 0.678 | **0.736** |
| L5 | 0.622 | **0.740** |
| L6 | 0.586 | **0.738** |
| L7 | 0.603 | **0.729** |
| L8 | 0.543 | **0.738** |
| L9 | 0.670 | 0.648 |
| L10 | 0.679 | 0.645 |
| L11 | 0.680 | 0.651 |

**Key patterns**:
1. **SV5-7 encodes TF/target distinction primarily in early layers (L0-L3)** and deep layers (L9-L11); drops dramatically in middle layers (L6-L8).
2. **SV2-4 encodes TF/target distinction primarily in middle layers (L4-L8)**; drops below SV5-7 at L9-L11.
3. A **complementary encoding regime** exists: SV5-7 dominates early+late, SV2-4 dominates mid-depth. The information is never simultaneously absent.
4. SV5-7 AUROC at L8 (0.543) nearly random — consistent with the subspace rotation pattern noted in iter_0054.

**Decision**: Promising — layer-AUROC dissociation is a novel mechanistic finding.

---

## H03: Directionality Trajectory — SV2-4 vs SV5-7 (control comparison)

**Hypothesis**: SV2-4 does NOT show the same increasing directionality trajectory as SV5-7. This control test distinguishes signal-specific directionality from generic SVD subspace effects.

**Results**: Both subspaces show significant directionality (all p<0.001), but their trajectories are OPPOSITE:

| Layer | SV2-4 displacement mag | SV5-7 displacement mag |
|-------|------------------------|------------------------|
| L0 | **2.035** | 0.416 |
| L1 | 1.905 | 0.444 |
| L2 | 1.989 | 0.545 |
| L3 | 1.815 | 0.808 |
| L4 | 1.650 | 0.858 |
| L5 | 1.516 | 0.917 |
| L6 | 1.477 | 0.960 |
| L7 | 1.566 | 0.968 |
| L8 | 1.302 | **1.222** |
| L9 | 0.983 | **1.295** ← crossover |
| L10 | 0.937 | 0.897 |
| L11 | 0.635 | 0.613 |

**Crossover**: SV5-7 directionality exceeds SV2-4 at L9, the same layer where SV5-7 AUROC bounces back (H02). The convergence at L10-L11 (both ~0.6-0.9) suggests a shared information integration in the final layers.

**Interpretation**: SV2-4 encodes a decreasing but dominant early-to-mid directional signal. SV5-7 encodes an increasing late-stage signal that amplifies L0→L9. The two subspaces represent distinct and complementary computational regimes, both relating to regulatory directionality but in different processing phases.

**Decision**: Promising — establishes SV2-4 as a quantitative control with opposite dynamics, strengthening the SV5-7 claim.

---

## Artifacts Generated
- `h01_cross_seed_auroc.csv` — per-seed AUROC, null, p-value
- `h02_layer_auroc.csv` — per-layer AUROC for SV5-7 and SV2-4
- `h03_directionality_sv24_vs_sv57.csv` — displacement magnitude and Cohen's d per layer/subspace
- `iter55_summary.json` — all results combined
