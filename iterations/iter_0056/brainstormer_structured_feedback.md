# Brainstormer Structured Feedback — iter_0056

**Date**: 2026-02-23
**Gate**: PASSED (`passed_min_research_gate: true`)

---

## Outcome Assessment

### H01: Joint SV2-7 6D Classifier — Promoted (Claim 50)
**Strong.** Mean AUROC=0.744, beats both individual subspaces at 11/12 layers, all p_perm=0.000 at every layer. The dual-subspace combination unifies Claims 42/48 into a single predictive object. Correctly promoted.

One concern: "layer-stable" is strong language when AUROC drops to 0.688 at L11. The claim says "most layers (9/12 ≥ 0.72)" which is accurate. No issue with the decision.

Cross-seed validation is the essential follow-up. We have seeds 43/44; the infrastructure is already in place.

### H02: L9 Crossover Cross-Seed — Inconclusive
**Biologically interesting failure.** Seed43 shows sv57_mag never exceeding sv24_mag across all 12 layers — this is a qualitative, not marginal, difference. Two possibilities:
1. Seed43's SVD basis is rotated such that the regulatory signal is distributed differently (different SV assignment to indices 2-4 vs 5-7).
2. Seed43 represents a genuinely different attractor of the same model (different local minimum in training) with different spectral structure.

Option 1 is more likely and testable: check cross-seed alignment of right singular vectors. If seed43's SV5-7 right singular vectors are closer to main's SV2-4 than its SV2-4, the labeling is swapped. This would rescue Claim 49.

**Do not retire Claim 49.** It holds in 2/3 seeds and the L9 transition is scientifically coherent. Condition it.

### H03: Subspace Rotation Angle → AUROC — Retire
**Clean negative.** rho=-0.27, p=0.42. The L2→L3 spike (69.7°, near-orthogonal rotation in SV5-7) is notable but does not linearly predict AUROC — in fact AUROC is high after this transition. This suggests near-orthogonal rotation is compatible with information preservation: the subspace re-orients without losing discriminative content. Interesting mechanistically but not a predictive signal.

**Key residual insight**: The L7→L8 AUROC minimum (SV5-7 AUROC drops to 0.543) coincides with a rotation of 46.3° AND the collapse of both displacement magnitudes in H02 (sv24: 2.03→0.46, sv57: 0.96→0.54). L8 may be a genuine computational bottleneck worth investigating via mechanism (attention pattern analysis) rather than geometry alone.

---

## Pattern Recognition Across Recent Iterations

The research has converged on a well-characterized 2-subspace picture (Claims 42-50):
- SV5-7: early regulatory geometry (L0-L3), amplifying directional asymmetry peaking L9
- SV2-4: late co-expression-confounded regulatory geometry (L4-L8)
- Joint 6D: layer-stable separator (Claim 50)
- L9: putative integration point

**What's missing**: mechanism, generalization, and topological depth.
- No cross-model test since iter_0003 (and that was inconclusive).
- No persistent homology test in the SV2-7 subspace specifically.
- No test of whether the seed43 anomaly is basis-alignment vs training-diversity.
- No test of whether the 6D classifier generalizes OOD (other cell types, other datasets).
- No spectral test linking graph Laplacian structure to embedding geometry.

---

## Direction Fitness Assessment

| Direction | Status | Assessment |
|-----------|--------|------------|
| Subspace rotation → AUROC | **retire_now** | Clean negative this iter |
| STRING PPI → proximity | retired (iter_0031) | — |
| Annotation → embedding distance | retired (iter_0032) | — |
| TDA loops on circuit genes | retired (iter_0052) | — |
| BCL6 divergence | retired (iter_0041) | — |
| Seed43 crossover investigation | **active** | One cheap iteration can resolve |
| Joint 6D cross-seed | **high priority** | Natural Claim 50 follow-up |
| SV8-14 secondary signal | **rescue_once_with_major_change** | Only tested at L8 briefly |
| Graph Laplacian → SV alignment | **new, high-reward** | Never tested |
| Persistent homology in SV2-7 | **new, targeted** | Different from prior PH on full space |
