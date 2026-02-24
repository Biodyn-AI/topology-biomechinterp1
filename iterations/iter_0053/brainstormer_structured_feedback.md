# Brainstormer Structured Feedback — iter_0053

**Date**: 2026-02-23
**Research gate**: PASSED (all 3 hypotheses tested, artifacts present, paper updated)

---

## Assessment of iter_0053 Results

### H01 — Bootstrap + Spectral Decay (POSITIVE)
Strong result. SV5-7 at L0 is now the project's best-validated finding: 100/100 bootstrap positive, mean rbc=0.147, CI=[0.084, 0.199]. The spectral decay table is a clean structural finding: early-layer encoding (SV5-7) vs deep-layer encoding (SV2-4) are orthogonal regulatory subspaces. This map of the spectral landscape is paper-ready as a standalone figure.

**What it opens**: Comparative geometry of SV5-7 vs SV2-4. Edge-level decomposition. Whether the two subspaces capture different regulatory event types.

### H02 — TF→Target Directional Asymmetry (POSITIVE, new finding)
This is the most important new structural discovery: the embedding geometry is directed. All three SV5-7 axes show systematic TF vs target displacement at L0 (Cohen's d≈0.12–0.14). The effect amplifies dramatically at L8 (p_combined<0.001). This is not a proximity-only story — the manifold has an intrinsic orientation aligned to regulatory flow.

**What it opens**: Classification of TF vs target role from geometry alone. Regulatory hierarchy depth in the SV5-7 space. Whether L8 SV2-4 encodes even stronger directionality. Whether the direction vector generalizes across seeds.

### H03 — SV1 Circuit Identity (NEUTRAL, confirms prior)
Confirms iter_0051 H02: circuit membership (binary) is the primary SV1 discriminant at L0. Continuous degree is a weak predictor (rho=−0.086). The finding is solid but not novel. Value: confirms SV1 axis is a regulatory-vs-background axis, not a hub-centrality axis. This makes SV1 + SV5-7 a natural multi-axis classifier system.

**What it opens**: Joint SV1 + SV5-7 classification of gene regulatory roles. SV1 might separate circuit membership while SV5-7 separates TF vs target.

---

## Directions Status

| Direction | Status | Verdict |
|-----------|--------|---------|
| SV5-7 regulatory proximity | FULLY VALIDATED | Keep, extend |
| TF→target directional asymmetry | NEW POSITIVE | Exploit aggressively |
| SV1 circuit identity | CONFIRMED (3rd time) | Integrate, don't re-test alone |
| SV2-4 at L8 regulatory proximity | PARTIALLY CONFIRMED | Test directionality at L8 |
| SV8-14 at L8 secondary signal | SUGGESTIVE | Low priority |
| Betti loops for regulatory circuits | RETIRED (iter_0052) | Do not revisit |
| TwoNN intrinsic dimensionality | RETIRED (iter_0047) | Do not revisit |
| kNN purity for cell-type markers | RETIRED (iter_0045) | Do not revisit |

---

## Priority Signals for Next Iteration

1. **TF role classification (AUROC)** — H02 directional asymmetry is significant but Cohen's d≈0.12 is small. The key test is whether AUROC > 0.6 from SV5-7 coordinates alone. This determines whether the geometry is predictively useful or just statistically significant.

2. **L8 SV2-4 directionality** — L8 directionality in SV5-7 (p_combined<0.001) is much stronger than L0. Test whether SV2-4 at L8 shows even stronger directional encoding — if yes, the deep-layer subspace may be the more useful predictor.

3. **Cross-seed replication of SV5-7** — bootstrap validates stability within one run. Seed replication validates independence. Essential before any paper claim about robustness.

4. **Multi-axis joint classifier (SV1 + SV5-7)** — natural synthesis: SV1 separates circuit vs non-circuit; SV5-7 separates TF vs target within circuit genes. A 4D logistic regression (SV1 + 3×SV5-7) might classify TF identity with high AUROC.
