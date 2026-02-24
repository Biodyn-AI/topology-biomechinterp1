# Brainstormer Structured Feedback — iter_0060
Date: 2026-02-23

## Gate Status
`passed_min_research_gate: true` — all three hypotheses ran to completion with artifacts.

---

## H01: Signed Displacement Projection AUROC
**Decision: NEGATIVE → RETIRE**

Directional projection onto mean TF→target displacement gives max AUROC=0.563 at L0, statistically real (perm_p=0.001) but effectively identical to scalar distance (0.565). The brainstormer prediction of >0.62 is falsified. Two successive iterations (H02 iter_0059, H01 iter_0060) confirm a hard ceiling at ~0.56-0.57 for SV5-7 pairwise displacement geometry on edge prediction. This method family is exhausted. Retire all directional-displacement-in-SV5-7 variants.

**What this tells us:** The ~0.56-0.57 AUROC ceiling is a property of the SV5-7 subspace geometry, not method choice. The TRRUST edge signal in pairwise geometry is weak. The strong signal (AUROC 0.73-0.76) found in earlier iterations for TF-vs-target *classification* operates in a different geometric regime than pairwise edge geometry.

---

## H02: TF Family Margin Trajectory
**Decision: NEUTRAL — Rescue with focused follow-up**

**What's real:**
- bHLH (HIF1A, n=1): Sign reversal at L8 — consistently negative margins L0-L7, strongly positive L8-L11 (margin +1.68 → +2.76). This is striking but single-gene.
- bZIP (BATF/FOS/JUN/JUNB, n=4): Consistently target-like across all layers (margin -0.77 to -0.98). Robust pattern.
- Overall AUROC on TRRUST-subset genes increases monotonically L0=0.570 → L11=0.739. This contradicts earlier reports on full gene set and needs explanation.

**What's unclear:**
- The AUROC increase may be an artifact of training+evaluating on the same TRRUST-subset. Not cross-validated, not a true generalization test.
- bHLH n=1 is insufficient for a family-level claim; HIF1A may be unique because it is both a regulated gene (oxygen-dependent) and a master hypoxia TF.

**Next action:** Follow HIF1A neighborhood composition change at L8 (nearest-neighbor identity before/after flip). Verify bZIP target-like pattern against GO/GSEA.

---

## H03: FFL Triangle Geometry
**Decision: INCONCLUSIVE — High rescue potential**

**What's real:**
- t_mean signal: significantly >0 at L0-L6 (p<0.0001 at L2). The intermediate TF B is displaced ~0.6 of the way from master regulator A toward target C on average. This is a novel geometric finding.
- Signal collapses sharply at L8-L11 (t_mean → 0.1, p>0.45). This aligns with the HIF1A flip boundary from H02.

**What's uncertain:**
- N=22 FFLs is small. Betweenness permutation test (binary, lower-powered) is NS at p=0.30.
- SV5-7 may not be the optimal subspace. SV2-4 or SV1-4 not yet tested for FFL geometry.

**Critical observation:** The same L8 boundary appears in both H02 (bHLH flip) and H03 (t_mean collapse). L8 is a structural transition in this model's representation space — this is worth characterizing directly.

**Next action:** Expand FFL set (target N≥50 via STRING TF-TF links), test in SV2-4, use t_mean (continuous) as primary metric rather than binary betweenness.

---

## Cross-Hypothesis Synthesis

Three converging signals point to **L8 as a structural boundary**:
1. HIF1A margin sign reversal at L8 (H02)
2. FFL t_mean collapse at L8 (H03)
3. The monotonic AUROC increase across layers (H02 TRRUST-subset)

This is the most promising emergent pattern from iter_0060. The next iteration should directly characterize what changes at L8 globally, not just within specific motifs.

**Persistent unresolved question:** What drives the strong AUROC ~0.73-0.76 TF/target boundary at L2-L3 (prior iterations)? We have ruled out: hub degree, effective rank, directional edge geometry. Candidates remaining: intrinsic dimension, local curvature, neighborhood composition, attention head assignment.
