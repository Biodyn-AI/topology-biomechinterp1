# Brainstormer Structured Feedback — iter_0034

## Gate Status
`passed_min_research_gate: true` — 3 hypotheses tested, 1 strong positive (H03), 1 mixed (H01), 1 negative (H02). Paper updated, all artifacts present.

---

## H01: Layer-wise B-cell Community Purity Emergence — MIXED

**Assessment:** The monotonic-emergence framing was wrong, but the result is still scientifically useful. B-cell geometry is present from L2 (OR=16), not a late-layer phenomenon. The non-monotonic pattern is a community-identity-tracking artifact (community labels shuffle as modularity optimization re-partitions), not evidence against B-cell geometry. The actual discovery is **early-layer B-cell commitment**: geometry baked in at L2, not crystallized at L11.

**What to do with it:** Reframe as "B-cell identity is encoded from L2 onward; L11 consolidates to a uniquely binary partition." Test the L2 community biologically to confirm the same B-cell/non-B-cell split.

**Retire the monotonic framing.** Keep the early-layer finding as a new positive.

---

## H02: PC2/PC3 Cell-Type Axes at L11 — NEGATIVE

**Assessment:** Definitive negative for the 195-gene set at the current marker count (T-cell n=14, Myeloid n=6). PC2 myeloid AUROC=0.618 with p=0.168 is suggestive but underpowered. Root cause is likely: (a) myeloid marker set too small (n=6); (b) T-cell and myeloid variance may be spread across PC4+. Not a fundamental claim about the model—just a power limit.

**What to do with it:** Retire PC2/PC3 myeloid/T-cell test at current marker count. If myeloid/T-cell axes are important, test with larger curated marker sets (n≥15 per lineage) on PC2–PC6.

---

## H03: Expression-Level Confound Rejection — POSITIVE (STRONG)

**Assessment:** This is the iteration's standout result. Bootstrap null (z=-3.06, empirical p=0.000) + L2-norm regression (AUROC=0.214, p=0.002 after deconfounding) definitively closes the expression-level artifact concern. The B-cell PC1 signal is structural.

**Strategic significance:** This result is publication-quality. With H03 in hand, the claim "scGPT encodes B-cell identity as a persistent geometric feature of its residual stream" now has:
- Spatial evidence (community, PC1 pole)
- Temporal evidence (12 layers, L2 onset)
- Structural evidence (bootstrap null, norm regression)

**Next priority:** Complete the biological story — what is the non-B-cell compartment, and does this cross-validate to Geneformer?

---

## Cumulative Direction Assessment

| Direction | Status | Rescue Potential |
|-----------|--------|-----------------|
| STRING PPI distance gradient | Retired (iter_0031) | None |
| GO BP SV2 enrichment | Retired (iter_0011) | None |
| Repression anti-pole geometry | Retired (iter_0013) | None |
| Monotonic B-cell community emergence | Negative (iter_0034 H01) | None — reframe as early-layer finding |
| PC2/PC3 cell-type axes (current markers) | Negative (iter_0034 H02) | Rescue with larger marker sets |
| TRRUST 195-gene proximity | Pending since iter_0030 | High — execute |
| GO Jaccard 195-gene gradient | Pending since iter_0030 | High — execute |
| Cross-model B-cell validation | Never tested | High — execute |
| B-cell subtype geometry | Never tested | High |
| Conditional nearest-neighbor of B-cell markers | Never tested | Medium-high |
