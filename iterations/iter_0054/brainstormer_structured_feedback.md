# Brainstormer Structured Feedback — iter_0054

**Date**: 2026-02-23
**Gate**: PASSED (all 3 hypotheses significant, 2 rated promising)

---

## Assessment of iter_0054 Results

### H01: TF vs Target-only Classification (AUROC=0.694)
**Verdict: Strong positive, ready for consolidation.**
- AUROC 0.694 vs null 0.494, p=0.000 — this is a clean, permutation-validated result.
- Early-layer specificity confirmed (L8 at chance). SV2-4 showed similar point estimate but high variance — instability suggests SV5-7 is the genuinely informative subspace.
- Immediate priority: cross-seed replication. If seed43/44 both exceed 0.65, this is a paper-anchor result: "the model's early-layer geometry functionally distinguishes transcription factors from their targets without any network features."
- Next upgrade path: test SV1-7 combined (7D). If AUROC breaks 0.70, that's the final claim. If not, 0.694 stands.

### H02: BFS Depth vs SV5-7 (rho~0.17)
**Verdict: Signal exists but weak; needs partial control before claiming.**
- rho=0.167-0.169 on SV5/SV7, absent at L8. Effect is real but small.
- Key question: is the depth-SV correlation an artifact of the TF/target binary signal? Partial Spearman controlling for TF/target label is the mandatory next test.
- Do not amplify this finding in the paper until confound is resolved.
- If partial correlation survives, this is a genuinely new claim about cascade position encoding.

### H03: Layer Directionality Trajectory (L0=0.0057 → L9=0.0200)
**Verdict: Best result of the iteration. Publication-quality trajectory.**
- 3.5× amplification from L0→L9, peak at L9, all 12 layers significant on ≥1 axis.
- NEW observation: axis rotation across layers (SV6 dominant L1-L4, SV5 dominant L5/L8-L9, SV7 dominant L6-L7/L10, SV6 at L11). This is not noise — it suggests the geometry is doing something structured through the network, not just uniformly amplifying.
- The L9 peak (not L11 = final layer) is conceptually interesting: the model reaches maximum regulatory discriminability one layer before the end.
- SV2-4 control trajectory is the essential complement to make the two-subspace story complete.

---

## Key Gaps Identified

1. **No cross-seed validation yet** — H01 and H03 are main-seed-only. This is the single biggest risk to paper credibility.
2. **Axis rotation phenomenon untested** — the layer-by-layer shift of dominant directionality axis is observed but not characterized. It could be a key mechanistic finding.
3. **No biological interpretation of directionality axis** — GO enrichment of the displacement axis (what genes are high-projection?) remains undone.
4. **No test at peak layer (L9) for classification** — H01 tested L0 and L8. L9 is now identified as peak directionality layer; classification test at L9 is missing.
5. **TRRUST edge subset analysis absent** — activator vs repressor edges have not been tested for differential directionality signature.

---

## Directions with Rescue Potential

| Direction | Status | Rescue Path |
|-----------|--------|-------------|
| BFS depth encoding (H02) | Neutral, weak | Partial Spearman → if survives, claim stands |
| Edge prediction from SV geometry | Proposed, untested | Medium cost; test at L9 (peak layer) rather than L0 |
| Topological persistence by gene class | Proposed, untested | Worth one cheap iteration using 0-dim PH |
| GO enrichment of displacement axis | Proposed, untested | Very cheap; should be done regardless |

## Directions to Retire

| Direction | Reason |
|-----------|--------|
| TwoNN intrinsic dimension | Retired iter_0047; no revival candidate |
| Betti loops (PH on small gene sets) | Consistently negative across multiple iterations |
| kNN cell-type purity | Wrong framing for this data |
| Continuous degree → SV1 correlation | Confirmed negative 3 times |
| SV8-14 secondary signal | rbc too small for dedicated iteration |
