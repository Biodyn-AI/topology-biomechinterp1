# Brainstormer Structured Feedback: iter_0052

## Research Gate
Passed. Three hypotheses tested, one positive (H03), one negative (H02), one inconclusive (H01).

---

## Per-Hypothesis Assessment

### H01 — Housekeeping gene enrichment in SV1-high: INCONCLUSIVE
- Root cause: 7/238 HK reference genes present in the 2039-gene nonzero set. Any signal is noise.
- At L0, 4 HK genes fall in SV1-**low** (not high), direction opposite to hypothesis.
- This does NOT refute the SV1 anti-TF axis; it just says the HK gene interpretation is untestable with this reference set.
- Action: **rescue_once_with_major_change** — use GO:0003723 / GO:0005737 constitutive process annotations or derive a "high-connectivity" node set from STRING degree within the 2039-gene graph as a functional proxy for housekeeping status.

### H02 — H1 Betti loops on circuit genes at L8: NEGATIVE (retire)
- Circuit genes form fewer (108 vs 119) and shorter (mean 0.00143 vs 0.00178) H1 loops than non-circuit genes in SV2-4 at L8.
- Direction is reversed: non-circuit genes are geometrically more loopy, likely because they are an incoherent mixture of functional classes.
- No rescue path for this specific framing. H1 topology does not detect regulatory feedback topology in this projection.
- Action: **retire_now** — documented in executor_hypothesis_screen.json.

### H03 — SV5-7 regulatory signal after co-expression residualization: POSITIVE
- rbc = 0.148, 0.119, 0.083 at L0/L1/L2 (all p<0.002); drops to non-significant at L3+.
- Contrast: SV2-4 retains signal only at L8. Together these imply a layer × spectral-subspace interaction in regulatory encoding.
- This is the most actionable positive result in recent iterations.
- Immediate priorities:
  1. Bootstrap validation (100 resamples, stratified).
  2. Extend spectral scan to SV8-10, SV11-20 to map full decay profile.
  3. Biological annotation of SV5-7 top-loading genes at L0 (GO enrichment).
  4. Probe whether the L0-L2 vs L8 dissociation maps onto known scGPT architecture (early attention vs late MLP layers).

---

## Thematic Signal Summary

The clearest emerging pattern across ~52 iterations:

| Subspace | Layer of peak signal | rbc | Interpretation |
|----------|---------------------|-----|----------------|
| SV2-4    | L8                  | ~0.08 (residualized) | Deep regulatory encoding |
| SV5-7    | L0-L2               | 0.083–0.148 (residualized) | Early regulatory encoding |
| SV1      | —                   | anti-TF axis (SV1-low = TF-annotated genes) | Dominant variance axis ≠ regulatory |

This layered spectral structure warrants systematic characterization. It is not a noise artifact — two independent spectral subspaces show orthogonal but complementary regulatory signals at different processing depths.
