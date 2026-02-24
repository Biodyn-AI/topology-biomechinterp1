# Brainstormer Structured Feedback: iter_0048

## Research Gate Status
`passed_min_research_gate: true`
Two of three hypotheses returned meaningful results; execution is on track.

---

## H01: Gene Annotation + SV1/SV2 Identity — NEUTRAL

**What worked**: Zero-norm gene biology validation is a clean, publishable supporting observation.
66 named zero-norm genes are systematically non-immune (TFAP2B, SNAI2, TYR, PGR, TBX5, CCL22, AMH).
This confirms the structural exclusion boundary between immune-active and immune-irrelevant gene programs.

**What failed**: SV1/SV2 identity remains unresolved. All cell-type markers score near zero on SV1
(B-cell=-0.017, T-cell=-0.016, Myeloid=-0.009). SV1 captures something other than lineage identity.
Top SV1 genes are unnamed — this is the critical gap.

**Critical next step**: Recover full gene index→name mapping for 2039 active genes. The SV1 axis
(18.6% of variance at L11) must be identified. Without this, the spectral decay story is incomplete.

**Status: Not retired. Rescue with gene name recovery immediately.**

---

## H02: Full Spectral Decay Curve — STRONG POSITIVE

This is the iteration's main result. Four independent metrics converge on the same signal:
- eff_rank: 236.9 → 48.7 (4.87×)
- SPR: 90.1 → 15.2 (5.93×)
- k90: 307 → 102 (3.01×)
- Frobenius norm: 558 → 204 (2.74×)

**Novel observation**: SV1 singular value peaks at L8 (143.9) then drops sharply to L11 (87.8, -39%),
while SV1 variance fraction grows monotonically (5.6% → 18.6%). This two-phase structure —
consolidation through L8, global contraction L8→L11 — is mechanistically distinct and not previously
characterized in this project.

**Critical next step**: Distinguish whether spectral compression is uniform across all gene programs
or biologically differentiated (B-cell subgroup compresses faster than T-cell or random).
Also: characterize the SV1 direction stability across layers — does the dominant axis rotate or persist?

**Status: Promising — extend biologically and mechanistically.**

---

## H03: Wasserstein/Transport Distance — INCONCLUSIVE

SWD range 0.29-0.55, with the largest step at L0→L1 and the smallest at L10→L11.
Spearman ρ(SWD, eff_rank_drop) = -0.46, p=0.15. The key problem: without a shuffle null,
absolute SWD values are uninterpretable. The centroid displacement (grows at late layers)
and SWD (shrinks at late layers) capture orthogonal processes — this is potentially interesting
but currently not validated.

**Next step**: Add shuffle null for SWD. If real SWD < shuffled SWD, transport is
constrained by structure. This is a low-cost resolution.

**Status: Inconclusive but salvageable in one cheap experiment.**

---

## Cumulative Direction Assessment

**Retiring now:**
- Lineage centroid cosine orthogonality (iter 46: z≤0.26, null)
- TCR circuit attractor (iter 44-45: n=2 circuit genes, never reaches significance)
- Cross-lineage ID compression law via 5-point TwoNN (iter 44: methodological incompatibility)
- SV1/SV2 cell-type identity approach without full gene vocabulary (blocked three times)

**Directions still open and productive:**
- Spectral decay characterization (two-phase dynamics, biological differentiation)
- SV1 identity (high priority — blocked only by missing gene names)
- Wasserstein/transport with null (one experiment away from resolution)
- PPI geometry replication in 2039-gene active set (SV2-SV4 axes)
- GC attractor characterization (well-supported, room for convergence geometry)
- Cross-model alignment (under-explored with corrected gene set)
