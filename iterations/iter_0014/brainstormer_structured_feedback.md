# Brainstormer Structured Feedback — iter_0014

## Gate Status
`passed_min_research_gate: true` — 3 hypotheses tested, 2 positive (H01, H03), 1 neutral (H02). All artifacts present and valid.

---

## What Worked

**H01 (GO enrichment SV3)**: Best result this iteration. Establishes SV3 as a biologically interpretable axis with a depth-progressive biological transition (kinase/signaling early → nuclear/immune regulatory late). 12/12 layers, all emp_p ≤ 0.017. This is a standalone publishable finding: scGPT's third singular geometry axis encodes a cell-biology function gradient across transformer depth.

**H03 (Hub-degree control)**: Definitive methodological closure on the hub confound. Non-hub mean z=3.11 (12/12 sig) exceeds hub z=2.60 (9/12 sig). The paper's PPI geometry claim now has a robust null control. This closes a critical reviewer attack surface.

## What Was Neutral

**H02 (Spearman layer trend)**: Low information gain. Marginal SV2 decline (rho=-0.56, p=0.059) and non-monotonic SV3 — neither finding sharpens the core argument. The result is worth one sentence in the paper but doesn't open a new direction.

---

## Direction Assessment

| Direction | Status | Rationale |
|---|---|---|
| SV1 as PPI axis | RETIRED (prior) | z=0.39, 2/12 layers |
| Repression anti-pole | RETIRED (prior) | 0/12 sig, z=-1.41 |
| Layer-depth monotonicity of SV2/SV3 | CLOSE — one more test at best | Marginal SV2, non-monotonic SV3. Not a productive direction. |
| SV2 PPI geometry (core claim) | CONFIRMED, consolidation phase | 12/12 layers, hub confound cleared |
| SV3 biological characterization | ACTIVE — high value remaining | GO done; drug/perturbation and community structure untested |
| Joint SV2+SV3 structure | ACTIVE — not yet tested | 2D geometry, community clustering, angle trajectory |
| GO co-annotation vs PPI prediction | ACTIVE — not yet tested | High-value discriminator |
| SV4+ axes | NEW — unexplored | Open question: are SV2/SV3 the only geometry axes? |

---

## Gaps in Current Evidence

1. **No perturbation/drug anchor for SV3**: The GO results identify kinase/immune biology but don't connect to any experimental readout. This is the biggest clinical relevance gap.
2. **No cross-axis geometry test**: SV2 and SV3 are treated independently. Their joint structure in 2D has not been analyzed.
3. **No higher SVD axes explored**: SV4/SV5 may carry additional structured biology.
4. **No test of GO co-annotation vs PPI as predictors**: We haven't asked whether scGPT geometry reflects functional vs physical proximity.
5. **No edge-weight gradient test**: Do STRING confidence scores correlate with SV2 proximity (higher-confidence PPI → closer in SV2 space)?
