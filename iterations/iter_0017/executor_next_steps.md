# Executor Next Steps — iter_0017

## Project Status: Near-Complete Results Narrative

The multi-axis PPI geometry story is now fully characterized:

| Axis | Confidence rho | Signed regulation | GO theme |
|------|---------------|-------------------|----------|
| SV2  | 1.000 | Yes (activation 12/12) | Innate immunity / exosome |
| SV3  | 0.900 | No (0/12) | Plasma membrane / cell surface |
| SV4  | 0.900 | No (2/12 borderline) | Receptor binding / apoptosis |

## Top Priorities for iter_0018

### Priority 1: Paper consolidation (mandatory)
- Write a complete Results section summarizing the 4-axis story
- Include quantitative table: axis, rho, GO theme, signed regulation
- Add Figure description for the quintile gradient plot
- Ensure Methods section has reproducible details for SVD, co-pole, quintile sweep

### Priority 2: SV2+SV3 joint 2D co-pole test (H-D from brainstormer)
- Test if STRING pairs appearing in BOTH SV2 and SV3 co-poles simultaneously have higher confidence than single-axis co-pole pairs
- This would strengthen the "multiple axes provide additive information" claim
- Estimated: ~10-line extension of existing scripts

### Priority 3: Cross-model (Geneformer) — if data surfaces
- Geneformer data not found in this iteration. If discovered in future, run SV2 co-pole test.
- Fallback: document as future work in paper Discussion.

### Priority 4: SV5+ axis screening
- Test if SV5 (index 4) also shows a confidence gradient
- If positive, extend the story to 4+ axes

## Retired Directions
- TRRUST signed regulation on SV3/SV4: retired (0/12, 2/12 layers)
- Single-gene module permutation tests: retired (prior iterations)
- TRRUST repression anti-pole (SV2): retired (iter_0013)
- GO BP enrichment single-axis SV2 poles: retired

## Key Evidence Locked In
1. SV2 quintile rho=1.000 (iter_0015 H02)
2. SV3 quintile rho=0.900 (iter_0016 H01)
3. SV4 quintile rho=0.900 (iter_0017 H01)
4. Axis independence: max inter-axis r=0.247 (iter_0016 H02)
5. Biological orthogonality: GO Jaccard ≤ 0.081 (iter_0017 H02)
6. Hub-degree confound ruled out (iter_0014 H03)
7. GO co-annotation confound ruled out (iter_0015 H03)
8. SV2 activation co-pole 12/12 layers (iter_0012 H02)
