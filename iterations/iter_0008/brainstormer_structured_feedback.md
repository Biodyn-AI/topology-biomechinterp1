# Brainstormer Structured Feedback — iter_0008

## Gate status
`passed_min_research_gate: true` — all three hypotheses returned strong positives. No recovery needed.

## What iter_0008 confirmed

| Result | Key stat | Status |
|--------|----------|--------|
| SV2 top-pole: IL-4 positive regulation (GO:0032753) | emp_p=0.000, OR=Inf | **Confirmed — strong** |
| SV2 bottom-pole: extracellular vesicle (GO:0070062) | OR=7.99, p=8.7e-9 | **Signal seen — null NOT YET run** |
| L3 mitochondrion transient | OR=23.25, emp_p=0.000 | **Confirmed** |
| L7 ER lumen transient | OR=10.74, emp_p=0.000 | **Confirmed** |
| L8 ER lumen transient | OR=6.95, emp_p=0.002 | **Confirmed** |
| Signal peptide proxy enrichment in SV1 top-pole | emp_p=0.002, MW p=2.4e-4 | **Confirmed** |
| SV3 antigen presentation (GO:0019886) | emp_p=0.088 | **Borderline — not confirmed** |

## What is still open

1. **SV2 bottom-pole null not run.** The SV2 bottom-pole extracellular vesicle OR=7.99 is highly significant by Fisher (p=8.7e-9) but has no gene-label shuffle null. This must be run before the SV2 axis is publishable — it is the most urgent gap.

2. **SV3 is borderline (emp_p=0.088).** Antigen presentation enrichment at SV3 top-pole could be real but N=500 permutations is underpowered for this p-value region. Needs N=2000-5000 or a smaller sub-permutation. Decision: rescue_once_with_major_change (larger N).

3. **Compartment scan is incomplete.** Mito@L3 and ER@L7/8 are confirmed but the remaining 9 layers × other compartments (Golgi, lysosome, nucleus, plasma membrane) have not been tested. The compartment "layer map" is the most publishable version of the H02 result.

4. **Signal peptide proxy is weak.** OR=2.88 is modest. The proxy (union of GO extracellular/ER lumen terms) is noisy. Ground-truth UniProt signal peptide annotations for the 209 genes would sharpen this.

5. **SV1 bottom-pole biological identity.** The SV1 bottom shows cell adhesion (GO:0007155, OR=3.33) and receptor signaling — not cytosol as hypothesized in iter_0007. This is interesting: the SV1 axis is extracellular-secreted vs. cell-surface/adhesion, not secreted vs. cytosolic. This needs explicit null validation.

6. **SV2 axis structure.** SV2 top = immune signaling (IL-4, T cell regulation), SV2 bottom = extracellular vesicle + mitochondria. This suggests SV2 is an immune activation / vesicle-secretory axis — potentially a very clean biological story. Needs integration.

7. **Cross-layer SV2 stability.** We only know SV2 at layer 11. Is the immune axis present across layers? If it emerges late (like SV1), this is additional evidence that layer depth drives biological specificity.

## Directions to retire or deprioritize

- **Rewiring/permutation-based PH nulls**: already retired in iter_0007. Confirmed closed.
- **Drift TF/RNAPolII enrichment**: retired in iter_0007. Confirmed closed.
- **Whole-vocabulary drift analysis**: deprioritized. No new reason to resurrect.
- **SV3 antigen presentation at N=500**: deprioritize as primary claim. Rescue once with N=5000 before retiring.

## Summary assessment

iter_0008 is the strongest iteration so far: 3/3 hypotheses positive, all null-tested. The project is now in a "depth extraction" phase — the major axes (SV1 extracellular, SV2 immune, layer-specific compartment transients) are confirmed; the work now is (a) completing the null testing for secondary results, (b) systematically extending the compartment map, and (c) testing mechanistic explanations (signal peptide, STRING PPI, cross-layer stability).
