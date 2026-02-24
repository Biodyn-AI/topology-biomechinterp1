# Brainstormer Structured Feedback — iter_0036

## Gate Status
`passed_min_research_gate: true` — all checks passed, 3 hypotheses tested, paper updated.

## What Worked
All three hypotheses returned positive results, completing the B-cell geometric story within scGPT:
- **H01** (specificity): B-cell clustering is genuinely cell-type specific (z=4.55 vs T-cell z=0.29, 15× difference).
- **H02** (permutation null): artifact hypothesis fully excluded — perm z-scores collapse to mean -0.01 while real z=4.35, emp_p=0.005.
- **H03** (neighborhood content): top-3 neighbors are BATF (7/7 B-cell genes), SPIB (5/7), BACH2 (3/7) — the three canonical B-cell identity TFs. Enrichment 2.98× at L2.

The B-cell story now has 14 converging quantitative claims (see paper). This is publishable as a standalone result.

## What This Iteration Closes
- The "is B-cell clustering real?" question is closed. It is real, specific, non-artifactual, and biologically interpretable.
- The 2-track structure of the paper is now clear:
  - Track A (iters 0002–0017): SV2/SV3/SV4 axes encode PPI/TRRUST co-regulation (lung immune gene set, 209 genes)
  - Track B (iters 0031–0036): B-cell geometric clustering in scGPT (195 in-vocab genes, 7 B-cell markers)
- Both tracks have reached internal completion within scGPT.

## Strategic Assessment
**Current state: ripe for cross-model validation and expansion.**

The highest-value move is no longer exploring new properties within scGPT — it is validating the two established findings (SV2/PPI geometry, B-cell kNN clustering) in Geneformer. Cross-model agreement would dramatically elevate the paper's claim from "scGPT has this property" to "foundation models for genomics encode this as a shared learned representation."

Secondary priority: expand the cell-type panel. T-cell z=0.29 and Myeloid z=1.17 are null results but with small marker sets. Testing NK, DC, and plasma cell markers would map the full immune cell-type geometric landscape.

## Directions to Retire
1. **TRRUST/GO co-annotation → embedding proximity** (repeated negatives at iter_0033, iter_0035): retire permanently. The signal is not there at this resolution.
2. **Intrinsic dimensionality (SVD AUROC for TF sets)**: retire. No signal in multiple attempts.
3. **GO BP enrichment in SV2 poles** (iter_0011 H03): already retired, stays retired.

## No False Positives Detected
The permutation null (H02) is exactly the right control. The neighborhood content (H03) is biologically interpretable without cherry-picking. No concerns about inflation.

## Next Iteration Priority Order
1. **Geneformer B-cell kNN precision@10** — cross-model validation (highest scientific value)
2. **SPIB/BACH2/BATF STRING PPI quantification** — convert "neighbors are B-cell TFs" to "neighbors are known B-cell TF regulators with STRING evidence" (low cost, high narrative value)
3. **Expand cell-type panel to NK/DC/Plasma cells** — map the immune geometric landscape

See roadmap for full hypothesis portfolio.
