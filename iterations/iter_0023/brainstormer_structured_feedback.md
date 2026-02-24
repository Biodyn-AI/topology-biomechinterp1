# Brainstormer Structured Feedback — iter_0023

## Gate Status
`passed_min_research_gate: true` — all checks green. No recovery needed.

---

## Iteration Quality Assessment

**Strongest result**: H03 (TRRUST directional asymmetry) is the most novel and scientifically significant finding to date. Activation AUROC=0.640 vs repression AUROC=0.459 (0/12 sig layers) is a clean, interpretable mechanistic signal: scGPT embedding proximity is a proxy for *co-activation* not *co-regulation* generically. The growing effect with depth (L0=-0.094 → L11=-0.215 distance difference) suggests progressive specialization in later layers.

**H02** establishes GO BP as 3rd biological anchor. Effect size (rho=-0.077) is modest but large-N guarantees high significance (p<1e-19). More useful as a calibration baseline than a discovery in isolation — needs deeper dissection (pathway-level, ontology comparison).

**H01** contamination control closes the validation loop on cell-type geometry. The AUROC=0.488 null result is decisive.

---

## Critical Follow-Up Flags

1. **TRRUST N is low (47 activation, 26 repression)**: the asymmetry needs validation in a larger dataset. Dorothea has hundreds of high-confidence TF-target pairs. This is the single highest-priority validation experiment.

2. **GO BP rho peak at L7**: the non-monotonic layer profile (rho deepens to L7, returns at L11) is unexplained. Could reflect attention head specialization near the middle layers. Compare with STRING and TRRUST layer profiles — if they all peak at L7 this becomes a unifying geometry claim.

3. **Repression AUROC ≈ 0.459 (below chance)**: not quite significant but directional. Repressor TF-target pairs may be *anti-clustered* — they are pushed apart. Test this directly with one-sided MW.

4. **Cell-type cross-model validation is overdue**: we have Geneformer word embeddings from iter_0019. CKA or direct pairwise correlation between scGPT and Geneformer cell-type geometry should be run before over-interpreting scGPT-specific structure.

---

## Directions to Retire or Deprioritize

- **Proxy persistence entropy (histogram entropy)**: confirmed negative, retired by executor.
- **Quintile-binned STRING effects**: superseded by continuous Spearman; retire.
- **Standalone SVD direction counting**: subsumed by co-polarity results; retire.
- **Expanding cell-type panel to macrophage/NK/myeloid/endothelial without named genes**: blocked by named-gene availability; do not pursue unless panel is first confirmed.

---

## Priority Signal for Next Loop

The activation/repression asymmetry is the narrative core of iter_0024. All three top candidates should extend this finding or triangulate it against independent datasets and models.
