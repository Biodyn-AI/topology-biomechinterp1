# Brainstormer Structured Feedback — iter_0047

## Gate status
`passed_min_research_gate: true`. Three hypotheses tested; critical validity correction propagated to paper.

---

## What was genuinely accomplished this iteration

1. **Zero-norm diagnosis closed**: 2902/4941 genes have zero norm at ALL 12 layers — same constant set. These are structurally absent from the embedding (out-of-vocabulary or excluded by the scGPT tokenizer). The valid embedding population is n=2039 throughout all analyses.

2. **TwoNN compression claim consolidated**: L0=34.7 → L11=20.0 (ratio=0.58) holds cleanly on nonzero-only genes. This is now the validated, artifact-free number for the paper.

3. **SVD claim corrected from iter_0046**: The 14× eff_rank collapse and SV1=93.4% variance were driven by zero vectors collapsing onto the origin. The real collapse is 4.9× (L0=236 → L11=48), SV1=18.6% of variance. Still a significant finding, but not the dramatic one previously claimed.

4. **ID compression shape characterized**: Linear (R²=0.986), slope=-1.30/layer. No phase transition. The process is smooth and progressive, not phasic.

---

## What is still missing / blocked

- **Gene names file not found at expected path**. This blocks ALL biological annotation: zero-norm gene identity, SV1/SV2 axes, TRRUST mapping to the nonzero subset. This is the single largest bottleneck. Must be resolved.

- **SVD spectrum shape unexplored**: We know SV1=18.6%, SV2=11.7%, but the full spectral decay curve across all 12 layers has not been computed. The shape of spectral decay (power-law vs exponential vs step) is diagnostically useful.

- **Biological meaning of the 4.9× eff_rank collapse**: We know it exists but not what drives it — is it a few large singular values growing or many small ones shrinking?

- **Cross-cycle replication not done**: All results are still cycle4_immune_main only.

---

## Iteration integrity notes

- The correction from iter_0046 is handled well — both corrected values and the contamination mechanism are clearly described.
- No overclaiming in the current iteration report.
- Biological annotation is genuinely blocked, not avoided; this is an appropriate partial result.

---

## Direction assessment

| Direction | Accumulated result | Recommendation |
|---|---|---|
| TwoNN intrinsic dimensionality | Positive, validated, corrected | Continue: characterize shape, biology |
| SVD eff_rank / spectral structure | Real but weaker than claimed; spectrum is distributed | Continue: spectral decay curve |
| ID phase transition / breakpoint | Negative — linear compression | Retire |
| SV1 biological identity | Blocked (no gene names) | Urgent: find gene names, then execute |
| Zero-norm gene biology | Critical unresolved | High priority |
| TRRUST co-target clustering (PCA-20) | Inconclusive (iter_0004) | Rescue with nonzero-only + larger gene set |
| Cross-model alignment (low-dim summaries) | Inconclusive (iter_0003) | Deprioritize |
| kNN clustering coefficient | Strongly positive (iter_0002) | Use as benchmark; test layer-specificity |
| H1 persistent homology | Positive (iter_0003) | Retest on nonzero-only clean set |
| Cross-layer CKA | Strongly positive (iter_0004) | Mine further: non-linear CKA, layer-pair structure |
