# Next Steps — iter_0028

## H01 (Hub Centrality) — follow-up
- **Immediate**: Characterize per-layer hub centrality effect size with confidence intervals (bootstrap).
- **Biological anchor**: Test whether the most-central hub genes (ETS1, IRF8, REL, ITGB2) are known master regulators (TF DB, BioGrid hub list). Compute fraction of top-10 central genes that are known TFs.
- **Cross-model**: Test same STRING-degree vs centrality correlation in Geneformer embeddings (if available). Would provide cross-model validation.
- **Novel extension**: Test hub centrality in Geneformer residuals or apply to alternative connectivity metrics (PageRank vs raw degree).

## H02 (PC1 axis) — retired
- L11 PC1 does not encode T-cell/APC axis. Both gene families collapse to same PC1 position due to dimensional collapse.
- Do not revisit unless radically different framing (e.g., higher PCs, or use of cell-type embedding rather than gene embedding).

## H03 (Spectral Gap) — follow-up
- **Priority**: Add random-Gaussian-embedding null to confirm spectral gap trend exceeds baseline.
  ```python
  null_embs = rng.standard_normal((N_NAMED, 512))
  gap_null, _, _, _ = compute_spectral_gap_and_entropy(null_embs, k=10)
  ```
- **Extension**: Compute spectral gap as a function of k (k=5, 10, 15, 20) to assess robustness.
- **Biological anchor**: Identify the connected components at L11 (n_components=2). Which genes form each component? Do they correspond to cell-type clusters?

## New directions for iter_0029
1. **Hub centrality bootstrap CI + TF enrichment** (follow-up H01)
2. **Spectral gap random-embedding null** (follow-up H03)
3. **New family: Ollivier-Ricci curvature** on kNN graph per layer — curvature captures local geometry beyond spectral properties.

## Retirements confirmed
- PC1 T-cell/APC axis: retired (H02 iter_0028 = negative)
- Dorothea same-pool regulatory proximity: already retired (iter_0027 H03)
- TRRUST activation/repression polarity at population level: already retired (iter_0027 H02)
