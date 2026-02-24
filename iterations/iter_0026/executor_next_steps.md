# Next Steps — iter_0026 → iter_0027

## Confirmed findings to build on

1. **Immune family clustering (H01):** HLA-I AUROC=1.0, AP1 AUROC=0.97, RUNX 0.93 at layer 8.
   - Next: Test larger panel of immune families; test across all 12 layers to identify family-specific peak layers.
   - Open question: Do structural (MHC) vs regulatory (AP1, RUNX) families peak at different layers?

2. **Intrinsic dimensionality collapse (H03):** PR drops from 21→1.7 monotonically.
   - Next: Bootstrap this estimate (resample gene subsets) to confirm it's not driven by outliers.
   - Cross-check: Do PR curves differ between functional gene subsets (e.g., TFs vs surface markers)?
   - Compare PR progression across biological families vs shuffled gene labels.

3. **Dorothea confidence gradient (H02):** High-conf AUROC=0.66 vs null, moderate high-vs-low gradient.
   - Need: Dorothea file with MOR (mode-of-regulation, ±1) to test activation vs repression direction polarity.
   - Alternative source: OmnipathR via Python, or use TRRUST activation/repression split (already shown in iter_0023).

## Top 3 hypotheses for iter_0027

### H-A: Family-specific layer peak analysis (new_method)
- For each immune family (AP1, RUNX, HLA-I/II, IL2_path, etc.), compute AUROC curve across all 12 layers.
- Test: Does each family peak at a different layer? (structural vs regulatory divergence)
- Expected: HLA peaks early (structural), TF families (AP1, RUNX) peak mid-to-late.
- Cost: very low (reuse existing infrastructure).

### H-B: PR bootstrap robustness (topology_stability, new_method)
- Bootstrap gene subsets (100 resamples, 80% of 208 genes) and recompute PR at each layer.
- Test: Does monotonic collapse survive? Compute 95% CI bands.
- Expected: Strong finding — confirms dimensional collapse is robust.

### H-C: Per-family intrinsic dimensionality (cross-family ID divergence)
- Compute PR separately for each immune family at each layer.
- Test: Do structurally similar families (HLA-I, HLA-II) have lower intrinsic dim than functionally diverse families (TNFSF, CCL)?
- Novel: links ID to biological coherence.

### H-D: TRRUST activation vs repression proximity (already validated families)
- Use TRRUST pairs (iter_0023 showed activation AUROC=0.68 in lung).
- Extend: compare activation vs repression TF-target distances at each layer.
- Already have TRRUST data in prior iter infrastructure.

## Retired directions (do not repeat)
- Null shuffle of gene labels: already validated (iter_0020)
- Chromosomal proximity: definitively negative (iter_0025 H03)
- Pure STRING univariate sweeps without new split regime
