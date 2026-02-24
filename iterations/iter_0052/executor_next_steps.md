# Next Steps: iter_0053

## Priority 1 (high-probability): SV5-7 Layer-Stratified Validation

H03 found strong co-expression-independent regulatory signal in SV5-7 at L0-L2 (rbc=0.148, 0.119, 0.083).
- **Next**: test robustness with bootstrap resampling of positive/negative pairs (100 bootstraps).
- **Also**: compare against SV8-10 and SV11-20 to determine if the signal decays monotonically with SV index, or if there are other non-trivial spectral axes.
- Family: intrinsic_dimensionality / module_structure
- Method: same OLS residualization framework, extended to higher SVs

## Priority 2 (new family): SV5-7 Biological Annotation

The SV5-7 signal at L0-L2 needs biological grounding:
- What gene families cluster in SV5-7 projection at L0?
- Are they enriched for specific GO biological processes or KEGG pathways?
- Use gene set over-representation test (Fisher exact against GO slim categories using named genes)
- Family: module_structure
- Method: GO enrichment on SV5-7 top-loading genes

## Priority 3 (housekeeping gene rescue): Broader Reference Gene Set

H01 failed due to only 7 HK genes in nonzero set.
- **Rescue**: Use a reference set matched to the dataset — e.g., the set of all 2039 nonzero genes' GO annotations (if obtainable), or use expression-level data to identify constitutively expressed genes from the scRNA-seq data itself.
- Alternative: characterize SV1-high genes by their network connectivity in STRING/TRRUST (nodes with 0 known interactions = unnamed/uncharacterized).
- Family: intrinsic_dimensionality
- Split: same edge_stratified

## Retired directions (do not revisit)
- H1 Betti loops for regulatory circuit detection: NEGATIVE (iter_0052 H02)
- SV2-4 regulatory proximity (co-expression confounded): RETIRED (iter_0051 H01)
- Graph topology kNN purity: NEGATIVE (iters 0045, 0046)
- TwoNN intrinsic dimensionality (multi-layer): inconclusive

## Summary of active leads
1. **SV5-7 regulatory signal at L0-L2** — most urgent, new positive result
2. **SV1 anti-TF axis** — needs housekeeping confirmation with broader gene coverage
3. **SV1-high vs SV1-low identity** — biological characterization open
