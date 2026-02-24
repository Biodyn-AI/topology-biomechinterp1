# Next Steps after iter_0047

## Critical corrections needed in paper

1. **SVD eff_rank finding from iter_0046 must be corrected**: The 14× collapse and SV1=93.4% claims were artifacts of including 2902 zero-norm genes. Correct claims: 4.9× collapse (236→48), SV1=18.6% in nonzero-only analysis.
2. **TwoNN claim stands**: L0=34.7→L11=20.0 (ratio=0.58) on nonzero-only genes confirmed.

## Immediate priorities for iter_0048

### H1 (NEW — HIGH PRIORITY): Find gene names file and do SV1 biological annotation
- Look for gene names in other output directories or implementation code
- Paths to check: `cycle4_immune_main/` directory listing, implementation scripts
- Goal: identify what SV1 (18.6% variance) and SV2 (11.7% variance) represent biologically
- Key question: does SV1 separate B/T/Myeloid cells, or reflect expression levels?

### H2 (NEW): Characterize zero-norm vs nonzero-norm genes biologically
- 2902 genes always zero = out-of-vocabulary in scGPT training
- 2039 genes always nonzero = within-vocabulary
- Test: are nonzero genes enriched for housekeeping genes? immune-relevant genes? high-expression?
- Could use GO annotation or HGNC gene family info

### H3 (NEW family — cross_model_alignment): Compare ID compression profile with Geneformer
- Geneformer should show similar or different compression if architecture/training drives it
- If both show ~L0→L11 compression of similar magnitude, it's a general LM phenomenon
- If different magnitudes, model-specific mechanism

## Retired directions
- SVD eff_rank as stand-alone claim (invalidated by zero-norm artifact; needs restatement)
- Lineage centroid orthogonality (H03, iter_0046: negative, z<0.3)
- Graph topology kNN purity (iter_0045 H01: negative)

## Active validated findings (to build on)
- TwoNN ID compression: L0=34.7 → L11=20.0 (ratio=0.58), cross-seed consistent, nonzero-only validated
- ID compression is linear across layers (R²=0.986), no phase transition
- Zero-norm structure: constant 2902 genes absent throughout all layers
