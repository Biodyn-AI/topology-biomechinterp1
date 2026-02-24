# Next Steps: iter_0051

## Key findings to build on

1. **H01 critical revision**: SV2-4 encodes co-expression proximity, not regulatory proximity specifically. Co-expression pairs (mean cosine=0.885) are closer than TRRUST pairs at 12/12 layers. This means iter_0049's strongest claim needs reframing.

2. **H02 SV1 sign-flip**: SV1 loading anti-correlates with embedding norm at L0-L2 (rho≈-0.33), goes neutral at L3-L8, then reverses to positive at L9-L11 (rho≈+0.22). This tracks the geometric rotation of SV1 documented in iter_0049.

3. **H03 spectral hierarchy**: SV2-4 > SV5-7 > SV8-10 >> SV1 for TF enrichment. SV5-7 is a new positive.

---

## Priority hypotheses for iter_0051

### H01_next (HIGH PRIORITY): Are co-expression pairs TF-biased?
- **Question**: Is the SV2-4 co-expression signal biologically meaningful because co-expressed genes include co-regulators?
- **Method**: Among top-300 co-expression proxy pairs, what fraction involve a TF? Compare TF-TF pairs vs non-TF pairs in SV2-4.
- **Why**: Determines whether the co-expression confound fully explains the regulatory signal.

### H02_next (MEDIUM): SV1 high-loading gene identity (early vs late layers)
- **Question**: Which genes have high |SV1| loading at L0-L2 (anti-expr) vs L9-L11 (pro-expr)?
- **Method**: Rank genes by SV1 loading at L1 and L11. Top 50 each. Check GO/TRRUST membership.
- **Why**: Mechanistic understanding of SV1 rotation — what genes "flip" in SV1 importance.

### H03_next (NOVEL): SV5-7 residual TF enrichment after SV2-4 projection
- **Question**: Is SV5-7 TF enrichment independent of SV2-4, or is it a downstream effect?
- **Method**: Project each gene's embedding onto SV2-4 subspace, subtract, then SVD residual → SV5-7 of residual. Test TF enrichment.
- **Why**: If SV5-7 retains TF enrichment in residuals, it's a genuinely independent biological axis.

### H04_next (NOVEL): Housekeeping gene membership of SV1
- **Question**: Do housekeeping genes (high constitutive expression) cluster at either extreme of SV1?
- **Method**: Use published housekeeping gene list (e.g. Eisenberg & Levanon 2013, ~3804 genes). Map to embedding vocab. Test Mann-Whitney: housekeeping vs non-housekeeping SV1 loading at L0-L2 and L9-L11.
- **Why**: Directly tests whether SV1 early-layer anti-expr signal = rare/specific genes vs late-layer pro-expr = housekeeping.

## Retired directions
- manifold_distance: 3+ negative iterations, no rescue rationale
- Pure graph_topology (kNN purity): multiple negatives

## Updated paper claims
- SV2-4 encodes co-expression proximity (stronger claim than regulatory proximity specifically)
- SV5-7 shows secondary TF enrichment (new finding)
- SV1 semantic sign-flip at L8-9 (mechanistic)
