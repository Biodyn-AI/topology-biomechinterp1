# Brainstormer Structured Feedback — iter_0007

## Gate Status
`passed_min_research_gate: true`

## What Worked

**Core positive confirmed:** SV1 extracellular/cytosol axis at layer 11 passes gene-label shuffle null (empirical p=0.004). OR=6.37 for extracellular space (GO:0005615), OR=2.96 for cytosol complement. This is the first finding in this project to survive an appropriate null control, making it the anchor result.

**Layer-wise SVD trajectory:** Technically strong result. SV1/SV2 ratio 4.07→7.70 across 12 layers; extracellular axis present (p<0.006) at all 12 layers; effective rank drops 19.54→1.63 (12x). The monotonic compression with SV1 biological axis strengthening is novel and concrete. Transient deviations (L3: mitochondrion, L7-8: ER lumen, L5: GPCR) are potentially very interesting — they suggest layer-specific biological encoding that gets overwritten by the dominant secretory axis.

**Null methodology resolved:** Identifying that feature-column shuffle is degenerate for L2-norm statistics is a genuine methodological contribution. Gene-label shuffle is now the validated null approach.

## What Failed

**Drift TF enrichment:** Nominally positive (OR=2.75, p=0.00198) but does not survive gene-label shuffle null (empirical p=0.124). Drop from primary evidence. It may still be a real biological trend at weaker effect size, but cannot be claimed as a confirmed finding without a cleaner signal.

## What Is Underexplored

1. **SV1 bottom-pole cytosol enrichment** (p=0.010, OR=2.96) has not been tested against a gene-label shuffle null. Quick validation needed.
2. **SV2 and SV3** have never been analyzed — only SV1 examined throughout this project.
3. **Transient layer deviations** (L3 mitochondrion, L7-8 ER lumen) are unexplained and potentially the most biologically interesting observation so far.
4. **STRING/PPI network geometry connection** — the SV1 axis likely reflects secretory pathway biology. Has never been tested against known interaction networks.
5. **Annotation completeness confounder** — well-annotated genes may be preferentially represented in GO enrichments regardless of geometry. Not yet controlled for.
6. **Within-layer structure of the transients** — do L3 genes actually form a spatially distinct cluster, or is it just a coincidental quartile-level enrichment?

## Confidence in Current Claims

| Claim | Confidence | Basis |
|-------|-----------|-------|
| SV1 top-pole = extracellular space | HIGH | gene-label shuffle p=0.004, OR=6.37 |
| SV1 bottom-pole = cytosol | MEDIUM | p=0.010, OR=2.96, not yet null-validated |
| SV1 axis traverses all 12 layers | MEDIUM | p<0.006 at all layers, no per-layer null |
| Effective rank compression 12x | HIGH | deterministic computation, no null needed |
| Drift TF enrichment | LOW | fails gene-label shuffle null (p=0.124) |
| Layer 3 transient = mitochondrion | LOW | single-layer observation, no null, may be noise |

## Strategic Assessment

The project has converged on a clean, validated geometric finding: the dominant axis of scGPT's layer-11 representation space encodes the extracellular/cytosolic localization contrast. The layer-wise trajectory shows this axis strengthens monotonically with progressive compression. The next iteration must: (1) validate the layer-wise specificity with per-layer nulls, (2) probe SV2/SV3 for additional axes, and (3) connect the finding to known biology (STRING PPI, signal peptide presence, cell type expression). The transient deviations are a high-reward direction if they can be confirmed as real.
