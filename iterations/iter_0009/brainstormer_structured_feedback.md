# Brainstormer Structured Feedback — iter_0009

**Date:** 2026-02-22
**Gate status:** PASSED (3/3 positive)

---

## Assessment of iter_0009 Results

### H01: SV2 Bottom-Pole EV Null (emp_p=0.000/1000)
Strong. Closes the missing control from iter_0008. The gene list (MIF, APP, VIM, LGALS1, CLU, HLA-I/II) is biologically coherent and matches known exosome cargo databases. This axis is now null-controlled and biologically anchored. **Ready to extend.**

### H02: 12-Layer × 8-Compartment Scan (33/96 cells significant)
Strong systematic result. Key finding: ER lumen and secreted are persistent across all 12 layers (implying the secretory axis is baked in at embedding initialization, not learned progressively). Mito signal is transient and resurges at L10 — this is a potentially interesting layer-specific signature. Only 3/8 compartments show significant signal on SV1; EV is confirmed to be off-SV1. **This scan design should be applied to SV2 and SV3 immediately.**

### H03: SV1 Cross-Layer Rank Stability (mean r=0.929, L10-L11 break)
Strong methodological contribution. The monotonically increasing variance explained (19.1% → 76.7%) paired with the L10-L11 discontinuity (r=0.663 vs r>0.94 elsewhere) is the most structurally interesting finding this iteration. This transition point coincides with the layer where mito resurges (H02) and ER lumen OR peaks (OR=18.45 at L11). **The L10→L11 transition is a priority mechanistic target.**

---

## Directions That Need Attention

### Rewiring-Null PH Branch
Ran 4 consecutive iterations (iter_0006–0008) without a positive result under any null formulation (degree-preserving, metric-matched, edge-length-constrained, bridge-conditioned). The failure is consistent and robust. No further investment warranted unless a genuinely new null design emerges (e.g., topology-aware spectral rewiring). **Retire.**

### GO BP Cosine Clustering
Two attempts (iter_0005 H01, iter_0005 H02 partial) show z-scores consistently ~-1.1 and non-significant. Under-powered by design at 209 genes. The positive signal lives in SVD projections, not in pairwise distances. **Do not revisit without full-vocabulary expansion.**

### TRRUST Co-Target Euclidean Clustering
iter_0004 H02: mean z=+0.222, 1% significant. Inconclusive at 209-gene scale. Rescue option exists with a better null design (projection-matched rather than raw Euclidean), but lower priority given the stronger SVD signal. **Deprioritize; rescue only with SVD-projection-space version.**

---

## Open Questions from iter_0009

1. What is SV3? The third singular axis is uncharacterized. It may encode a distinct biological axis.
2. What drives the L10→L11 rank break in SV1? Is it attention-layer-specific or residual stream restructuring?
3. Does SV2's EV enrichment emerge at all 12 layers or only in late layers?
4. Is annotation density a confounder for SV1 position? (genes with more GO annotations may systematically cluster in pole)
5. Are SV1/SV2 axes orthogonal biological programs or partially correlated?
6. Does the mito transience at L2-4 then L10 reflect two distinct mito-regulatory mechanisms?
