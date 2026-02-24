# Brainstormer Structured Feedback — iter_0026

**Date:** 2026-02-22

---

## Research Gate Status

`passed_min_research_gate: true` — 3 hypotheses tested, all positive. Clean iteration.

---

## Per-Hypothesis Assessment

### H01: Immune Gene Family Clustering AUROC
**Verdict: Strong positive — deepen immediately.**

Peak AUROC=0.754 is solid. The per-family breakdown is the real story:
- HLA-I (1.000) and AP1 (0.969) are near-perfect. These are structurally and functionally tight families.
- KLF (0.310) is actively anti-clustered — KLF14 is a likely outlier (liver-enriched, not immune-relevant). This is diagnostic information, not just noise.
- BCL2fam (0.537) and TNFSF (0.480) are weak/negative — functionally defined "families" with diverse regulatory contexts in the training corpus.

The AUROC-vs-layer curve peaks at L8 and declines at L11 (0.694), consistent with information concentration into the dominant PC after dimensional collapse (see H03). The flattish plateau from L0–L8 (0.729–0.754) suggests family structure is encoded early and persists, not a late-layer phenomenon.

**Actionable next step:** Per-family layer curves — does HLA-I peak differently than AP1? HLA-I is structural (protein sequence similarity) while AP1 is functional (co-regulation). If they peak at different layers, this maps protein-vs-regulatory family to the geometric hierarchy.

---

### H02: Dorothea Confidence-Tier Proximity Gradient
**Verdict: Confirmatory positive — retire as standalone, use as refinement target.**

AUROC_high_vs_null=0.658 replicates the core Dorothea signal from iter_0024 H01 (AUROC was 0.66 at L7 there too — exact match). The high-vs-low gradient (0.522) is small but consistently positive.

The missing MOR column blocks the directional polarity test. TRRUST has the A/R field and is already available at the external networks path. This is the correct pivot.

---

### H03: Intrinsic Dimensionality Progression (Participation Ratio)
**Verdict: Major finding — this is a centerpiece result for the paper.**

PR 21→1.68 with rho=-1.000 is striking. Key structural observation from the raw data:
- sv_max grows monotonically L0→L10 (78→106), then drops slightly at L11 (96.5).
- sv_2nd drops sharply at L10→L11 (34→27), creating a gap-widening phase transition.
- The top-10 PC variance fraction at L11 is 0.91 — the representation is essentially 1-2 dimensional.

This is the geometric compression / low-rank routing pattern. The question is whether this compression is uniform (all genes collapse together) or selective (different gene classes collapse at different rates). That's the immediate follow-up.

Also note: the PR elbow occurs around L7–L8 (where the PR curve has the steepest relative drop from 5.3 to 2.1 over L7→L10). This is exactly where Dorothea and GO BP signals peak. This temporal alignment between geometric compression and biological signal peak is a meaningful mechanistic observation.

---

## Cross-Hypothesis Synthesis

Three findings now align into a coherent story:
1. **Biological structure is pre-encoded by L0–L5** (GO CC, STRING peaks early-mid).
2. **Regulatory/functional identity is sharpened L5–L8** (GO BP, Dorothea, family AUROC peak at L7–L8).
3. **Dimensional collapse concentrates representation into 1–2 axes at L9–L11** (PR drops from 5.3 to 1.7).

The model appears to: (a) scaffold structural/cellular context early, (b) refine biological function mid-network, (c) compress into low-dimensional output representation at the end.

The next critical question: **what does the dominant axis (PC1 at L11) encode biologically?** This is the mechanistic keystone.

---

## What's Still Not Executed from iter_0025 Portfolio

| Hypothesis | Status | Priority |
|-----------|--------|----------|
| H-A: Geneformer replication | Not started | HIGH |
| H-B: CORUM protein complex | CORUM artifact exists (h01_corum_complex_auroc.json) but was replaced by immune family | Inspect artifact |
| H-C: Attention head attribution | Not started | HIGH-cost |
| H-D: Dorothea direction polarity | Blocked by missing MOR | → pivot to TRRUST |
| H-E: Disease gene clustering | Not started | MEDIUM |
| H-F: PH of GO CC clusters | Not started | MEDIUM |
| H-G: Intrinsic dim per bio subset | Not started | HIGH — directly extends H03 |
| H-H: Geneformer layer timeline | Not started | HIGH |
| H-I: STRING threshold sensitivity | Not started | LOW cost |
| H-J: Hub vs. peripheral density | Not started | LOW cost |
| H-K: Pathway-specific timeline | Not started | MEDIUM |
| H-L: Co-expression mediation | Not started | HIGH-cost |
| H-M: Oncogene/tumor suppressor | Not started | MEDIUM |
