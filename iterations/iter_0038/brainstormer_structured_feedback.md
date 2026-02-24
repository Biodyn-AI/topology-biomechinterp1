# Brainstormer Structured Feedback — iter_0038

## Research Gate: PASSED

---

## Assessment of This Iteration

### H01: GC-TF vs Plasma-TF Centroid Proximity Divergence
**Grade: Strong positive.** The plasma-TF z-score trajectory (−0.62 → +0.55) is the cleanest quantitative result in the B-cell geometry line so far. GC-TFs at stable z≈−1.5 while plasma-TFs cross null at L8 and become *farther than null* by L11 is a genuine layer-resolved geometric signal. This finding is ready for consolidation and cross-model testing.

**Key risk**: Plasma panel is only 2 genes (IRF4, PRDM1); PRDM1 is also in B-cell centroid panel, creating partial overlap. This contaminates the distance computation. Must flag this artifact in reporting and re-test with plasma-exclusive genes when extended vocab allows.

### H02: Directional Drift
**Grade: Mixed, but geometrically interesting.** The 94° GC-plasma angle at L11 is a publishable geometric signature — it says the two differentiation programs are encoded in near-perpendicular directions in the final layer. The drift magnitude (26.4 units) with zero alignment to plasma axis is puzzling but informative: the B-cell centroid moves strongly, but NOT toward its differentiated successor. This demands investigation of what attractor the drift IS heading toward.

**Executable follow-up**: Find the gene closest to the drift endpoint (B-cell centroid at L11 + displacement vector from L0). This will identify what semantic direction the embedding is being pulled in across layers.

### H03: Blocked
No data. Must switch to cycle2_maxgenes1024 for PAX5/EBF1/BCL6 screen. This is a 0-cost vocab change that unblocks multiple hypotheses.

---

## Most Productive Open Threads

1. **GC-plasma orthogonality at L11** — confirm with principal angle (SVD), extend to Geneformer
2. **Plasma-TF divergence z-trajectory** — replicate with larger plasma set (cycle2 vocab: XBP1, MZB1, JCHAIN)
3. **B-cell centroid drift target** — identify what gene/direction the massive 26.4-unit drift points toward
4. **Cross-model** — Geneformer precision@10 and GC-TF proximity: required to generalize claims

---

## Directions to Retire

- `NK/myeloid specificity screen` (current form): retire_now. Repeatedly blocked by vocab; no in-vocab representation. Only rescue: use a completely different embedding (cycle2 or Geneformer with broad vocab).
- `directional drift cosine alignment` (current formulation): deprioritize standalone. Drift does not align with any tested direction; it is more useful as a *background characterization* than as a primary hypothesis.
