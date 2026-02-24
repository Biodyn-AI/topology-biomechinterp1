# Brainstormer Structured Feedback: iter_0045

## Gate Status
`passed_min_research_gate: true` — 3 hypotheses tested, machine artifacts present, paper updated.

---

## Result Assessment

### H02: Full-Manifold TwoNN ID Compression — PRIMARY POSITIVE
The main finding: scGPT gene manifold intrinsic dimension L0=32.57 → L11=18.05 (−44.6%), monotone
across all 12 layers, versus Gaussian null ~123D. This is 4× below the null. Every single layer
compresses — no reversals. This is a clean, replicable, numerically large signal.

**What makes it strong:**
- Monotone across all 12 layers — no noise reversals
- Gaussian null is 4× higher; feature-shuffle null matches real (expected — shuffle doesn't
  change pairwise distances, confirming the signal is about geometric structure)
- Numerically large: 44.6% compression across 12 layers; L0 already far below Gaussian null
- Consistent with B-cell-specific ID finding in iter_0042/0043

**What it needs:**
1. Cross-seed/cross-cycle replication (cycle1, cycle4_seed43, cycle4_seed44)
2. Bootstrap CI for ID estimates (confirm n_valid=846 subsampling is stable)
3. Functional subset analysis: does B-cell submanifold compress faster than random?
4. The rate-of-compression-per-layer as a quantity to compare across contexts

### H01: kNN Lineage Purity — RETIRED (correct decision)
With n=1 for Myeloid and NK, same-lineage purity is impossible by construction (0/9 neighbors
can be same-lineage when set size=1). This is a design flaw, not a biological negative.
Do not retry without ≥4 in-vocab genes per lineage.

### H03: TCR Circuit Convergence — INCONCLUSIVE but CD28 trajectory is the real signal
Key observation: CD28 shows a V-shape (595→86 at L6→149 at L11), not monotone convergence.
LCK is constitutively embedded (rank 5-8, always). The V-shape is mechanistically interesting —
it suggests transient convergence at intermediate layers followed by divergence, potentially
reflecting a distinct computational role for CD28 (co-stimulation priming?) at L6.

**What's needed**: more circuit genes (FOS/JUN if in vocab), cross-seed verification of V-shape,
alternative centroid definitions. The single-gene CD28 story is suggestive but not publishable
without corroboration.

---

## Stale Directions to Exit

| Direction | Status | Reason |
|---|---|---|
| kNN lineage purity (n<3 per lineage) | RETIRE | By-design underpowered |
| BCL6 trajectory tests | RETIRED (iter_0044) | 3+ negatives, confirmed dead |
| STRING/GO proximity → embedding distance | RETIRED (iter_0031/0035) | AUROC≈0.5, dead |
| 5-point TwoNN isolated subsets | RETIRED (iter_0044) | Invalid method |

---

## Navigation Assessment

The ID compression finding is the strongest geometric result so far. It is:
- A claim about the overall representational architecture of scGPT
- Robust (monotone, large effect size, clear null comparison)
- Extensible (can be cross-validated, decomposed by gene subset, compared cross-model)

The next 2-3 iterations should primarily exploit this finding: replicate, validate, and decompose it
into mechanistically informative sub-questions. Secondary: rescue the CD28 V-shape signal with
more genes. Tertiary: explore new directions (Ricci curvature, cross-model alignment) that could
produce independent positive signals.
