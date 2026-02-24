# Brainstormer Structured Feedback — iter_0005

## Gate Status
`passed_min_research_gate: true` — two promising results, one inconclusive.

---

## What worked

**H02 (residual drift + GO enrichment)** is the most significant finding to date. It establishes the first functional biological anchor: high-drift genes (transcription regulation + immune) vs low-drift genes (metabolic/structural). OR=3.58 for GO:0006357 (p=0.004) is noteworthy even at n=209. The biology is coherent with scGPT's training objective. This is a direction worth doubling down on.

**H03 (effective rank)** independently confirms the progressive dimensionality compression seen by TwoNN in iter_0004. Pearson r=0.793 is strong. The 6× ER compression (7.89→1.28) vs 1.2× TwoNN compression is mechanistically interesting: ER is collapsing into a near-rank-1 subspace while the manifold geometry (TwoNN) barely moves. This dissociation is a real finding worth characterizing.

**H01 (GO clustering)** is directionally correct (all 12 layers z < 0) but noise-limited at 209 genes. The signal is real but needs power — not a failure, just a sample-size problem. Retire the 209-gene version, promote the 4803-gene version.

---

## What is stale / needs retirement

- **TRRUST co-target clustering**: inconclusive x2 across iters 0004–0005. Both attempts at small gene pool. Retire until full vocab is available, then try once more.
- **Cross-model scGPT/Geneformer alignment**: infeasible without Geneformer residual embeddings. Full retirement.
- **Rewiring-null PH survival**: uniformly negative across 4 iterations (0005–0008 in master log). No improvement with constrained/quantile variants. Retire this null variant.
- **Distance-permutation null PH**: uniformly negative, over-adversarial. Retired.
- **Cross-layer linear CKA** (as a primary question): answered — CKA≈1.0 across all layer pairs. Demote to supporting evidence only.

---

## Key dissociation to exploit next

The ER / TwoNN dissociation is the most mechanistically interesting signal not yet exploited:
- ER: 7.89 → 1.28 (6× collapse, monotone)
- TwoNN: 10.33 → 8.52 (1.2× change, irregular)

Interpretation: scGPT is collapsing the gene representation into a rank-1-like spectral subspace by layer 11, while the manifold geometry (local neighborhood structure) stays roughly constant. This means the dominant direction is absorbing energy at the expense of minor directions, but local topology is preserved. The dominant singular vector at layer 11 likely encodes a biologically meaningful axis.

---

## Critical data dependency for iter_0006

The single biggest unlock is recovering the **full 4803-gene vocabulary** from scGPT. This unblocks:
- GO clustering at full power
- Residual drift enrichment at FDR level
- STRING analysis
- Betti curves on GO-defined gene subsets
- Singular vector biology

Search path: `subproject_38/.../scgpt_vocab.json` or the scGPT model checkpoint directory for `vocab.json`.

---

## Confidence in current findings

| Finding | Confidence | Notes |
|---------|-----------|-------|
| kNN CC elevation | High | z>9 in all 36 tests, layer-invariant |
| PH H1 positive (feature-shuffle) | Medium | Fails stronger nulls; may reflect kNN geometry artifact |
| TwoNN + ER progressive ID reduction | High | Two independent methods, r=0.79 |
| Residual drift + GO enrichment | Medium | Coherent biology, needs FDR replication |
| CKA residual stream stability | High | Near-identity updates confirmed |
| GO term clustering direction | Low-Medium | Directional but underpowered |
