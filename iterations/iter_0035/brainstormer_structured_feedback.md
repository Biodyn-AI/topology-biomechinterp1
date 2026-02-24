# Brainstormer Structured Feedback — iter_0035

## Gate Status
`passed_min_research_gate: true` — 3 hypotheses tested, 2 positive, 1 definitive null.

---

## Outcome Assessment

### H01 (TRRUST + GO proximity) — NEGATIVE, RETIRE NOW
This is the second failure for annotation→proximity (STRING iter_0031, TRRUST/GO now). AUROC=0.506 for co-regulated pairs and 0.473 for GO pairs at L11. The GO trend across layers (rho=+0.720) is a noise artifact, not actionable — GO AUROC never clears 0.5. The full functional-annotation-to-geometric-proximity hypothesis is dead in this embedding space. Any future annotation-based test must use a fundamentally different framing (e.g., functional annotation of the *neighborhood* rather than pairwise co-annotation as a *predictor* of distance).

### H02 (B-cell kNN precision@10) — STRONG POSITIVE
7x enrichment at L2 (z=5.20, emp_p=0.000), sustained through L11 (z=2.25, emp_p=0.044). This is a clean, interpretable result: B-cell markers reliably appear in each other's top-10 nearest neighbors at all tested layers. The early-layer peak (L2) is mechanistically interesting. This is now the strongest operational measure of B-cell geometric coherence — stronger than PC1 because it is not susceptible to singular vector orientation artifacts.

### H03 (B-cell centroid separation, normalized) — POSITIVE WITH CAVEAT
Normalized rho=+0.972 (p<0.0001) confirms monotonic B-cell amplification. The bootstrap null caveat (raw rho z=0.82, p=0.088) is important: raw centroid distance is not distinguishable from null. Only the normalized metric clears significance. This means the claim must be framed carefully: "B-cell separates at an increasing fraction of the overall manifold scale" — not "B-cell centroid moves further away in absolute terms."

---

## Stale Direction Audit

| Direction | Status | Iterations negative | Decision |
|-----------|--------|---------------------|----------|
| STRING PPI → embedding distance | Confirmed null (iter_0031) | 2 | retire_now |
| TRRUST/GO co-annotation → proximity | Confirmed null (iter_0035 H01) | 2 | retire_now |
| PC2/PC3 T-cell/myeloid axes (small marker set) | Inconclusive (n<14 markers) | 1 | rescue_once_with_major_change — expand to ≥20 markers |
| Monotonic community emergence (Spearman rho) | Null (iter_0034 H01) | 1 | retire_now — subsumed by kNN precision |
| Chromosomal proximity | Null (iter_0025 H03) | 1 | retire_now (control use only) |
| GO enrichment at SV2 poles | Null (iter_0011 H03) | 1 | retire_now |
| Hub-uncorrected Dorothea proximity | Null after population null (iter_0027) | 1 | retire_now |

---

## Confirmed Positive Signals (Cumulative)

1. **kNN clustering coefficient** elevated vs feature-shuffle null — all 12 layers, all seeds (iter_0002)
2. **H1 persistent homology** elevated vs feature-shuffle null — 11/12 layers (iter_0003)
3. **TwoNN intrinsic dimensionality** lower than noise — mild compression across depth (iter_0004)
4. **Cross-layer CKA ~1.0** — residual stream near-identity across all 12 layers (iter_0004)
5. **B-cell PC1 pole** — strongest AUROC at L11 (iter_0029+, multiple replications)
6. **B-cell kNN precision@10** — 7x enrichment at L2, persists L11 (iter_0035 H02)
7. **B-cell normalized centroid separation** — rho=+0.972 monotonic increase (iter_0035 H03)
8. **Dorothea activation proximity** — AUROC ~0.60 at early layers (OOV-corrected, iter_0032)

---

## What is Missing to Build a Publishable B-cell Geometry Claim

1. Cross-model validation (Geneformer) — not yet done
2. Expanded marker set kNN precision (5 → 15-20 markers)
3. Characterization of B-cell geometric neighborhood (what non-B-cell genes are near B-cell markers?)
4. Mechanistic probe: what changes at L2 (the peak enrichment layer)?
5. Intrinsic dimensionality of B-cell sub-manifold vs full manifold
6. Negative control validation: does the B-cell signal vanish with random gene-name permutation?
