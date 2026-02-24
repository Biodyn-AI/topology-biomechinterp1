# Brainstormer Structured Feedback — iter_0042

## Gate Status
`passed_min_research_gate: true`. All three hypotheses returned positive or mixed-positive results. Paper updated. Full forward pass.

---

## Assessment of iter_0042 Findings

### H01: TwoNN Changepoint at L3
**Quality: Good signal, weak null.**
- B-cell slope reversal at L3 (-4.832 change) is visually strong and cross-lineage contrast is compelling: B-cell compresses, T-cell expands post-L7.
- Null p=0.386 (n=12 data points, 500 permutations) is too weak to publish standalone. The permutation is operating on only 12 layer-points — sensitivity is fundamentally limited by the small sequence length.
- **Fix needed**: Use segment R² ratio as test statistic rather than slope change. Bootstrap the changepoint selection over sub-sampled gene sets. The cross-lineage divergence (B vs T vs Myeloid) is the publishable claim; the null for B-cell alone is secondary.
- Myeloid flat profile (ID ~23-26, slope change ~0.13) is the most stable null lineage — strong control.

### H02: BCL6 Metabolic Isolation Specificity
**Quality: Clean, high-impact.**
- 90x enrichment over random, 0 overlap for any B-cell master regulator. Zero-vs-nine is a crisp separation.
- The 9→4/20 trajectory from L0→L11 is unexplored: BCL6 is shedding metabolic neighbors at late layers. Where is it moving? This is a strong follow-up signal.
- GO annotation of the metabolic cluster (NAMPT, GLUL, PFKFB3, ACSL1...) is the immediate priority for biological interpretation. These genes are canonical Warburg/aerobic glycolysis regulators.

### H03: BATF/BACH2 Convergence
**Quality: Strongest result of the iteration.**
- rho=-0.972 and -0.844 with p<0.0001 and p=0.0006. Both convergence trajectories are highly significant.
- BATF–BACH2 co-convergence (18.66→4.53) is geometrically tight.
- PAX5 pre-wired vs BATF/BACH2 recruited is a coherent mechanistic narrative.
- Critical gap: PRDM1 and IRF4 are the downstream targets repressed by BCL6 during GC reaction. Their trajectories are untested. Do they anti-converge (move away from the B-cell centroid at the same rate BATF/BACH2 approach)?

---

## Emerging Coherent Model (iter_0042 synthesis)

> PAX5 encodes B-cell identity from L0 (anchor). At L3, B-cell manifold begins compressing (ID decrease) as GC-TFs (BATF, BACH2) are recruited into the PAX5 neighborhood. BCL6 is stably isolated in a metabolic stress module (NAMPT, GLUL, PFKFB3) — separate from the B-cell regulatory axis. The T-cell manifold, by contrast, expands after L7, suggesting a qualitatively different geometric regime. Myeloid is flat throughout.

This model is publishable as-is with two additions: (1) BCL6 biological annotation, (2) PRDM1/IRF4 anti-convergence test to close the repression side.

---

## Directions Assessment

| Direction | Status | Rationale |
|-----------|--------|-----------|
| B-cell manifold compression (TwoNN, ID) | **Continue** | Extend to stronger null, cross-lineage comparative |
| GC-TF convergence trajectories | **Continue / Expand** | Add PRDM1/IRF4; compute full pairwise matrix |
| BCL6 metabolic isolation | **Continue** | GO annotate; trace layer-resolved cluster membership changes |
| Cross-model alignment (scGPT vs Geneformer) | **Retire** | Blocked by vocab mismatch; no rescue path without new data |
| GC-plasma subspace angles | **Retire** | Inconclusive (iter_0040); superseded by trajectory approach |
| Chromosomal proximity | **Retire** | Negative |
| PC2/PC3 axes | **Retire** | Superseded by B-cell-specific PC1 findings |
| T-cell manifold expansion (post-L7) | **New priority** | Opposite geometry from B-cell; understudied; cheap with existing data |
| PRDM1/IRF4 anti-convergence | **New priority** | Tests repression geometry; closes GC biology narrative |
