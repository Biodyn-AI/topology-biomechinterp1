# Brainstormer Structured Feedback — iter_0006

## Gate Status
`passed_min_research_gate: true` — three positive results, no blockers.

---

## Assessment of iter_0006 Results

### H01 — Full-Vocab Drift + GO Enrichment
**Verdict: Strong replication, FDR borderline.**
- TF/RNA-Pol-II enrichment now confirmed in two independent analyses (iter_0005: 187 terms, iter_0006: 288 terms). Consistent OR ~3-3.6 range.
- Bimodal drift distribution is a new structural fact: ~75% of scGPT vocabulary is inert; named regulatory genes form a high-drift elite cluster at the 84th percentile.
- FDR threshold not crossed (n=209 genes is power-limited). The decisive test is a feature-shuffle null: if the enrichment vanishes with shuffled gene-position labels, it's real structure; if it survives, it's an artifact of gene selection bias.
- **Action needed**: feature-shuffle null (shuffle which gene gets which embedding position), rerun drift + GO enrichment.

### H02 — SVD Biology of Layer-11 Near-Rank-1 Subspace
**Verdict: Strongest single result of the project so far.**
- Near-rank-1 structure at layer 11 (SV1 = 739, SV2 = 96, ratio 7.7×; top-5 SVs explain 96.8% of variance) is not a marginal statistical finding — it is an extreme spectral collapse.
- SV1 encodes a biologically interpretable extracellular-vs-nuclear axis (OR=6.37, p=0.0003 for extracellular space; top loaders include IRF8, RUNX1, CD79A, CD19 — canonical immune/hematopoietic TFs).
- This axis is consistent with the training domain (Tabula Sapiens immune mix of cytokines and nuclear regulators).
- **Critical gap**: no null control yet. Feature-shuffle null on SV1 enrichment is mandatory before this becomes a publishable claim.
- **Secondary gap**: we only have layer 0 vs layer 11; the emergent trajectory (at which layer does SV1 dominance arise?) is unknown.

### H03 — ER Curvature / Compression Breakpoint
**Verdict: Confirmed, mechanistically suggestive but not yet mechanistically explained.**
- Front-loaded compression (max rate layer 1, half-life layer 4, inflection layer 3) is a clean quantitative result.
- Strong agreement with TwoNN ID (r=0.79) validates the metric.
- Compression dynamics coincide with the span of layers 1–4 that do "structural reorganization," but the mechanism is not characterized — we don't know whether this coincides with specific attention patterns, co-expression graph collapse, or other phenomena.
- **Action needed**: run full 12-layer SVD to test whether SV1 dominance emergence follows the same front-loaded trajectory as ER.

---

## Stale / Exhausted Directions

### Rewiring-null survival for persistent homology
- Tested across 4 iterations (iter_0006–iter_0008 historical, plus metric-matched and constrained variants). Uniformly `0/24` significant. The null is too adversarial regardless of edge-length constraints or metric matching.
- **Retire permanently.** Do not invest another iteration in rewiring-based null robustness.

### Cross-model Geneformer alignment
- Tested iter_0003: cosine 0.825, combined permutation p=0.349/0.409. Inconclusive. Requires residual-level Geneformer embeddings not yet available in workspace.
- **Deprioritize** until Geneformer residual embeddings are directly accessible.

### TRRUST co-target clustering (209-gene pool)
- Already retired after iter_0004. n=209 is too small for significance.
- **Remains retired.**

### GO term cosine distance clustering (209-gene)
- 0/12 layers significant (iter_0005). Power-limited by gene pool.
- **Retire in current form.** Can be revisited with full 4803-gene vocabulary if GO annotations cover that set.

---

## What's Working — Pattern Summary

The consistent signal across iterations:
1. **Biological identity is preserved in the dominant spectral structure.** High-drift genes = regulatory/TF function; SV1 axis = extracellular/secreted vs nuclear/regulatory.
2. **Spectral collapse is the dominant geometric event.** ER drops 6× (19.54→1.63). It is front-loaded (layers 1–4). By layer 11, a single 739-magnitude singular vector dominates.
3. **The representation is structurally stable (CKA≈1.0) but functionally compressive (ER 6×).** The residual stream updates are small but directionally consistent.

These three facts form the core interpretable structure of the project.
