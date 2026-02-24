# Brainstormer Structured Feedback — iter_0028

## Gate status
`passed_min_research_gate: true` — two promising results, one clean negative. Strong iteration.

---

## H01 — Hub Gene Embedding Centrality: KEEP AND EXTEND

Signal is solid: rho=+0.225 at L10, null_p=0.000/500 permutations, monotonic increase L0→L10. The five top-central genes at L11 (ETS1, IRF8, TGFBR2, REL, ITGB2) are all transcriptional hubs with broad immune regulatory roles.

**Next actions already identified by executor:**
- Bootstrap CI per layer (quantify uncertainty in the rho trajectory)
- TF enrichment in top-k most central genes vs least central
- Cross-model test (Geneformer)

**Brainstormer addition:** The centrality-connectivity correlation tested STRING degree. Next natural step is to test TRRUST regulatory target count as a second, independent "hubness" metric. If rho holds with TRRUST count at L10, the finding generalizes from PPI degree to transcriptional regulatory degree — more mechanistically interpretable.

---

## H02 — PC1 T-cell vs APC Axis: RETIRE

Both T-cell and APC markers collapse to PC1 ≈ +2.2 to +2.4 at L11. This was tested in iter_0027 (TF vs non-TF, HLA, AP1, hub degree: all negative) and again here with cell-type marker curation. Two iterations, same null result. The dimensional collapse to PC1_var=76.7% means the dominant axis is not separating any tested biological binary. Retire this direction.

**Retirement verdict:** No rescue value. The collapse-to-1D is a real geometric phenomenon but the 1D axis does not encode interpretable binary biology. The interesting question is not "what does PC1 separate?" but "what drives the collapse?" — which is covered by the hub centrality and spectral gap families.

---

## H03 — kNN Spectral Gap: KEEP AND EXTEND

rho=−0.993, p=1.3×10⁻¹⁰ for layer vs spectral gap — near-perfect monotonic signal. The identified null limitation (permuted gene labels → same gap, which is expected) is not a flaw; it confirms geometry not gene identity drives the signal. Missing: random-Gaussian baseline and identification of the two connected components.

**Critical next step:** The data already shows n_components=2 at every layer. These components are not gene-label artifacts — they emerge from the actual embedding geometry. Identifying which genes belong to each component, and what biological function they share, is the highest-value next test in this family.

---

## Cumulative pattern assessment

Across iters 0026–0028, a coherent picture is emerging:
- Deep scGPT layers show near-1D collapse (PR 1.68, PC1_var 76.7%)
- Within this collapsed space, hub genes (STRING degree) are geometrically privileged (more central)
- The kNN graph becomes more modular (decreasing spectral gap), consistent with genes clustering into tighter communities
- Immune gene families HLA-I, AP1, RUNX show strong within-family clustering (AUROC 0.97–1.0) even in collapsed space
- PC1-based binary separation tests all fail — the collapse does not encode simple categorical biology

This suggests the model is performing **continuous re-weighting of biological importance** (hub genes central, families cohesive) rather than discrete categorization. The next productive question: do the spectral components at deep layers correspond to the immune families already known to cluster strongly?

---

## Directions to retire now

| Direction | Retirement reason |
|-----------|------------------|
| PC1 biological axis (TF/non-TF, HLA, AP1, cell-type) | Tested twice across iter_0027 and iter_0028, all negative; 1D collapse makes all binary separation tests weak |
| TRRUST activation/repression polarity | Tested iter_0027 H02, negative; null construction confound confirmed |
| Dorothea confidence-tier distance | iter_0026 H02, mixed/inconclusive; repeated refinement has not improved signal |
