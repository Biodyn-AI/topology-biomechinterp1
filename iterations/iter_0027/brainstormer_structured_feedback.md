# Brainstormer Structured Feedback — iter_0027

**Date:** 2026-02-22

---

## Research Gate Status

`passed_min_research_gate: true` — 3 hypotheses tested, all artifacts present.

---

## Assessment of iter_0027 Results

### H01: PC1 Biological Identity at L11 — NEGATIVE
PC1 captures 67.6% of variance at L11 but does not separate any simple binary biological feature (TF/non-TF, HLA, AP1, hub degree). **However:** the gene extremes are substantively interesting — JUN at top, FOS at bottom (same heterodimer partner!); LCK at top, HLA-A/HLA-DPB1 at bottom. This is not a random pattern. The JUN/LCK vs FOS/HLA polarity looks like a T-cell effector vs antigen-presentation axis, or alternatively an expression-level axis in a particular cell type context. The binary tests were the wrong test for this question; continuous annotation is needed.

**Key insight to carry forward:** PC1 is not null — it's just not separable by simple binary labels. Test it with continuous scores (expression specificity, cell-type marker weights, single-cell variance).

### H02: TRRUST Activation/Repression Polarity — NEGATIVE
Clean null across all layers and distance metrics. The exception — ETS1 delta = +10.85 — has only 2 repression targets and is not actionable. Direction of regulatory effect is not encoded in the population-level geometry of these embeddings.

**Retired direction:** Population-level polarity via PC1 or L2 distance. ETS1-specific within-TF analysis needs ≥5 repression targets to revisit.

### H03: Dorothea Threshold Sensitivity — CRITICAL METHODOLOGICAL FINDING
**This is the most important result in iter_0027.** All AUROC values < 0.5 with population-level null directly contradict iter_0026 AUROC=0.658. The explanation is confirmed: prior positive results used a same-TF-gene-pool null, which controls for hub centrality. With population-level null, regulatory pairs are *farther* than random — i.e., hub TFs pull their target genes toward the embedding center, and so TF-target pairs look closer only because they're selected from a biased subset.

**Implication:** The entire manifold_distance family of results (iter_0022–0026) is now under suspicion. They likely measure "hub genes cluster centrally" rather than "regulatory relationship = geometric proximity." This is still a finding — hub centrality IS a real signal — but it's not the specificity claim the prior results implied.

**What survives this critique:**
- Dimensional collapse (PR profile across layers) — purely geometric, null-independent.
- Gene family clustering (within-family vs cross-family null) — uses family membership as both signal and null basis; less affected by hub centrality IF family members don't share disproportionate hub status.
- PH-based results (if any) using within-subset filtration.

**What is compromised:**
- All Dorothea/TRRUST regulatory proximity AUROC from iter_0022–0026 with same-gene-pool null.
- STRING Spearman ρ if computed with same-pool null.

---

## Stale Direction Audit

| Direction | Recommended Action | Reason |
|-----------|-------------------|--------|
| Dorothea regulatory proximity (same-pool null) | **RETIRE** | H03 confirmed null-sensitivity. Misleading specificity claim. |
| TRRUST activation/repression polarity (population-level) | **RETIRE** | H02 confirmed null across all layers and metrics. |
| PC1 identity via binary features (TF/HLA/AP1) | **RETIRE** | Three binary tests all AUROC ~0.43–0.49. No rescue with same approach. |
| STRING Spearman ρ without matched null | **DEPRIORITIZE** | Likely same hub confound. Must retest with degree-matched null before reporting. |
| Hub-degree/connectivity as predictor (per se) | **RESCUE** | Hub centrality is *the signal* — characterize it properly instead of treating it as confound. |

---

## Key Productive Threads from iter_0027

1. **Hub centrality as a first-class finding:** Hub TFs occupy central positions in embedding space. This is mechanistically interesting (scGPT probably saw these genes co-occur across many cell types in training data). Characterize it directly rather than controlling it away.

2. **PC1 continuous annotation:** The JUN-top / FOS-bottom / LCK-top / HLA-bottom pattern at L11 is biologically coherent at the cell-type level. A T-cell activation axis interpretation is testable with publicly available PBMC cell-type marker gene lists.

3. **Gene family clustering validity check:** Need to confirm that prior immune family AUROC results are not driven by hub-degree concentration within families. HLA family members are not hubs; AP1 members (JUN, FOS) are hubs. This is a concrete confound to dissect.

4. **Matched-degree null design:** The correct approach for regulatory proximity going forward is to match null pairs by (gene degree, embedding centrality) rather than drawing from all genes or all TF-target genes. This turns an embarrassment into a proper methods contribution.
