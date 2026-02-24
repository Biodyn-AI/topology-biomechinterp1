# Brainstormer Structured Feedback: iter_0032

**Date:** 2026-02-23
**Gate status:** PASSED (passed_min_research_gate = true)

---

## Assessment of iter_0032 Results

### H01: Community Detection at L11
**Status: Promising — needs biological identity resolution**

The geometry is unambiguously real (z=34 vs null). The two communities (sizes 107/88) are stable but their biological meaning is unclear: B-cell markers trend into community 1 (p=0.09), IL2-pathway genes trend into community 0 (p=0.163), T-cell markers also into community 0 (p=0.251). The core problem is that the 10 curated gene families tested are all too small (n=3–7 in vocabulary) to reach significance for a ~50/50 community split. The communities are likely real immune functional partitions but we tested the wrong labels. This is a tractable problem.

**Next move:** Test the 2-community membership against continuous PCA scores and functional gene annotations beyond curated families (e.g., secreted vs membrane-bound, TF vs non-TF, cytokine vs receptor, STRING subgraph connectivity within community vs across).

### H02: Dorothea Regulatory Proximity Decays with Layer
**Status: Positive — important revision of prior finding**

This is the cleanest new finding in several iterations. AUROC 0.57 at L0–L2 → 0.50 at L8–L11, rho=-0.853, p=0.0004. The OOV-corrected result fundamentally revises iter_0024 (stable AUROC ~0.67 across all layers), showing that late-layer regulatory signal was an OOV artifact.

The mechanistic interpretation is now coherent: early layers preserve regulatory co-membership as geometric proximity; late layers undergo PR collapse (PR 58→9.5) which reorganizes the manifold around other organizing principles (the 2-community structure, PC1 dominance). The question is: **what does the collapsed late-layer manifold encode positively?**

**Immediate follow-on:** Test whether activation vs. repression pairs show different decay rates (repression was never tested on the OOV-corrected set). Test whether early-layer geometry (L0–L2) can predict TF-target relationships better than late layers for any other annotation databases.

### H03: PC1 Polarity (B-cell vs TF)
**Status: Inconclusive — correct direction, insufficient power**

3 B-cell markers in vocabulary is too few. The signal is directionally consistent at multiple layers. This is worth rescuing with expanded gene sets, not retiring.

---

## Revised State of Cumulative Claims (OOV-corrected)

| Claim | Status | Confidence |
|-------|--------|------------|
| kNN spectral gap monotonic decrease | Confirmed (195 genes) | High |
| PR collapse L0→L11 | Confirmed (195 genes) | High |
| 2-community structure at L11 (z=34) | Confirmed, biology unclear | Medium |
| Regulatory proximity early-layer only | New, strong (rho=-0.853) | High |
| Cell-type marker clustering (AUROC=0.851) | Needs 195-gene retest | Pending |
| STRING PPI proximity | RETIRED on 195-gene set | — |
| GO CC/BP Jaccard proximity | Needs 195-gene retest | Pending |
| TRRUST activation proximity (AUROC=0.640) | Needs 195-gene retest | Pending |
| PC1 = B-cell/TF axis | Inconclusive (n=3 B-cell) | Low |

---

## Key Open Scientific Questions Going Forward

1. **What does the late-layer manifold encode positively?** PR collapse suggests variance concentrates in PC1 (26%). If regulatory geometry is gone, and STRING is gone, what drives the 2-community split and PC1 structure?

2. **Is the early-layer regulatory signal specific to activation-type regulation?** On the full 209-gene set, TRRUST activation showed strong proximity (AUROC=0.640), repression did not. Does this asymmetry also appear in the Dorothea decay curve?

3. **Does Geneformer show an analogous early-layer regulatory geometry?** The cross-model claim (iter_0019) used static embeddings. Geneformer transformer layers may show the same decay pattern.

4. **What is the 2-community split actually encoding?** This is the most actionable open question given the strong geometry.

5. **Do GO CC/BP and cell-type marker proximity hold on the 195-gene set?** These were the three strongest prior claims and need OOV-corrected replication.
