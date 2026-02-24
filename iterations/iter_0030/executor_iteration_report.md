# Executor Iteration Report — iter_0030

## Summary

Three hypotheses tested. Most significant outcome: **the 14 "diverging" genes (FOS, JUNB, TNF, HLA-A, etc.) are OOV (out-of-vocabulary) tokens in scGPT** with zero input embeddings. All previous "bifurcation" findings (iter_0028–0029) reflect an OOV artifact, not learned geometric structure. This resolves the structural mystery and pivots the research to the 195 in-vocabulary named genes.

---

## Command Trace

```bash
# Setup
mkdir -p /Volumes/Crucial\ X6/.../iterations/iter_0030

# H01 + H02 (first pass, H03 failed on STRING key format):
conda run -n subproject40-topology python run_iter0030_screen.py

# H03 standalone (corrected STRING format):
conda run -n subproject40-topology python run_h03.py

# H04 (OOV discovery — inline quick check):
conda run -n subproject40-topology python3 -c "... OOV analysis ..."
```

---

## H01: GO + TRRUST Enrichment of 14 Diverging Genes

**Family:** module_structure
**Method:** Mann-Whitney slope test (div vs conv 209-gene trajectories); mygene API GO BP query for the 14 diverging genes; TRRUST annotation (cache not available).
**Artifact:** `h01_diverging_annotation.json`

**Results:**
- Slope divergence confirmed: div mean slope = +0.9790, conv mean slope = −0.9542
- Mann-Whitney AUROC (div > conv) = **1.000**, p = 2.34e−11 (perfect separation)
- Top GO BP terms for 14 diverging genes (from mygene API):
  - "positive regulation of transcription by RNA polymerase II" (20 hits)
  - "regulation of DNA-templated transcription" (8 hits)
  - **"positive regulation of canonical NF-kappaB signal transduction" (8 hits)**
  - **"immune response" (7 hits)**
  - **"inflammatory response" (4 hits)**
  - **"extrinsic apoptotic signaling pathway via death domain receptors" (4 hits)**
- Unique GO BP terms across 14 genes: 379
- TRRUST cache unavailable (separate annotation not fetched)

**Decision:** positive — diverging genes are enriched for immune/stress/NF-κB terms. **BUT see H04 below for critical reinterpretation.**

---

## H02: Layer-of-Bifurcation Anatomy

**Family:** topology_stability
**Method:** Compute inter-group distance (centroid-centroid) and intra-group variance for diverging (n=14) vs converging (n=195) at each of 12 layers. Compute separation ratio = inter / mean(intra).
**Artifact:** `h02_bifurcation_anatomy.json`

**Results:**
| Layer | Inter-group dist | Intra-conv | Ratio |
|-------|-----------------|------------|-------|
| L00 | 17.738 | 11.339 | 3.129 |
| L05 | 19.785 | 8.004 | 4.944 |
| L10 | 24.075 | 5.165 | 9.323 |
| L11 | 21.752 | 3.803 | 11.438 |

- **Intra-diverging group = 0.000 at ALL 12 layers** — all 14 genes are IDENTICAL vectors at every layer
- Spearman rho(layer, ratio) = **1.000**, p = 0.0 (perfect monotonic increase in separation)
- Spearman rho(layer, inter-group distance) = +0.979, p = 3.09e−8

**Critical interpretation:** Intra-div = 0 is not rounding error — all 14 genes have exactly equal embeddings at every layer (confirmed by pairwise L2 check: max dist = 0.000000 at all layers). This triggers the OOV investigation in H04.

**Decision:** positive (finding is real), but **reinterpreted as OOV artifact** — not learned geometry.

---

## H03: Trajectory Slope vs STRING Degree (Continuous)

**Family:** graph_topology
**Method:** Per-gene trajectory slope = Spearman rho(layer, dist_to_centroid). Correlate with STRING degree and weighted degree across 209 named genes.
**Artifact:** `h03_slope_vs_degree.json`

**Results:**
- Named genes with STRING edge: 201/209
- rho(degree, slope) all: **+0.018**, p = 0.793 — no correlation
- rho(weighted-degree, slope) all: +0.012, p = 0.863 — no correlation
- High-degree (Q75) slope mean: −0.803; Low-degree (Q25) slope mean: −0.882 — MW AUROC = 0.525, p = 0.683
- Diverging genes degree mean: **37.14**, converging: **29.05** — div genes are NOT degree-depleted (AUROC = 0.575, p = 0.825)

**Decision:** negative — STRING degree does not predict trajectory slope. No structure-function gradient.

---

## H04: OOV Gene Discovery (Bonus — decisive)

**Family:** null_sensitivity
**Method:** Check L0 embedding norm for all 209 named genes; identify zero-norm genes.
**Artifact:** `h04_oov_analysis.json`

**Results:**
- **14/209 named genes have zero L0 embedding vectors**: FOS, HLA-A, HLA-DPB1, JUNB, KLF6, LDHA, LGALS1, NCAM1, NCOA3, NR4A3, PAX5, PTGS2, TBXAS1, TNF
- These 14 are EXACTLY the "diverging" cluster from iter_0029
- All 14 have unique vocabulary indices (no index collision) — the zero embedding is genuine OOV
- Non-zero genes: L0 norm mean = 19.68, max = 21.90

**Implication:** All prior "bifurcation" findings (iter_0026 H01, iter_0028 H03, iter_0029 H01/H02/H03) reflect scGPT's zero-embedding of OOV genes, not learned biological geometry. The transformer passes all OOV tokens through identically, producing a collapsed single point that diverges from the centroid of in-vocabulary genes.

**Decision:** decisive negative/reinterpretation — retires the "diverging cluster" narrative as an artifact.

---

## Quantitative Summary

| Hypothesis | Family | Primary Metric | Value | Direction | Decision |
|---|---|---|---|---|---|
| H01 | module_structure | AUROC(div>conv slope) | 1.000 | positive | positive (artifact confirmed) |
| H02 | topology_stability | rho(layer, ratio) | 1.000 | positive | reinterpreted as OOV |
| H03 | graph_topology | rho(degree, slope) | +0.018, p=0.79 | negative | negative |
| H04 | null_sensitivity | n_zero_L0 genes | 14/209 = 6.7% | decisive | retires bifurcation narrative |

---

## Next Priority

Rerun core geometric analyses (kNN graph topology, spectral gap, manifold distance, STRING-distance correlation) on the **195 in-vocabulary named genes only** (excluding the 14 OOV genes). This is the correct analysis for testing learned geometry in scGPT.
