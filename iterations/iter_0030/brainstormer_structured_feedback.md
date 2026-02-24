# Brainstormer Structured Feedback — iter_0030

## Overall Assessment

**Gate:** PASSED. iter_0030 delivered a decisive methodological finding that resets the research direction cleanly.

The OOV discovery (H04) is the most important result in the last 5 iterations. It retroactively invalidates the entire "bifurcation" narrative (iter_0026–0029) and provides a hard filter: all future analyses must use the 195 in-vocabulary genes exclusively.

---

## Per-Hypothesis Evaluation

### H01 (GO enrichment of 14 diverging genes)
- Result is real and correctly reported (NF-κB, immune, apoptosis enrichment).
- The biological interpretation is now: scGPT's OOV token set happens to be enriched for immune/stress/context-dependent genes. This is an interesting model-vocabulary observation but not a geometry finding.
- **Status:** Published as methodology note. Retire as a research direction.

### H02 (Layer bifurcation anatomy)
- rho=1.000 for layer vs inter/intra ratio is not a discovery — it measures the trivial property that a zero-vector set (OOV tokens) diverges from a non-zero set more strongly as the non-zero set converges. No learned geometry involved.
- intra-div = 0.000 at ALL layers was the critical tell that should have stopped analysis immediately; executor correctly identified this.
- **Status:** Retire. The finding is: all 14 OOV tokens are processed identically by scGPT (as expected). Not publishable as geometry.

### H03 (Trajectory slope vs STRING degree)
- rho=0.018, p=0.793 — genuinely null. However, the slope calculation was run on all 209 genes including 14 OOV genes whose slope is an artifact. The null result should be verified on 195 in-vocab genes only before final retirement.
- **Status:** Rescue-once: rerun on 195 in-vocab genes. If still null, retire permanently.

### H04 (OOV discovery)
- Decisive. The most valuable finding of the iteration.
- Identifies 14 specific genes (FOS, HLA-A, HLA-DPB1, JUNB, KLF6, LDHA, LGALS1, NCAM1, NCOA3, NR4A3, PAX5, PTGS2, TBXAS1, TNF) as scGPT OOV.
- Establishes that the correct analysis universe is 195 genes.
- **Status:** Include in paper as methodology finding. Retire as active research direction (it's resolved).

---

## Strategic Pivot

The OOV discovery creates a clean opportunity: we now know the valid dataset (195 in-vocab genes) and can run all core geometric hypotheses with correct inputs. Previous results from the STRING/spectral/kNN/manifold analyses on 209 genes are suspect and need clean replications on 195 genes.

Priority order for iter_0031:
1. Establish clean geometric baseline on 195 in-vocab genes (kNN topology, spectral gap, centroid trajectories)
2. Test STRING → embedding distance correlation on 195 genes (the core structure-function hypothesis)
3. Probe manifold geometry (intrinsic dimension, PCA per layer, curvature) for the 195-gene set

---

## What to Not Repeat

- Do not run any analysis that includes the 14 OOV genes unless the purpose is explicitly to characterize OOV behavior.
- Do not retest STRING degree vs slope until the clean 195-gene slope data is established.
- Do not run more cluster stability or bifurcation anatomy tests — that narrative is fully closed.
