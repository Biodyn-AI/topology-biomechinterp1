# Brainstormer Structured Feedback — iter_0029

## Overall Assessment

Three-for-three positive. iter_0029 is the most productive iteration to date. The combination of H01+H02 crystallizes a concrete, testable biological story: scGPT partitions 209 immune genes into two trajectory classes with near-perfect geometric separation (silhouette=0.905), and the 14-gene outlier class is identifiable and thematically coherent (stress/inflammatory response + antigen presentation). H03 closes the k-robustness question decisively.

The project now has a **core narrative**: monotonic spectral compression + bifurcated gene trajectories + stress/AP1/HLA outlier identity. The next iterations must (a) characterize *why* these 14 genes diverge, (b) test whether this is model-specific or general, and (c) probe the internal geometry of the converging 195-gene manifold.

---

## Iter_0029 Hypothesis-by-Hypothesis Assessment

### H01 — kNN Connected Components Biological Enrichment
- **Status:** Strong positive. AP1 enrichment p=0.023 in the 14-gene outlier component.
- **What's new:** Gene identity of the outlier component is now fixed and biologically interpretable.
- **Gaps remaining:** Fisher exact was only computed for 3 family categories (AP1, KLF, HLA-I). Broader GO enrichment not yet done. Cell-type-specific expression variance of the 14 genes not characterized.
- **Decision:** Continue. High rescue value in downstream characterization.

### H02 — Gene Trajectory Clustering
- **Status:** Best result so far. Silhouette=0.905 is unusually high for biological embeddings. Two-cluster solution exactly matches H01 connected components — this is not coincidental and deserves further characterization.
- **What's new:** Establishes the *temporal* nature of the separation: the 14 genes start farther from centroid (L0=16.6 vs 11.4) and continue to move away. The separation is not a deep-layer artifact but a trajectory difference apparent from early layers.
- **Gaps remaining:** Is the separation already present at L0? What is the mechanism — is it input expression pattern or model architecture? Is the trajectory sigmoidal or linear? Bootstrap stability of the 14-gene cluster membership not yet quantified.
- **Decision:** Continue aggressively.

### H03 — Spectral Gap k-Robustness
- **Status:** Confirmed and closed. k-robustness is complete; this finding needs no further testing.
- **Decision:** Mark as confirmed. Do not re-test; focus on interpreting it in the context of H01+H02.

---

## Stale / Overworked Directions

- **PC1 binary axes**: Retired (iter_0028 roadmap). Still retired.
- **TRRUST activation/repression polarity**: Retired (iter_0027). Still retired.
- **Dorothea confidence-tier pairwise distance**: Retired. Still retired.
- **Fiedler vector biological partition (H-A iter_0028)**: Superseded — the 2-component connected component finding (H01) is a cleaner and stronger result. The Fiedler vector idea is now a secondary analysis at best.
- **Spectral gap k-robustness (H-C iter_0028)**: COMPLETED by H03 iter_0029. Retire from queue.
- **STRING score gradient (H-J iter_0028)**: Still pending, cheap, worth one more attempt.

---

## What the Data Is Telling Us to Do Next

1. **Characterize the 14-gene outlier set biologically** — GO enrichment, expression variance across cell types, literature support for "stress/context-dependent" expression.
2. **Identify the layer of bifurcation** — at what layer does the trajectory slope sign change for the 14 diverging genes vs 195 converging genes? Is it abrupt (single-layer break) or gradual?
3. **Cross-model test** — do the same 14 genes diverge in Geneformer embeddings? If yes, this is a universal property of the training gene set, not scGPT-specific.
4. **Bootstrap stability** — is the 14-gene identity stable under bootstrap resampling of the distance-to-centroid matrix?
5. **Intrinsic dimension of the L11 converging manifold** — what is the effective dimensionality of the 195-gene converging cluster at L11? Is it truly near-1D or a low-D hyperplane?
