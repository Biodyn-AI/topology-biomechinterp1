# Brainstormer Structured Feedback — iter_0004

## Research Gate: PASSED

## Signal Assessment

### H01: TwoNN Intrinsic Dimensionality
**Verdict: Genuine signal, follow up.**
ID range 8.5–10.3, decreasing across depth, mean z = -2.29 vs shuffle null. Consistent with manifold compression hypothesis. Effect is moderate but real. Note: z-scores are non-monotone (layers 7–8 are near-null), suggesting depth-dependent compression is not smooth — worth characterizing. The PCA-20 projection may be losing some ID signal at higher layers where the manifold is more curved.

Next natural move: (a) multi-seed replication to confirm trend, (b) MADA estimator as cross-check (different geometric assumption), (c) local ID per gene neighborhood to identify subpopulation structure.

### H02: TRRUST Co-Target Clustering
**Verdict: Retire this specific implementation.**
Mean z = +0.222, significant clustering only 1.0% of tests. The gene pool of 209 genes is far too small — same-pool random groups are not an independent null, they are highly correlated with the test group. Result is uninterpretable without expanding to full gene coverage. This is not a biological negative; it's a power/design issue. Rescue requires: larger gene pool (all annotated genes in scGPT vocabulary, ~4000+) or switching to GO/STRING biological annotations where coverage is higher. Do not retry TRRUST on the same 209-gene pool.

### H03: Cross-Layer Linear CKA
**Verdict: Strongly confirmed, but raises a sharper question.**
CKA ≈ 1.000 for all 12×12 layer pairs (z ≈ 174 vs null) is unambiguous. scGPT transformer blocks make near-identity linear perturbations to the residual stream. This is structurally important: the model is dominated by skip connections, not block-wise transformations.

But this uniformity actually makes the *deviations* interesting. The CKA can't distinguish which 1% of genes drive any residual variation. Two follow-ups are immediately valuable:
1. Per-gene layer-to-layer displacement (L2 norm of the change vector across layers) — identifies which genes are processed vs passed through.
2. Kernel CKA (RBF) to check if the near-identity also holds in nonlinear geometry.

## Cumulative Pattern Synthesis (iter_0001 → iter_0004)

| Family | Best result | Status |
|--------|------------|--------|
| kNN graph topology (CC) | z ≈ 27.64, all layers | Confirmed strong |
| kNN transitivity | z ≈ -34.78 | Retired (negative) |
| H1 persistence (Ripser) | Fisher p = 0.0056 | Confirmed positive |
| Intrinsic dimensionality (TwoNN) | mean z = -2.29, 6/12 sig | Promising |
| Cross-layer CKA | CKA ≈ 1.000, z ≈ 174 | Confirmed strong |
| TRRUST biological anchoring | mean z = +0.22, 1% sig | Inconclusive / underpowered |
| Cross-model feature alignment | Spearman r = 0.83, Fisher p = 0.41 | Inconclusive |

**Emerging picture:** scGPT gene embeddings sit on a structured, lower-dimensional manifold with robust local topology (elevated clustering coefficient, significant H1 loops), and this geometry is essentially layer-stable (CKA ≈ 1.0). The biological question — does this structure map to known gene function? — remains unanswered. That is the priority gap.

## Key Open Questions for Next Iteration

1. **Which genes change the most across layers?** CKA uniformity hides per-gene variability.
2. **Does biological grouping (GO, KEGG, STRING) predict embedding proximity?** TRRUST failed due to gene pool size, not biological concept.
3. **Is the manifold structure nonlinear?** Kernel CKA and geodesic distance tests are cheap.
4. **Does the ID decrease monotonically with better methodology?** TwoNN z-scores are noisy — MADA would confirm or refute.
5. **Is the H1 signal (loops) biologically interpretable?** Which gene subsets drive the H1 classes?
