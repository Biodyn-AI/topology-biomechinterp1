# Executor Next Steps — iter_0002 → iter_0003

## Portfolio for iter_0003 (Broad Screening)

Target 2–3 new hypotheses, prioritizing biological grounding and cross-model validation.

---

## Priority 1: TRRUST Biological Anchoring (NEW FAMILY)

**Hypothesis ID:** N02_trrust_geodesic_coregulation

**Rationale:**
- iter_0002 confirmed robust CC structure in kNN embeddings
- **Next question**: Is this structure *biologically* meaningful?
- TRRUST co-regulatory pairs should be closer in embedding space than random pairs
- If geodesic distance correlates with TRRUST co-regulation, CC signal is anchored to real biology

**Method:**
1. Download TRRUST human transcription factor network (co-regulation edges)
2. Measure embedding geodesic distance between all TRRUST pairs (top 500–1000 edges by confidence)
3. Measure geodesic distance for 1000 random gene pairs
4. Statistical test: Wilcoxon rank-sum (TRRUST vs random)
5. Compute Spearman ρ(TRRUST_confidence, 1/geodesic_distance)

**Expected outcome:**
- Positive: TRRUST pairs cluster in embedding space (validates CC findings)
- Negative: No distance difference (CC is spurious to biology)

**Effort:** 2–3 hours; data already available

---

## Priority 2: Geneformer Cross-Validation (BREADTH)

**Hypothesis ID:** N03_geneformer_knn_cc

**Rationale:**
- If CC elevation is model-agnostic, it indicates genuine residual space property
- Geneformer is orthogonal architecture to scGPT
- Expect similar effect size if signal is robust

**Method:**
1. Locate Geneformer residual embeddings from subproject_38 (if available)
2. Run identical kNN CC screening (k=10, PCA-20, feature-shuffle null)
3. Compare z-score distributions vs scGPT

**Expected outcome:**
- Positive: CC elevation similar magnitude and pattern as scGPT
- Negative/Inconclusive: CC is architecture-specific artifact

**Effort:** 1 hour (reuse script)

---

## Priority 3: H0/H2 Persistence Homology Fast Screen (CHEAP BREADTH)

**Hypothesis ID:** N04_betti_profile_stratification

**Rationale:**
- iter_0002 tested graph metrics (0-dimensional CC/TR)
- Can quickly extract Betti-0 (connected components) and Betti-2 (voids/cycles) from existing Ripser outputs
- Tests complementary topological fingerprint

**Method:**
1. Check for existing Ripser persistence diagrams in iter_0001 or subproject_38 cache
2. If available: stratify Betti counts by layer and seed
3. Compute H0 lifecycle (mean edge where components merge) and H2 lifetime (cycle birth–death)
4. Compare to feature-shuffle null (generate 5 shuffled replicates, compute Ripser)

**Expected outcome:**
- Positive: H0/H2 Betti profiles significantly differ from null
- Inconclusive: Ripser data unavailable (fallback to skip)

**Effort:** 30 min–1 hour (if data exists)

---

## Retired/Archived Directions

- ❌ **Transitivity metrics**: H02 showed negative signal; transitivity is lower in real embeddings
- ❌ **Degree-preserving rewiring null**: Historical (iter_0006–0008) showed 0/24 significance; not rescuing

---

## Schedule for iter_0003

1. **N02 (TRRUST)**: Execute first — highest biological stakes
2. **N03 (Geneformer CC)**: Execute parallel if time permits (independent data)
3. **N04 (H0/H2)**: Light execution; add to results if quick, else defer to iter_0004

**Success gate**: At least one of N02/N03 executed with machine-readable results and decision (promising/negative/inconclusive)

---

## Cumulative Hypothesis Portfolio

| ID | Family | Name | Status | Decision |
|----|----|----|----|-----------|
| H01 | graph_topology | CC vs feature-shuffle | tested | promising |
| H02 | graph_topology | Transitivity vs feature-shuffle | tested | negative |
| N02 | biology_grounding | TRRUST coregulation vs geodesic | pending | — |
| N03 | cross_model | Geneformer kNN CC | pending | — |
| N04 | homology | H0/H2 Betti profile | pending | — |

---

## Brainstormer Guidance Addressed

✅ **N01 (Recover iter_0001 gate)**: Executed graph topology screen → H01 promising, H02 negative
✅ **New hypothesis portfolio**: N02–N04 designed to test biological anchoring, model robustness, and higher-order topology

**Next frontier**: Biological grounding (N02 TRRUST) is the highest-value next move, as confirmed by brainstormer feedback.

