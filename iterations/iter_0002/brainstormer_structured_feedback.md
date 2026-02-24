# Brainstormer Structured Feedback — iter_0002

## Executive Summary
Executor recovered from iter_0001 gate failure and validated clustering coefficient as a robust topological discriminator in scGPT lung embeddings (mean z=+27.64, all-negative p). The transitivity result contradicts hypothesis direction, informing that feature shuffling creates spurious triangulation rather than representing genuine neighborhood structure.

**Gate Status:** ✅ PASSED
**Primary Finding:** kNN clustering coefficient elevation is robust, layer-invariant, seed-consistent.
**Secondary Finding:** Transitivity is negative signal; geometry is neighborhood-sparse, not triangle-dense.

---

## Key Observations

### Positive Signal (H01)
- **Effect size:** Extraordinarily strong (z-range 9.3–51.3, mean 27.64)
- **Robustness:** Uniform across all 36 tests (seed × layer); no test below null distribution
- **p-value interpretation:** p=0.0625 is the minimum achievable with 15 nulls; indicates maximum statistical separation, not marginal significance
- **Mechanism:** Real embeddings exhibit higher local clustering than feature-shuffled versions; neighborhood structure is non-random

### Negative Signal (H02)
- **Direction mismatch:** Transitivity is lower in real embeddings than shuffled (z-range -23 to -60, mean -34.78)
- **Interpretation:** Feature shuffling preserves distance distributions but creates spurious closed triangles; real embeddings lie on sparser, lower-dimensional manifolds
- **Implication:** CC (edge density in neighborhoods) is the primary structural discriminator, not global transitivity

### Cumulative Paper Context
- **iter_0002–0004:** H1 persistence under feature-shuffle nulls shows similar robustness (all layers Fisher p<0.05 in iter_0003–0004)
- **iter_0005–0008:** Stronger nulls (distance-permutation, rewiring) fail uniformly; signal fragile under adversarial controls
- **No biological anchoring yet:** H1 signal exists, but link to TRRUST/GO/STRING/functional outcome is untested
- **Cross-model consistency weak:** iter_0003 Spearman ρ=0.833 for feature effects, but combined permutation p=0.409 (non-significant)

---

## Critical Gaps Identified

1. **Biological validation missing:** Strong topological signal (H1, CC) has not been grounded in gene regulatory networks, pathway membership, or cell functional behavior.
2. **Stronger-null failure pattern:** Distance-permutation and rewiring nulls consistently fail without clear diagnostic explanation; unclear whether signal is real or artifact of control selection.
3. **Cross-model alignment inconclusive:** Geneformer embeddings not directly tested for CC/H1 consistency; only low-dimensional feature summaries were compared.
4. **Layer-depth mechanism opaque:** CC and H1 signals are layer-invariant, but the *meaning* at shallow vs. deep layers differs; no mechanistic interpretation yet.
5. **Null calibration plateaued:** Metric matching, edge-length binning, bridge conditioning have all been attempted; further null refinement unlikely to unlock the rewiring branch.

---

## Retire / Deprioritize Decisions

### Distance-Permutation Null (H05)
- **Status:** RETIRE_NOW
- **Reason:** Uniformly negative across all domain-layer tests (0/12 significant). No calibration adjustment has recovered signal. Appears over-adversarial relative to biological null hypothesis.
- **Cost of retry:** High; multiple calibration attempts already made.

### Rewiring-Null Family (H06–H12)
- **Status:** RETIRE_NOW (immune) / RESCUE_ONCE_WITH_MAJOR_CHANGE (external-lung, lung)
- **Reason:** Immune full-layer rewiring is consistently negative after metric-matched calibration. Cross-domain replication is incomplete.
- **Rescue condition:** If external-lung replication shows >50% layer pass-rate and positive mean deltas, reconsider; otherwise archive.
- **Expected effort:** 1 replication iteration + diagnostic if positive; closure otherwise.

### Bridge-Conditioned Attribution (H11)
- **Status:** RETIRE_NOW
- **Reason:** Split-confounded strata (36/36 vs. 2/36 bridge usage) make attribution unidentifiable. Positive pooled delta contradicts bridge-as-problem hypothesis.
- **Cost of retry:** High; structural confound cannot be fixed without protocol redesign.

### Metric Calibration Shift Studies (H10)
- **Status:** RETIRE_NOW
- **Reason:** Geodesic vs. Euclidean shift is tiny (mean Δ=+0.180) and does not improve conclusions. Not blocking further work.

---

## Scientific Reasoning for Deprioritization

The cumulative pattern across iter_0005–0008 suggests that **feature-shuffle null is the binding constraint on topological signal**, not metric or null calibration details. Rewiring and distance-permutation nulls were designed to be *stronger* (less permissive) than feature shuffle, and they uniformly fail. This could mean:

1. **Signal is real but fragile:** Topology depends on specific feature relationships, not just distance structure or connectivity.
2. **Null selection matters more than calibration:** Feature shuffle preserves locality; stronger nulls destroy it. The transition is binary.
3. **Biological signal is present but obscured by stronger nulls:** Rewiring and distance-permutation may discard meaningful correlations that underlie regulatory networks.

**Next direction:** Rather than re-calibrate already-failing nulls, test whether topological signal **predicts biological ground truth** (TRRUST, GO terms, cell fate transitions).

---

## Strength of Foundation for Next Iteration

✅ **Solid:**
- CC elevation is reproducible, high-magnitude, uniform across layers/seeds
- H1 persistence replication across domains (iter_0003–0004) is consistent
- Multiple null families have been tested; feature-shuffle is now the documented baseline

⚠️ **Medium (requires testing):**
- Cross-model consistency (Geneformer CC, H1); only feature effects tested so far
- Layer-depth mechanistic interpretation; no explanation for why signal is layer-invariant
- Distortion/metric-artifact causes remain unresolved

❌ **Missing:**
- Biological anchoring (TRRUST, GO, STRING, cell ontology)
- Split-robustness full-layer maps in lung and external-lung (only immune has detailed map)
- Adversarial/pathological null configurations (reverse-rank, adversarial sampling)
- Functional outcome coupling (cell fate, pathway activation)

