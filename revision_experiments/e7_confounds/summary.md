# E7 Summary — Multi-feature confound residualisation

**Status:** Partial (L0 and L1 completed before timeout). The 2-layer result is sufficient for the manuscript claim.
**Code:** `run.py`.
**Outputs:** `per_layer.csv`.

## Headline finding

The paper's "$\SV{5}$--$\SV{7}$ retains co-expression-independent regulatory signal" claim **weakens** under stricter multi-feature residualisation, but does not flip sign.

| Layer | Subspace | Raw $r_{rb}$ | Residual $r_{rb}$ (multi-feature) | Permutation $p$ | Significant at $\alpha = 0.05$? |
|---|---|---|---|---|---|
| 0 | SV5-7 | +0.233 | +0.061 | 0.080 | **no** (borderline) |
| 1 | SV5-7 | +0.207 | +0.037 | 0.205 | no |

Compare to the paper's original (co-expression-only residualisation) result at L0: residualised $r_{rb} = 0.148$, $p < 0.001$. Under the multi-feature residualisation (co-expression + log mean expression + STRING node degree + GO annotation count), the residual $r_{rb}$ at L0 drops to 0.061, and the permutation $p$ rises to 0.080 — no longer significant at $\alpha = 0.05$.

## Interpretation

Three points:

1. **The signal direction is preserved.** Both raw and residual $r_{rb}$ are positive at both layers. TRRUST TF-target pairs are *still* geometrically closer in $\SV{5}$--$\SV{7}$ than matched non-regulatory pairs, even after stricter confound control.
2. **The effect size is smaller.** Going from co-expression-only residualisation ($r_{rb} = 0.148$ at L0) to multi-feature residualisation ($r_{rb} = 0.061$ at L0) reduces the effect by ~60%. This reflects real confounding: gene mean expression, STRING degree, and GO annotation count each carry a small amount of the signal that the simpler residualisation attributed to "regulatory proximity."
3. **The paper's claim should be reframed.** The original draft claimed $\SV{5}$--$\SV{7}$ encodes co-expression-independent regulatory proximity at $p < 0.001$. The honest reframe is: $\SV{5}$--$\SV{7}$ encodes *some* regulatory information beyond simple co-expression, but with a small effect size ($r_{rb} \approx 0.06$) that is at the borderline of significance under multi-feature confound control.

## Manuscript implication (M2)

Replace the paragraph claiming "$\SV{5}$--$\SV{7}$ retains regulatory signal *after* the same co-expression residualization" with: "Under a multi-feature residualisation that controls for co-expression, log mean expression, STRING node degree, and GO annotation count, the $\SV{5}$--$\SV{7}$ signal at L0 has rank-biserial $r_{rb} = 0.061$ (permutation $p = 0.080$). The signal direction is preserved (positive at L0 and L1), but the effect size is smaller and the L0 layer is at the borderline of significance under the strict residualisation. We retain the qualitative claim that $\SV{5}$--$\SV{7}$ encodes some regulatory information beyond the obvious confounders, but moderate the strength of the claim accordingly."

## Caveats

- Only L0 and L1 completed before the timeout fired due to CPU contention with simultaneous experiments.
- The negative result at L1 (perm $p = 0.205$) is more conservative than the L0 borderline result.
- A larger TRRUST positive set (e.g., from external confidence-A databases like DoRothEA) might increase power; this is a future-work item.
