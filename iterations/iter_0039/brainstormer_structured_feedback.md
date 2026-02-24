# Brainstormer Structured Feedback — iter_0039

## Gate Status
PASSED. 3/3 hypotheses tested; machine-readable artifacts present; paper updated.

---

## Assessment of Each Result

### H01 — Drift Target (POSITIVE)
Strong result. BATF rank-2, SPIB rank-5 among 195 genes closest to the L11 B-cell centroid endpoint. This closes the open question from iter_0038: the drift endpoint IS the GC-TF neighborhood. The near-zero cosine alignment (0.008–0.064) correctly rules out directional drift, confirming the drift is a global geometric compression toward this attractor.

**Action**: Verified. Do NOT re-test. Build on it — next step is to determine whether the L11 GC-TF attractor is present across all layers (i.e., when does it first emerge?).

### H02 — Extended Plasma Panel + Clean Anchor (MIXED)
The primary value here is the corrected framing: plasma divergence is real (rho=+1.0, p<0.001) but relative — z goes from −1.35 to −0.81, not to positive. BCL6 and PAX5 added to the GC-TF core (now 4 genes). The PRDM1 artifact from iter_0038 is now explained.

**Correction urgency**: The paper likely still contains the "plasma goes z>0" claim from iter_0038. This must be corrected in the next iteration's paper update pass.

**Action**: Retire the "plasma goes positive" framing. The relative divergence (monotonically increasing z, rho=1.0) is the publishable finding. Add BCL6/PAX5 to GC-TF cluster description.

### H03 — LOO Ablation (POSITIVE)
CD19 is the critical anchor gene (LOO drops precision@10 by 50% at both L2 and L11). BLK matters at L11 but not L2, suggesting it encodes a layer-depth-specific aspect of B-cell geometry. MS4A1 and CD79A are redundant. This mechanistically dissects WHY the precision@10 signal is strong — it is not distributed equally across markers.

**Action**: CD19+BLK is now the minimal 2-gene anchor. Use this in future precision screens to reduce noise. The differential (BLK matters at L11 but not L2) deserves its own follow-up: does BLK uniqueness increase with layer depth?

---

## Open Questions Generated This Iteration

1. **When does the GC-TF attractor first emerge?** The drift endpoint is known at L11. At what layer does BATF first enter top-20? Full 12-layer scan of centroid-to-BATF distance needed.

2. **Is CD19 also the critical anchor in cycle4_immune?** LOO was done only in cycle1. Replication in cycle4_immune tests generalization across cell-type composition.

3. **Is IRF4 available in cycle4_immune for plasma-side enrichment?** H02 next-action notes this as open. IRF4 bridges the GC→plasma transition; its geometric behavior is mechanistically important.

4. **Does the GC-TF cluster replicate in Geneformer?** This is the highest-priority cross-model validation and has not yet been tested.

5. **What is the principal angle between GC-TF and plasma-TF subspaces using the extended 4-gene GC cluster?** Prior iter_0038 found ~94°. Should retest with BCL6/PAX5 included.

---

## Direction Health Summary

| Direction | Status | Action |
|-----------|--------|--------|
| B-cell centroid geometry / GC-TF proximity | **Active, strong** | Extend (layer-scan, Geneformer) |
| Plasma divergence (relative) | **Active, refined** | Report correctly; seek IRF4 |
| LOO ablation / anchor sensitivity | **Active, new** | Extend to cycle4, all-layer scan |
| GC-TF subspace orthogonality (principal angle) | **Pending** | Test now |
| Geneformer cross-model validation | **Pending, HIGH priority** | Test now |
| NK/myeloid screen | **Retired** | vocab limitation |
| Directional drift cosine alignment | **Retired** | confirmed null |
| STRING PPI Euclidean | **Retired** | Prior iters confirmed null |
| intrinsic_dimensionality SVD AUROC | **Retired** | Prior iter |
| Plasma z goes positive (absolute) | **Retired** | PRDM1 artifact |
