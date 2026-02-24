# Executor Next Steps — iter_0004

## Top Priority for iter_0005

### 1. Cross-seed validation of H01 (TwoNN ID) and H03 (CKA)
Both H01 and H03 tested on seed42 only. Need to replicate on seeds 43 and 44 to confirm:
- TwoNN ID decreasing trend across layers (positive if consistent across seeds)
- CKA ≈ 1.0 is not an artifact of a single seed's embedding structure

**Commands:**
```bash
# Add seeds 43, 44 to the analysis loop
```

### 2. H01 follow-up: MADA and local PCA estimators
TwoNN is sensitive to noise; cross-check with:
- MADA (maximum likelihood intrinsic dimensionality)
- Local PCA scree analysis per neighborhood
- Test whether ID decreases monotonically with layer or shows non-monotonic structure

### 3. H03 follow-up: Examine layer delta magnitudes
If CKA ≈ 1.0 everywhere, the representation barely changes. Quantify:
- L2 norm of layer-to-layer residual changes
- Whether CKA is exactly 1.0 (up to numerical precision) or has meaningful variation
- Compare to permutation-invariant CKA (uncentered) to separate mean-shift from structure

### 4. H02 rescue: Expand gene pool and use cosine similarity
Current TRRUST test used only 209 genes (too small pool for meaningful null).
- Extract all gene names from a larger source (e.g., scGPT vocabulary file)
- Test co-target clustering using cosine similarity in original 512-dim space
- Compare to GO-term annotation-based groupings

### 5. New hypothesis: Persistent homology (H1 Betti) on scGPT gene subspace
Prior iterations tested H1 on immune embeddings with rewiring null (negative). Try:
- H1 on scGPT lung gene embeddings (different domain)
- Feature-shuffle null (shown to work for CC in iter_0002)
- Cross-layer H1 comparison (does topology change across depth?)

## Retired Directions

- Distance-permutation null (prior iters): retired (0/12 sig)
- Degree-preserving rewiring null (immune): retired (0/24 sig)
- kNN transitivity (iter_0002 H02): negative, do not revisit
- TRRUST co-target clustering with 209-gene pool: inconclusive; revisit only with larger gene pool

## Hypothesis Roadmap (top candidates for iter_0005)

| Priority | Hypothesis | Family | Rationale |
|----------|-----------|--------|-----------|
| 1 | Cross-seed validation H01+H03 | intrinsic_dim / cross_model | Confirm promises with replication |
| 2 | H1 Betti on lung gene embeddings + feature-shuffle null | persistent_homology | New domain, compatible null |
| 3 | Expanded TRRUST/GO test with full gene vocabulary | module_structure | Fix gene pool limitation |
| 4 | Layer-delta analysis (residual norm vs CKA) | intrinsic_dim | Explain CKA=1 mechanistically |
