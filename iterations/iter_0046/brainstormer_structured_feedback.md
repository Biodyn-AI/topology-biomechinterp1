# Brainstormer Structured Feedback — iter_0046

## Research Gate Status
PASSED. Three hypotheses tested, two strong positives (H01, H02), one decisive negative (H03).

---

## Finding Assessment

### H01: Cross-seed ID compression (POSITIVE — primary validated finding)
- TwoNN ID: L0~33 → L11~19, ratio≈0.57–0.59 across three independent seeds
- Bootstrap CI L0=[27,37] vs L11=[16,20] — non-overlapping
- This is the strongest replicated structural finding of the project
- Next: generalize to other tissue contexts (cycle1, cycle6); characterize the compression rate curve

### H02: SVD spectral collapse (POSITIVE — novel, independent corroboration)
- Effective rank: L0=23.6 → L11=1.64 (14× collapse)
- Top-1 variance fraction: 53.7% → 93.4%
- Null (feature-shuffled): eff_rank=28.86 — real L11 is 17.6× below null
- Critical insight: the two measures tell complementary but compatible stories — global linear structure is nearly 1D, but local nonlinear structure retains ~19 dimensions. This means the representation lies on a curved (nonlinear) manifold with a single dominant linear direction.
- The dominant SV1 direction is unidentified biologically. Identifying it is now high priority.

### H03: Lineage centroid orthogonality (NEGATIVE — retired)
- All z-scores < 0.3 at all layers; indistinguishable from random gene sets
- The dominant shared direction (SV1 at 93.4%) causes all centroid cosine similarities to approach 1.0 equally
- Side finding: 2902/4941 (58.7%) embedding rows have zero L2-norm at L11; only ~2039 genes have active embeddings
- This zero-norm sparsification is a structural feature that has not yet been analyzed

---

## Core Tension to Resolve

The TwoNN (local, nonlinear) says ID≈19; SVD eff_rank (global, linear) says ~1.64. This is consistent with a curved manifold having one dominant direction that explains ~93% of global variance, but with genuine nonlinear variation in the remaining directions. The key question: **are the ~18 residual nonlinear dimensions biologically structured, or noise?**

Testing this requires projecting out SV1 and re-measuring ID and biological enrichment in the residual subspace.

---

## Retired Directions

- **Lineage centroid orthogonality via cosine similarity**: killed by dominant shared direction. Any cosine-based centroid comparison is uninformative until SV1 is removed.
- **Raw centroid comparisons without dominant direction subtraction**: generally retire all such tests unless SV1 is first projected out.

---

## Critical Unanalyzed Feature

**Zero-norm embeddings**: 58.7% of genes are zero at L11. This is not noise — it's a structural property. Questions:
1. Which genes go to zero? Are they biologically distinguishable from nonzero genes?
2. At which layer do genes first go to zero? Is zero-onset layer-coordinated within GO terms?
3. Does zero-fraction increase monotonically with layer depth?
4. Is ID compression partly driven by sparsification (more zeros → lower effective dimensionality)?

This is cheap to test and could reframe all prior geometry results.
