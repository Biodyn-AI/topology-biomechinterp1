# Brainstormer Structured Feedback — iter_0055

**Date**: 2026-02-23
**Gate**: PASSED (all 3 hypotheses positive, cross-seed replication confirmed)

---

## Iteration Assessment

Three out of three hypotheses yielded clean, significant, interpretable results. This is the most consistent positive batch since iter_0053. The findings converge to a coherent mechanistic picture:

1. **H01 (Cross-seed replication)**: SV5-7 TF/target-only AUROC=0.71–0.74 at L0 is fully reproducible. This closes the replication question definitively.

2. **H02 (Layer-resolved AUROC)**: The complementary encoding pattern — SV5-7 dominant at early+late, SV2-4 dominant at mid — is a structurally novel finding. The near-chance SV5-7 AUROC at L8 (0.543) coinciding with the directionality crossover at L9 is suspicious and mechanistically important.

3. **H03 (Opposite directionality dynamics)**: The crossover at L9 is the most striking quantitative signature in the current dataset. Two subspaces with opposite trajectories that equalize precisely at L9, then both decay at L10-L11, suggests L9 is a genuine phase transition layer.

## What the Data Is Telling Us

The current findings establish a two-subspace regulatory architecture:
- **SV5-7**: carries early-layer (L0-L3) regulatory identity (TF vs. target) encoded before the transformer has processed context deeply; signal decays through L4-L8; recovers at L9-L11.
- **SV2-4**: carries mid-layer (L4-L8) regulatory identity; dominant when SV5-7 is at minimum; then fades at L9+.
- **L9 crossover**: the point where SV5-7 magnitude exceeds SV2-4 for the first time AND SV5-7 AUROC rebounds. This is a second-phase activation of early-layer information or an information routing event.

The natural next frontier is: **what is L9 doing, and can we dissect the joint SV2-7 space to produce a layer-stable classifier?**

## Directions That Need to Be Explored

### Immediate consolidation (high-probability)
- Joint SV2-7 (6D) logistic regression classifier across all 12 layers: expect higher, flatter AUROC if the two subspaces are complementary as claimed.
- Directionality trajectory replication on seed43/seed44: verifies the crossover at L9 is not a main-seed artifact.
- Label-shuffle null for displacement magnitude: validates the absolute direction signal (both subspaces p<0.001 per layer, but need null for magnitude trajectory).

### Mechanistic dissection (high-reward)
- GO/biological enrichment of top-displacement genes at L9: which genes move most in SV5-7 at the crossover layer?
- Biological identity of the SV5-7 directionality vector: what direction does the mean TF→target displacement point in SV5-7 space at L0 and at L9?
- L8 minimum anatomy: why does SV5-7 AUROC drop to 0.543 at L8? Is this a specific set of genes losing or gaining TF-like coordinates?

### New directions
- Persistent homology on SV5-7 projections of circuit genes at L0 vs L9: is there different loop/cluster topology at the two AUROC peaks?
- SV5-7 early-layer specificity: does the early-layer SV5-7 signal survive a within-TF-family null (testing not just TF vs non-TF, but whether SV5-7 separates gene classes beyond co-expression)?
- Cross-model test: does Geneformer show any analogous early-layer subspace regulatory encoding in SV5-7-equivalent directions?
- Information routing hypothesis: Mutual information between SV5-7 coordinates at L0-L2 and SV2-4 coordinates at L5-L8 — are early SV5-7 signals preserved in mid-layer SV2-4 (rotation) or truly independent?

## Stale Directions

- H1 Betti loops on circuit genes: retired (iter_0052 H02).
- OOV bifurcation narrative: retired (iter_0030).
- TCR circuit attractor: insufficient vocabulary, retired after iter_0045.
- Lineage centroid cosine orthogonality: retired (iter_0046).
- BCL6 metabolic isolation as divergence narrative: retired (iter_0044).
- Dorothea proximity (hub-corrected): indefinitely deferred.

## Meta-observation

The paper is now rich in mechanistic claims but missing two things that would make it publication-ready:
1. A null-robust, joint-subspace classifier showing the combined SV2-7 space is a stable, high-AUROC regulatory encoder.
2. A biological identity for the SV5-7 directionality vector — what does the direction TFs-to-targets actually mean in gene-space terms?
