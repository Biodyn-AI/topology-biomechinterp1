# Brainstormer Hypothesis Roadmap: iter_0045 → iter_0046+

---

## Retire / Deprioritize

| Direction | Decision | Justification |
|---|---|---|
| kNN lineage purity (n<3/lineage) | **retire_now** | Impossible by construction; same result guaranteed |
| BCL6 B-cell divergence variants | **retire_now** | 4 tests across iter_0042–0044 all negative |
| 5-point TwoNN isolated-subset ID | **retire_now** | Invalid method (simplex ≠ manifold) |
| STRING PPI proximity → distance | **retire_now** | AUROC≈0.5, confirmed dead iter_0031 |
| GO/co-reg annotation → distance | **retire_now** | Two negatives iter_0035, confirmed dead |
| TCR circuit (single CD28 gene, p=0.32) | **rescue_once_with_major_change** | Need ≥3 circuit genes + V-shape reframe |

---

## New Hypothesis Portfolio

### H-A: Cross-Seed ID Compression Replication
**Hypothesis**: The monotone ID compression (L0=32.57→L11=18.05) observed in cycle4_immune is
reproducible across independent training seeds (cycle4_seed43, cycle4_seed44) and across training
cycles (cycle1), confirming it is a stable architectural property of scGPT and not a single-run artifact.
**Test**: Compute TwoNN on cycle1 [12,4803,512] and cycle4_seed43/44 at all 12 layers.
Compare ID trajectories: same delta direction, similar absolute values?
**Expected signal**: All seeds show monotone compression; |delta| within ±5D of iter_0045 result.
**Null**: Different seeds produce ID trajectories that cross (non-monotone or different direction).
**Value**: high | **Cost**: low (same code, different input files)

---

### H-B: Bootstrap CI for ID Compression Stability
**Hypothesis**: The TwoNN ID estimates at L0 and L11 are stable across subsampling regimes
(n_samples=500..3000), and the L0–L11 gap persists at all subsampling levels.
**Test**: At L0 and L11, compute TwoNN ID for n_samples in {500, 1000, 1500, 2000, 3000}.
10 bootstrap replicates each. Report 95% CI. Confirm CI bands do not overlap at L0 vs L11.
**Expected signal**: Non-overlapping CIs at L0 vs L11 for all n_samples ≥ 500.
**Null**: ID estimates are unstable (CIs overlap L0 vs L11 at small n).
**Value**: high | **Cost**: low

---

### H-C: Functional Gene Subset Local ID (B-cell vs T-cell vs Random)
**Hypothesis**: B-cell markers occupy a lower-dimensional local submanifold than random gene sets
at deep layers (L8+), and B-cell local ID decreases faster across layers than T-cell local ID or
random gene sets.
**Test**: For B-cell markers ({MS4A1, CD79A, BLK, PAX5, BATF, BACH2, JCHAIN} n~7),
T-cell markers ({CD8A, PRF1, CD247, RUNX3, CD28, LCK} n~6), and 50 random 7-gene sets:
compute TwoNN ID on the 50-NN neighborhood of each group's centroid (local ID in ambient space).
Track across L0–L11. Compare compression rates (Spearman rho of ID vs layer).
**Expected signal**: B-cell local ID compresses faster (rho < -0.85) than random (rho ~ -0.5).
**Null**: All groups show identical local ID trajectories.
**Value**: high | **Cost**: medium

---

### H-D: ID Compression Rate Correlates with Attractor Formation Layer
**Hypothesis**: The per-layer ID compression is largest at L2–L4, coinciding with the established
GC-TF attractor onset (L3), indicating that ID compression and functional attractor formation are
mechanistically linked.
**Test**: Compute delta ID between consecutive layers (ID_L_k − ID_L_{k+1}). Plot per-layer
compression delta. Test if max delta occurs at L2–L4 window. Compare to random permutation of
layer order.
**Expected signal**: Max compression delta at L3 or L4; significant peak vs flat permutation null.
**Null**: Compression is uniform across all layers (no peak).
**Value**: high | **Cost**: low (derivative of H-A/H-B data, no new computation needed)

---

### H-E: CD28 V-Shape Replication + AP1 Circuit Expansion
**Hypothesis**: The CD28 transient convergence trajectory (rank 595→86 at L6→149 at L11)
is reproducible across cycles and extends to AP1 family members (FOS, JUN, JUNB if in vocab),
revealing a layer-6-specific convergence zone for T-cell co-stimulatory circuit genes.
**Test**: In cycle1 and cycle4_seed43, recompute CD28 rank toward T-cell centroid {CD247, CD8A, RUNX3}.
Also check FOS, JUN, JUNB, FOSL2 if present in vocab. Report layer of minimum rank for each gene.
**Expected signal**: CD28 V-shape reproduced in ≥2 of 3 independent cycles; AP1 members show L6 minimum.
**Null**: V-shape is noise; CD28 rank trend not reproducible across cycles.
**Value**: high | **Cost**: low

---

### H-F: Singular Value Spectrum Concentration Across Layers
**Hypothesis**: The singular value spectrum of the gene embedding matrix concentrates (top-k SVs
capture larger fraction of variance) from L0 to L11, consistent with the ID compression finding
and providing a complementary characterization of representational geometry.
**Test**: SVD of [4941 × 512] embedding matrix at each layer. Compute fraction of variance in top
k=10, 20, 32, 50 singular values. Track cumulative variance explained. Null: shuffled-row matrix.
**Expected signal**: Cumulative variance at top-32 SVs increases from ~40% at L0 to ~65%+ at L11.
**Null**: Shuffled-row matrix shows flat cumulative variance (each gene contributes equally).
**Value**: high | **Cost**: low (standard SVD, fast on 4941×512)

---

### H-G: B-Cell vs T-Cell Subspace Orthogonality
**Hypothesis**: The lineage identity directions (B-cell centroid − background mean) and
(T-cell centroid − background mean) become more orthogonal (lower cosine similarity) at deep
layers, quantifying geometric separation of lineage-specific representation axes.
**Test**: At each layer, compute cosine similarity between unit vectors:
  v_B = normalize(B-cell centroid − all-gene mean)
  v_T = normalize(T-cell centroid − all-gene mean)
Track cosine similarity L0–L11. Null: random same-size gene sets (cosine sim between random pairs).
**Expected signal**: cosine similarity decreases from ~0.7 at L0 to ~0.2 at L11.
**Null**: cosine similarity flat or indistinguishable from random pairs.
**Value**: medium | **Cost**: low

---

### H-H: Cross-Lineage Mutual Repulsion Test
**Hypothesis**: B-cell markers are pushed away from the T-cell centroid across layers (rank
increases), and T-cell markers are pushed away from B-cell centroid — active mutual repulsion
as an additional lineage separation mechanism beyond attraction.
**Test**: Track rank of B-cell markers {MS4A1, CD79A, BLK, PAX5} toward T-cell centroid, and
T-cell markers {CD8A, PRF1, CD247} toward B-cell centroid, at L0–L11. Spearman rho of rank vs layer.
Null: same-size random gene sets tracked toward opposing centroid.
**Expected signal**: Rank increases monotonically (rho > +0.8 for real markers, ~ 0 for random).
**Null**: Rank flat or decreasing (no repulsion signal).
**Value**: medium | **Cost**: low

---

### H-I: Plasma Cell TF Trajectory — Late Divergence from GC Program
**Hypothesis**: Plasma cell TFs (IRF4, XBP1) track B-cell centroid proximity early (L0–L5) then
diverge at L6+, whereas GC-TFs (BATF, BACH2, PAX5) maintain or increase B-cell proximity across
all layers — reflecting the GC-to-plasma fate decision geometry.
**Test**: Track rank of {IRF4, XBP1} (if in vocab) vs {BATF, BACH2, PAX5} toward B-cell centroid.
Test: do the two groups diverge in rank trajectory at L6 (rank divergence = different Spearman rho)?
Null: random same-size gene sets show flat ranks.
**Expected signal**: GC-TF ranks decrease L0→L11 (rho < -0.8); IRF4/XBP1 ranks decrease then
reverse at L6 (rho near 0 overall; negative at L0–L5, positive at L5–L11).
**Null**: Both groups co-move (like BCL6 co-move from iter_0044).
**Value**: high | **Cost**: low

---

### H-J: Ollivier-Ricci Curvature at B-Cell Edges Across Layers
**Hypothesis**: kNN graph edges between B-cell markers show a sign change or inflection in
Ollivier-Ricci curvature at Layer 3 (established attractor onset), reflecting a topological
phase transition from negative (divergent) to positive (convergent) curvature.
**Test**: Build k=10 kNN graph on cycle4_immune at each layer. Compute edge-level Ollivier-Ricci
curvature (GraphRicciCurvature). Extract curvature of edges with ≥1 endpoint in B-cell markers.
Track mean curvature L0–L11. Null: random same-size edge sets.
**Expected signal**: B-cell edge curvature negative at L0–L2, inflects positive at L3+.
**Null**: Flat curvature trajectory; no inflection at L3.
**Value**: high | **Cost**: medium (library installation risk; curvature computation is slow)

---

### H-K: Wasserstein Distance Between Adjacent-Layer Persistence Diagrams
**Hypothesis**: The Wasserstein/bottleneck distance between H1 persistence diagrams of adjacent
layers peaks at L3 (GC-TF attractor onset), indicating topological reorganization concentrated
at that transition.
**Test**: Compute H1 persistence diagram on cycle4_immune embeddings (195 or 361 in-vocab gene
subset) at each of L0–L11. Compute Wasserstein distance between adjacent pairs (L0,L1), (L1,L2)...(L10,L11).
Null: random permutation of layer order.
**Expected signal**: Wasserstein(L2,L3) or Wasserstein(L3,L4) > all other adjacent pairs.
**Null**: Distances are uniform across adjacent layers (no peak at L3).
**Value**: medium | **Cost**: medium (Ripser + Wasserstein, moderate compute)

---

### H-L: Per-Layer ID Decomposition by Biological Gene Group
**Hypothesis**: The 18D total manifold at L11 is partitioned into distinct geometric subspaces
by lineage — B-cell, T-cell, myeloid gene sets occupy non-overlapping regions of the 18D space,
and this partition is absent at L0.
**Test**: At L0 and L11, perform PCA to 32D. Project B-cell, T-cell, myeloid marker gene sets.
Compute: (1) principal angle between subspaces spanned by each group, (2) variance explained
by first 3 PCs of each group separately. Test if principal angles are larger at L11 than L0.
Null: random same-size gene sets.
**Expected signal**: Principal angles between lineage subspaces larger at L11 vs L0; each lineage
occupies more orthogonal directions.
**Null**: Lineage subspaces are not orthogonalized across layers.
**Value**: high | **Cost**: medium

---

### H-M: Cross-Model ID Comparison (scGPT vs Geneformer)
**Hypothesis**: If another transformer model (Geneformer) trained on similar scRNA data shows
a similar ID compression trajectory (monotone decrease across layers), the finding is model-class
general. If Geneformer shows flat or different ID, it is scGPT-specific architecture.
**Test**: Extract Geneformer residual embeddings for same gene set. Compute TwoNN at each layer.
Compare L0→L11 ID trajectory to scGPT results.
**Expected signal**: Geneformer also shows monotone compression (lower ID at final vs initial layer).
**Null**: Geneformer ID trajectory is flat or non-monotone.
**Value**: high | **Cost**: high (requires Geneformer embedding extraction)

---

### H-N: Embedding Norm as Secondary Geometry Signal
**Hypothesis**: GC-TFs (BATF, BACH2, PAX5) show an embedding L2 norm increase at L3+
coinciding with attractor formation, distinct from background gene norms.
**Test**: Record L2 norm of BATF, BACH2, PAX5 at each layer vs mean background norm and 200
random same-size gene sets. Spearman rho of norm vs layer.
**Expected signal**: GC-TF norms increase at L3+ (rho > 0.7); background is flat.
**Null**: Norms indistinguishable from random gene sets.
**Value**: low | **Cost**: low (trivial computation)

---

## Top 3 for Immediate Execution

### Slot 1: High-Probability Discovery Candidate
**H-A + H-B + H-D combined: ID Compression Replication + Bootstrap + Per-Layer Rate**

Rationale: The iter_0045 ID compression finding (L0=32.57→L11=18.05, monotone) is the strongest
result in the project. Before building narrative on it, it needs cross-seed replication (H-A) and
bootstrap CI (H-B). H-D (per-layer compression delta) is a free derivative requiring no new data.
These three are a single coherent validation package that either confirms the finding (expected,
likely publishable) or reveals instability (critical to know). All are low-cost on existing data.

**Implementation**: One script computing TwoNN at all layers for cycle1 + cycle4_seed43 + cycle4_seed44,
plus bootstrap loop (n_samples sensitivity), plus per-layer delta computation. Expected runtime:
<30 min total.

---

### Slot 2: High-Risk / High-Reward Candidate
**H-F: Singular Value Spectrum Concentration**

Rationale: SVD of the embedding matrix is a complementary geometric characterization that captures
the same phenomenon as TwoNN ID (dimensional collapse) but from a different angle (covariance
structure). If SVD cumulative variance at top-32 SVs increases from L0 to L11, it independently
corroborates the ID finding and adds spectral analysis as a methodological tool. If it diverges
from the TwoNN result, it will constrain the interpretation. This is fast (single SVD per layer),
directly interpretable, and could reveal the specific axes along which collapse occurs.
High reward: if top SVs align with lineage axes, this becomes a story about spectral lineage encoding.

**Implementation**: np.linalg.svd (full=False) on [4941×512] at each layer; cumulative variance
at k=10/20/32/50. ~10 min total. Null: row-shuffled matrix.

---

### Slot 3: Cheap Broad-Screen Candidate
**H-G + H-H combined: Lineage Subspace Orthogonality + Cross-Lineage Repulsion**

Rationale: Both are cheap (cosine similarity and rank computations using existing embeddings and
gene sets, ≤50 lines of new code each). H-G tests the directionality hypothesis — do lineage
identity axes become orthogonal at deep layers? H-H tests repulsion — do B-cell genes actively
move away from T-cell centroid? Together they provide a geometric picture of lineage separation
complementing the attractor (convergence) story. Low cost, interpretable outcomes, either result
advances the narrative.

**Implementation**: Two independent rank/cosine-tracking loops on cycle4_immune. 200 random null
sets for each. ~15 min total.
