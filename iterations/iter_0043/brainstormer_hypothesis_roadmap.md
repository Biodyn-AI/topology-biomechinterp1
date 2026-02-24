# Brainstormer Hypothesis Roadmap — iter_0043

---

## Retire / Deprioritize

| ID | Hypothesis family | Reason |
|----|-------------------|--------|
| GC-repressor-anti-convergence | Repressors diverge from B-cell centroid | PRDM1 rho=-1.000 (iter_0043 H01), BACH2/BATF also rho=-1.000; no residual signal |
| BCL6-B-cell-divergence | BCL6 separates from B-cell centroid | rho=-0.993 two independent tests (iter_0042 H03, iter_0043 H02); falsified |
| Cross-layer CKA stability variants | Incremental CKA tests | Fully established (rho~1.000 all pairs); no new information expected |
| TRRUST co-target basic clustering | Co-targets cluster without stratification | Inconclusive (iter_0004); only worth revisiting with full gene vocab and layered stratification |

---

## New Hypothesis Portfolio

### H-A: Cross-lineage functional-specialization → ID compression law
**Hypothesis**: Functionally specialized terminal-effector gene sets compress in ID across layers; pan-lineage gene sets expand. This is a general law independent of lineage.
**Test**: Measure TwoNN ID trajectory for: (1) effector myeloid (ELANE, MPO, CTSG — neutrophil effectors, 3 genes), (2) pan-myeloid (CD14, LYZ, ITGAM — 3 genes), (3) NK effector (GNLY, NKG7, KLRD1, GZMB — 4 genes), (4) plasma cell (JCHAIN, SDC1 — 2 genes, if vocab). Compare ΔID (post-L7 vs pre-L7) and full-trajectory rho for each. Null: 20 random gene sets of matched size.
**Expected signal if true**: effector subsets rho < 0 (compression); pan-lineage rho > 0 or flat.
**Null/control**: Random matched-size gene sets; myeloid pan-lineage (already measured flat in iter_0042).
**Value**: HIGH | **Cost**: LOW (direct TwoNN computation, 4-5 gene sets, same pipeline)

### H-B: T-cell activation circuit attractor (analog of B-cell GC attractor)
**Hypothesis**: T-cell activation circuit genes (CD28, LAT, LCK, ZAP70, CD3G) form a progressive attractor toward T-cell centroid analogous to the B-cell GC-TF attractor, with layer-specific onset.
**Test**: Compute rank of CD28/LAT/LCK/ZAP70/CD3G near T-cell centroid (CD3E, CD3D, TRAC, TRBC1) at each layer. Measure Spearman rho of rank vs layer. Compare to: (a) null random genes, (b) B-cell GC attractor convergence curve. Check for onset layer.
**Expected signal if true**: rho < -0.8, onset layer 3-6, ranks improve monotonically.
**Null/control**: Random non-T-cell genes; the B-cell attractor is the positive benchmark.
**Value**: HIGH | **Cost**: LOW

### H-C: Naive vs effector T-cell ID (intermediate differentiation test)
**Hypothesis**: Naive T-cell markers (CCR7, LEF1, TCF7, SELL) have intermediate ID between effector CD8 (compresses) and pan-T (expands), consistent with naive T cells being less functionally committed than effectors but more committed than undifferentiated.
**Test**: TwoNN ID trajectory for naive T-cell markers (check vocab: CCR7, LEF1, TCF7, SELL). Compare to effector CD8 (ΔID=-45.8) and general T-cell (ΔID=+26.2) from iter_0043. If naive has ΔID between -45 and +26, hypothesis supported.
**Expected signal if true**: naive T-cell ΔID in range [-10, +10]; ordered effector < naive < pan-T.
**Null/control**: Matched-size random gene sets.
**Value**: MEDIUM | **Cost**: LOW

### H-D: GC circuit unity internal structure (sub-cluster analysis)
**Hypothesis**: Although activators (BATF, BACH2) and repressors (PRDM1) co-locate at L11, they are not uniformly distributed — PCA of the 5-member cluster reveals at least 2 sub-clusters (activator pole and repressor/plasma pole) explaining >30% of within-cluster variance.
**Test**: PCA of the 5×512 matrix {PAX5, BATF, BACH2, BCL6, PRDM1} at each layer. Compute variance explained by PC1. Project each gene onto PC1. Compare BATF/BACH2/PAX5 vs PRDM1/BCL6 PC1 scores. Also compute pairwise distances within activator vs activator-repressor subsets.
**Expected signal if true**: PC1 explains >30% variance; PRDM1 and BCL6 score on opposite side from BATF/BACH2.
**Null/control**: Random 5-gene sets; expect PC1 ~20% for isotropic cloud.
**Value**: MEDIUM | **Cost**: LOW

### H-E: PRDM1/plasma paradox — plasma markers test
**Hypothesis**: Despite PRDM1 co-localizing with GC TFs at L11, true plasma cell identity markers (JCHAIN, SDC1, MZB1) do NOT converge toward B-cell centroid — the GC circuit unity finding does not extend to downstream effector plasma cell genes.
**Test**: Track Euclidean distance of JCHAIN, SDC1, MZB1 (check vocab) to B-cell centroid across L0-L11. Compare rho to GC circuit members (rho~-1.000). Also compare pairwise distances PRDM1-JCHAIN vs BATF-JCHAIN at L0 and L11.
**Expected signal if true**: Plasma markers rho > 0 (diverge) or rho ~0 (flat); PRDM1-JCHAIN distance increases relative to BATF-JCHAIN.
**Null/control**: Null random genes; B-cell convergence curve is positive benchmark.
**Value**: HIGH | **Cost**: LOW

### H-F: BCL6 rank anomaly tracks metabolic cluster density
**Hypothesis**: BCL6 rank (relative distance to B-cell centroid) at each layer is inversely correlated with BCL6 metabolic cluster size (overlap count), i.e., when BCL6 is most metabolically isolated (layers 0-4), it ranks farthest from B-cell centroid in relative terms.
**Test**: Compute layer-by-layer Spearman correlation between BCL6 rank (from iter_0043 H02: [654, 772, 1311, 1077, 1164, 1096, 929, 1011, 790, 750, 770, 927]) and BCL6 metabolic overlap (from iter_0042: [9, 9, 9, 10, 10, 9, 9, 7, 6, 5, 4, 4]). Null: 20 random genes' rank vs BCL6 metabolic overlap.
**Expected signal if true**: rho > +0.5 (high metabolic isolation ↔ high B-cell rank).
**Null/control**: The two series are already measured; this is a pure correlation test.
**Value**: MEDIUM | **Cost**: VERY LOW (no new computation needed — use existing artifacts)

### H-G: Effector CD8 compression breakpoint layer
**Hypothesis**: The effector CD8 ID compression has a structural breakpoint (analogous to B-cell L3 breakpoint) between L0-L6 that can be precisely located, and this breakpoint does NOT coincide with the pan-T expansion breakpoint at L7.
**Test**: Apply exhaustive piecewise linear fit to effector CD8 ID trajectory [176.7, 197.2, 160.9, 165.5, 105.5, 80.0, 89.2, 87.8, 123.1, 98.3, 84.5, 73.9]. Find best breakpoint layer. 500 permutation null. Also fit general T-cell trajectory. Test if breakpoints are at different layers.
**Expected signal if true**: Effector CD8 breakpoint at L4-L5 (inflection from ~165 to ~80); pan-T breakpoint at L7 (already measured).
**Null/control**: 500 layer-permutation null; cross-lineage comparison.
**Value**: MEDIUM | **Cost**: LOW

### H-H: Myeloid functional specialization — neutrophil effector vs pan-myeloid
**Hypothesis**: Neutrophil effector genes (ELANE, MPO, CTSG, S100A8, S100A9) compress in ID across layers, while pan-myeloid genes (CD14, LYZ, ITGAM, FCGR3A) remain flat or expand, mirroring the effector CD8 vs general T-cell dissociation.
**Test**: TwoNN ID trajectory for neutrophil effector vs pan-myeloid. Check vocab: MPO, ELANE, CTSG, S100A8, S100A9, CD14, LYZ, ITGAM. 500 permutation null. Compute ΔID (post-L7 mean - pre-L7 mean) for each subset.
**Expected signal if true**: Neutrophil effector ΔID < -10; pan-myeloid ΔID > 0.
**Null/control**: Random matched-size gene sets; pan-myeloid is the in-lineage non-effector control.
**Value**: HIGH | **Cost**: LOW

### H-I: Persistent homology (H1 loops) in lineage-specific subspaces
**Hypothesis**: B-cell marker gene subspaces have systematically different H1 persistence (loop count/lifetime) than T-cell subspaces, reflecting their geometric attractor vs expansion geometry.
**Test**: Run Ripser on gene embedding subspaces per layer: B-cell markers (4), T-cell markers (9), effector CD8 (5), pan-T (4). Count H1 births, compute total persistence. Null: matched-size random gene subsets (20 replicates). Compare across lineages and layers.
**Expected signal if true**: B-cell subspace has different loop structure (potentially fewer/simpler) than T-cell expansion subspace.
**Null/control**: Random matched-size subsets; null PH from shuffled distances.
**Value**: MEDIUM | **Cost**: MEDIUM (Ripser on small matrices — fast)

### H-J: Geneformer input embedding — T-cell analog of B-cell null result
**Hypothesis**: Geneformer input token embeddings also show no T-cell activation circuit proximity signal (as B-cell was null in iter_0040), confirming that attractor-like organization is specific to scGPT learned contextual representations.
**Test**: Compute T-cell activation circuit (CD28/LAT/LCK/ZAP70) precision@10 near T-cell centroid in Geneformer input embeddings. Compare to scGPT contextual result (H-B above).
**Expected signal if true**: Geneformer input shows no T-cell activation circuit proximity; scGPT contextual does show proximity.
**Null/control**: Random gene sets in Geneformer space; scGPT B-cell null (iter_0040 H02) as precedent.
**Value**: MEDIUM | **Cost**: MEDIUM

### H-K: Angle between effector CD8 subspace and pan-T subspace across layers
**Hypothesis**: The principal angle between the effector CD8 gene subspace and the pan-T gene subspace increases across layers (they diverge geometrically), mirroring the GC-plasma angle increase (77°→94°) found in iter_0038.
**Test**: Compute principal angle between span({CD8A, GZMB, PRF1, NKG7, GNLY}) and span({CD3E, CD3D, TRAC, TRBC1}) at each layer L0-L11. Spearman rho of angle vs layer. Null: random gene subsets of same sizes.
**Expected signal if true**: rho > +0.6; angle increases from L0 to L11.
**Null/control**: Random matched-size subspaces; B-cell/plasma angle trajectory as comparison benchmark.
**Value**: MEDIUM | **Cost**: LOW

### H-L: B-cell centroid drift terminal destination with larger gene set (cycle4_immune)
**Hypothesis**: Using the full 4941-gene cycle4_immune embedding, the B-cell centroid drift vector (L0→L11) points toward GC circuit genes at its terminus, consistent with the GC circuit unity finding that all GC members converge there. Specifically, the top-10 genes nearest the L11 B-cell centroid in cycle4 should be enriched for GC circuit members vs the cycle1 result.
**Test**: In cycle4_immune (4941 genes), identify top-20 nearest neighbors of B-cell centroid (MS4A1, CD79A, BLK, PAX5) at L11. Count GC circuit member representation (BATF, BACH2, PRDM1, BCL6, PAX5). Compare cosine alignment of drift vector with each GC member's L11 position.
**Expected signal if true**: ≥3/5 GC circuit members in top-20; PRDM1 now included (absent from cycle1 analysis).
**Null/control**: Random 4-gene centroid; cycle1 result (top-5: FAM162A, BATF, EOMES, HDAC9, SPIB).
**Value**: MEDIUM | **Cost**: LOW

---

## Top 3 for Immediate Execution

### Rank 1 — High-probability discovery candidate
**H-A: Cross-lineage functional-specialization → ID compression law**

Rationale: iter_0043 H03 gave one clean positive (effector CD8 ΔID=-45.8 vs general T ΔID=+26.2). This is a single data point for a potentially general law. Testing NK effector, neutrophil effector (H-H is subsumed here), and plasma cell extends the same exact method to 3 additional lineages. If 2+ additional lineages confirm the pattern, this becomes a major cross-lineage finding about how scGPT encodes functional specialization. Cost is minimal — same TwoNN pipeline, 4-5 small gene sets. Directly executable with current artifacts.

**Concrete test plan**:
- Gene sets to test (all TwoNN ID across L0-L11):
  - NK effector: GNLY, NKG7, KLRD1, GZMB (4 genes — check vocab)
  - Plasma cell: JCHAIN, SDC1, MZB1 (3 genes — check vocab)
  - Neutrophil effector: MPO, ELANE, S100A8, S100A9 (4 genes — check vocab)
  - Pan-myeloid (control): CD14, LYZ, ITGAM (3 genes — already measured)
- Output: ΔID (post-L7 - pre-L7), full rho, trajectory plot per lineage
- Expected result: NK and neutrophil compress; plasma cell compresses; pan-myeloid flat

---

### Rank 2 — High-risk/high-reward candidate
**H-B: T-cell activation circuit attractor**

Rationale: The B-cell GC attractor is the project's strongest mechanistic finding (established over 10+ iterations). If the T-cell activation circuit (CD28/LAT/LCK/ZAP70) forms an analogous attractor — with layer-specific onset, monotone rank improvement, and co-convergence — this would extend the attractor mechanism beyond B-cells and suggest it is a general mechanism for encoding key regulatory circuits. The risk is that T-cell markers showed no precision@10 signal in earlier broad screens (iter_0037: T-cell z=-1.37), but the activation circuit is more specific than pan-T markers. Worth one focused test.

**Concrete test plan**:
- Centroid: CD3E, CD3D, TRAC, TRBC1
- Target genes: CD28, LAT, LCK, ZAP70, CD3G (check vocab for each)
- Metric: rank near T-cell centroid at each layer L0-L11, Spearman rho
- Null: 20 random non-T genes
- Comparison: B-cell GC-TF attractor curve (BATF/BACH2: rho=-0.972/-0.844, onset L3)

---

### Rank 3 — Cheap broad-screen candidate
**H-F: BCL6 rank anomaly tracks metabolic cluster density (pure correlation)**

Rationale: Both data series already exist in iter_0042 and iter_0043 artifacts. Zero new computation. The hypothesis is testable by a single Spearman correlation between two length-12 vectors. If BCL6 rank peaks (relative distance from B-cell centroid) systematically coincide with BCL6 metabolic isolation peaks, this connects two independent findings (metabolic isolation and B-cell contextual rank) into a single coherent mechanistic picture. Very cheap; provides interpretive value regardless of direction.

**Concrete test plan**:
- BCL6 rank: [654, 772, 1311, 1077, 1164, 1096, 929, 1011, 790, 750, 770, 927] (from iter_0043 H02)
- BCL6 metabolic overlap: [9, 9, 9, 10, 10, 9, 9, 7, 6, 5, 4, 4] (from iter_0042 H02)
- Compute Spearman rho between these two series
- Null: shuffle one series 10000 times
- Expected: rho > +0.5 (high metabolic isolation → farther from B-cell in rank)
- This can be a standalone mini-analysis appended to H-A or H-B execution

---

## Portfolio coverage summary

| Family | Hypotheses assigned |
|--------|-------------------|
| Intrinsic dimensionality | H-A, H-C, H-G, H-H |
| Manifold distance / attractor | H-B, H-D, H-E, H-F, H-L |
| Subspace angle / geometry | H-K |
| Persistent homology | H-I |
| Cross-model | H-J |
| Mechanistic integration | H-D, H-F |
