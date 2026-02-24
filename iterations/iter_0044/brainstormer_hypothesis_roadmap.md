# Brainstormer Hypothesis Roadmap: iter_0044 → iter_0045+

---

## Retire / Deprioritize

| ID | Direction | Decision | Justification |
|---|---|---|---|
| BCL6-metabolic-isolation | BCL6 selectively drifts toward metabolic cluster | **RETIRE** | H03 iter_0044: co-moves, rho=0.5; third negative overall |
| BCL6-B-cell-divergence | BCL6 diverges from B-cell centroid | **RETIRE** | iter_0042 H01/H02 + iter_0043 H01 + iter_0044 H03 all negative |
| 5-point-TwoNN-ID | ID compression law via isolated-subset TwoNN | **RETIRE** | Methodologically invalid; 5-point simplex ≠ local manifold ID |
| STRING-PPI-proximity | PPI score predicts embedding distance | **RETIRED** (iter_0031) | AUROC≈0.5, confirmed dead |
| Gene-annotation-proximity | GO/co-reg annotation predicts distance | **RETIRED** (iter_0035) | Two negatives, confirmed dead |

---

## New Hypothesis Portfolio

### H-A: TCR Circuit Attractor (Clean Retest)
**Hypothesis**: TCR signaling circuit genes (CD28, LCK) converge toward the T-cell centroid (CD3E/CD3D/TRAC/TRBC1) across layers, with convergence only at deep layers (L7+).
**Test**: Remove CD247 from both centroid and circuit definitions. Use centroid = {CD3E, CD3D, TRAC, TRBC1}; circuit = {CD28, LCK}. Track mean rank and null p-value (200 random 2-gene sets). Additionally try PTPRC (CD45) as a third circuit gene.
**Expected signal if true**: rank decreases L0→L11, null p<0.05
**Null**: rank change indistinguishable from random 2-gene sets
**Value**: high | **Cost**: low
**Lineage**: iter_0044 H02 (rescue)

---

### H-B: Full-Manifold TwoNN ID Trajectory
**Hypothesis**: The global intrinsic dimensionality of the ~4941-gene cycle4_immune embedding manifold decreases monotonically from L0 to L11, consistent with progressive dimensional collapse (confirmed qualitatively for B-cell subset in iter_0042).
**Test**: Compute TwoNN ID on the full 4941-gene manifold at each layer L0–L11. Report global ID and local ID within 50-NN of B-cell centroid, T-cell centroid, myeloid centroid. Plot trajectories.
**Expected signal if true**: Global ID decreases L0→L11; B-cell local ID decreases with breakpoint ~L3
**Null**: Permuted embedding (feature-shuffled) ID trajectory is flat
**Value**: high | **Cost**: medium
**Note**: This also validates/replicates the iter_0042 B-cell ID breakpoint with correct method.

---

### H-C: Plasma Cell Fate Divergence from GC Attractor
**Hypothesis**: Plasma cell TFs (IRF4, XBP1) and GC-TFs (BATF, BACH2, PAX5) share B-cell-proximal positions at early layers but diverge geometrically at late layers, reflecting the GC-to-plasma differentiation decision.
**Test**: Track rank of IRF4 and XBP1 (if in vocab) toward B-cell centroid AND toward plasma centroid ({JCHAIN, SDC1, IGHG1}) at each layer. Compare rank trajectory to BATF/BACH2/PAX5 reference trajectories. Test for rank divergence at L7+.
**Expected signal if true**: IRF4/XBP1 and GC-TFs move together through L6, then IRF4/XBP1 diverges toward plasma centroid
**Null**: Both tracks co-move throughout (like BCL6 co-move finding)
**Value**: high | **Cost**: low

---

### H-D: kNN Community Structure at Full Scale (4941 genes)
**Hypothesis**: At the full-gene scale (4941 genes), kNN communities (k=10, Louvain) at L2 and L11 contain lineage-sorted B-cell, T-cell, and myeloid modules; community separation increases from L2 to L11.
**Test**: Build kNN graph on 4941 genes at L2 and L11. Run Louvain community detection. Compute Fisher OR for B-cell, T-cell, myeloid markers in communities. Compare modularity L2 vs L11.
**Expected signal if true**: Clear lineage-sorted communities, modularity increases from L2 to L11
**Null**: Permuted gene-label Fisher OR
**Value**: high | **Cost**: medium
**Note**: iter_0036 was positive on 195-gene set; scaling to 4941 is the natural extension.

---

### H-E: Cross-Lineage Circuit Unity Test
**Hypothesis**: The GC circuit unity (activator BATF/BACH2 + repressor PRDM1 all converge to same B-cell neighborhood) is NOT universal — T-cell activation circuits (AP1 activators JUN/FOS vs inhibitors NFATC1/CTLA4) show geometric separation of activators vs repressors.
**Test**: Track mean rank of T-cell AP1 activators (JUN, JUNB, FOS, FOSL2) and T-cell inhibitors/regulators (CTLA4, PDCD1, NFATC1) toward T-cell centroid. Compare convergence rates. Null: random gene sets.
**Expected signal if true**: Activators and inhibitors diverge; contrast with GC unity
**Null**: All converge similarly (unity is universal)
**Value**: high | **Cost**: low
**Note**: Would determine whether GC circuit unity is a discovery-worthy exception or a generic model property.

---

### H-F: T-Cell Exhaustion Sub-Lineage Attractor
**Hypothesis**: Exhaustion TFs TOX and NR4A1 (if in vocab) converge toward the T-cell centroid specifically at late layers, forming an exhaustion-specific sub-attractor within the T-cell neighborhood.
**Test**: Track TOX and NR4A1 rank toward T-cell centroid. Compare to naive T-cell TFs (TCF7, LEF1) and to GC attractor rank trajectory as positive control.
**Expected signal if true**: TOX/NR4A1 decrease rank at L7+ (late activation), distinct from TCF7/LEF1
**Null**: ranks flat or random
**Value**: medium | **Cost**: low

---

### H-G: Subspace Angle Between B-Cell and T-Cell Attractor Directions
**Hypothesis**: The geometric direction from background centroid to B-cell centroid and from background centroid to T-cell centroid becomes more orthogonal at deep layers, reflecting increasing lineage separation.
**Test**: At each layer, compute unit vector (B-cell centroid − background mean) and (T-cell centroid − background mean). Compute cosine similarity between these vectors. Track across L0–L11.
**Expected signal if true**: Cosine similarity decreases (more orthogonal) from L0 to L11
**Null**: Random subsets of same sizes show flat cosine similarity trajectory
**Value**: medium | **Cost**: low

---

### H-H: Ollivier-Ricci Curvature Transition at Layer 3
**Hypothesis**: The kNN graph Ollivier-Ricci curvature of edges between B-cell markers undergoes a sign change or inflection at Layer 3 (the established GC-TF attractor onset layer), reflecting a topological phase transition.
**Test**: Compute edge-level Ollivier-Ricci curvature (GraphRicciCurvature library) on k=10 kNN graph built on cycle4_immune embeddings at each layer. Extract mean curvature of edges involving B-cell markers. Track L0–L11. Compare to random edges.
**Expected signal if true**: B-cell edge curvature shifts from negative to positive (or shows inflection) at L3
**Null**: T-cell edge curvature and random edge curvature show no L3 transition
**Value**: high | **Cost**: medium
**Risk**: GraphRicciCurvature may not be installed; Ricci computation is expensive.

---

### H-I: Cross-Model B-Cell Axis (scGPT → Geneformer)
**Hypothesis**: The B-cell identity axis (PC1 negative pole) is a model-invariant biological feature, present in both scGPT and Geneformer gene embeddings.
**Test**: Extract Geneformer residual embeddings for the same 195-gene set. Compute PC1. Test if B-cell markers are enriched at the negative PC1 pole (same protocol as iter_0033/0034). Compute CKA similarity between scGPT and Geneformer PC1 directions.
**Expected signal if true**: Geneformer PC1 also has B-cell negative pole; CKA > 0.5
**Null**: Geneformer PC1 has no B-cell enrichment; CKA ~ random
**Value**: high | **Cost**: high
**Note**: High cost because Geneformer embeddings must be extracted/loaded.

---

### H-J: Persistent Homology on Immune Cycle4 B-Cell Neighborhood
**Hypothesis**: The B-cell neighborhood (top 50 nearest genes to B-cell centroid at L11) contains H1 loops that are absent from equally-sized random neighborhoods, and these loops are more persistent (longer bars) at L11 than L0.
**Test**: At L0 and L11, extract 50-NN of B-cell centroid. Run Ripser H1. Null: 50 random genes from the 4941-gene set (100 replicates). Compare sum of persistence lengths (total H1 persistence).
**Expected signal if true**: B-cell neighborhood H1 persistence > null; increases from L0 to L11
**Null**: B-cell neighborhood H1 indistinguishable from random 50-gene sets
**Value**: medium | **Cost**: medium

---

### H-K: Cross-Lineage Repulsion Test
**Hypothesis**: B-cell markers are actively pushed AWAY from the T-cell centroid as layers progress, and T-cell markers are pushed away from the B-cell centroid — mutual repulsion increases with depth.
**Test**: Track rank of B-cell markers (MS4A1, CD79A, BLK, PAX5) toward T-cell centroid, and rank of T-cell markers toward B-cell centroid, at each layer. Test if ranks increase (repulsion) with Spearman rho.
**Expected signal if true**: Rank increases monotonically (rho > 0.7); T-cell/B-cell markers are increasingly distant from the opposing centroid
**Null**: Random gene sets show flat rank trajectories toward centroids
**Value**: medium | **Cost**: low

---

### H-L: Bottleneck Distance Between Adjacent Layer Persistence Diagrams
**Hypothesis**: The Wasserstein/bottleneck distance between persistence diagrams of adjacent layers has a peak at Layer 3 (the GC-TF onset layer), indicating a topological phase transition coinciding with the structural breakpoints already measured.
**Test**: At each of L0–L11, compute H1 persistence diagram on 195-gene (or 4941-gene) embeddings. Compute pairwise bottleneck distance between adjacent layers (L_k, L_{k+1}). Plot and test if distance peaks at L3.
**Expected signal if true**: Bottleneck distance(L3, L4) > bottleneck distance of other adjacent pairs
**Null**: Random permutation of layer order produces equal bottleneck distances
**Value**: medium | **Cost**: medium

---

### H-M: GC Attractor Radius Contraction Across Layers
**Hypothesis**: The B-cell attractor "radius" (mean distance from GC-TFs to B-cell centroid) contracts monotonically from L0 to L11, and this contraction is faster than for T-cell TF distances to T-cell centroid.
**Test**: At each layer, compute mean Euclidean distance from {BATF, BACH2, PAX5} to B-cell centroid, and {RUNX3} to T-cell centroid. Track and compare rates of contraction.
**Expected signal if true**: GC-TF to B-cell centroid distance decreases steeply; T-cell TF shows flatter profile (pre-wired from L0)
**Null**: Both distance profiles are flat
**Value**: medium | **Cost**: low

---

### H-N: Embedding Norm Trajectory of GC-TFs vs Background
**Hypothesis**: GC-TFs (BATF, BACH2, PAX5) show a distinct embedding L2 norm trajectory across layers compared to background genes — either consistently high norm (indicating persistent representation strength) or increasing norm at L3 coinciding with attractor onset.
**Test**: At each layer L0–L11, record L2 norm of BATF, BACH2, PAX5 embeddings. Compare to mean norm of 195-gene background and to random same-size subsets.
**Expected signal if true**: GC-TF norms increase at L3+; background norms are flat
**Null**: GC-TF norms are indistinguishable from random same-size gene sets
**Value**: low | **Cost**: low

---

## Top 3 for Immediate Execution

### Slot 1: High-Probability Discovery Candidate
**H-A: TCR Circuit Attractor (Clean Retest)**

Rationale: The rank improvement signal (L0=92 → L11=41) already exists — it failed significance only because of the CD247 confound. This is a one-fix retest. The B-cell GC attractor is our positive control and we know what a real convergence signal looks like (rho ~ -0.95, p < 0.0001). If TCR circuit shows even a fraction of this, it's significant and extends the "delayed convergence" finding to a second lineage. Low cost, clear protocol.

**Implementation**: centroid = {CD3E, CD3D, TRAC, TRBC1}; circuit = {CD28, LCK, PTPRC}. Rank + null comparison at each layer. Same code pattern as iter_0040/iter_0041.

---

### Slot 2: High-Risk / High-Reward Candidate
**H-E: Cross-Lineage Circuit Unity Test**

Rationale: The GC circuit unity finding (iter_0043) — that PRDM1 repressor converges with BATF/BACH2 activators — could be either (a) a discovery showing scGPT doesn't encode regulatory directionality, or (b) a generic property of all circuits. Testing the T-cell AP1 circuit (JUN/JUNB/FOS as activators vs inhibitors CTLA4/PDCD1/NFATC1) will resolve this. If T-cell circuit shows separation (activators ≠ inhibitors), GC unity becomes a specific finding about the GC program. If T-cell also shows unity, it's a broader model property worth documenting. Either outcome is publishable.

**Implementation**: Track JUN, JUNB, FOS, FOSL2, CTLA4, PDCD1 ranks toward T-cell centroid {CD3D, CD3E, TRAC, TRBC1}. Same rank/null protocol.

---

### Slot 3: Cheap Broad-Screen Candidate
**H-C + H-K combined screen: Plasma Cell Divergence + Cross-Lineage Repulsion**

Rationale: Both tests are rank-tracking experiments using the existing cycle4_immune embeddings and established gene sets. They can be implemented in a single script with ~30 lines of new code each. H-C tests whether plasma TFs (IRF4, XBP1) share early trajectory with GC-TFs then diverge — directly extending the GC circuit unity story. H-K tests whether lineage separation involves active repulsion (B-cell genes moving away from T-cell centroid, and vice versa). Together they cost one script and one execution pass.

**Implementation**: Two rank tracking loops (one for plasma divergence, one for cross-lineage repulsion). 200 random null sets for each. Report rho and p-values.
