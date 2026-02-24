# Brainstormer Hypothesis Roadmap — iter_0040

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| GC-plasma subspace convergence | rho=-0.098, p=0.762 — null in 2 tests | **retire_now** |
| BCL6 as GC attractor member | rank 750-1500 throughout all 12 layers | **retire_now** (revise paper claims) |
| Geneformer INPUT token cross-model | floor-level signal, CD19 absent, olfactory receptor neighbors | **retire_now** (reopen only with layer-wise engineering) |
| Generic cycle1 vocab limitations | superseded by cycle4_immune with 2039 genes | **deprioritize** |

---

## New Hypothesis Portfolio

### H-A: L3 Attention Structure vs L2 — Norm and Rank-Topology Change
**Hypothesis**: The L3 transition marks a qualitative change in embedding geometry, detectable by (a) embedding norm inflection and (b) kNN graph topology (connectivity or clique structure) at L3 vs L2.
**Test**: For cycle4_immune, compute per-gene L2 norm at each layer; compute kNN graph (k=10) at L2 and L3; compare number of connected components, average clustering coefficient, and modularity. Permutation null for kNN metrics.
**Expected signal if true**: Norm inflection at L3; kNN graph shows sharply increased clustering of GC-TFs at L3 vs L2.
**Null/control**: Gradual change across all layers, no L3 discontinuity.
**Value**: high | **Cost**: low

---

### H-B: BCL6 Neighborhood Identification — What Program Does BCL6 Belong To?
**Hypothesis**: BCL6 (rank ~900 near B-cell centroid) actually clusters with memory B-cell or GC B-cell markers (not TFs), meaning it represents B-cell state rather than TF regulatory identity.
**Test**: cycle4_immune L11. Find top-20 nearest genes to BCL6 embedding. Annotate top-20 with cell-type markers (TRRUST/manual). Compare BCL6 neighborhood to neighborhoods of BATF, BACH2, PAX5.
**Expected signal if true**: BCL6 neighbors include memory B or activated B markers (e.g., AICDA, CD38, MKI67) — different from GC-TF program.
**Null/control**: BCL6 neighbors are random/incoherent.
**Value**: high | **Cost**: low

---

### H-C: T-Cell Lineage Attractor (parallel test)
**Hypothesis**: scGPT encodes a T-cell TF attractor analogous to the GC-TF attractor: FOXP3, GATA3, TBX21 converge toward a T-cell centroid (CD3D, CD3E, TRAC) across layers with a similar L3-range onset.
**Test**: cycle4_immune. T-cell anchor = CD3D, CD3E, TRAC. TF panel = FOXP3, GATA3, TBX21, RORC. Rank scan L0-L11, Spearman rho. Onset layer detection.
**Expected signal if true**: TF ranks decrease monotonically; onset at L3-L5 (possibly same or slightly later than GC attractor).
**Null/control**: TF ranks flat or noisy across layers.
**Value**: high | **Cost**: low

---

### H-D: PAX5 Pre-Wiring at L0 — Neighborhood Composition
**Hypothesis**: PAX5's anomalously low rank (36) at L0 reflects scGPT input encoding of PAX5's known role as master B-cell TF — its L0 neighbors should include other B-cell regulators (EBF1, IKZF1, POU2F2) rather than random genes.
**Test**: cycle4_immune L0. Find top-50 nearest genes to PAX5. Annotate with B-cell TF/marker ontology. Compare to L11 top-50. Fisher test for B-cell TF enrichment at L0 vs null.
**Expected signal if true**: L0 PAX5 neighbors enriched for B-cell TFs (EBF1, IKZF1, CD19, etc.) — pre-existing regulatory co-embedding.
**Null/control**: Random gene set of size 50 from cycle4_immune.
**Value**: high | **Cost**: low

---

### H-E: Persistent Homology on Gene Embedding Cloud — L3 Topological Transition
**Hypothesis**: The gene embedding point cloud undergoes a topological change at L3 detectable via persistent homology (H0 connected components, H1 loops). The GC-TF cluster becomes a persistent feature at L3.
**Test**: Use ripser or gudhi on cycle4_immune embeddings at each layer (subsample to 500 genes for tractability). Compute Betti-0 and Betti-1 persistence diagrams. Track lifetime of GC-TF cluster as a persistent feature. Compare bottleneck distance between consecutive layers.
**Expected signal if true**: Bottleneck distance spike at L2→L3 transition; GC-TF subcloud appears as a distinct Betti-0 component at L3.
**Null/control**: Random permutation of gene labels within each layer.
**Value**: high | **Cost**: medium

---

### H-F: Intrinsic Dimensionality Change at L3
**Hypothesis**: The intrinsic dimensionality of the gene embedding manifold decreases at L3, reflecting specialization of the representation into lineage-specific subspaces.
**Test**: cycle4_immune. Estimate local intrinsic dimension using TwoNN or MLE estimator at each layer. Also compute explained variance (PCA) for the top-10 PCs per layer. Spearman rho of dimensionality vs layer.
**Expected signal if true**: Dimension drops or inflects at L3; variance more concentrated in fewer PCs.
**Null/control**: Shuffled embeddings at each layer.
**Value**: medium | **Cost**: low

---

### H-G: IRF4 Plasma-Bridge Test in Cycle1
**Hypothesis**: IRF4 (OOV in cycle4) is present in cycle1 and shows a distinct trajectory from GC-TFs — converging toward the B-cell centroid less strongly than BATF/BACH2, consistent with its dual GC→plasma bridge role.
**Test**: cycle1 embeddings (4803 genes). Check IRF4 in-vocab. If present: rank IRF4 vs BATF, BACH2, SPIB at each layer relative to B-cell centroid (MS4A1, CD19, CD79A, BLK, PRDM1). Compare rank trajectories.
**Expected signal if true**: IRF4 rank trajectory diverges from BATF/SPIB — stays more intermediate, not converging as tightly.
**Null/control**: Genes from plasma panel (JCHAIN, SDC1) show flat or increasing ranks.
**Value**: medium | **Cost**: low

---

### H-H: Myeloid/Monocyte TF Attractor Screen
**Hypothesis**: Myeloid-specific TFs (SPI1/PU.1, CEBPA, CEBPB) converge toward a monocyte centroid (CD14, FCGR3A, LYZ) across layers — analogous to GC attractor but for myeloid lineage.
**Test**: cycle4_immune. Anchor = CD14, FCGR3A, LYZ. TF panel = SPI1, CEBPA, CEBPB, IRF8. Rank scan L0-L11, Spearman rho. p@10 with null.
**Expected signal if true**: SPI1 (master myeloid TF) shows monotonic rank decrease; onset possibly L3-L5.
**Null/control**: Randomly chosen TFs not associated with myeloid lineage.
**Value**: high | **Cost**: low

---

### H-I: GC-TF Inter-Pair Distance Contraction — Convergence Geometry
**Hypothesis**: The three confirmed GC-TFs (BATF, BACH2, PAX5) not only converge toward the B-cell centroid but also converge toward each other — pairwise L2 distances decrease across layers, with steepest contraction at L3-L6.
**Test**: cycle4_immune. Compute pairwise L2 distances among {BATF, BACH2, PAX5} at each layer. Spearman rho vs layer. Compare to pairwise distances among random TF triplets (null).
**Expected signal if true**: All 3 pairwise distances decrease; inter-TF contraction onset matches L3.
**Null/control**: Random gene triplets show no systematic distance decrease.
**Value**: medium | **Cost**: low

---

### H-J: GC Attractor Replication in Cycle1 — Cross-Dataset Robustness
**Hypothesis**: The GC attractor (onset L3, BATF/BACH2/PAX5 convergence) replicates in cycle1 embeddings, demonstrating the finding is not cycle4-specific.
**Test**: cycle1 (layers available? — check shape). B-cell anchor = MS4A1, CD19, CD79A, BLK. GC-TF panel = BATF, BACH2, SPIB (PAX5 was OOV in cycle1 — use SPIB instead). Rank scan L0-L11. Onset detection.
**Expected signal if true**: Similar rank trajectory and L3 onset in cycle1.
**Null/control**: Different onset layer (> L5) or flat trajectory would indicate cycle4-specific effect.
**Value**: high | **Cost**: low

---

### H-K: Precision@K Sweep (K=5,10,20,50) for Attractor Robustness
**Hypothesis**: The GC-TF attractor signal (p@10 = 0.1 at L11) is robust across different K thresholds, indicating a genuine cluster rather than an edge-case proximity.
**Test**: cycle4_immune at L11. For GC-TFs {BATF, BACH2, PAX5}, compute precision@K for K ∈ {5,10,20,50} against CD19+BLK anchor. Report fraction of GC-TFs in top-K at each K. Null: 200 random 3-gene panels.
**Expected signal if true**: GC-TFs appear in top-20 at significant enrichment; p@50 > p@10 proportionally.
**Null/control**: Random 3-gene panels show flat p@K across K values.
**Value**: medium | **Cost**: low

---

### H-L: NK Cell Attractor Screen (EOMES, TBX21, GZMB)
**Hypothesis**: NK cell TFs (EOMES, TBX21) and effector markers (GZMB) converge toward NK centroid (NCAM1/CD56, KLRD1, GNLY) across layers.
**Test**: cycle4_immune. NK anchor = NCAM1, KLRD1, GNLY. TF/effector panel = EOMES, TBX21, GZMB, PRF1. Rank scan L0-L11. Spearman rho.
**Expected signal if true**: EOMES and TBX21 show rank decrease; onset possibly L3-L5 (mirroring GC attractor).
**Null/control**: Random gene sets vs NK centroid.
**Value**: medium | **Cost**: low

---

### H-M: Layer-Wise UMAP Topology Visualization (Qualitative Validation)
**Hypothesis**: UMAP embedding of gene vectors at each layer, colored by cell-type annotation, shows progressively clearer lineage clustering from L0 to L11, with visible GC-TF / B-cell co-localization emerging by L3.
**Test**: cycle4_immune, all 2039 genes, UMAP at L0, L3, L6, L11. Highlight GC-TFs, B-cell markers, plasma markers, T-cell markers. Qualitative assessment plus silhouette score for lineage clusters.
**Expected signal if true**: GC-TFs and B-cell markers co-localize at L3+; silhouette score for lineage clusters increases from L0 to L11.
**Null/control**: Random gene color assignment.
**Value**: medium | **Cost**: low (visualization only)

---

## Top 3 for Immediate Execution

### #1 — High-Probability Discovery Candidate
**H-C: T-Cell and H-H: Myeloid Lineage Attractor Screen (bundled)**
- Run both in a single script: T-cell (FOXP3/GATA3/TBX21 → CD3D/CD3E/TRAC) and myeloid (SPI1/CEBPA/CEBPB → CD14/FCGR3A/LYZ) rank scans over all 12 layers.
- If either shows rho < -0.7 with L3-ish onset, we have a generalizable principle (not B-cell-specific).
- This would be the most impactful next finding: "scGPT encodes lineage-specific TF attractors for multiple immune cell types."
- Cost: low (same pipeline as H01).

### #2 — High-Risk/High-Reward Candidate
**H-E: Persistent Homology Topological Transition at L3**
- First use of topological data analysis on the scGPT gene embedding cloud.
- If the GC-TF cluster emerges as a distinct topological feature (Betti-0 component) precisely at L3, this provides a new mathematical characterization of the attractor onset.
- Risk: PH computation on 512-dim embeddings may require careful distance threshold selection; subsample may miss signal.
- Use ripser, subsample 500 genes, run all 12 layers.

### #3 — Cheap Broad-Screen Candidate
**H-B + H-D bundled: BCL6 neighborhood + PAX5 L0 neighborhood (single query)**
- Two nearest-neighbor lookups at L0 and L11 in the already-loaded cycle4_immune embeddings.
- BCL6: identify its actual biological neighborhood (resolves the mystery of its non-convergence).
- PAX5: confirm whether L0 pre-wiring neighbors are B-cell regulators.
- Cost: near-zero compute, high biological interpretability payoff.

---

## Strategic Notes

1. **Generalization is the next frontier.** We have a highly characterized GC-TF attractor in B cells. The single highest-value next experiment is showing this pattern is a general principle across immune lineages (T-cell, myeloid). If confirmed, this reframes the paper as a general discovery about how scGPT encodes cell-type regulatory identity.

2. **L3 mechanism is unresolved.** We know WHEN but not WHY. The attention structure / norm / topology probes (H-A, H-E, H-F) target this gap. Even one positive result here would substantially deepen the paper.

3. **BCL6 exclusion creates a mystery worth solving.** BCL6 is a canonical GC B-cell TF in biology yet appears at rank ~900 throughout. Resolving this (H-B) either reveals a limitation of the model or a nuanced biological finding (BCL6 encodes B-cell state identity differently from regulatory TF identity).
