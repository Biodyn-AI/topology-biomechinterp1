# Brainstormer Hypothesis Roadmap — iter_0039

---

## Retire / Deprioritize

| Direction | Reason | Action |
|-----------|--------|--------|
| Plasma z goes absolutely positive | PRDM1 overlap artifact confirmed; replace with relative divergence framing | **retire_now** |
| NK/myeloid cell type screens (cycle1) | Zero in-vocab genes across 3+ attempts | **retire_now** (rescue only with Geneformer or cycle3/cycle4 vocab) |
| Directional drift cosine alignment (standalone) | Cosines ~0 across all top genes; confirmed non-directional | **retire_now** |
| STRING PPI Euclidean proximity | Retired in prior iterations; null result confirmed multiple times | **retired** |
| intrinsic_dimensionality SVD AUROC | Retired iter_0037 | **retired** |
| TRRUST/GO Euclidean proximity | Retired iter_0037 | **retired** |

---

## New Hypothesis Portfolio

### H-A: Full 12-Layer Scan — BATF/BACH2/BCL6/PAX5 Distance Trajectory
**Hypothesis**: The GC-TF cluster's proximity to B-cell centroid is not constant — it increases progressively, with a defined onset layer where the attractor first becomes detectable.
**Test**: In cycle4_immune, compute mean distance of {BATF, BACH2, BCL6, PAX5} to B-cell centroid at all 12 layers. Find earliest layer where z < −1.0. Compare to random 4-gene panels.
**Expected signal**: GC-TF z becomes significantly negative (< −1.0) starting at L4–L6 and deepens to L11; monotonic trend with an inflection point.
**Null/control**: 200 random 4-gene panels; plasma-TF {JCHAIN, SDC1} as divergent reference.
**Value**: high | **Cost**: low

---

### H-B: CD19+BLK Minimal 2-Gene Anchor — All-Layer Precision@10 Scan
**Hypothesis**: The CD19+BLK pair achieves precision@10 ≥ 0.2 for GC-TFs across all 12 layers with z consistently at ≥97th percentile, establishing a minimal sufficient anchor.
**Test**: cycle1 or cycle4_immune. Compute precision@10 for GC-TFs using {CD19, BLK} at all 12 layers. Compare to full 5-gene anchor and null. Find the best layer.
**Expected signal**: CD19+BLK yields equal or nearly equal precision@10 to full 5-gene panel; confirms these 2 genes are the geometric core of B-cell identity.
**Null/control**: 200 random 2-gene anchors; L0 as within-sample baseline.
**Value**: high | **Cost**: low

---

### H-C: IRF4 Plasma-Side Trajectory (cycle4_immune)
**Hypothesis**: IRF4 (plasma-TF bridging GC→plasma transition) is available in cycle4_immune and shows an intermediate z-score trajectory — less divergent than JCHAIN/SDC1 but more divergent than core GC-TFs — placing it geometrically between the two programs.
**Test**: Identify IRF4 in cycle4_immune vocab. Compute z-score distance to B-cell centroid at L0,L2,L5,L8,L11. Compare trajectory to BATF/BACH2 (GC reference) and JCHAIN/SDC1 (plasma reference).
**Expected signal**: IRF4 z-trajectory rises (more distal) at late layers but less steeply than JCHAIN; sits between GC and plasma bands.
**Null/control**: Random same-sized gene panels. GC-TF and plasma-TF trajectories as reference bands.
**Value**: high | **Cost**: low

---

### H-D: Geneformer B-Cell Precision@10 Cross-Model Replication
**Hypothesis**: Geneformer token embeddings show B-cell marker geometric clustering (precision@10 z > 3) at at least one transformer layer, confirming model-class generality.
**Test**: Load Geneformer embeddings; identify MS4A1, CD19, CD79A in its token vocabulary; compute precision@10 for GC-TFs (BATF, SPIB) at each layer using 200-permutation null. Compare peak layer to scGPT.
**Expected signal**: z > 3 at mid-to-late layers; T-cell markers show z < 1 (specificity).
**Null/control**: T-cell markers as negative; random gene sets (n=200).
**Value**: high | **Cost**: medium

---

### H-E: Geneformer GC-TF Proximity Cross-Model
**Hypothesis**: Geneformer reproduces BATF/BACH2 proximity to B-cell centroid (≥85th percentile) at equivalent transformer depths, confirming GC-TF geometry is not scGPT-specific.
**Test**: Geneformer embeddings, compute distance of BATF, BACH2, SPIB to B-cell centroid at each available layer. Rank vs all in-vocab genes.
**Expected signal**: ≥2/3 GC-TFs at ≥80th percentile. Peak at late layers.
**Null/control**: IRF4 as plasma reference; random gene percentiles.
**Value**: high | **Cost**: medium

---

### H-F: GC-Plasma Subspace Principal Angle (Extended 4-Gene GC Cluster)
**Hypothesis**: The GC-TF program (BATF, BACH2, BCL6, PAX5) and plasma-TF program (JCHAIN, SDC1, IRF4) occupy near-orthogonal subspaces at L11 in cycle4_immune, with principal angle > 80°.
**Test**: At L11 (cycle4_immune): compute SVD of GC-TF embedding matrix [4, 512] and plasma-TF [3, 512]. Compute principal angles between the two subspaces. Compare to 500 random gene-pair subspaces of same sizes.
**Expected signal**: Principal angle > 80° (consistent with 94° centroid-angle finding from iter_0038); z > 2 vs random pairs.
**Null/control**: 500 random gene-set pairs of same sizes.
**Value**: high | **Cost**: low

---

### H-G: LOO Ablation Replicated in cycle4_immune
**Hypothesis**: CD19 is also the critical anchor gene in cycle4_immune (different cell type composition than cycle1), confirming cross-data generality of the LOO finding.
**Test**: cycle4_immune at L2 and L11. Full anchor {MS4A1, CD79A, BLK} precision@10 for 4-gene GC-TF panel. LOO each anchor gene. Compare delta precision@10 per gene.
**Expected signal**: CD19 (if in cycle4 vocab) or BLK shows largest drop; confirms which gene is truly critical across cell-type compositions.
**Null/control**: 200 random 2-gene anchor panels.
**Value**: medium | **Cost**: low

---

### H-H: BLK Layer-Specificity Analysis — Why Does BLK Matter Only at L11?
**Hypothesis**: BLK encodes a layer-depth-specific aspect of B-cell identity that is progressive — its importance for GC-TF proximity increases monotonically with layer depth, reflecting increasing specialization.
**Test**: cycle1. Compute precision@10 for GC-TFs with and without BLK at all 12 layers (not just L2 and L11). Track delta precision@10 (BLK contribution) as a function of layer.
**Expected signal**: BLK delta precision@10 increases from ~0 at L2 toward +0.1 at L11; there is a transition layer where BLK becomes critical.
**Null/control**: MS4A1 and CD79A as redundant-gene controls (expect delta ≈ 0 at all layers).
**Value**: medium | **Cost**: low

---

### H-I: Non-B-Cell Immune Centroid Drift Targets
**Hypothesis**: Other immune cell centroids (T-cell, monocyte) also drift toward functionally meaningful gene neighborhoods at L11 — each cell type has a distinct attractor that matches its key regulators.
**Test**: cycle4_immune (295 in-vocab genes). Define T-cell anchor (CD3D, CD3E, CD8A), monocyte anchor (CD14, CSF1R). Compute centroid drift L0→L11. Find top-20 genes nearest to L11 centroid. Annotate.
**Expected signal**: T-cell centroid drifts toward TCR-related TFs (TBX21/T-bet, GATA3); monocyte toward myeloid TFs (SPI1/PU.1, CEBPB).
**Null/control**: Random same-size gene sets as anchor controls.
**Value**: medium | **Cost**: low

---

### H-J: Persistent Homology of GC-TF vs. Plasma-TF Cluster (L2 vs L11)
**Hypothesis**: The GC-TF cluster (BATF, BACH2, BCL6, PAX5) shows decreasing H0 persistence (tighter clustering) from L2→L11, while the plasma-TF cluster (JCHAIN, SDC1) shows increasing H0 persistence (dispersal) — PH as a fingerprint of differentiation trajectory geometry.
**Test**: cycle4_immune. Compute Vietoris-Rips H0 (connected components) for GC-TF 4-gene cloud and plasma-TF cloud at L2 and L11. Compare H0 death times (max pairwise distance within each cluster). Compare to 100 random 4-gene permutation clouds.
**Expected signal**: GC-TF max pairwise distance decreases L2→L11; plasma max pairwise distance increases.
**Null/control**: 100 random 4-gene gene sets; layer-matched.
**Value**: medium | **Cost**: medium

---

### H-K: Seed Robustness of CD19+BLK Critical-Anchor Finding
**Hypothesis**: The LOO result (CD19 critical at both layers) replicates across cycle4_immune seed replicates (seed43, seed44), ruling out seed-specific geometry as the explanation.
**Test**: If seed replicates exist for cycle4_immune, run same LOO protocol as H03 at L2 and L11. Compute precision@10 delta per removed gene across seeds.
**Expected signal**: CD19 removal shows ≥−0.05 delta at both layers across seeds; MS4A1/CD79A removal shows ~0 delta.
**Null/control**: Within-seed permutation nulls.
**Value**: medium | **Cost**: low (if data exists)

---

### H-L: Gene Embedding Velocity Clustering — Displacement Direction Communities
**Hypothesis**: The 195 in-vocab genes fall into biologically interpretable displacement-direction clusters (L0→L11), with at least one cluster enriched for GC-TFs and a separate one for plasma markers.
**Test**: cycle1. For all 195 in-vocab genes, compute L0→L11 displacement vector, unit-normalize. Cluster by cosine similarity of displacement direction (k-means, k=5–8). Annotate each cluster by GO/cell-type enrichment.
**Expected signal**: ≥1 cluster enriched for B-cell/GC markers; ≥1 cluster showing different trajectory (differentiation markers); displacement clustering is non-random (cluster silhouette > null).
**Null/control**: 200 displacement vectors from permuted coordinates; cluster silhouette vs null.
**Value**: medium | **Cost**: medium

---

### H-M: BATF/BACH2 Co-Regulation Subspace Alignment with SV2 Axis
**Hypothesis**: The known BATF-BACH2 heterodimer co-regulatory relationship (they form a functional complex) is reflected in their geometric proximity being even closer than the broader GC-TF cluster — their pairwise distance is at ≤1st percentile within the GC-TF cluster.
**Test**: cycle4_immune. Compute all pairwise distances within GC-TF cluster {BATF, BACH2, BCL6, PAX5} at L11. Test whether BATF-BACH2 distance is smaller than BATF-BCL6, BATF-PAX5, etc. Compare to known protein interaction partners vs non-partners.
**Expected signal**: BATF-BACH2 pairwise distance < BATF-BCL6 distance; functional complex partners are geometrically closer than non-complex regulatory partners.
**Null/control**: All pairwise distances within the 4-gene GC-TF set; 500 random 4-gene cluster distances.
**Value**: medium | **Cost**: low

---

## Top 3 for Immediate Execution

### 1. HIGH-PROBABILITY DISCOVERY: H-A + H-C (Full 12-Layer GC Attractor Scan + IRF4 Trajectory)
**Rationale**: These two experiments complete the picture of the GC-plasma geometry across layers and across the full differentiation axis. H-A determines when the GC-TF attractor first emerges (finding the "onset layer" is a mechanistically important result). H-C adds IRF4 as a bridge gene, testing whether the geometry encodes the GC→plasma transition order. Both run entirely on existing cycle4_immune data with trivial code changes to the H02 script. Expected result: strong positives that extend the current paper narrative with new mechanistic depth. Cost: very low.

### 2. HIGH-RISK/HIGH-REWARD: H-D + H-E (Geneformer Cross-Model Precision@10 + GC-TF Proximity)
**Rationale**: If BATF/BACH2/BCL6 are geometrically proximal to B-cell centroid in Geneformer, the paper moves from "scGPT observation" to "universal property of transformer gene FMs." This is the single highest-impact experiment remaining. Risk: Geneformer embeddings may require setup or may not be precomputed. Payoff if positive: the central paper claim becomes model-class general and publishable in stronger venues. Run both H-D and H-E together since they require the same data access step.

### 3. CHEAP BROAD-SCREEN: H-F + H-B (GC-Plasma Subspace Principal Angle + CD19+BLK Minimal Anchor)
**Rationale**: H-F (principal angle) requires only a 4-line SVD computation on already-loaded cycle4_immune embeddings and directly tests the orthogonality claim (extending the 94° finding from iter_0038 to the correct 4-gene GC cluster). H-B (CD19+BLK precision@10 all-layer) confirms the minimal anchor finding and identifies the "peak precision layer" — useful for identifying where the geometry is sharpest. Both have near-zero engineering overhead and will close two open questions from this and prior iterations. Cost: low.
