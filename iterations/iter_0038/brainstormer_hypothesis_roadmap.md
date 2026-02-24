# Brainstormer Hypothesis Roadmap — iter_0038

---

## Retire / Deprioritize

| Direction | Reason | Action |
|-----------|--------|--------|
| NK/myeloid specificity screen (cycle1 vocab) | Zero in-vocab representation across 3 attempts | **retire_now** (rescue only with cycle2/Geneformer vocab) |
| Directional drift cosine alignment (standalone) | No alignment found; drift is directionally random w.r.t. tested axes | **deprioritize** — fold into drift-target investigation |
| intrinsic_dimensionality SVD AUROC | Retired in iter_0037 | **retired** |
| TRRUST/GO proximity | Retired in iter_0037 | **retired** |
| STRING PPI Euclidean | Retired in iter_0037 | **retired** |
| T-cell as primary clustering target | Retired in iter_0037 | **retired** |

---

## New Hypothesis Portfolio

### H-A: GC-Plasma Subspace Principal Angle (SVD Orthogonality Test)
**Hypothesis**: The GC-TF program and plasma-TF program occupy near-orthogonal subspaces in scGPT layer 11, quantifiable as a principal angle >80° between the 2D PCA planes of each gene set.
**Test**: At L11, compute 2D PCA of GC-TFs (BATF, SPIB, BACH2) and plasma-TFs (IRF4, PRDM1). Compute the principal angle between the two planes. Compare to 500 permuted random 2-gene/3-gene subsets.
**Expected signal**: Principal angle >80° (consistent with 94° centroid angle); z > 2 vs random pairs.
**Null/control**: 500 random gene-set pairs of same sizes.
**Value**: high | **Cost**: low

---

### H-B: B-Cell Centroid Drift Target Identification
**Hypothesis**: The B-cell centroid drift direction (L0→L11, 26.4 units) points toward a semantically meaningful gene or gene cluster, identifiable as the nearest named gene to the drift endpoint.
**Test**: Compute drift endpoint = L0 centroid + displacement vector. Find top-10 named genes nearest to that endpoint in L11 embedding space. Annotate with cell-type/function.
**Expected signal**: Drift endpoint neighbors are enriched for antigen-presentation (HLA-DR, CD74) or memory B-cell markers (CR2/CD21, BANK1) — reflecting the "default B-cell resting program" the model encodes as layer depth increases.
**Null/control**: Drift endpoint distance to random gene set centroids.
**Value**: high | **Cost**: low

---

### H-C: Plasma-TF z-Trajectory Replication with Extended Plasma Panel
**Hypothesis**: The plasma-TF z-score shift (L0: −0.62 → L11: +0.55) replicates and strengthens when additional plasma markers (XBP1, MZB1, JCHAIN, SDC1) are included from cycle2_maxgenes1024 vocabulary.
**Test**: Load cycle2_maxgenes1024 embeddings. Identify XBP1, MZB1, JCHAIN, SDC1 in vocab. Compute mean distance to B-cell centroid (using plasma-exclusive anchor: MS4A1, CD19, CD79A, BLK only) at L0,L2,L5,L8,L11. z-score vs null.
**Expected signal**: Plasma z-trajectory still shifts from negative to positive; larger panel resolves the PRDM1 overlap artifact.
**Null/control**: All in-vocab named genes as null. GC-TF trajectory as positive reference.
**Value**: high | **Cost**: low

---

### H-D: Geneformer B-Cell Precision@10 Cross-Model Replication
**Hypothesis**: Geneformer token embeddings show B-cell marker geometric clustering (precision@10 z>3) at at least one transformer layer, confirming the signal is model-class general.
**Test**: Load Geneformer embeddings; identify B-cell markers (MS4A1, CD19, CD79A) in its token vocabulary; compute precision@10 and z-score vs all-gene null at each layer. Compare peak layer to scGPT.
**Expected signal**: z>3 in middle layers; T-cell markers show z<1 (specificity preserved).
**Null/control**: T-cell markers as negative; random gene sets (n=200).
**Value**: high | **Cost**: medium

---

### H-E: Geneformer GC-TF Proximity Cross-Model
**Hypothesis**: Geneformer reproduces BATF/SPIB/BACH2 proximity to B-cell centroid (≥85th pctile) at equivalent transformer depths, confirming GC-TF geometry is not scGPT-specific.
**Test**: In Geneformer, compute distance of BATF, SPIB, BACH2 to B-cell centroid at each layer. Rank vs all in-vocab null. Compare pctile ranks to scGPT (86-96th).
**Expected signal**: At least 2/3 GC-TFs at ≥80th pctile in Geneformer.
**Null/control**: IRF4 as cross-model plasma reference; random gene percentiles.
**Value**: high | **Cost**: medium

---

### H-F: PAX5 / EBF1 / BCL6 Proximity in Extended Vocabulary
**Hypothesis**: PAX5 and EBF1 (master B-cell TFs absent in cycle1 vocab) are present in cycle2_maxgenes1024 and cluster at ≥90th pctile proximity to B-cell centroid, extending the GC-TF cluster finding.
**Test**: Load cycle2_maxgenes1024 embeddings; scan for PAX5, EBF1, BCL6, IKZF1, IKZF3 in vocab; compute centroid distances; rank vs all-gene null.
**Expected signal**: PAX5 ≥90th, EBF1 ≥90th; BCL6 ≥80th. Together with BATF/SPIB/BACH2, this would constitute a 5-6 gene GC-TF cluster.
**Null/control**: All in-vocab named genes; IRF4 and PRDM1 as plasma references.
**Value**: high | **Cost**: low

---

### H-G: Layer-Resolved Principal Angle Scan Across Cell-Type Pairs
**Hypothesis**: The near-orthogonal GC-plasma relationship at L11 (94°) is specific to that differentiation pair, not a generic property of any two cell-type-associated gene sets.
**Test**: Compute centroid-angle between 5+ cell-type gene pairs (e.g., B/T, B/NK, GC/memory, naive/plasma) at L11. Test whether GC/plasma angle is an outlier vs other pairs.
**Expected signal**: Most pairs show <70° angles; GC/plasma is ≥85°, making it a specific geometric feature.
**Null/control**: 200 random gene-set pairs of matched size; inter-pair angle distribution.
**Value**: medium | **Cost**: low

---

### H-H: Gene-Level Embedding Velocity Clustering
**Hypothesis**: Genes can be clustered by their displacement *direction* from L0→L11, and the resulting direction clusters will correspond to biological function modules (GC program, plasma program, antigen presentation, etc.).
**Test**: For all 195 in-vocab genes, compute displacement vector L0→L11 (unit-normalized). Cluster by cosine similarity of displacement vectors (k-means or spectral, k=5-8). Annotate clusters by gene function.
**Expected signal**: At least one cluster is enriched for B-cell identity markers; another for differentiation markers; gene set annotations should align with known biology.
**Null/control**: Displacement direction distributions from random 195-gene permutations.
**Value**: medium | **Cost**: medium

---

### H-I: B-Cell Cluster Ablation Stability (Leave-One-Out)
**Hypothesis**: The B-cell precision@10 signal (z=7.55) is robust to removal of any single marker gene, confirming it is not driven by one dominant gene.
**Test**: For each of 5 B-cell markers (MS4A1, CD19, CD79A, BLK, PRDM1), compute precision@10 for the 4-gene panel. Compare z-scores to full 5-gene panel. Identify most/least critical gene.
**Expected signal**: All 4-gene LOO panels maintain z>4; no single gene is solely responsible.
**Null/control**: Random 4-gene sets; T-cell 4-gene LOO as negative.
**Value**: medium | **Cost**: low

---

### H-J: Plasma-TF Trajectory Comparison Across Differentiation-Proximal vs Distal Genes
**Hypothesis**: Genes co-expressed with plasma cells in scRNA-seq (secretory pathway: SEC61A1, HSPA5, SSR1) show the same z-score divergence pattern (becoming farther from B-cell centroid across layers) as IRF4/PRDM1, while B-cell identity genes remain close.
**Test**: Identify 3-5 secretory pathway genes in cycle1 or cycle2 vocab. Track their mean z-score vs B-cell centroid across layers. Compare trajectory shape to GC-TF (stable) and plasma-TF (diverging) profiles.
**Expected signal**: Secretory genes show positive-drifting z, confirming the plasma-TF trajectory reflects a general differentiation-axis geometry rather than TF-specific effect.
**Null/control**: GC-TF trajectory as stable reference; random gene trajectories.
**Value**: medium | **Cost**: low

---

### H-K: Persistent Homology of Full 195-Gene Embedding at L2 vs L11
**Hypothesis**: The H0/H1 persistent homology of the full 195-gene cloud changes from L2 to L11, with fewer connected components (fewer H0 generators dying early) and possible new H1 loops at L11, reflecting global geometric reorganization.
**Test**: Compute Vietoris-Rips PH (H0, H1) of all 195-gene embeddings at L2 and L11. Compare Betti numbers, max persistence, and total persistence across layers.
**Expected signal**: L11 shows lower max H0 persistence (genes pulled into tighter clusters) but possibly higher H1 persistence (loop structures from differentiation axes).
**Null/control**: 50 random 195-gene permutation clouds; layer-matched.
**Value**: medium | **Cost**: medium

---

### H-L: Memory B-Cell Marker Proximity at Late Layers
**Hypothesis**: Memory B-cell markers (CR2/CD21, BANK1, IGHM) occupy an intermediate position in L11 embedding space — closer to B-cell centroid than plasma-TFs but farther than GC-TFs — reflecting memory as a partial-differentiation state.
**Test**: Identify CR2, BANK1, IGHM in cycle1 or cycle2 vocab. Compute distance to B-cell centroid at L11. Rank vs null; compare z-score to GC-TF (z≈−1.5) and plasma-TF (z≈+0.55).
**Expected signal**: Memory markers at z between −0.5 and −1.3 — intermediate positioning.
**Null/control**: All 195 in-vocab genes as null.
**Value**: medium | **Cost**: low

---

## Top 3 for Immediate Execution

### 1. HIGH-PROBABILITY DISCOVERY: H-B + H-C (Drift Target + Extended Plasma Panel)
**Rationale**: Both are cheap (low-cost), directly extend strong positive findings from iter_0038, and require no new data sources — only computation on existing embeddings (cycle1 for H-B, cycle2 for H-C). H-B answers the open question of "what direction is the B-cell centroid actually drifting?" which is needed to complete the H02 story. H-C removes the PRDM1 overlap artifact from H01. Can be bundled into a single script run.

### 2. HIGH-RISK/HIGH-REWARD: H-D + H-E (Geneformer Cross-Model)
**Rationale**: If the B-cell specificity (z=7.55) and GC-TF geometry replicate in Geneformer, this elevates the paper from "scGPT-specific observation" to "universal property of transformer gene foundation models." This is the single highest-impact pending experiment. Risk: Geneformer embeddings may not be precomputed; requires setup. But the payoff — confirming model-class generality — is critical for the central paper claim.

### 3. CHEAP BROAD-SCREEN: H-F + H-I (Extended Vocab + LOO Ablation)
**Rationale**: H-F (PAX5/EBF1 in cycle2) is a direct vocab-switch away from a meaningful biological extension; H-I (LOO ablation) is a 5-iteration loop on existing code with near-zero engineering overhead. Together they provide: (a) a broader GC-TF cluster covering master regulators, and (b) robustness confirmation for the precision@10 metric. Both can be run in the same iteration as H-B/H-C.
