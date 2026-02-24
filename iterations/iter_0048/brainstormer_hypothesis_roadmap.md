# Brainstormer Hypothesis Roadmap: iter_0048 → iter_0049+

---

## Retire / Deprioritize

| Direction | Reason | Verdict |
|-----------|--------|---------|
| Lineage centroid cosine orthogonality | z≤0.26 vs null at all layers; SV1 explains >90% pre-correction, 18.6% post-correction. Both regimes make it uninformative. | `retire_now` |
| TCR circuit attractor | Tested twice (iter 44-45). n=2 in-vocab circuit genes after confound removal; power insufficient. Not rescuable without new data. | `retire_now` |
| Cross-lineage ID compression via 5-point TwoNN | Methodological incompatibility confirmed (iter 44). Results ~100× smaller than global TwoNN. | `retire_now` |
| SV1/SV2 cell-type identity (annotation-only approach) | Cell-type markers all score near zero on SV1 (three attempts). Dead end until full gene vocabulary is recovered. | `rescue_once_with_major_change` (needs full gene names) |
| H1 persistence on 209-gene set | Tested thoroughly (iters 3, 20-21). Signal confirmed. No new signal expected at this scale. | `deprioritize` (extend to 2039-gene set instead) |
| TRRUST co-regulation annotation → embedding proximity | Definitively retired (iter 35: AUROC=0.506, p=0.867). | `retire_now` |

---

## New Hypothesis Portfolio

### H-A: SV1 Axis Biological Identity (Full Gene Vocabulary)
**Hypothesis**: SV1 at L11, which explains 18.6% of active-gene variance, encodes a biologically coherent
program distinct from cell-type identity — likely housekeeping function, expression activity, or metabolic state.
**Test**: Recover gene index→name mapping for all 2039 active genes. Extract top-50 positive and top-50
negative SV1 loadings. Run GO term Fisher exact test (BP, CC, MF) and test correlation of SV1 loading
with STRING degree, expression level proxy (embedding norm), and known housekeeping gene lists (HPA).
**Expected signal**: Significant GO enrichment in top SV1 poles (e.g., metabolic, ribosomal, stress response).
**Null**: Random 50-gene sets from same 2039 active universe.
**Value**: HIGH | **Cost**: LOW
**Why now**: SV1 explains 18.6% of variance and is completely unidentified. This is the single largest
unexplained signal in the embedding.

---

### H-B: Cell-Type Subgroup Spectral Decay
**Hypothesis**: B-cell-specific genes undergo faster spectral compression than T-cell-specific or random
genes, indicating that cell-type attractor formation is visible in per-group spectral statistics.
**Test**: For B-cell markers (n~20 in-vocab), T-cell markers (n~25), myeloid markers (n~15), and
matched random sets (n~20 each), compute SVD eff_rank, SPR, and k90 per layer using the same
centered-submatrix approach as H02 iter_0048. Track compression ratios and test group differences.
**Expected signal**: B-cell subgroup eff_rank compression ratio > T-cell > random.
**Null**: Same-size random gene subsets (bootstrap distribution of compression ratios).
**Value**: HIGH | **Cost**: LOW
**Why now**: Direct extension of the spectral decay result with biological differentiation. Low implementation
cost since all machinery exists.

---

### H-C: Wasserstein Distance with Shuffle Null
**Hypothesis**: SWD between real consecutive layers is smaller than SWD between shuffled-embedding
consecutive layers, indicating that representation geometry is not random and constrains inter-layer transport.
**Test**: For each layer transition L_i→L_{i+1}, compute SWD for real embeddings. Also compute SWD
for shuffled embeddings (row-permute genes within each layer). Compare real vs shuffled SWD distributions.
Also: compare early-layer SWD (L0→L1: 0.552) vs late-layer (L10→L11: 0.291) after null normalization.
**Expected signal**: Real SWD < shuffled SWD; late-layer gap larger (more constrained structure).
**Null**: Row-permuted embedding null.
**Value**: MEDIUM | **Cost**: LOW
**Why now**: H03 iter_0048 is one experiment away from interpretability.

---

### H-D: SV1 Direction Stability Across Layers
**Hypothesis**: The dominant SV1 direction is established early and preserved across layers, such that
cosine similarity between SV1 vectors at successive layers is high from L3 onward.
**Test**: Compute SV1 unit vector at each of 12 layers. Compute pairwise cosine similarities (12×12 matrix).
Specifically: cosine(SV1_L0, SV1_Lk) for k=1..11, and consecutive-layer cosine(SV1_Li, SV1_{i+1}).
Track when SV1 direction "locks in" relative to L11.
**Expected signal**: SV1 direction becomes stable after L8 (when singular value peaks and drops).
**Null**: Random 512-dim unit vectors have expected cosine ≈ 0.
**Value**: MEDIUM | **Cost**: VERY LOW
**Why now**: The two-phase SV1 dynamics (peak at L8, drop to L11 while variance fraction grows)
suggest a direction-lock followed by amplitude contraction. This is mechanistically testable cheaply.

---

### H-E: SV2-SV4 PPI Geometry Replication in 2039-Gene Active Set
**Hypothesis**: SV2-SV4 axes in the 2039-gene embedding encode STRING PPI confidence as co-polarity
enrichment, replicating the iter_0015-017 findings in the larger, corrected gene set.
**Test**: SVD of 2039-gene centered embeddings at L8 (established as peak layer). For 295 named genes
with STRING annotations, compute SV2-SV4 scores and test whether STRING quintile gradient (Q1→Q5)
shows monotonic enrichment. Spearman ρ(quintile, mean z-score) per axis.
**Expected signal**: SV2 quintile gradient ρ ≈ 0.9-1.0 as in prior results; SV3/SV4 weaker.
**Null**: Shuffled gene labels on same SVD.
**Value**: HIGH | **Cost**: LOW
**Why now**: Prior PPI geometry work (iters 15-17) used 209 named genes. Replication in the corrected
2039-gene active set is needed for validity. If confirmed, substantially strengthens the primary PPI claim.

---

### H-F: kNN Graph Topology Transitions at L3 and L8
**Hypothesis**: The 2039-gene kNN graph undergoes discrete topological changes at L3 (GC attractor onset)
and L8 (SV1 singular value peak), detectable as changepoints in Fiedler value, community count, or
Betti number trajectories.
**Test**: Compute k=10 kNN graph on 2039-gene embeddings at each layer. Track: Fiedler value,
number of communities (greedy modularity), H0 count (connected components). Test for piecewise
linear breakpoints at each of 12 layers. Compare to prior 195-gene results.
**Expected signal**: Slope change in Fiedler value at L3 and/or L8.
**Null**: Smooth monotone trend (no breakpoints).
**Value**: MEDIUM | **Cost**: LOW

---

### H-G: Persistent Homology on Full 2039-Gene Active Set
**Hypothesis**: H1 persistence on 2039 active genes reveals more complex loop structure (more H1 features)
than the 209-gene subset, with loop density decaying monotonically at a different rate.
**Test**: Run Ripser (H0+H1) on PCA(30) projections of 2039 active genes at L0, L5, L8, L11.
Feature-shuffle null (20 replicates per layer). Track H1 count, mean lifetime, and persistence diagram.
Compare lifetime decay trajectory to 209-gene results (iter_0021: ρ=-0.916).
**Expected signal**: Higher H1 count, same monotone decay direction; possibly different inflection points.
**Null**: Feature-shuffle replicates.
**Value**: MEDIUM | **Cost**: MEDIUM

---

### H-H: Zero-Norm Gene Biology Characterization
**Hypothesis**: The 2902 zero-norm (excluded) genes are systematically enriched for tissue-specific,
developmental, and non-immune biological programs, forming a structurally coherent out-of-vocabulary set.
**Test**: Map zero-norm gene indices to names using scGPT vocabulary. Run GO enrichment (Fisher exact)
for BP, CC, MF on the full 2902-gene zero-norm set vs 2039 active genes. Test STRING connectivity:
are zero-norm genes less connected to each other in STRING than active genes?
**Expected signal**: Zero-norm genes enriched for developmental/tissue-specific GO terms; lower STRING
density than active set.
**Null**: Random 2902-gene sets from combined vocabulary.
**Value**: MEDIUM | **Cost**: MEDIUM

---

### H-I: SV1 Anisotropy and Two-Phase Dynamics Mechanism
**Hypothesis**: The SV1/SV2 singular value ratio (anisotropy) peaks at L8, coinciding with maximum
spectral concentration before global embedding contraction; this reflects a mechanistic phase transition
in representation organization.
**Test**: Compute per-layer SV1/SV2 ratio, SV1/SV3 ratio, and SV1/mean(SVs) anisotropy. Identify
the layer that maximizes each metric. Test whether anisotropy layer peaks are concentrated around L8.
Also compute Frobenius-normalized SV1 (SV1/Frob) to separate amplitude from anisotropy.
**Expected signal**: Anisotropy peaks at L8; Frobenius-normalized SV1 increases monotonically (direction
sharpens even as amplitude drops L8→L11).
**Null**: Feature-shuffled embeddings show no peak.
**Value**: HIGH | **Cost**: LOW

---

### H-J: GC Attractor Spatial Statistics
**Hypothesis**: The GC-TF triangle (PAX5, BATF, BACH2) and B-cell markers occupy a statistically
compact, geometrically isolated cluster in 2039-gene embedding space, measurable as a local
density anomaly relative to the global manifold.
**Test**: For the GC-TF triangle {PAX5, BATF, BACH2} + B-cell markers, compute: mean within-group
distance, mean distance to 50 nearest non-group neighbors, local density ratio (within vs background).
Test at each layer; compare to 1000 random size-matched groups.
**Expected signal**: GC attractor density ratio > null at all layers; increases with depth.
**Null**: Random size-matched groups from 2039-gene universe.
**Value**: MEDIUM | **Cost**: LOW

---

### H-K: Cross-Model SWD Comparison (Geneformer vs scGPT)
**Hypothesis**: For 295 named genes shared between Geneformer and scGPT, the SWD between
cross-model embedding spaces is smaller for STRING-interacting pairs than random pairs,
reflecting shared inductive bias toward biologically proximate representations.
**Test**: Compute pairwise distances between Geneformer token embeddings and scGPT L8 embeddings
for 295 annotated genes. Compute AUROC(STRING pairs have smaller cross-model distance than random).
**Expected signal**: AUROC > 0.55 for cross-model STRING proximity.
**Null**: Shuffled cross-model gene assignments.
**Value**: HIGH-RISK/HIGH-REWARD | **Cost**: MEDIUM

---

### H-L: Hub Gene Centrality Evolution in 2039-Gene Set
**Hypothesis**: Hub genes (high STRING degree among annotated genes) move toward the global centroid
with increasing layer depth, replicating iter_0028 finding in the corrected active-gene set.
**Test**: For 295 named genes with STRING edges: compute distance to 2039-gene centroid at each layer.
Spearman ρ(STRING degree, -dist_to_centroid) per layer, with 500-permutation null.
**Expected signal**: ρ increases L0→L11 (hub genes converge toward centroid); peak ρ > 0.2.
**Null**: Random gene degree permutations.
**Value**: MEDIUM | **Cost**: LOW

---

### H-M: Attention Head Biological Specialization
**Hypothesis**: Individual attention heads specialize for different biological information types
(PPI proximity vs. regulatory co-occurrence vs. cell-type identity), with specialization
increasing at deeper layers.
**Test**: For each of 8 heads at L0, L5, L11: compute co-attendance AUROC for STRING pairs,
TRRUST pairs, and B-cell marker pairs. Measure head specialization entropy (variance of per-head AUROC
across databases).
**Expected signal**: Head AUROC variance increases with depth; specific heads achieve AUROC > 0.65.
**Null**: Random head assignments.
**Value**: HIGH | **Cost**: HIGH (needs attention weight extraction)

---

## Top 3 for Immediate Execution

### 1. High-Probability Discovery Candidate: H-A (SV1 Axis Biological Identity)
SV1 at L11 explains 18.6% of active-gene variance and is completely unidentified — the largest open
gap in the current findings. The fix is cheap: recover gene index→name mapping for 2039 active genes,
extract top-50 SV1 loadings, run GO enrichment and correlation with gene-level properties. This has
high probability of identifying a biologically coherent axis (metabolic, housekeeping, or cell-cycle)
that completes the spectral decay story.
**Implementation**: Load gene vocabulary from scGPT tokenizer/vocab file; SVD of 2039-gene embeddings
at L11; Fisher exact GO test on top-50 loading poles; Spearman ρ(SV1 loading, STRING degree/norm/etc.)

---

### 2. High-Risk / High-Reward Candidate: H-E (SV2-SV4 PPI Geometry in 2039-Gene Set)
The core PPI geometry finding (SV2-SV4 encode STRING interaction confidence as co-polarity enrichment)
was established on 209 named genes. If it replicates in the 2039-gene active set with 295 annotated
genes, it dramatically increases the power and credibility of the central PPI-encoding claim.
If it fails to replicate (lower enrichment or no gradient), this refines what the axes encode.
Either outcome is publication-informative.
**Implementation**: SVD of 2039-gene centered embeddings at L8; compute SV2-SV4 scores for 295 named
genes; STRING quintile gradient and Spearman ρ per axis; shuffle null.

---

### 3. Cheap Broad-Screen Candidate: H-D + H-C Combined (Direction Stability + SWD Null)
Both can be computed from existing embeddings in a single script:
- SV1 direction stability: cosine(SV1_Li, SV1_L11) per layer — tells when the dominant axis locks in
- SWD with shuffle null: resolves whether transport distances reflect structure

Together these characterize the two-phase spectral dynamics (consolidation L0→L8, contraction L8→L11)
with minimal new computation. Strong expected signal for direction stability (should lock near L8);
SWD null comparison should immediately disambiguate whether the transport distances are structured.
**Implementation**: Extend iter_0048 script with: (a) cosine similarity between successive SV1 vectors;
(b) SWD between shuffled-row embeddings at same layer pairs.
