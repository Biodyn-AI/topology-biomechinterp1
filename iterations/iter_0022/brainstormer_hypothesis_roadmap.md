# Brainstormer Hypothesis Roadmap — iter_0022

---

## Retire / Deprioritize

| Direction | Reason | Disposition |
|-----------|--------|-------------|
| Quintile-binned STRING gradient | Replaced by continuous Spearman; no additional information | `retire_now` |
| Standalone SVD direction counting | Absorbed into co-polarity result; publish-ready | `retire_now` |
| Transitivity (global) negative result | Well-characterized negative from iter_0002; no rescue signal | `retire_now` |
| NULL sensitivity / permutation tests as standalone hypotheses | Covered in every test; not a standalone hypothesis | `retire_now` |
| Proxy persistence entropy (histogram entropy) | p=0.16; proxy is invalid; blocked on ripser data | `rescue_once_with_major_change` — only if full ripser rerun is scheduled |
| TwoNN ID as primary test | Signal found (iter_0004), mild; not worth more standalone tests unless conditioned on cell type | `deprioritize` |

---

## New Hypothesis Portfolio

### H-A: Cell-type marker expansion (5–6 cell types)
**Hypothesis**: The AUROC=0.853 cell-type separation signal extends to a 5-6 cell-type vocabulary including macrophage, NK, myeloid, and endothelial markers.
**Test**: Add macrophage (ALOX5, PTPRC, CSF1R), NK (NCAM1, KLRD1), myeloid (CD14, LYZ, S100A8), endothelial (PECAM1, CDH5) markers if present among 209 named genes. Recompute within- vs cross-type AUROC for expanded vocabulary. Permutation null: 500 shuffles.
**Expected signal**: AUROC remains ≥0.80; expanding the test increases statistical power and confirms generality.
**Null/control**: Non-marker named genes shuffled into same-sized groups should give AUROC≈0.50.
**Value**: high | **Cost**: low

### H-B: Cell-type contamination control (non-marker genes)
**Hypothesis**: The cell-type separation is marker-specific, not a general property of all 209 named genes.
**Test**: Take remaining ~195 named genes (non-cell-type-markers). Sample 3 random groups of sizes matching T/B/fibroblast (7/4/3). Compute within-group vs cross-group AUROC × 100 bootstrap samples. Compare distribution to the observed marker AUROC=0.853.
**Expected signal**: Non-marker AUROC distribution centers around 0.50; marker AUROC is an outlier (>5 SD above).
**Null/control**: The non-marker bootstraps are the null.
**Value**: high | **Cost**: low

### H-C: Cross-model cell-type replication (Geneformer word embeddings)
**Hypothesis**: Geneformer token/word embeddings show the same cell-type marker cluster structure as scGPT, indicating architecture-agnostic biological encoding.
**Test**: Load Geneformer static gene embeddings (token embedding matrix). Map same 14 marker genes. Compute within- vs cross-type AUROC. Compare to scGPT layer-11 AUROC=0.853.
**Expected signal**: Geneformer AUROC ≥0.75 → geometry is model-agnostic; if lower, scGPT-specific.
**Null/control**: Permutation shuffle of cell-type labels.
**Value**: high | **Cost**: medium

### H-D: Layer-depth specialization curve — functional form
**Hypothesis**: The monotonic deepening of cell-type effect (L0: −0.164 → L11: −0.276) follows a sigmoid or linear ramp, not a step function, revealing continuous progressive specialization.
**Test**: Fit linear and sigmoid models to the 12-point within/cross distance-effect vs layer-index curve. Report R², test if late-layer acceleration exists (layers 8-11 vs 0-7 slope ratio). Compare to co-polarity ratio trend (no trend, rho=0.084) to contrast the two signals.
**Expected signal**: Linear or sigmoid fit explains >90% variance; if sigmoid, inflection point identifies the "specialization layer."
**Null/control**: Flat model (constant effect).
**Value**: medium | **Cost**: low

### H-E: TRRUST activation vs repression directional proximity
**Hypothesis**: Repressor TF-target pairs and activator TF-target pairs occupy geometrically distinct regions of the embedding space.
**Test**: Split 141 TRRUST-exclusive pairs (and 185 STRING-overlapping) by annotation (Activation / Repression / Unknown). Mann-Whitney test: activation pairs vs repression pairs distances. Effect vs both-direction TRRUST background.
**Expected signal**: Repressor pairs at systematically different distances than activator pairs (either closer or further); or no difference — both informative.
**Null/control**: All TRRUST pairs pooled; permutation of direction labels.
**Value**: medium | **Cost**: low

### H-F: TF hub structure — do TFs with many targets cluster together?
**Hypothesis**: TFs with ≥5 TRRUST targets occupy a geometrically distinct, clustered subspace relative to their own targets and to non-TF genes.
**Test**: From 209 named genes, identify TF-only subsets (from TRRUST edge table, source nodes). Compute TF–TF pairwise distances vs TF–target pairwise distances vs target–target pairwise distances. Kruskal-Wallis test. Null: equal distance across all three groups.
**Expected signal**: TF–TF distances < TF–target < target–target, or TF–TF < all others. Hub TFs cluster.
**Null/control**: Random subsets of 209 genes matched by size.
**Value**: medium | **Cost**: low

### H-G: STRING confidence decile nonlinearity
**Hypothesis**: The STRING confidence–distance relationship has a nonlinear threshold: high-confidence edges (>0.7) show much stronger proximity than the linear rho=−0.093 suggests.
**Test**: Bin 3092 STRING pairs into 10 deciles by confidence score. Compute mean pairwise distance per decile per layer. Fit linear vs piecewise-linear model (threshold at decile 7 vs 8 vs 9). Report threshold where distance curve bends.
**Expected signal**: Strong nonlinearity above score 0.7; bottom 50% of STRING scores contribute noise.
**Null/control**: Linear fit across all deciles.
**Value**: medium | **Cost**: low

### H-H: GO biological process term proximity
**Hypothesis**: Genes sharing Gene Ontology (GO) biological process annotations are closer in embedding space than genes without shared GO terms, independent of STRING/TRRUST overlap.
**Test**: For the 209 named genes, pull GO BP annotations (e.g., from mygene.info or a local OBO file). Build shared-GO-term pairs (threshold: ≥2 shared terms) vs non-shared pairs. Mann-Whitney distance test. Subtract STRING/TRRUST-annotated pairs (overlap correction).
**Expected signal**: GO-sharing pairs show AUROC >0.55 at most layers; if positive, this is a third independent biological anchor.
**Null/control**: Non-GO-sharing pairs; STRING/TRRUST-overlap-corrected baseline.
**Value**: high | **Cost**: medium

### H-I: Cell-type-conditional intrinsic dimensionality
**Hypothesis**: Within-cell-type marker gene subspaces have lower intrinsic dimensionality than cross-cell-type subsets, reflecting tighter functional manifold structure.
**Test**: TwoNN ID estimator applied to: (a) T-cell marker genes only (N=7), (b) B-cell markers (N=4), (c) random sets of 7/4 from remaining genes. Bootstrap 200 random-set samples for null CI. Report layer profile of conditional ID.
**Expected signal**: T-cell/B-cell marker ID significantly below random-set CI at layers 8–11.
**Null/control**: Random same-size gene samples from 209-gene pool.
**Value**: medium | **Cost**: low

### H-J: Cell-type centroid geometry — inter-cluster metric space
**Hypothesis**: The pairwise distances between cell-type cluster centroids (T/B/fibroblast) form a stable, interpretable metric space that changes predictably across transformer layers.
**Test**: Compute per-layer centroids for T, B, fibroblast marker genes. Track 3 centroid–centroid distances (T–B, T–fibro, B–fibro) across 12 layers. Test: do centroid distances increase with depth (progressive differentiation) or stabilize? Compare to cell-type AUROC deepening curve.
**Expected signal**: Centroid distances increase monotonically with layer, tracking the within/cross separation effect.
**Null/control**: Centroid distances from 3 random gene subsets matched by size.
**Value**: medium | **Cost**: low

### H-K: True persistence entropy per layer (ripser rerun)
**Hypothesis**: Shannon entropy of H1 lifetime distribution decreases with transformer depth, consistent with progressive topological simplification.
**Test**: Re-run ripser at each of 12 layers, storing full dgms[1] lifetime arrays. Compute entropy: −Σ p_i log p_i where p_i = lifetime_i / sum(lifetimes). Compare entropy vs layer Spearman. Null: feature-shuffle dgms[1].
**Expected signal**: Spearman rho < −0.6, p < 0.01 (prior H1 compaction rho=−0.916 suggests this should hold).
**Null/control**: Feature-shuffled ripser at each layer.
**Value**: medium | **Cost**: high (ripser is slow; budget 30–60 min)

### H-L: Cell-type boundary sharpness — do marker genes cluster near centroid or periphery?
**Hypothesis**: Cell-type marker genes lie near the centroid of their within-type cluster (not at the periphery), indicating tight functional identity encoding.
**Test**: For each cell type, compute distance of each marker gene to its own centroid vs distance to nearest cross-type centroid. Ratio: centroid_distance / nearest_cross_distance. Compare to random gene groups.
**Expected signal**: Marker genes have ratio <<1 (much closer to own centroid than to other types' centroids). Non-marker genes have ratio ≈1.
**Null/control**: Random same-size groups from 209-gene pool.
**Value**: medium | **Cost**: low

### H-M: Geodesic (angular) distance vs Euclidean — which better predicts biological proximity?
**Hypothesis**: On the unit hypersphere, angular (geodesic) distance is the geometrically natural metric and outperforms L2 distance for predicting STRING/TRRUST/cell-type structure.
**Test**: Re-run H01 (STRING Spearman, TRRUST AUROC) and H02 (cell-type AUROC) using 1−cosine_similarity as the distance metric instead of L2. Compare Spearman rho and AUROC under each metric.
**Expected signal**: Angular distance gives higher AUROC for cell-type and higher |rho| for STRING; or equivalent (informative either way).
**Null/control**: L2 baseline from iter_0022.
**Value**: medium | **Cost**: low

### H-N: Biological-program co-membership (STRING + TRRUST union) vs purely geometric clustering
**Hypothesis**: A combined biological-program score (union of STRING edge + TRRUST edge + shared GO terms) is a better predictor of pairwise embedding distance than any single database alone.
**Test**: Score each pair with a 3-bit binary vector (STRING∈{0,1}, TRRUST∈{0,1}, shared_GO∈{0,1}). Compute mean distance per score level (0,1,2,3). Test monotonic trend. Null: random pair-score assignment.
**Expected signal**: Distance decreases monotonically with co-membership count (score 3 << score 0). Each added annotation layer increases geometric proximity.
**Null/control**: Random permutation of biological annotations across pairs.
**Value**: high | **Cost**: medium

---

## Top 3 for Immediate Execution

### Slot 1 — High-probability discovery candidate
**H-A + H-B (Cell-type expansion + contamination control, bundled)**
Run together in one script. Direct extension of the strongest result. H-B is required to validate H-A's specificity claim. Collectively takes the headline finding from "3 cell types" to "5-6 cell types with demonstrated marker-specificity." Cost: low. Expected: AUROC remains ≥0.80 for expanded vocabulary; non-marker AUROC distribution centers at 0.50.

### Slot 2 — High-risk/high-reward candidate
**H-H (GO biological process proximity)**
If positive, establishes a third independent biological anchor beyond STRING and TRRUST. Opens GO term enrichment as a geometric interpretability tool. Requires external data (GO annotations for 209 genes). Risk: many of the 209 named genes may have sparse GO coverage. If positive at AUROC >0.55 (overlap-corrected), this is a major result.

### Slot 3 — Cheap broad-screen candidate
**H-D + H-E + H-G (layer curve fit + TRRUST direction + STRING decile, bundled)**
All three re-use existing iter_0022 result arrays (no new embedding loads needed). H-D: fit functional form to the 12-point specialization curve (30 lines of code). H-E: split TRRUST pairs by direction annotation (already in edge table). H-G: bin STRING pairs by decile. Combined runtime: <5 minutes. Maximum information per compute unit.
