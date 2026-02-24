# Brainstormer Hypothesis Roadmap — iter_0021

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| STRING quintile Spearman (5 quintile means) | Fundamentally underpowered (n=5 test points). Design flaw, not a negative. Replace with continuous. | `retire_now` |
| Permutation null sensitivity tests | Established: z=4.88 in iter_0020. Adding more permutation tests yields no new information. | `retire_now` |
| SVD axis counting / co-polarity direction detection | Well-characterized across 6+ iterations. Refinements yield diminishing returns unless tied to a specific biological question. | `deprioritize` |
| Persistence diagram vectorization (PD images) without biological anchoring | Good idea in principle but prior unanchored topology tests were inconclusive. Only pursue if paired with biological group labels. | `rescue_once_with_major_change` |

---

## New Hypothesis Portfolio

### GROUP A — BIOLOGICAL INTERACTION BREADTH (extend the STRING/TRRUST finding)

**A1. TRRUST overlap-corrected proximity**
- *Hypothesis*: TRRUST TF-target pairs show significant geometric proximity in scGPT embeddings even after removing pairs that are also STRING edges, indicating regulatory structure contributes independently of PPI co-membership.
- *Test*: Identify TRRUST pairs not in STRING (any score threshold). Mann-Whitney(TRRUST_only vs non-TRRUST_non-STRING) at all 12 layers.
- *Expected signal*: Effect remains negative and p<0.05 in ≥8/12 layers.
- *Null/control*: Randomly matched gene pairs not in TRRUST or STRING.
- *Value*: **high** | *Cost*: **low** (data already loaded, 5-line filter)

**A2. Reactome pathway co-membership proximity**
- *Hypothesis*: Gene pairs sharing a Reactome pathway are geometrically closer in scGPT embedding space than random pairs, independent of STRING score.
- *Test*: Download Reactome gene-sets (GMT), compute intra-pathway pairs, Mann-Whitney vs non-pathway pairs at each of 12 layers.
- *Expected signal*: Effect ~−0.03 to −0.06 comparable to STRING; possibly stronger for tight complexes.
- *Null/control*: Pairs in no shared pathway.
- *Value*: **high** | *Cost*: **medium** (need GMT download and parsing)

**A3. Protein complex membership (CORUM)**
- *Hypothesis*: Subunits of the same protein complex (CORUM database) are geometrically closer than STRING neighbors, reflecting stoichiometric co-expression requirements.
- *Test*: CORUM complex gene lists → intra-complex pairs. Mann-Whitney vs STRING-matched pairs controlling for complex size.
- *Expected signal*: Effect larger than STRING (~−0.06 to −0.08), since complex members must truly co-express.
- *Null/control*: STRING pairs with score 0.4–0.6.
- *Value*: **high** | *Cost*: **medium**

**A4. Activating vs repressive TF-target directionality in TRRUST**
- *Hypothesis*: TRRUST contains directionality annotations (activation/repression). Activating pairs should be geometrically closer (co-regulated) while repressive pairs should not show the same proximity, or may even show greater distance.
- *Test*: Split 288 TRRUST pairs by interaction type (activation, repression, unknown). Mann-Whitney per type vs non-TRRUST pairs at layer 8.
- *Expected signal*: Activating pairs effect < repressive pairs effect (more negative = closer for activating).
- *Null/control*: Repressive pairs vs random.
- *Value*: **high** | *Cost*: **low** (TRRUST already loaded, just need to parse the type column)

**A5. GO Biological Process co-annotation proximity**
- *Hypothesis*: Gene pairs sharing at least one GO Biological Process term are geometrically closer than random pairs, and the effect is proportional to the number of shared GO terms.
- *Test*: From GO annotation file, compute shared-term count for all gene pairs in the named-gene set. Spearman(shared_GO_count, pairwise_distance) at each layer.
- *Expected signal*: Negative Spearman rho, stronger than STRING quintile gradient due to continuous predictor.
- *Null/control*: Pairs with zero shared GO terms.
- *Value*: **medium** | *Cost*: **medium**

---

### GROUP B — LAYER STRUCTURE & REPRESENTATION GEOMETRY

**B1. Per-layer bootstrap CIs for co-polarity enrichment (complete the H03 result)**
- *Hypothesis*: Co-polarity enrichment is significantly higher in early layers (1–4) than late layers (9–11), confirming a layer-dependent transition from broad co-expression to specialized encoding.
- *Test*: Bootstrap N=500 CIs for enrichment at each of 12 layers. Test: do CI[layer 1–4] and CI[layer 9–11] overlap?
- *Expected signal*: CI for early layers ≥1.45x, CI for late layers ≤1.35x, non-overlapping.
- *Null/control*: If CIs overlap, the gradient claim is weakened.
- *Value*: **high** | *Cost*: **low** (uses pre-computed data, just run bootstrap 12x)

**B2. Intrinsic dimensionality trajectory by layer**
- *Hypothesis*: The intrinsic dimension of the gene embedding manifold decreases monotonically with layer depth, paralleling the H1 topological compaction.
- *Test*: TwoNN estimator (or PCA 95% variance cutoff) on named-gene embeddings at all 12 layers. Spearman(layer, intrinsic_dim).
- *Expected signal*: Negative correlation matching H1 trajectory (rho ~ −0.9).
- *Null/control*: Intrinsic dimension of shuffled embeddings (should be constant or show no trend).
- *Value*: **medium** | *Cost*: **low**

**B3. STRING effect magnitude trajectory by layer (directional depth analysis)**
- *Hypothesis*: The TRRUST/STRING effect ratio changes with layer: deeper layers encode regulatory (TF-target) relationships preferentially over PPI, reflecting scGPT's training signal shifting from co-expression to regulatory context.
- *Test*: Compute TRRUST_effect / STRING_effect ratio at each layer. Spearman(layer, ratio). Compare to constant-ratio null.
- *Expected signal*: Ratio trends toward 1.0 (convergence) or increases in deep layers if regulation is emphasized.
- *Null/control*: Permuted gene-pair labels; flat ratio = no layer-dependent specialization.
- *Value*: **medium** | *Cost*: **low** (data already available in h02 JSON)

**B4. Geodesic distance vs Euclidean distance for biological proximity**
- *Hypothesis*: Geodesic distance along the manifold is a stronger predictor of biological interaction than Euclidean distance, because the manifold geometry captures non-linear structure.
- *Test*: Compute k-NN graph (k=15) and approximate geodesic distances via Dijkstra. Mann-Whitney(STRING_pairs, non-STRING_pairs) using geodesic vs Euclidean. Compare AUROCs.
- *Expected signal*: Geodesic AUROC > Euclidean AUROC by ≥0.02 units.
- *Null/control*: Euclidean AUROC as baseline.
- *Value*: **medium** | *Cost*: **medium**

---

### GROUP C — CROSS-MODEL & TRANSFERABILITY

**C1. Geneformer proximity replication**
- *Hypothesis*: STRING pairs show the same geometric proximity effect in Geneformer word embeddings as in scGPT, indicating this property is model-agnostic and reflects the underlying gene co-occurrence statistics in single-cell data.
- *Test*: Load Geneformer gene embeddings (already present from iter_0019). Compute pairwise L2-normalized distances. Mann-Whitney(STRING_pairs, non-STRING_pairs). Compare effect size to scGPT layer 0.
- *Expected signal*: Geneformer effect ~−0.03 to −0.05 (similar magnitude to scGPT).
- *Null/control*: Shuffled gene-label embeddings.
- *Value*: **high** | *Cost*: **low** (embeddings available)

**C2. Cross-model aligned distance correlation**
- *Hypothesis*: After Procrustes alignment, gene pairwise distances in scGPT and Geneformer are correlated, and the correlation is higher for STRING pairs than random pairs.
- *Test*: Align Geneformer to scGPT layer 8 via Procrustes. Compute Pearson(scGPT_distance_vector, Geneformer_distance_vector) for (a) all pairs, (b) STRING pairs, (c) random pairs. Fisher z-test for difference.
- *Expected signal*: r_STRING > r_random.
- *Null/control*: Random-pair correlation.
- *Value*: **medium** | *Cost*: **medium**

---

### GROUP D — ALGORITHMIC / MECHANISTIC

**D1. STRING continuous AUROC (upgrade of failed H01)**
- *Hypothesis*: STRING interaction score is a continuous predictor of geometric proximity (Spearman across 3092 pairs), and AUROC for predicting proximity from score ≥ 0.5/0.6/0.7/0.8 vs 0.4 baseline increases monotonically with threshold.
- *Test*: Spearman(STRING_score, pairwise_distance) for all pairs at each layer. AUROC at 0.5, 0.6, 0.7, 0.8 vs 0.4 baseline using distance percentile rank as predictor.
- *Expected signal*: Negative Spearman rho across all layers; AUROC increases with threshold.
- *Null/control*: Shuffled STRING scores.
- *Value*: **medium** | *Cost*: **low** (direct upgrade of H01)

**D2. Cell-type marker gene cluster separation**
- *Hypothesis*: Known cell-type marker genes (e.g., T-cell markers CD3D/CD3E/CD3G, hepatocyte markers ALB/APOE/FGB) cluster by cell-type identity in scGPT embedding space more than expected by chance.
- *Test*: Select 5–8 cell types with known markers (≥5 markers each). Compute within-cell-type vs cross-cell-type pairwise distances at layer 11. Mann-Whitney + effect size. Silhouette score.
- *Expected signal*: Within-cell-type distance significantly lower; silhouette >0.15.
- *Null/control*: Random gene group assignment of same sizes.
- *Value*: **high** | *Cost*: **medium** (need marker gene list curation)

**D3. Persistence entropy trajectory**
- *Hypothesis*: Persistence entropy (Shannon entropy of normalized lifetime distribution) of H0 diagrams decreases monotonically with layer, measuring the collapse of topological complexity independently of mean lifetime.
- *Test*: Compute H0 persistence entropy at each layer. Spearman(layer, entropy). Compare to H1 entropy trend.
- *Expected signal*: Monotonic decline similar to mean lifetime trajectory.
- *Null/control*: Entropy of random distance matrix.
- *Value*: **medium** | *Cost*: **low** (entropy is 3 lines from existing Ripser output)

---

## Top 3 for Immediate Execution

### #1 — High-probability discovery: A1 + B1 combined (TRRUST overlap-corrected + layer bootstrap CIs)
- **Rationale**: Both are near-zero-cost (filter existing data), build directly on two confirmed positives, and close the two most critical gaps in the current paper: (a) independence of regulatory vs PPI signal, and (b) layer-stratified enrichment CIs. Combine into one script.
- **Expected outcome**: Either confirms regulatory independence (strong new claim) or shows it's mediated by STRING membership (informative negative). Either way, the paper is strengthened.
- **Key metric**: TRRUST-only effect across layers (should be negative and significant in ≥8/12); per-layer bootstrap CIs for co-polarity.

### #2 — High-risk/high-reward: D2 (Cell-type marker gene cluster separation)
- **Rationale**: If cell-type markers cluster by cell-type identity in embedding space, this directly connects geometric structure to biological function in a cell-biology-interpretable way. This is the most impactful downstream interpretability claim: "scGPT encodes cell-type identity geometrically." It requires curating marker lists but the computational test is standard (Mann-Whitney + silhouette).
- **Expected outcome**: Strong positive would transform the paper narrative from "interaction proximity" to "cell-type geometric encoding." Negative is also informative.
- **Key metric**: Silhouette >0.15 and within-type < cross-type distance, Mann-Whitney p<0.001 at layer 11.

### #3 — Cheap broad-screen: D1 + D3 (STRING continuous AUROC + persistence entropy)
- **Rationale**: D1 fixes the failed H01 with a 5-line change (use all 3092 pairs instead of 5 quintile means). D3 adds persistence entropy from existing Ripser output. Together they close two open quantitative questions in ≤30 min of runtime.
- **Expected outcome**: D1 will be positive (direction confirmed in H01); D3 will likely confirm monotonic decline. Both are paper-table ready.
- **Key metric**: Spearman rho negative across 12 layers; entropy declines monotonically.
