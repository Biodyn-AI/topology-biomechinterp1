# Brainstormer Hypothesis Roadmap — iter_0008 → iter_0009

---

## Retire / Deprioritize

| Direction | Reason | Disposition |
|-----------|--------|-------------|
| Rewiring/permutation PH nulls | Degenerate null confirmed in multiple prior iterations | **retire_now** — fully closed |
| Drift TF/RNAPolII GO enrichment | Failed gene-label shuffle (emp_p=0.124) | **retire_now** |
| Whole-vocabulary drift analysis | 4803 gene set without label-shuffle null; computationally expensive; no confirmed signal | **deprioritize** |
| SV3 antigen presentation at N=500 | emp_p=0.088 not significant at current permutation depth | **rescue_once_with_major_change**: re-run with N=5000; if still >0.05, retire |
| Feature-column shuffle null | Documented as degenerate; replaced by gene-label shuffle | **retire_now** — already closed |

---

## New Hypothesis Portfolio

### H-A: SV2 Bottom-Pole Extracellular Vesicle Gene-Label Null
**Hypothesis:** The SV2 bottom-pole enrichment for extracellular vesicle genes (GO:0070062, OR=7.99, Fisher p=8.7e-9) is specific to actual gene-to-embedding assignment, not a consequence of annotation density.
**Test:** Gene-label shuffle (N=1000) at layer 11, SV2 bottom-25%. Empirical p for GO:0070062. Also test GO:0005856 (cytoskeleton, OR=20.35) and GO:0005739 (mitochondrion, OR=6.36) in SV2 bottom.
**Expected signal:** emp_p < 0.05 for extracellular vesicle; confirms SV2 as a bipolar immune-vesicle axis.
**Null:** emp_p ≥ 0.05; SV2 bottom is enrichment noise.
**Value:** HIGH | **Cost:** LOW (existing infrastructure, <10 min)

### H-B: Complete Compartment-Layer Map (All 12 Layers × 8 Compartments)
**Hypothesis:** Different GO:CC compartments are selectively enriched at specific scGPT layers (SV1 top-pole), forming a reproducible "layer compartment map" beyond the confirmed mito@L3 and ER@L7/8.
**Test:** For each of 12 layers × compartments {mitochondrion (GO:0005739), ER lumen (GO:0005788), Golgi apparatus (GO:0005794), lysosome (GO:0005764), nucleus (GO:0005634), plasma membrane (GO:0005886), cytosol (GO:0005829), extracellular space (GO:0005615)}: compute SV1 top-25% Fisher enrichment. Run gene-label shuffle (N=500) for any nominally significant hit (raw Fisher p < 0.01).
**Expected signal:** At least 3-4 additional layer-compartment pairs survive null; reveals full trajectory of organelle encoding across transformer layers.
**Null:** Only mito@L3 and ER@L7/8 survive; no further structure.
**Value:** HIGH | **Cost:** MEDIUM (96 Fisher tests, ~20-30 selective shuffle runs, ~30-60 min)

### H-C: SV2 Cross-Layer Stability (Is Immune Axis Layer-11-Specific or Persistent?)
**Hypothesis:** The SV2 immune signaling axis (IL-4, T cell regulation) is present across multiple layers (not just layer 11), indicating a persistent secondary organizational principle throughout scGPT's transformer stack.
**Test:** Run SV2 GO enrichment at layers 0, 3, 5, 7, 9, 11. Gene-label shuffle (N=300) at each layer for top GO term. Compute stability score: number of layers where emp_p < 0.05 for immune-related terms.
**Expected signal:** SV2 immune enrichment present at ≥6/12 layers; complements the SV1 cross-layer analysis.
**Null:** Immune signal only at layer 11; SV2 content random at earlier layers.
**Value:** HIGH | **Cost:** MEDIUM (6 layers × 2 poles × shuffle, ~30 min)

### H-D: SV3 Antigen Presentation Re-Test at N=5000
**Hypothesis:** The SV3 antigen-presentation MHC-II enrichment (emp_p=0.088 at N=500) crosses significance threshold with adequate statistical power.
**Test:** Gene-label shuffle N=5000 at layer-11 SV3 top-pole. Record empirical p for GO:0019886 and all MHC-II related terms.
**Expected signal:** emp_p < 0.05 (rescue); 500-rep noise was hiding a real ~0.03-0.04 effect.
**Null:** emp_p ≥ 0.05 → retire SV3 antigen presentation.
**Value:** MEDIUM | **Cost:** LOW (10× the current permutation count, same code, ~10-15 min)

### H-E: Ground-Truth Signal Peptide vs. SV1 (UniProt)
**Hypothesis:** Genes with experimentally confirmed signal peptides (UniProt annotation "Signal peptide") are enriched in SV1 top-quartile more strongly than the GO-proxy (H03, OR=2.88).
**Test:** Retrieve UniProt "Signal peptide" field for 209 named genes via UniProt REST API. Fisher test + Mann-Whitney on SV1 scores. Gene-label shuffle (N=500). Also compare OR vs. GO-proxy OR to see if ground-truth signal is stronger.
**Expected signal:** UniProt signal peptide OR > 2.88; confirms and sharpens the H03 proxy result.
**Null:** OR similar or lower than proxy; proxy was already capturing the effect well.
**Value:** HIGH | **Cost:** LOW (UniProt REST API for 209 genes, fast)

### H-F: Cross-Layer SV1 Gene Rank Stability
**Hypothesis:** Individual genes maintain consistent relative SV1 rank across layers (genes at SV1 top-pole at L11 were already high-ranking at L0), indicating the biological gradient is encoded from early layers.
**Test:** Extract SV1 projection for 209 named genes at all 12 layers. Compute Spearman rank correlation between L0 and L11 SV1 ranks. Compute all consecutive-layer correlations (L0-L1, L1-L2, ..., L10-L11). Plot as stability profile.
**Expected signal:** Spearman r(L0, L11) > 0.6; consecutive-layer stability r > 0.9. Rank crystallizes early.
**Null:** r(L0, L11) < 0.3; axis identity shifts substantially between early and late layers.
**Value:** HIGH | **Cost:** LOW (pure computation on existing data, <5 min)

### H-G: SV1 vs SV2 Orthogonality Test (Biological Axis Independence)
**Hypothesis:** The SV1 (extracellular secretion) and SV2 (immune signaling) axes are biologically independent — genes high on SV1 are not systematically high or low on SV2.
**Test:** Compute Spearman correlation between SV1 and SV2 projections for 209 genes at layer 11. Test whether extracellular-annotated genes are preferentially SV1-high/SV2-neutral, and immune-annotated genes are SV2-high/SV1-neutral. Visualize 2D scatterplot SV1 vs SV2 colored by compartment.
**Expected signal:** |Spearman(SV1, SV2)| < 0.15 globally; extracellular genes cluster near SV2=0; immune genes cluster near SV1=0.
**Null:** SV1 and SV2 are correlated (|r| > 0.3); axes are not truly orthogonal biologically.
**Value:** HIGH | **Cost:** LOW (<5 min pure computation)

### H-H: SV2 Immune Axis Gene Identity and Mechanistic Interpretation
**Hypothesis:** The 8 genes in SV2 top-quartile annotated to GO:0032753 (IL-4 positive regulation) are known immune cell-type marker or cytokine genes, not random inclusions, confirming the mechanistic interpretation.
**Test:** Extract the 8 specific gene IDs driving GO:0032753 enrichment in SV2 top-pole. Cross-reference with: (1) IL-4 signaling pathway databases (KEGG, Reactome), (2) known cell-type marker gene lists (T helper cells, mast cells, eosinophils), (3) TRRUST transcription factor-target database for IL-4 responsive genes. Document biological coherence qualitatively.
**Expected signal:** All 8 genes are canonical IL-4 pathway members or known IL-4-responsive genes; mechanistic interpretation fully supported.
**Null:** Genes are incidental inclusions with indirect GO annotation inheritance.
**Value:** HIGH | **Cost:** LOW (manual/database lookup, no permutation needed)

### H-I: SV2 Top vs. Bottom Pole Contrast (Immune vs. Vesicle Axis)
**Hypothesis:** SV2 encodes a single bipolar axis where immune signaling genes are at the top and extracellular vesicle / mitochondrial genes are at the bottom, forming a "intracellular organelle vs. immune activation" contrast.
**Test:** After running H-A null for SV2 bottom, create a joint gene list: SV2 top-quartile genes vs. SV2 bottom-quartile. Compare GO:CC annotations between poles. Test whether this axis has a biological interpretation as "intracellular energy/vesicle" vs. "extracellular immune signaling."
**Expected signal:** Clear biological contrast — one pole = immune cell surface signaling, other pole = intracellular organelles/vesicles; axis describes a cell state or cell-type contrast.
**Null:** No coherent contrast; poles encode unrelated biology.
**Value:** HIGH | **Cost:** LOW (integration/annotation of existing data)

### H-J: CKA Cross-Model Alignment vs. SV1/SV2 Ratio Correlation
**Hypothesis:** The per-layer CKA alignment score (from iter_0004) correlates with the SV1/SV2 singular value ratio, suggesting biological axis emergence and model-convergence are coupled phenomena.
**Test:** Load iter_0004 CKA alignment scores per layer. Load SV singular value ratio (SV1/SV2) per layer from existing data. Compute Pearson/Spearman r across 12 layers. Test significance via permuting layer order (N=1000).
**Expected signal:** r > 0.6; both curves peak at layers 9-11.
**Null:** r < 0.3; CKA alignment and biological axis strength are unrelated.
**Value:** HIGH | **Cost:** LOW (<5 min, data already exists)

### H-K: Annotation Density Confounder Control
**Hypothesis:** The GO enrichment results (especially SV1 and SV2) are not confounded by genes with more GO annotations being over-represented in top quartiles.
**Test:** For 209 named genes, count total GO:CC annotations per gene. Compute Spearman correlation between annotation count and SV1 projection. If confounded (r > 0.3), re-run enrichment with annotation-count-stratified random controls.
**Expected signal (clean):** |r| < 0.2 between annotation count and SV1/SV2 scores; enrichment results stand without correction.
**Expected signal (confound present):** r > 0.3; correction needed.
**Value:** MEDIUM | **Cost:** LOW (<5 min)

### H-L: Persistent Homology on 209 Named Genes at Layer 11
**Hypothesis:** The 209 named gene embeddings at layer 11 form a topologically non-trivial point cloud (H1 loops) that exceeds what is expected from gene-label shuffles.
**Test:** Compute Vietoris-Rips persistent homology (H0, H1) on 209-gene layer-11 submatrix using Euclidean distances. Compare observed total H1 persistence mass to N=200 gene-label shuffles. Focus on small gene set (209 genes) which is computationally tractable.
**Expected signal:** H1 persistence mass exceeds 95th percentile of shuffles; indicates non-trivial loop structure in biological gene space.
**Null:** H1 features not distinguishable from shuffled assignments.
**Value:** MEDIUM | **Cost:** MEDIUM (PH on 209×512 matrix, feasible; ~20 min including setup)

### H-M: STRING PPI Network Distance vs. SV1 Axis Proximity
**Hypothesis:** Pairs of genes close on the SV1 axis have higher STRING protein-protein interaction scores than distant pairs, linking geometric proximity to functional co-regulation.
**Test:** Compute all pairwise SV1 axis distances for 209 genes (|SV1_i - SV1_j|, 21736 pairs). Retrieve STRING scores for these pairs. Spearman correlation. Gene-label shuffle null (N=500).
**Expected signal:** Spearman r < -0.1 (closer = higher PPI score); survives null.
**Null:** No significant correlation.
**Value:** HIGH | **Cost:** MEDIUM (STRING data retrieval for gene pairs)

### H-N: SV4/SV5 Axes Scan at Layer 11
**Hypothesis:** SVD dimensions 4 and 5 of the layer-11 embedding encode further distinct biological axes (e.g., metabolic, proliferative, developmental) that are not co-linear with SV1-SV3.
**Test:** Project 209 genes onto SV4 and SV5 at layer 11. GO enrichment (top/bottom 25%) with gene-label shuffle null (N=500). Check for non-redundant GO terms vs. SV1-SV3.
**Expected signal:** At least one of SV4/SV5 encodes a significant axis (emp_p < 0.05); the gene embedding space has at least 4 interpretable dimensions.
**Null:** No significant enrichment; biological information is concentrated in SV1-SV3.
**Value:** MEDIUM | **Cost:** LOW (trivial extension of existing SVD infrastructure)

### H-O: Per-Gene SV1 Trajectory Cluster Analysis (Cross-Layer)
**Hypothesis:** The 209 named genes cluster into distinct groups based on their SV1 score trajectory across 12 layers (how their relative position evolves), revealing co-regulated modules with coherent layer-dependent encoding.
**Test:** Build 209×12 matrix of SV1 projections across layers. Apply hierarchical clustering (Ward linkage, cosine distance on z-scored trajectories). Identify 3-6 clusters. Run GO enrichment on each cluster. Validate cluster coherence with silhouette score.
**Expected signal:** 2-3 clusters with distinct trajectory shapes; each cluster enriched for coherent GO terms (e.g., secreted proteins crystallize early, immune genes crystallize late).
**Null:** No cluster structure; genes have idiosyncratic trajectories.
**Value:** HIGH | **Cost:** MEDIUM (requires per-layer SV1 extraction for all 12 layers, ~15-20 min)

---

## Top 3 for Immediate Execution

### Execution Candidate 1: HIGH-PROBABILITY DISCOVERY
**H-A (SV2 bottom-pole null) + H-I (SV2 bipolar axis interpretation)**

Rationale: The SV2 bottom-pole extracellular vesicle enrichment (OR=7.99, Fisher p=8.7e-9) is the most significant untested result in the project. It likely survives the gene-label null (OR is very large). If confirmed, SV2 becomes a fully validated bipolar axis: IL-4 immune signaling (top) vs. extracellular vesicle/organelle (bottom). This transforms the SV2 finding from "one pole confirmed" to a complete second biological axis of the scGPT embedding — a major interpretability result.

Test design:
1. Gene-label shuffle N=1000 at layer-11 SV2 bottom-25% for GO:0070062, GO:0005856, GO:0005739
2. Compute empirical p for each
3. If confirmed: document the full SV2 axis as immune-activation vs. vesicle/organelle biology
4. Generate SV1 vs. SV2 2D plot colored by compartment (H-G) — visualization of the two-axis structure

Expected result: emp_p < 0.05 for extracellular vesicle; SV2 axis confirmed as bipolar. The 2D SV1×SV2 plot would clearly show the two independent biological gradients.

---

### Execution Candidate 2: HIGH-RISK / HIGH-REWARD
**H-B: Complete Compartment-Layer Map**

Rationale: If mito@L3, ER@L7/L8, and extracellular@L11 are nodes in a reproducible trajectory, there may be 3-4 more compartments enriched at specific layers. This would be a "biological computation map" of the transformer — showing which organelle/compartment information is being processed at each layer. This is the kind of result that would constitute the central figure of a paper. Risk: the remaining compartments may not show clean transients.

Test design:
1. For each of 12 layers, test SV1 top-25% enrichment for 8 GO:CC compartments (list above)
2. Only run gene-label shuffle (N=500) for entries where raw Fisher p < 0.01
3. Build a 12×8 heatmap of -log10(p_value) where available
4. Record which compartments show confirmed transients (emp_p < 0.05)

Expected result: 4-6 compartment-layer pairs confirmed; reveals a structured "processing sequence" hypothesis across scGPT layers.

---

### Execution Candidate 3: CHEAP BROAD-SCREEN
**H-F (cross-layer SV1 rank stability) + H-J (CKA correlation) + H-D (SV3 at N=5000) + H-K (annotation density confounder)**

Rationale: All four use only existing data and code with minimal new computation. H-F tests whether the extracellular axis is crystallizing gradually or abruptly. H-J tests whether CKA and biological axis emergence are linked (mechanistic claim). H-D either confirms or retires SV3. H-K validates that no GO enrichment confound exists. Together these four close the most important open questions about data quality and interpretation, at essentially zero additional cost.

Test design:
1. H-F: Extract SV1 projections per layer (existing infrastructure), compute Spearman(L0, L11) and consecutive pairs. Runtime: <5 min
2. H-J: Load iter_0004 CKA data, correlate with SV1/SV2 ratio per layer. Runtime: <5 min
3. H-D: Re-run SV3 shuffle at N=5000. Runtime: ~10-15 min
4. H-K: Count GO annotations per gene, correlate with SV1/SV2 scores. Runtime: <5 min
Total: ~30 min, high information density per runtime cost.
