# Brainstormer Hypothesis Roadmap — iter_0042 → iter_0043

---

## Retire / Deprioritize

| ID | Direction | Reason |
|----|-----------|--------|
| R01 | Cross-model scGPT/Geneformer alignment | Blocked by vocab mismatch; no new data path |
| R02 | GC-plasma subspace angles | Inconclusive in iter_0040; trajectory approach (H03) supersedes it |
| R03 | Chromosomal proximity | Negative result; no rescue path |
| R04 | PC2/PC3 axes | Retired; B-cell PC1 result is the canonical finding |
| R05 | GO BP enrichment in SV poles | Retired in iter_0011; 0/591 terms |

---

## New Hypothesis Portfolio

### H-A: Full GC-TF Pairwise Distance Matrix (PAX5 / BATF / BACH2 / BCL6 / PRDM1 / IRF4)
**Hypothesis**: The 6-gene GC-TF regulatory network (activation: PAX5→BATF/BACH2; repression: BCL6→PRDM1/IRF4) is geometrically encoded in scGPT embeddings as convergence (activators) and divergence (repressed targets), with the repression axis visible as anti-convergence.
**Test**: Compute full pairwise Euclidean distance matrix among {PAX5, BATF, BACH2, BCL6, PRDM1, IRF4} at each layer L0–L11. Track PRDM1/IRF4 rank near B-cell centroid across layers (same method as H03_iter0042).
**Expected signal**: PAX5–BATF, PAX5–BACH2, BATF–BACH2 distances decrease (confirmed for 3/3); BCL6–PAX5, BCL6–BATF, BCL6–BACH2 remain large or decrease slowly; PRDM1 and IRF4 diverge from B-cell centroid (rank increases at same layers BATF/BACH2 converge — anti-convergence).
**Null/control**: 10 random non-B-cell genes as in H03_iter0042; shuffle null for rank trajectories.
**Value**: high | **Cost**: low (reuse iter_0042 infrastructure)

### H-B: BCL6 Metabolic Cluster Layer-Resolved Composition
**Hypothesis**: BCL6's metabolic neighborhood (9/20 at L0, 4/20 at L11) dissolves progressively — specific genes are lost at specific layers, revealing when and how the metabolic identity is reorganized.
**Test**: At each layer L0–L11, recompute BCL6's k=20 NN. Track which genes are in the neighborhood at each layer. Identify genes that join/leave at L3 (the B-cell ID inflection layer).
**Expected signal**: Some metabolic genes (NAMPT, PFKFB3) persist throughout; others leave after L3 — replaced by genes closer to B-cell or immune-function modules.
**Null**: No genes should change neighborhood membership more than expected by distance diffusion (compare to random gene neighborhood stability).
**Value**: high | **Cost**: low

### H-C: GO Annotation of BCL6 Metabolic Cluster
**Hypothesis**: The 10-gene BCL6 metabolic cluster (NAMPT, GLUL, PFKFB3, ACSL1, NIBAN1, FNDC3B, VMP1, STAT3, CEBPD, TRIB1) is enriched for aerobic glycolysis / Warburg effect gene ontology terms, matching BCL6's known role in metabolic reprogramming of GC B-cells.
**Test**: Fetch GO annotations for all 10 genes via MyGene.info API. Test Fisher enrichment for metabolic GO terms (glycolysis, fatty acid synthesis, oxidative phosphorylation, metabolic stress). Compare to 100 random gene sets of size 10.
**Expected signal**: Glycolysis (PFKFB3, NAMPT), fatty acid (ACSL1), stress response (CEBPD, TRIB1) enriched. OR > 3 for relevant metabolic terms.
**Null**: 100 random gene sets of size 10 from vocab.
**Value**: medium | **Cost**: low

### H-D: T-cell Manifold Expansion — Post-L7 Driver Analysis
**Hypothesis**: T-cell ID increases after L7 (slope change +4.965 in iter_0042 H01), reflecting T-cell activation or differentiation geometry. Specific T-cell marker subsets (cytotoxic vs regulatory vs exhaustion) may drive this expansion differently.
**Test**: Stratify T-cell markers by subtype: effector (GZMB, PRF1), regulatory (FOXP3, IL2RA), exhaustion (TIGIT, PDCD1), helper (CD4, IFNG). Compute TwoNN ID per subtype across layers. Identify which subtype has the strongest L7 breakpoint. Also compare T-cell ID breakpoint at L7 vs random gene subsets of size 9.
**Expected signal**: Effector or exhaustion markers show steepest L7 breakpoint (expansion); regulatory markers are flat. Null p < 0.05 for at least one subtype.
**Null**: 500-permutation null as in H01_iter0042; compare to Myeloid flat profile as reference.
**Value**: high | **Cost**: medium

### H-E: SPIB Trajectory — Pre-wired or Recruited?
**Hypothesis**: SPIB (ETS-family TF, B-cell identity, found in B-cell neighborhood at L2) is pre-wired at L0 like PAX5 — not recruited like BATF/BACH2. If true, this divides GC TFs into two geometric classes: anchor (PAX5, SPIB) and recruited (BATF, BACH2).
**Test**: Apply same convergence trajectory analysis as H03_iter0042 to SPIB: rank near B-cell centroid at L0–L11, Spearman rho, distance to PAX5 across layers. Compare trajectory shape to PAX5 (near-flat, pre-wired) vs BATF (steeply converging).
**Expected signal**: SPIB rank near B-cell centroid starts low (< 200) at L0 and stays low (flat profile), unlike BATF (1510→189). Distance to PAX5 starts small at L0 for SPIB vs large for BATF.
**Null**: 10 random non-B-cell genes.
**Value**: medium | **Cost**: low

### H-F: BCL6 vs B-cell Centroid — Does BCL6 Diverge?
**Hypothesis**: BCL6 starts close to the B-cell centroid at L0 (it is a B-cell gene) but diverges progressively as its metabolic cluster identity strengthens — moving away from the B-cell cluster in opposite direction from BATF/BACH2.
**Test**: Track BCL6 rank near B-cell centroid across L0–L11 (same method as H03_iter0042). Compute Spearman rho for BCL6 vs layer. Compare sign of rho to BATF/BACH2 (converging) vs BCL6 (potentially diverging or flat).
**Expected signal**: BCL6 rank near B-cell centroid stays high (i.e., stays far) or increases (diverges further). Rho near 0 or positive (opposite sign from BATF rho=-0.97).
**Null**: Same 10 random gene null as H03.
**Value**: high | **Cost**: low (trivial addition to H03 infrastructure)

### H-G: Myeloid Flat ID — Housekeeping Gene Baseline
**Hypothesis**: Myeloid ID being flat (23-26 across all layers) reflects broad metabolic/housekeeping gene expression, not lineage-specific geometric identity. Test: housekeeping gene sets (ribosomal, mitochondrial) should show similarly flat ID profiles.
**Test**: Select ribosomal protein genes (RPS, RPL series) and mitochondrial complex genes (NDUFA, ATP5 series) if in vocab. Compute TwoNN ID across layers. Compare flatness (variance of ID values across layers) to B-cell (high variance, inflection) vs Myeloid (low variance, flat).
**Expected signal**: Housekeeping genes show flat profiles similar to myeloid, confirming myeloid ID flatness reflects generic gene behavior, not myeloid-specific geometry.
**Null**: B-cell ID profile as the high-variance reference.
**Value**: medium | **Cost**: low

### H-H: Persistent Homology on Lineage-Specific Gene Subsets
**Hypothesis**: PH (H1 loops) computed on B-cell-only gene subsets shows a topological transition at L3, coinciding with the TwoNN changepoint. PH may detect higher-order structure (cycles, voids) that TwoNN misses.
**Test**: Run Ripser H0+H1 on B-cell marker genes (n=7) + their k=20 NN (n~140 unique) at each layer L0–L11. Track H1 sum persistence and number of loops across layers. Compare to T-cell and Myeloid equivalent subsets.
**Expected signal**: B-cell H1 persistence peaks or changes character at L3; T-cell changes at L7; Myeloid is flat.
**Null**: Feature-shuffle null per layer (as in iter_0003).
**Value**: high | **Cost**: medium

### H-I: B-cell Community Phase Transition — First Layer Detection
**Hypothesis**: The B-cell kNN community emerges not gradually but as a phase transition — precision@10 is near-null before a specific layer and jumps sharply. The L3 ID inflection may coincide with community emergence.
**Test**: Compute B-cell precision@10 (same as iter_0035) at all 12 layers L0–L11 with 500-draw bootstrap null. Fit a step function to the precision@10 vs layer trajectory. Identify the first layer where precision@10 exceeds null p < 0.01. Compare to L3 (ID inflection) and earlier work.
**Expected signal**: Sharp transition at L2–L3 with precision@10 jumping from ~null to ~7x null. Continuous increase is less interesting than a step.
**Null**: Bootstrap null per layer (existing infrastructure).
**Value**: medium | **Cost**: low

### H-J: Geodesic vs Euclidean Distance Comparison for GC-TF Trajectories
**Hypothesis**: GC-TF (BATF, BACH2) convergence trajectories measured by geodesic graph distance differ from Euclidean — manifold structure may mean Euclidean underestimates actual convergence on the embedding manifold.
**Test**: Build kNN graph (k=10) on full gene embedding per layer. Compute shortest-path graph distance between {PAX5, BATF, BACH2, BCL6} at each layer. Compare trajectory shape and rate of convergence vs Euclidean from H03_iter0042.
**Expected signal**: Geodesic distances decrease faster than Euclidean, suggesting stronger actual convergence on the manifold. Alternatively, if geodesic > Euclidean, suggests the path is topologically constrained.
**Null**: Expected path length for random gene pairs at same embedding distances.
**Value**: medium | **Cost**: medium

### H-K: BCL6 Metabolic Module — Network Centrality Test
**Hypothesis**: BCL6 is the topological hub (highest degree centrality) of the metabolic cluster kNN graph, not just a member. If BCL6 is the hub, it is the core organizer of this functional module in the scGPT geometry.
**Test**: Build k=5 NN graph among the BCL6 metabolic cluster members (BCL6 + 10 neighbors). Compute degree centrality and betweenness centrality for each node at each layer. Test: is BCL6 consistently the highest-centrality node?
**Expected signal**: BCL6 degree centrality > all 10 neighbors at L0–L5. Centrality decreases at late layers (consistent with neighborhood dissolution from 9→4/20).
**Null**: Random re-labeling of node identities within the cluster.
**Value**: medium | **Cost**: low

### H-L: Cross-Dataset Validation — Lung scGPT Embeddings for B-cell Markers
**Hypothesis**: The B-cell geometry (PAX5 pre-wired, BATF/BACH2 convergence) is an immune-context-specific finding, not a general scGPT property. Test: do the same trajectories appear in the scGPT lung embeddings (different tissue context)?
**Test**: Run BATF, BACH2, PAX5, BCL6 convergence trajectory analysis (same method as H03_iter0042) on cycle3_lung or scGPT lung embeddings. Compute Spearman rho for rank near B-cell centroid vs layer.
**Expected signal**: Lung embeddings should NOT show the same convergence (no GC biology expected in lung dataset). Weak or absent convergence in lung would confirm the finding is immune-context-specific.
**Null**: If rho is equally strong in lung, the finding may be a general scGPT property rather than biology-specific.
**Value**: high | **Cost**: medium

### H-M: Layer-Stratified Intrinsic Dimension Ratio: B-cell vs Non-B-cell Genes
**Hypothesis**: The B-cell ID compression (after L3) is specific to B-cell-related genes — non-B-cell genes should not show the same compression. Test ID separately for B-cell gene subsets vs randomly sampled non-B-cell genes across layers.
**Test**: Split cycle4_immune vocab into B-cell markers (n~20) and non-B-cell genes (n~200 random). Compute TwoNN ID for each group across layers. Breakpoint analysis for each group.
**Expected signal**: B-cell genes show L3 breakpoint; non-B-cell genes show no clear breakpoint or show a different layer. Cross-group comparison tests specificity of L3 inflection.
**Null**: Permutation of gene-to-group assignment (500 replicates).
**Value**: high | **Cost**: low

---

## Top 3 for Immediate Execution

### Candidate 1 — HIGH PROBABILITY DISCOVERY
**H-A: Full GC-TF pairwise distance matrix + PRDM1/IRF4 anti-convergence**
- Directly extends the strongest result of iter_0042 (H03, rho=-0.97).
- Adds repression geometry: PRDM1 and IRF4 are the BCL6-repressed plasma cell TFs. If they anti-converge while BATF/BACH2 converge, the model is encoding the GC reaction regulatory logic geometrically.
- Cost: low — reuse iter_0042 infrastructure exactly.
- Output: a "GC-TF pairwise distance triangle + anti-convergence" figure that anchors the germinal center biology chapter of the paper.

### Candidate 2 — HIGH RISK / HIGH REWARD
**H-D: T-cell manifold expansion driver analysis (post-L7 subtype stratification)**
- The T-cell ID expansion is the counter-narrative to B-cell compression: two lineages show geometrically opposite behavior. This is an unexpected and publishable finding.
- If T-cell effector/exhaustion markers drive L7 expansion, this reveals that scGPT encodes different functional programs at different transformer depths: B-cell identity consolidates at L3 while T-cell effector function consolidates at L7.
- Risk: T-cell subtype markers may be OOV. Requires vocab check before full run.
- Reward: positions the paper as discovering lineage-stratified geometric phase transitions, not just B-cell biology.

### Candidate 3 — CHEAP BROAD SCREEN
**H-F + H-B combo: BCL6 divergence from B-cell centroid + layer-resolved cluster composition**
- H-F is 10 lines of code: track BCL6 rank near B-cell centroid across layers using iter_0042 infrastructure. Expected sign is opposite to BATF/BACH2 — if confirmed, BCL6 is geometrically anti-correlated with B-cell identity while BATF/BACH2 are correlated. That's the full regulatory story.
- H-B adds which genes leave BCL6's neighborhood at L3+ — biological specificity of the metabolic module dissolution.
- Combined cost: low. Combined payoff: completes the BCL6 narrative (divergence + cluster dissolution).

---

## Summary Table

| ID | Hypothesis | Value | Cost | Priority |
|----|-----------|-------|------|----------|
| H-A | GC-TF full pairwise + PRDM1/IRF4 anti-convergence | high | low | **Top-1** |
| H-D | T-cell expansion driver (subtype stratification) | high | medium | **Top-2** |
| H-F+H-B | BCL6 divergence + cluster dissolution | high | low | **Top-3** |
| H-C | BCL6 GO annotation | medium | low | Next |
| H-E | SPIB pre-wired vs recruited | medium | low | Next |
| H-M | B-cell vs non-B-cell ID specificity | high | low | Next |
| H-I | B-cell community phase transition | medium | low | Next |
| H-H | PH on lineage-specific subsets | high | medium | Queued |
| H-K | BCL6 hub centrality | medium | low | Queued |
| H-G | Myeloid flat ID housekeeping baseline | medium | low | Queued |
| H-J | Geodesic vs Euclidean convergence | medium | medium | Later |
| H-L | Cross-dataset lung validation | high | medium | Later |
