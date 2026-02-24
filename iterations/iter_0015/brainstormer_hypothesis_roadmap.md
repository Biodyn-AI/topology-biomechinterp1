# Brainstormer Hypothesis Roadmap — iter_0015

---

## Retire / Deprioritize

The following are confirmed dead or saturated — do not revisit without a structural change to method or data:

| Direction | Reason | Status |
|---|---|---|
| SV1 PPI co-pole | z≈0 consistently | retire_now |
| Repression anti-pole in SV2 | Negative, iter_0013 H02 | retire_now |
| GO BP enrichment in SV2 poles | Negative, iter_0011 H03 | retire_now |
| TRRUST bootstrap CI on act-rep differential | Underpowered, not rescuable without new data | retire_now |
| Layer-depth Spearman monotonicity | Neutral, iter_0014 H02 | retire_now |
| SV5 PPI co-pole | z=1.65, non-monotonic layers, low consistency | deprioritize (low rescue potential without new geometry angle) |

**Near-saturation (diminishing returns without new angle):**
- SV2 co-pole statistics: core claim is proven; only run if testing new stratification logic
- Basic GO enrichment of SVD poles without hypothesis-driven stratification

---

## New Hypothesis Portfolio

### Consolidation / High-confidence follow-ups (proven pipeline, no new infrastructure)

**H-A: SV4 GO biology characterization**
- Hypothesis: SV4 encodes a biologically distinct function from SV2 and SV3, with early-layer enrichment for a specific pathway
- Test: At each of 12 layers, extract top-K=52 and bottom-K=52 genes on SV4; run Fisher exact vs GO BP+CC using gene2go_all.pkl. Compare enriched terms to SV3 profile.
- Expected signal: Early layers (L0-L6) should show a specific biological theme distinct from kinase/signaling (SV3) and immune (SV3 late)
- Null: Random gene sets of same size
- Value: high | Cost: low

**H-B: STRING confidence gradient on SV3 axis**
- Hypothesis: SV3 exhibits the same confidence-to-geometry gradient as SV2 (rho≈1.0)
- Test: Same 5-quintile split on STRING scores for 3092 pairs; compute mean SV3 co-pole z-score per quintile; Spearman correlation
- Expected signal: rho≥0.9 (if SV3 is a parallel PPI encoding axis); potentially higher gradient slope given SV3 is stronger in some layers
- Null: Same gene-shuffle null; compare gradient slope to SV2 as internal calibration
- Value: high | Cost: low

**H-C: GO-PPI factorial with extreme quintile separation**
- Hypothesis: With top-20% vs bottom-20% extreme split, the GO-only signal vanishes and the PPI-only signal strengthens
- Test: Rerun H03 using top quintile GO Jaccard vs bottom quintile GO Jaccard (instead of median split). Also test: top-quintile STRING only (score>0.7 approx) vs PPI median split.
- Expected signal: GO-only z drops toward null; PPI-only z stays robust; cleaner four-quadrant separation
- Null: Gene-label shuffle
- Value: medium | Cost: low

### New geometry: spatial structure of the SV code

**H-D: SV2/SV3/SV4 joint gene clustering — are the axes encoding independent subsets?**
- Hypothesis: Genes that are "co-pole" on SV2, SV3, and SV4 simultaneously form a tight, spatially coherent submanifold in residual stream; genes are non-overlapping across axes if they encode orthogonal biology
- Test: For each layer, find the intersection of top-K=52 genes on SV2, SV3, SV4. Compute overlap fractions. Then in the full embedding space, compute pairwise cosine distances for the top-K genes on each axis; visualize as 3D scatter colored by STRING cluster membership.
- Expected signal: Low inter-axis overlap (orthogonal biological processes) but within-axis high co-localization. STRING complexes should cluster to single axes.
- Null: Random axis permutation preserving K
- Value: high | Cost: medium

**H-E: Known protein complex members collapse to single SV pole**
- Hypothesis: All subunits of a known protein complex (e.g., ribosome, proteasome, SWI/SNF, MCM) are co-located at the same pole of SV2/SV3/SV4 in early layers
- Test: Select 5-10 manually chosen protein complexes with ≥5 members in the 209-gene set. For each complex at each layer, test whether all members fall in the same top or bottom K=52 of SV2/SV3/SV4. Compute empirical p vs gene-label shuffle.
- Expected signal: Ribosomal/proteasomal members all in same pole; complexes tile different axes
- Null: Random gene sets matched on complex size; gene-label shuffle
- Value: high | Cost: medium

**H-F: Layer trajectory of individual genes in SV space**
- Hypothesis: Individual genes follow consistent directional trajectories through SV2/SV3 space across layers, and trajectory direction encodes functional class
- Test: For each gene in 209-gene set, extract its SV2 and SV3 coordinates at layers L0-L11. Plot 2D trajectories. Cluster genes by trajectory shape (DTW or Euclidean on flattened 2D vector). Test whether trajectory clusters align with GO slim categories or STRING cluster membership.
- Expected signal: Kinases follow one trajectory type, immune genes another; gene trajectory clusters reproduce STRING interaction communities
- Null: Permuted layer ordering within each gene's trajectory; random gene trajectory assignment
- Value: high | Cost: medium

### Cross-model alignment (high-risk/high-reward)

**H-G: Geneformer SV2 co-pole test on same 209 genes**
- Hypothesis: Geneformer, trained on a different corpus with a different architecture, encodes STRING PPI proximity in its top singular vectors
- Test: Load Geneformer residual stream embeddings for the same 209 genes (check data availability in reference root). Compute SVD of mean-centered gene embeddings per layer. Test SV2 co-pole z-score for STRING pairs. Compare z-score and layer profile to scGPT.
- Expected signal: STRING co-pole z>1.96 in at least some Geneformer layers; possibly different layer depth profile
- Null: Gene-label shuffle on Geneformer embeddings; comparison to scGPT null z-distribution
- Value: high | Cost: high (data availability gated)

**H-H: Procrustes alignment between scGPT and Geneformer SV2 axes**
- Hypothesis: The SV2 axis in scGPT and Geneformer is structurally aligned (Procrustes distance near zero after rotation) because both models encode the same underlying biological geometry
- Test: After confirming Geneformer embeddings available, extract SV2/SV3 coordinates for 209 genes in both models per layer. Compute orthogonal Procrustes distance. Test if Procrustes residual < shuffle null.
- Expected signal: Low Procrustes residual (high alignment) in layers where both models show PPI co-pole signal
- Null: Layer-permuted Geneformer coordinates; gene-label shuffle
- Value: high | Cost: high (data availability gated)

### Topology / PH variants

**H-I: Persistent homology of PPI-stratified gene subsets**
- Hypothesis: Genes with high-confidence STRING neighbors (top quintile) form a topologically richer manifold (higher H1 persistence, more loops) than genes with no high-confidence STRING partners
- Test: Partition 209 genes into "hub" (≥3 STRING partners score>0.7) vs "periphery" (0 or 1 partners). At each layer, compute H1 persistence of each subset after PCA(10). Compare total persistence and max persistence.
- Expected signal: Hub genes have higher H1 persistence (more complex topology) — they encode network neighborhoods, not isolated points
- Null: Random gene subset of same size
- Value: medium | Cost: medium

**H-J: Topological stability of SV2 co-pole across STRING score threshold**
- Hypothesis: The co-pole test result is robust to STRING threshold choice, with a monotonic relationship between threshold and z-score (extending H02)
- Test: Run SV2 co-pole test at 8 STRING thresholds from 0.3 to 0.9 in 0.1 increments. Report z-score at each threshold. Fit monotonic regression.
- Expected signal: Monotonic increase in z with threshold; breaks down only below 0.3 (non-informative noise)
- Null: Within-test gene-label shuffle
- Value: medium | Cost: low

**H-K: SV2 signal in GO-annotated gene subsets by ontology depth**
- Hypothesis: The GO-only geometry (z=1.89 from H03) is driven by deep/specific GO terms (small, coherent pathways), not shallow/broad terms (cell, metabolism)
- Test: Partition the GO-only group by average GO term depth (max depth across shared terms). Low depth (1-3) vs high depth (5+). Compute SV2 co-pole z per group.
- Expected signal: High-depth GO pairs show stronger SV2 signal; low-depth GO pairs near null
- Null: Gene-label shuffle within each depth group
- Value: medium | Cost: low

### Mechanistic / algorithmic motifs

**H-L: Attention-head co-activation pattern for STRING pairs**
- Hypothesis: Gene pairs with high STRING score have correlated attention patterns in scGPT transformer layers — they attend to similar context genes
- Test: For each of 12 layers, extract attention weights for the 209 genes over context. For each gene pair in STRING high-confidence quintile vs low-confidence quintile: compute cosine similarity of attention weight vectors. Compare distributions across quintiles.
- Expected signal: High STRING pairs have higher attention vector cosine similarity; mirrors SV2 co-pole structure
- Null: Random gene pair matching on STRING score; shuffle attention vectors
- Value: high | Cost: high (requires attention weight extraction infrastructure)

**H-M: SV2 geometry preserved under cell-type perturbation**
- Hypothesis: The SV2 PPI co-pole signal is cell-type invariant — it appears in lung, immune, and external domains
- Test: Run SV2 co-pole test separately on embeddings from each available cell type / tissue. Compare z-scores and sig-layer counts across domains.
- Expected signal: Consistent co-pole signal across all tested cell types; may be stronger in cell types where the 209 genes are highly expressed
- Null: Gene-label shuffle per domain
- Value: high | Cost: medium (if embeddings already available for multiple domains)

---

## Top 3 for Immediate Execution

### Slot 1: High-probability discovery candidate
**H-B: STRING confidence gradient on SV3 axis**

Rationale: H02 showed rho=1.000 for SV2. Testing SV3 with identical infrastructure will either replicate (strengthening the paper claim to "multiple orthogonal axes encode interaction confidence as geometry") or reveal axis-specific differences (also interesting). Zero new infrastructure. Expected runtime: <5 min. This directly strengthens the most important quantitative finding in the project.

### Slot 2: High-risk/high-reward candidate
**H-D: SV2/SV3/SV4 joint gene clustering — axis independence test**

Rationale: We now have confirmed signal on 3 axes. The obvious next question is: are these three axes encoding the same genes (redundant) or different gene subsets (truly orthogonal biological code)? If orthogonal, the model is doing dimensionality-efficient biological encoding — a mechanistically significant finding. If overlapping, the axes are correlated projections. Either result is informative and paper-ready. Cost: medium (requires multi-axis extraction + overlap analysis + STRING cluster assignment), but all infrastructure exists.

### Slot 3: Cheap broad-screen candidate
**H-A: SV4 GO biology characterization** (combined with **H-J: SV2 threshold sweep**)

Rationale: SV4 GO profiling is a direct copy of the SV3 GO pipeline from iter_0014. Run both H-A and H-J together: SV4 GO profile characterizes the biological content of the fourth axis, while H-J (threshold sweep at 8 values) is a <2min computation that fully characterizes the score-geometry relationship. Bundle both into one script execution.
