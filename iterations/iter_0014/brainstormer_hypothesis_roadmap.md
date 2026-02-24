# Brainstormer Hypothesis Roadmap — iter_0014

---

## Retire / Deprioritize

| Hypothesis | Reason | Decision |
|---|---|---|
| SV1 PPI axis | z=0.39, 2/12 layers — noise | **RETIRED** |
| Repression anti-pole | 0/12 sig, z=-1.41 mean | **RETIRED** |
| Layer-depth monotonicity (Spearman) | Marginal SV2 (p=0.059), non-monotonic SV3; no new insight from repeating | **DEPRIORITIZE** — one consolidation mention in paper only |
| Direct TRRUST TF→target geometry | Tested iter_0010-0012 with mixed/inconclusive results; superseded by PPI direction | **RETIRE unless rescued with orthogonal angle** |
| Bootstrap act/rep z-score comparison | Negative, superseded | **RETIRED** |

---

## New Hypothesis Portfolio

### Group A: SV3 Downstream Anchoring

**A1: STRING confidence gradient correlates with SV2 Euclidean distance**
- Hypothesis: For interacting gene pairs, higher STRING confidence scores predict smaller Euclidean distance in SV2 (or SV2+SV3) space — i.e., scGPT encodes interaction strength as geometric proximity.
- Test: For all 1022 STRING pairs (score≥0.7), compute SV2 distance per layer. Spearman correlation between STRING score and SV2 distance. Repeat at score≥0.4 (3364 pairs) for power.
- Signal if true: Spearman rho < −0.15 with p<0.05 in ≥6/12 layers (negative = higher confidence → closer).
- Null: rho ≈ 0 (no gradient).
- Value: **high** | Cost: **low** (all data already cached)

**A2: SV3 pole membership predicts known kinase inhibitor target membership**
- Hypothesis: Genes in the SV3 TOP pole (kinase/signaling early layers, DNA-binding late) are enriched for known drug targets of kinase inhibitors catalogued in DrugBank or ChEMBL.
- Test: Extract TOP pole genes at L3 (peak SV3 layer). Map to DrugBank kinase inhibitor targets. Fisher exact test vs background 209 genes. Shuffle null N=500.
- Signal if true: Fisher p<0.05, emp_p<0.05, OR>2 for drug-target enrichment.
- Null: No enrichment (OR≈1).
- Value: **high** | Cost: **medium** (requires DrugBank mapping)

**A3: SV3 projection at L3 separates cancer driver genes from non-drivers**
- Hypothesis: Known cancer driver genes (from COSMIC Cancer Gene Census) project to the TOP pole of SV3 at L3 (kinase/signaling cluster), distinguishing them from non-driver genes.
- Test: Intersect 209 genes with COSMIC CGC. SV3 projection distribution for CGC genes vs non-CGC. Mann-Whitney U test + effect size. Shuffle null N=500.
- Signal if true: CGC genes significantly higher SV3 projection (p<0.05, large effect size).
- Null: No difference in SV3 distribution.
- Value: **high** | Cost: **low** (COSMIC list is public/static)

### Group B: Joint SV2+SV3 Geometry

**B1: Joint SV2+SV3 2D space reveals community clustering exceeding either axis alone**
- Hypothesis: Genes belonging to the same STRING Louvain community are more spatially clustered in the SV2×SV3 2D projection than in SV2 or SV3 alone, indicating synergistic multi-axis PPI geometry.
- Test: Run Louvain community detection on STRING (score≥0.7). For each community, compute average pairwise distance in SV2 only, SV3 only, and SV2+SV3 jointly. Compare to community-label shuffle null. Test at L3 (peak SV3) and L1 (peak SV2).
- Signal if true: Intra-community z-score in 2D > max(SV2-only z, SV3-only z) by >0.5 in ≥4/12 layers.
- Null: 2D distance ratio matches independent axis combination (no synergy).
- Value: **high** | Cost: **medium**

**B2: SV2–SV3 axis angle trajectory across layers is non-random**
- Hypothesis: The angle between SV2 and SV3 (as gene-rank vectors) changes systematically across layers 0→11, reflecting functional geometry evolution with transformer depth.
- Test: Compute Spearman correlation between gene SV2-rank and SV3-rank vectors at each layer (as proxy for axis alignment). Test Spearman rho across layers 0-11 for monotonicity.
- Signal if true: Spearman rho (layer vs axis-angle) > 0.6 or < −0.6, p<0.05.
- Null: Axis angle is constant or random across layers.
- Value: **medium** | Cost: **low**

**B3: Intrinsic dimensionality of PPI gene subgraph embedding is lower than background**
- Hypothesis: The subgraph of high-confidence PPI partners (score≥0.7) occupies lower intrinsic dimension in embedding space (all 12 layers) than a degree-matched random gene set, reflecting geometric compression by scGPT.
- Test: Use TwoNN or PCA variance-explained for PPI partner pairs vs random same-sized subsets (N=500 permutations). Report dimensionality ratio per layer.
- Signal if true: PPI gene subset has significantly lower intrinsic dimension (ratio < 0.8, p<0.05).
- Null: PPI subset has same intrinsic dimension as random.
- Value: **high** | Cost: **medium**

### Group C: GO Co-annotation vs PPI as Geometry Predictors

**C1: GO biological process co-annotation predicts SV2 co-pole as well as STRING PPI**
- Hypothesis: Genes sharing ≥2 GO BP terms are co-localized in SV2 at rates comparable to STRING PPI pairs, suggesting scGPT geometry encodes functional similarity beyond physical interaction.
- Test: Define GO-coannot pairs (≥2 shared GO BP terms, size 3-60). Compute SV2 co-pole Fisher z-score for GO pairs. Compare to STRING PPI z-scores per layer.
- Signal if true: GO-coannot z > 2.5 in ≥6/12 layers, within 1 SD of STRING z.
- Null: GO co-annotation z < 2 (PPI is unique predictor).
- Value: **high** | Cost: **medium** (requires GO term matrix computation)

**C2: Combining GO co-annotation + STRING PPI as predictor additively improves geometry signal**
- Hypothesis: Gene pairs that are BOTH GO co-annotated AND STRING PPI partners show the strongest SV2 co-pole enrichment, indicating independent biological signals compound in scGPT geometry.
- Test: Split 209-gene pairs into four groups: neither, GO only, PPI only, both. Compare SV2 co-pole z-scores across the four groups per layer.
- Signal if true: "Both" group z > max(GO-only z, PPI-only z) by >0.5, consistent across layers.
- Null: "Both" z ≈ max of two single predictors (no additivity).
- Value: **high** | Cost: **medium** (depends on C1 data)

### Group D: Higher SVD Axes

**D1: SV4 and SV5 show significant PPI co-pole enrichment**
- Hypothesis: Beyond SV2 and SV3, scGPT embedding encodes additional PPI axes in SV4/SV5 with distinct biological identities.
- Test: Compute SV4, SV5 co-pole Fisher z-scores (K=52) for STRING pairs (1022 pairs, score≥0.7) at all 12 layers. Run GO enrichment on SV4/SV5 poles.
- Signal if true: SV4 or SV5 z > 3 in ≥6/12 layers; distinct GO terms from SV2/SV3.
- Null: SV4/SV5 z < 2 (geometry confined to first 3 axes).
- Value: **high** | Cost: **low** (same infrastructure as H01/iter_0013)

**D2: SV4+ axes capture metabolic vs signaling biology not in SV2/SV3**
- Hypothesis: If SV4+ axes show PPI enrichment, their poles will map to metabolic processes (oxidative phosphorylation, glycolysis) orthogonal to the kinase/immune axes already found in SV2/SV3.
- Test: Contingent on D1 being positive. GO enrichment on SV4/SV5 poles (same method as H01 iter_0014).
- Signal if true: Metabolic GO terms appear as top hits for SV4 or SV5.
- Null: Same GO terms as SV2/SV3 or no enrichment.
- Value: **high** | Cost: **low** (contingent on D1)

### Group E: Topological / Manifold

**E1: Persistent homology H0 barcodes of PPI gene subsets are longer than random**
- Hypothesis: The VR persistence diagram H0 for the subset of high-confidence PPI partners (score≥0.7) at SV2+SV3 distances has longer birth-to-death lifetimes than random subsets, reflecting geometric clustering.
- Test: Compute pairwise Euclidean distances in SV2+SV3 for PPI genes vs matched-size random gene subsets (N=200 permutations). Compare mean H0 persistence lifetime via permutation test.
- Signal if true: Mean H0 persistence of PPI subset significantly greater than random (z>2).
- Null: Persistence matches random (z<1.96).
- Value: **medium** | Cost: **medium** (requires ripser/gudhi)

**E2: Geodesic distances in kNN graph of 209-gene SV2+SV3 embedding correlate with STRING confidence**
- Hypothesis: In a kNN graph built from SV2+SV3 distances, shortest-path (geodesic) distances between gene pairs are anti-correlated with STRING interaction confidence, beyond Euclidean distance.
- Test: Build kNN graph (k=5) on SV2+SV3 space. Compute all-pairs shortest paths for 209 genes. Spearman correlation between geodesic distance and STRING score. Compare rho to Euclidean-based rho (test if geodesic is a better predictor).
- Signal if true: Geodesic rho has larger magnitude than Euclidean rho in ≥6/12 layers.
- Null: Geodesic rho ≤ Euclidean rho (manifold structure doesn't add information).
- Value: **medium** | Cost: **medium**

**E3: Layer-to-layer geometric continuity — SV2 subspace angle between consecutive layers**
- Hypothesis: The SV2 subspace (top singular vector) rotates gradually across layers 0→11, with larger rotations at certain depth transitions, identifying "phase transitions" in scGPT's representational geometry.
- Test: Compute principal angles between the rank-1 SV2 subspaces at layers i and i+1. Plot rotation magnitude per layer transition. Test if rotation magnitude is non-uniform (variance test, Bartlett).
- Signal if true: Non-uniform rotation profile with identifiable jump points.
- Null: Uniform rotation (constant principal angle across transitions).
- Value: **medium** | Cost: **low**

---

## Top 3 for Immediate Execution

### Candidate 1: High-Probability Discovery — D1 (SV4/SV5 axes)
**Rationale**: Uses exactly the same pipeline as iter_0013 H01 and iter_0014 H01. Zero infrastructure cost. We confirmed SV2 and SV3 are PPI axes with distinct biology. SV4+ is the natural next question. If positive, substantially expands the paper's geometric characterization. High prior probability given SV2/SV3 both hit.

### Candidate 2: High-Risk/High-Reward — A1 (STRING confidence gradient)
**Rationale**: Tests a continuous prediction: does scGPT encode interaction *strength* not just *presence*? All data already cached (1022 pairs + 3364 pairs with scores). If the confidence gradient correlates with Euclidean distance, this transforms the story from "binary PPI co-pole" to "continuous interaction geometry" — a more powerful mechanistic claim.

### Candidate 3: Cheap Broad-Screen — C1 (GO co-annotation vs PPI)
**Rationale**: Defines GO-coannot pairs from existing GO term data, then runs the same co-pole test. Answers the fundamental question: is scGPT geometry specific to physical interactions or general to functional similarity? This directly addresses the most obvious interpretability question about what the geometry *means*.

---

## Full Hypothesis Priority Table

| ID | Name | Value | Cost | Priority |
|---|---|---|---|---|
| D1 | SV4/SV5 PPI co-pole test | high | low | **P1** |
| A1 | STRING confidence gradient | high | low | **P2** |
| C1 | GO co-annotation vs PPI | high | medium | **P3** |
| A3 | SV3 cancer driver enrichment | high | low | P4 |
| B1 | Joint SV2+SV3 community clustering | high | medium | P5 |
| C2 | GO+PPI additivity | high | medium | P6 (needs C1) |
| D2 | SV4+ metabolic biology | high | low | P7 (needs D1) |
| B3 | Intrinsic dimensionality PPI subset | high | medium | P8 |
| A2 | SV3→kinase inhibitor drug targets | high | medium | P9 |
| B2 | SV2–SV3 angle trajectory | medium | low | P10 |
| E3 | Layer-to-layer SV2 subspace rotation | medium | low | P11 |
| E1 | PH H0 barcodes PPI subset | medium | medium | P12 |
| E2 | Geodesic vs Euclidean predictor | medium | medium | P13 |
