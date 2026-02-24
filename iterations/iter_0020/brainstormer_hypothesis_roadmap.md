# Brainstormer Hypothesis Roadmap — iter_0020

---

## Retire / Deprioritize

| Direction | Reason | Status |
|---|---|---|
| TRRUST signed SVD | Negative ≥3 tries, no rescue candidate | `retire_now` |
| 1D single-SV TRRUST co-polarity | Dominated by multi-axis composite | `retire_now` |
| GO term membership sweeps | Negative ≥2 tries | `retire_now` |
| Attention for STRING prediction | AUROC=0.482, below chance, confirmed iter_0020 | `retire_now` |
| Joint attention-SVD additive predictor | No synergy confirmed H02, AUROC does not improve over best single feature | `retire_now` |
| Random Gaussian null experiments | Confirmed negative, not needed further | `retire_now` |
| SV3/SV4 signed pole direction (1D) | Negative iter_0016 | `retire_now` |
| H0 trivial connectivity (all-connected PH) | Uninformative — every filtration ends connected; replaced by H0 co-clustering at low threshold | `retire_now` |

---

## New Hypothesis Portfolio

### Cluster T: Topological Structure (H1 Follow-up)

**T1 — H1 loop membership gene-set enrichment**
- *Hypothesis*: Genes appearing in high-persistence H1 cycles (persistent loops in VR filtration of scGPT embeddings) are enriched for specific biological processes or protein complexes beyond random expectation.
- *Test*: Extract per-layer H1 birth/death pairs from stored ripser diagrams. Select top-10 most persistent H1 bars. For each, identify the cycle-generating simplices (gene pairs at birth filtration). Collect genes in those simplices. Run Fisher exact / hypergeometric test against STRING, GO-BP, CORUM. Repeat at layers 0, 6, 11.
- *Expected signal*: Top H1 cycles at mid-deep layers (L6–L8) contain genes from known complexes or shared-pathway clusters.
- *Null*: Random gene sets of same size, 1000 permutations.
- *Value*: high | *Cost*: medium

**T2 — H1 Betti number comparison: STRING-dense vs STRING-sparse subgraphs**
- *Hypothesis*: Genes with high STRING degree (top quartile) form a subgraph with more H1 persistent loops than genes with low STRING degree (bottom quartile) when analyzed as an embedding subcloud.
- *Test*: Split 209 named genes by STRING degree quartile. For top-Q4 and bottom-Q1 sets, run ripser independently at layers 0, 6, 11. Compare H1 Betti number (number of bars with lifetime > 0.01) and total H1 persistence.
- *Expected signal*: Q4 (hub genes) show significantly more H1 loops (higher connectivity in embedding space reflects cyclic geometric structure).
- *Null*: Randomly subsampled gene sets of matching size, 200 permutations.
- *Value*: high | *Cost*: low

**T3 — H0 co-clustering threshold for CORUM complex members**
- *Hypothesis*: Genes belonging to the same CORUM protein complex merge into the same H0 connected component at a significantly lower filtration threshold than random gene pairs of equal distance distribution.
- *Test*: Build pairwise distance matrix at each layer. For each CORUM complex with ≥3 named genes in the 209-gene set: record filtration threshold at which all complex members are in the same H0 component. Compare to null distribution from random gene sets of same size.
- *Expected signal*: CORUM complexes merge at threshold δ significantly below random expectation; effect strengthens at deeper layers.
- *Null*: 1000 random same-size gene subsets.
- *Value*: high | *Cost*: medium

**T4 — Layer trajectory of H1 total persistence vs STRING enrichment**
- *Hypothesis*: The layer-wise trajectory of total H1 persistence (sum of H1 bar lifetimes) is anti-correlated with the layer-wise STRING pair distance effect, indicating that as PPI geometry tightens, loop structure in the global embedding decreases.
- *Test*: Compute total H1 persistence per layer from stored PH data. Correlate with per-layer STRING distance effect size (|d|). Additionally correlate with H0 mean lifetime trajectory.
- *Expected signal*: Spearman < −0.4 between total H1 persistence and |STRING d-effect|; H0 mean lifetime and H1 persistence are positively correlated (both decrease together = progressive geometric compaction).
- *Null*: Layer-shuffled trajectories (permute layer order).
- *Value*: medium | *Cost*: low (data already in h03 artifact)

---

### Cluster G: Geometric Extensions of H03

**G1 — STRING confidence quintile × unit-sphere distance**
- *Hypothesis*: The STRING pair proximity effect (d=−0.237) scales monotonically with STRING confidence score — high-confidence interactions are geometrically closer than low-confidence ones.
- *Test*: Stratify STRING pairs by confidence quintile (Q1=lowest, Q5=highest). Compute Cohen-d effect size for each quintile vs non-STRING pairs at each layer. Compare d-values across quintiles.
- *Expected signal*: Monotone d-effect from Q1 (~−0.1) to Q5 (~−0.4); Spearman(quintile, |d|) > 0.9.
- *Null*: Randomly reassigned confidence scores, 1000 permutations.
- *Value*: high | *Cost*: low

**G2 — Distance decay profile: where does STRING proximity emerge?**
- *Hypothesis*: The STRING pair proximity effect is NOT present at layer 0 (it emerges during transformer processing), meaning the model learns PPI geometry rather than encoding it from token initialization.
- *Test*: Test whether layer 0 d-effect (already: d=−0.261) is statistically weaker than mid-layers (L5–L8, d≈−0.245). Layer 0 signal present means it's pre-existing; peak at mid layers means active learning. Use bootstrap CI at each layer to compare effect sizes.
- *Expected signal*: Layer 0 effect is already present (−0.261) but bootstrap CI overlaps with L6–L8; this would indicate the baseline model already encodes PPI. Alternatively, L0 may show weaker d than L6–L8 — cleaner narrative.
- *Null*: Permutation null already available from H03.
- *Value*: high | *Cost*: low

**G3 — TRRUST TF-target pairs: unit-sphere distance test**
- *Hypothesis*: TRRUST TF-target gene pairs are NOT closer in embedding space (unlike STRING pairs), consistent with the dissociation finding that attention (not SVD geometry) encodes regulatory signal.
- *Test*: Mann-Whitney test of pairwise distance: TRRUST pairs vs non-TRRUST non-STRING pairs at each layer. Compute Cohen-d.
- *Expected signal*: d-effect for TRRUST ≈ 0 or small negative (contrast with −0.237 for STRING); dissociation is the result.
- *Null*: Permutation null (label shuffle).
- *Value*: medium | *Cost*: low

**G4 — Ollivier-Ricci curvature on STRING PPI subgraph across layers**
- *Hypothesis*: STRING edges among the 209 named genes acquire increasingly positive Ollivier-Ricci curvature across scGPT layers, reflecting geometric organization of PPI neighborhood.
- *Test*: For each of 12 layers: compute Ollivier-Ricci curvature for STRING edges using embedding-derived distances (GraphRicciCurvature library or manual transport). Track mean curvature across layers. Compare to random graph of same density.
- *Expected signal*: Curvature increases from near-zero at L0 to positive at L8–L11; STRING edges more positive than random.
- *Null*: Random graph of same density; shuffled embeddings.
- *Value*: medium | *Cost*: medium

---

### Cluster D: Intrinsic Dimension and Manifold

**D1 — Local intrinsic dimension by STRING degree (TwoNN)**
- *Hypothesis*: Hub genes (top-quartile STRING degree) occupy regions of lower local intrinsic dimension in the scGPT embedding space than peripheral genes.
- *Test*: For each of 209 genes, compute local intrinsic dimension using TwoNN estimator in its k=20 nearest neighbors at each layer. Compute Spearman(STRING degree, local ID) across genes.
- *Expected signal*: Spearman < −0.2 (anti-correlation: high degree → low local ID); effect increases at mid-deep layers.
- *Null*: Shuffled gene-to-embedding assignment.
- *Value*: medium | *Cost*: low

**D2 — Mapper topological skeleton with STRING degree coloring**
- *Hypothesis*: KeplerMapper applied to 209-gene embeddings at mid-layer produces a topological graph with identifiable modules; nodes colored by STRING connectivity degree show spatial clustering of hubs.
- *Test*: Apply KeplerMapper at L6–L8 using first 2 PCA/SVD axes as filter. Color nodes by mean STRING degree. Identify whether hub-enriched nodes cluster in same region vs scattered.
- *Expected signal*: Hub genes (high degree) cluster in a central node or branch of the Mapper graph; peripheral genes scatter into outer branches.
- *Null*: Mapper on shuffled embeddings; comparison of node STRING enrichment to random.
- *Value*: medium | *Cost*: medium

---

### Cluster V: Validation and Publication-Quality Figures

**V1 — Layer-stratified bootstrap CI for multi-axis composite enrichment**
- *Hypothesis*: The multi-axis composite enrichment (mean 1.45x) has narrow bootstrap CIs at each layer, confirming the signal is stable and not driven by outlier layers.
- *Test*: 10,000 bootstrap samples of gene pairs at each of 12 layers. Report enrichment ± 95% CI at count=3. Plot enrichment trajectory with CI bands.
- *Expected signal*: Non-overlapping CI at count=3 vs count=0 for all layers; CI width < 0.3x.
- *Null*: Bootstrap null from permuted labels.
- *Value*: high | *Cost*: low

**V2 — Combine multi-axis composite + distance geometry into joint predictor**
- *Hypothesis*: A predictor combining multi-axis co-pole count (SVD) and unit-sphere distance achieves AUROC > 0.60 for STRING prediction — higher than either feature alone (SVD count AUROC=0.548).
- *Test*: At each layer: (i) co-pole count, (ii) unit-sphere distance, (iii) joint logistic regression. AUROC for STRING labels. Compare to SVD count baseline.
- *Expected signal*: Joint AUROC > 0.60 at mid-deep layers; distance feature improves over count alone.
- *Null*: Random feature assignment.
- *Value*: high | *Cost*: low

---

### Cluster C: Cross-Model Alignment (one new probe)

**C1 — Cycle12 GF feature metrics correlation with STRING distance**
- *Hypothesis*: The cycle12 Geneformer feature metrics (centered_cosine) for named-gene pairs correlate with scGPT unit-sphere distance, indicating cross-model agreement on PPI geometry.
- *Test*: Load cycle12 edge data; filter to 209-gene pairs; correlate GF centered_cosine with scGPT pairwise distance at each layer. Compute Spearman; test if STRING-labeled pairs show stronger negative correlation.
- *Expected signal*: Spearman(GF centered_cosine, scGPT distance) > 0.2 overall; stronger for STRING pairs.
- *Null*: Shuffled gene identities in cycle12 data.
- *Value*: medium | *Cost*: low

---

## Top 3 for Immediate Execution

### Rank 1 — High-Probability Discovery Candidate
**G1 + V2 combined: STRING confidence × distance gradient + joint SVD-distance predictor**

G1 directly extends the confirmed H03 result with a natural confidence stratification — if high-confidence STRING pairs are geometrically closer than low-confidence ones, this converts a binary proximity claim into a graded quantitative result. V2 then tests whether combining distance with multi-axis SVD count breaks the 0.60 AUROC barrier. Both tests use existing data (H03 embeddings + SVD results), take <2 hours, and yield publication-quality quantitative claims. This is the highest expected-value experiment given the H03 platform.

- Primary metrics: Spearman(confidence quintile, d-effect) > 0.8; joint AUROC > 0.60
- Cost: low | Confidence: high

### Rank 2 — High-Risk/High-Reward Candidate
**T1 + T3 combined: H1 loop gene-set enrichment + CORUM H0 co-clustering**

T1 uses the H1 diagrams already in h03_persistent_homology.json to ask whether persistent loops contain biologically meaningful gene cycles. T3 tests whether protein complex members merge into shared H0 components at lower filtration thresholds — the project's original topological mandate. Both require loading the stored PH diagrams and adding CORUM data lookup. T3 requires sourcing CORUM gene sets (downloadable) but otherwise uses existing distance matrices. If either finds enrichment, this is the strongest novel topological finding to date.

- Primary metrics: Fisher p < 0.01 for H1 cycle gene enrichment; CORUM merge threshold significantly below null
- Cost: medium | Confidence: medium

### Rank 3 — Cheap Broad-Screen Candidate
**T4 + G3: H1 total-persistence layer trajectory + TRRUST distance test**

T4 requires no new computation — compute total H1 persistence from the stored per-layer values and correlate with the stored STRING distance effects. Takes <30 minutes. G3 tests whether TRRUST pairs have a distance effect like STRING pairs (predicted: no), directly confirming the dissociation between geometry (STRING) and attention (TRRUST). Both are fast, use existing data, and either confirm or deny specific mechanistic predictions.

- Primary metrics: H1 persistence anti-correlated with |STRING d-effect| (Spearman < −0.3); TRRUST d-effect ≈ 0 vs STRING d-effect −0.237
- Cost: low | Confidence: high
