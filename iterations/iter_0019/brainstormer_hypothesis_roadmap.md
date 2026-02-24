# Brainstormer Hypothesis Roadmap — iter_0019

---

## Retire / Deprioritize

| Direction | Reason | Status |
|---|---|---|
| SV3/SV4 signed pole direction tests | Negative in iter_0016; no biological interpretation found | `retire_now` |
| GO term membership sweeps | Negative in iter_0011; GO is not encoded in SVD axes | `retire_now` |
| TRRUST signed direction in SVD geometry | Negative in iter_0013 and iter_0017 | `retire_now` |
| SV2 1D distance ranking (single-axis P@k) | 1.2x in iter_0018 — strictly dominated by multi-axis composite (2.18x) | `retire_now` |
| Random Gaussian null experiments | Confirmed negative in iter_0018 | `retire_now` |

---

## New Hypothesis Portfolio

### Cluster A: Multi-Axis Composite Validation and Extension

**A1 — Shuffle null for composite co-pole score**
- *Hypothesis*: The monotonic 0.86x→2.18x enrichment gradient is not an artifact of pole definition or multiple testing across axes.
- *Test*: Permute gene labels independently within each axis (SV2, SV3, SV4 separately), recompute co-polarity counts, regenerate enrichment table. Compare observed gradient to null distribution of gradients.
- *Expected signal if true*: Null gradient is flat (~1.0x at all counts); observed gradient significantly exceeds null at counts 2 and 3.
- *Null/control*: 1000 permutations of gene-to-pole assignments per axis.
- *Value*: high | *Cost*: low

**A2 — Magnitude-weighted composite score within 3-axis stratum**
- *Hypothesis*: Within the 3-axis co-polar stratum, pairs with larger total projection magnitude (|SV2| + |SV3| + |SV4|) are even more enriched for STRING edges than a binary co-pole cutoff captures.
- *Test*: For 3-axis co-polar pairs, rank by sum of absolute projections. Quintile-stratify. Report STRING enrichment per quintile.
- *Expected signal if true*: Monotonic enrichment within stratum, potentially >3x at top quintile.
- *Null/control*: Random rank assignment within the 3-axis stratum.
- *Value*: high | *Cost*: low

**A3 — STRING confidence quintile × multi-axis composite**
- *Hypothesis*: The STRING confidence quintile gradient (established as z~5 for SV2 alone in iter_0016/0017) is steeper for the 3-axis composite score than for single-axis.
- *Test*: Stratify STRING edges by confidence score quintile (Q1–Q5). For each quintile, compare multi-axis composite co-pole enrichment vs SV2-only enrichment.
- *Expected signal if true*: Q5 (high-confidence) edges show >3x enrichment for 3-axis composite; flatter curve for SV2 alone.
- *Null/control*: Random string confidence assignment.
- *Value*: high | *Cost*: low

---

### Cluster B: Attention-SVD Integration

**B1 — Attention + SVD joint predictor for TF-target prediction**
- *Hypothesis*: Combining multi-axis SVD co-pole count and attention score gives additive predictive power for TRRUST TF-target pairs, outperforming either alone.
- *Test*: For each named-gene pair: joint_score = (co_pole_count/3) + (attention > median_attention). ROC-AUC for TRRUST prediction using: (i) attention alone, (ii) composite co-pole alone, (iii) joint score.
- *Expected signal if true*: Joint AUC > max(attention AUC, SVD AUC) by >0.02.
- *Null/control*: Random score assignment.
- *Value*: high | *Cost*: low

**B2 — Head-specific attention decomposition for regulatory signal**
- *Hypothesis*: The ~2x TF regulatory co-attention signal is concentrated in a small subset of attention heads, not distributed uniformly.
- *Test*: Load per-head attention matrices (if available from scGPT). For each head, compute TRRUST enrichment separately. Identify head(s) with highest enrichment.
- *Expected signal if true*: ≤3 heads account for >80% of the regulatory enrichment signal.
- *Null/control*: Compare enrichment across heads; shuffle gene labels.
- *Value*: high | *Cost*: medium

**B3 — Cross-dataset attention transfer (lung → other tissue)**
- *Hypothesis*: The attention-based regulatory co-occurrence signal (TRRUST pairs ~2x) generalizes to scGPT embeddings from a different tissue dataset.
- *Test*: If a second scGPT processed dataset exists in the reference root, compute attention enrichment for TRRUST pairs and compare to lung result.
- *Expected signal if true*: Enrichment >1.5x in second dataset; Pearson correlation of attention scores across datasets > 0.3.
- *Null/control*: Per-dataset shuffled gene identity.
- *Value*: high | *Cost*: low (if data available)

---

### Cluster C: Cross-Model Alignment Extension

**C1 — Contextualized Geneformer embeddings (critical)**
- *Hypothesis*: Geneformer contextualized gene embeddings (post-transformer, averaged across cells) show stronger cross-model alignment with scGPT SV2 than static token embeddings (Spearman_abs > 0.446).
- *Test*: Run Geneformer inference on lung dataset; extract per-layer gene embeddings averaged across cells; compute SVD; correlate with scGPT SV2 layer-by-layer.
- *Expected signal if true*: Peak cross-model Spearman_abs > 0.6; alignment stronger at deeper GF layers.
- *Null/control*: Static GF embedding result (Spearman_abs=0.446) as baseline.
- *Value*: high | *Cost*: high (requires GF inference)

**C2 — Cycle12 GF feature metrics as contextualized GF proxy**
- *Hypothesis*: The cycle12 geneformer_feature_metrics.csv centered_cosine scores for named-gene pairs correlate with scGPT multi-axis composite co-pole membership.
- *Test*: Load cycle12 edge data; filter to named-gene pairs; compute Spearman correlation of GF centered_cosine with scGPT co-pole count (0–3); compare correlation at each scGPT layer.
- *Expected signal if true*: Spearman > 0.3 at peak layer; stronger than static GF result for the same pairs.
- *Null/control*: Shuffled gene identity in cycle12 data.
- *Value*: high | *Cost*: low (no inference needed)

**C3 — Cross-model attention comparison (GF vs scGPT)**
- *Hypothesis*: scGPT and Geneformer attention patterns are correlated for TF-target gene pairs — regulatory co-attention is model-invariant.
- *Test*: Extract Geneformer attention matrices (if accessible); compute symmetric attention for TRRUST pairs; correlate per-pair attention scores with scGPT attention scores.
- *Expected signal if true*: Spearman > 0.2 for TRRUST pairs; enrichment pattern consistent across models.
- *Null/control*: Non-TRRUST, non-STRING pairs.
- *Value*: high | *Cost*: medium

---

### Cluster D: Topological Methods (Original Mandate — Underexplored)

**D1 — Persistent homology H0/H1 of gene embedding cloud**
- *Hypothesis*: Vietoris-Rips persistent homology on pairwise cosine distances reveals topological features (H0: connected components, H1: loops) that correspond to known protein complex membership.
- *Test*: Compute pairwise cosine distance matrix for 209 named genes at each of 12 scGPT layers. Run VR persistent homology (gudhi or ripser). Record H0 birth/death times and H1 persistence. Test: do genes within the same CORUM complex cluster into same H0 component at low filtration threshold?
- *Expected signal if true*: Known complexes co-cluster at lower filtration thresholds than random gene sets; H1 loops contain functionally related gene cycles.
- *Null/control*: Random gene set of same size from CORUM universe.
- *Value*: high | *Cost*: medium

**D2 — Persistence diagram layer evolution**
- *Hypothesis*: The persistence diagram of gene embeddings changes systematically across scGPT layers, with H1 features (loops) emerging at specific layers corresponding to when PPI geometry is established.
- *Test*: Compute PH at each of 12 layers. Track total H1 persistence sum and number of bars with persistence > threshold. Correlate with layer-wise SVD PPI enrichment signal.
- *Expected signal if true*: H1 persistence peaks at layers where SV2 PPI enrichment is strongest (approximately L5–L8 based on prior iterations).
- *Null/control*: Shuffled-coordinates embeddings at each layer.
- *Value*: medium | *Cost*: medium

**D3 — Mapper topological skeleton of gene manifold**
- *Hypothesis*: The Mapper algorithm applied to scGPT gene embeddings reveals a topological skeleton with identifiable modules corresponding to protein complexes or pathway clusters.
- *Test*: Apply KeplerMapper (or scikit-tda) to 209-gene embeddings at mid-layer (L5–L8). Use SVD projections as filter function. Color nodes by STRING connectivity degree. Identify flares/branches.
- *Expected signal if true*: Mapper graph has discernible branches; genes within branches are enriched for common GO biological processes or CORUM complexes.
- *Null/control*: Mapper on shuffled gene embeddings.
- *Value*: medium | *Cost*: medium

---

### Cluster E: Local Geometry and Intrinsic Dimension

**E1 — Hub gene local intrinsic dimension**
- *Hypothesis*: Genes with high STRING degree (hub genes) occupy regions of lower local intrinsic dimension in the scGPT embedding space — they are more constrained by their many interaction partners.
- *Test*: For each gene, compute local intrinsic dimension in its k-NN neighborhood (k=20) using TWO-NN or PCA eigenvalue decay. Compute Spearman correlation between local ID and STRING degree.
- *Expected signal if true*: Spearman < -0.2 (anti-correlation: high degree → low local ID); effect increases at deeper layers.
- *Null/control*: Shuffled gene-to-embedding assignment.
- *Value*: medium | *Cost*: low

**E2 — Layer-by-layer manifold curvature (Ollivier-Ricci)**
- *Hypothesis*: scGPT gene embeddings develop increasing Ollivier-Ricci curvature on the STRING PPI graph across layers, reflecting geometric convergence toward a structured manifold.
- *Test*: For STRING edges among named genes, compute Ollivier-Ricci curvature using embedding-derived distances at each of 12 layers. Track mean curvature across layers.
- *Expected signal if true*: Curvature increases (becomes more positive or less negative) at mid-to-deep layers where PPI enrichment is strongest.
- *Null/control*: Curvature on random graph of same density.
- *Value*: medium | *Cost*: medium

---

## Top 3 for Immediate Execution

### Rank 1 — High-Probability Discovery Candidate
**A1 + A2 combined: Multi-axis shuffle null + magnitude stratification**

The 2.18x enrichment at 3-axis co-polarity is the strongest geometric result to date, but lacks a null control. Adding shuffle null is mandatory before publication. The magnitude-stratification within 3-axis stratum could push enrichment above 3x. Both are cheap (no new data, ~1 hour compute). Execute together in one script.

- Primary metric: shuffle null gradient (should be flat); quintile 5 enrichment within 3-axis stratum (expected >2.5x)
- Cost: low | Confidence: high

### Rank 2 — High-Risk/High-Reward Candidate
**D1: Persistent homology H0 of gene embedding cloud**

Topology is the project's original mandate. iter_0019 confirms the embedding manifold is multi-dimensional — exactly the regime where PH features are informative. CORUM complex co-clustering in H0 components is a concrete, testable prediction. This could be the first direct topological result in the project. Risk: requires gudhi/ripser setup and CORUM data; may show null result.

- Primary metric: CORUM complex H0 co-clustering at filtration threshold significantly below baseline
- Cost: medium | Confidence: medium-low

### Rank 3 — Cheap Broad-Screen Candidate
**B1 + C2 combined: Attention-SVD joint predictor + cycle12 GF proxy**

B1 uses existing data (attention + SVD, already computed) to test joint prediction. C2 uses cycle12 data already on disk to proxy contextualized GF alignment. Both are fast (<30 min combined). Together they advance two open threads (attention-SVD integration, cross-model alignment) at minimal cost.

- Primary metric: joint AUC vs individual AUC (B1); Spearman of GF cycle12 scores with co-pole count (C2)
- Cost: low | Confidence: medium-high
