# Brainstormer Hypothesis Roadmap: iter_0052 → iter_0053+

---

## Retire / Deprioritize

| Direction | Reason | Disposition |
|-----------|--------|-------------|
| H1 Betti loops for regulatory circuit detection | Negative at L8 SV2-4; non-circuit genes have more loops | retire_now |
| SV2-4 regulatory proximity at L0-L7 | Co-expression confounded; only L8 survives residualization | retired (iter_0051) |
| Graph topology kNN purity | Negative at iters 0045-0046 | retired |
| TwoNN intrinsic dimensionality multi-layer | Inconclusive across layers | retired |
| Housekeeping gene enrichment (Eisenberg reference set) | 7/238 gene overlap; untestable | rescue_once_with_major_change — switch reference |

---

## New Hypothesis Portfolio

### SV5-7 Spectral Axis Group (build on positive H03)

**A1 — SV5-7 Bootstrap Robustness**
- Hypothesis: The SV5-7 co-expression-independent regulatory signal at L0-L2 is stable across subsamples of TRRUST pairs.
- Test: 100 stratified bootstrap resamples of TRRUST pos/neg pairs; recompute residualized rbc at L0-L2 for SV5-7; report CI and stability.
- Expected signal: rbc 95% CI excludes 0 at L0 and L1; shrinks at L2.
- Null: Bootstrap CIs straddle 0 for all layers.
- Value: high | Cost: low (same pipeline, just loop).

**A2 — Full Spectral Decay Profile (SV8-10, SV11-20)**
- Hypothesis: Regulatory signal in successive spectral bands decays monotonically with SV index at L0, or there exist non-trivial secondary peaks.
- Test: Same OLS-residualization framework applied to SV8-10 and SV11-20 projections at L0-L2. Compare rbc values across SV2-4, SV5-7, SV8-10, SV11-20.
- Expected signal: monotonic decay; or a secondary bump indicating a second regulatory axis.
- Null: All higher bands below significance threshold.
- Value: high | Cost: low.

**A3 — SV5-7 Top-Loading Gene Biological Annotation**
- Hypothesis: Genes with extreme SV5-7 loadings at L0 are enriched for a coherent biological function (e.g., immune signaling, cell-cycle regulation) that explains why this subspace encodes regulatory proximity.
- Test: Extract top/bottom 5% SV5-7 loadings at L0 (index 4,5,6 each). Fisher exact GO-slim enrichment against all 2039 nonzero genes. Also check TRRUST TF overlap.
- Expected signal: Specific GO terms (p<0.05 FDR) with recognizable immune/regulatory biology.
- Null: No GO term survives FDR correction.
- Value: high | Cost: low.

**A4 — Layer × Subspace Interaction: Is L0-L2 Signal from Embedding vs Attention?**
- Hypothesis: The SV5-7 regulatory signal at L0-L2 corresponds to scGPT input embedding geometry (before transformer attention), while SV2-4 L8 signal emerges from deep attention layers.
- Test: If scGPT layer 0 is the raw gene embedding (pre-attention), compare SV5-7 rbc at L0 vs L1 vs L2 — if the signal decays over layers, it originates in the embedding. Use the existing H03 CSV data.
- Expected signal: rbc peaks at L0 and decays, consistent with embedding-layer origin.
- Null: Signal is flat across L0-L2, suggesting it is created (not destroyed) by early transformer operations.
- Value: high | Cost: zero (re-analysis of existing H03 CSV).

---

### Manifold Geometry Group

**B1 — Local Linearity of TRRUST Pairs vs. Non-Pairs in SV5-7**
- Hypothesis: TRRUST gene pairs (TF→target) lie on locally linear patches of the SV5-7 manifold at L0, whereas non-pairs span regions of higher local curvature.
- Test: For each TRRUST positive pair, compute the geodesic distance on a kNN graph (k=15) in SV5-7 L0 space vs straight-line Euclidean distance; the ratio is a local nonlinearity proxy. Compare pos vs neg pairs (MW test).
- Expected signal: pos pairs have lower nonlinearity ratio (more linear local neighborhood).
- Null: No difference in ratio.
- Value: medium | Cost: medium.

**B2 — Curvature Anisotropy Along Regulatory vs Non-Regulatory Directions**
- Hypothesis: The spectral axes encoding regulatory information (SV2-4 at L8, SV5-7 at L0-L2) have higher "functional anisotropy" — variance explained per unit regulatory information — than non-regulatory axes.
- Test: For each SVD axis group, compute: (rbc_regulatory_signal) / (fraction_variance_explained). Compare ratio across SV1, SV2-4, SV5-7, SV8-10.
- Expected signal: SV5-7 and SV2-4 have highest ratio; SV1 lowest despite explaining most variance.
- Null: ratio correlates with variance explained.
- Value: medium | Cost: low (re-uses existing data).

**B3 — Geodesic Distance in Full 512-dim Embedding vs SV Projections**
- Hypothesis: Regulatory gene pairs have systematically shorter shortest-path distances on the kNN graph of full 512-dim embeddings at L8 than non-regulatory pairs, independent of Euclidean distance.
- Test: Build kNN graph (k=20) on L8 full embeddings; compute shortest-path geodesic distance for 589 pos and 1523 neg TRRUST pairs. MW test.
- Expected signal: pos pairs have shorter geodesic (more connected) paths, rbc>0.
- Null: geodesic distances track Euclidean (no difference).
- Value: medium | Cost: medium.

---

### Persistent Homology (new variants, not H1 loops)

**C1 — H0 Persistence (Cluster Merging) of TRRUST Subgraph at L0**
- Hypothesis: The TRRUST-positive subgraph (genes with known regulatory interactions) merges into a single H0 connected component at a shorter Vietoris-Rips filtration radius than random subgraphs of equal size, indicating tighter clustering in SV5-7 space at L0.
- Test: Ripser H0 on TRRUST-annotated genes projected to SV5-7 at L0; compare death of last H0 bar (full connectivity radius) vs 1000 random gene sets of equal size.
- Expected signal: TRRUST annotated genes merge earlier (lower max H0 bar).
- Null: No difference in max H0 bar.
- Value: medium | Cost: low.

**C2 — Persistent Entropy of SV5-7 Point Clouds Across Layers**
- Hypothesis: The topological complexity (persistent entropy of H0 bars) of the SV5-7 gene cloud decreases from L0 to L3 as the regulatory structure is transformed, then stabilizes.
- Test: Ripser H0 on all 2039 nonzero genes in SV5-7 at each layer; compute persistent entropy = -Σ(p_i log p_i) where p_i = lifetime_i / total_lifetime. Plot across layers.
- Expected signal: monotonic decrease L0→L3, plateau L3-L11.
- Null: Flat or increasing entropy across layers.
- Value: medium | Cost: low.

---

### Cross-Model / Biological Anchoring Group

**D1 — SV5-7 Signal Replication in Different Gene Subsets**
- Hypothesis: The SV5-7 regulatory proximity signal at L0 is not specific to the 2039 nonzero-embedding gene set but holds for the full TSP14 dataset (4941 genes), confirming it is a property of scGPT's embedding geometry, not a sampling artifact.
- Test: Run H03 protocol on all 4941 genes with non-NaN embeddings; compare rbc_residual at L0 for TRRUST pairs.
- Expected signal: similar or stronger rbc (more power with more genes).
- Null: rbc drops or loses significance, suggesting 2039-gene result was a sampling artifact.
- Value: high | Cost: low.

**D2 — SV1 Gene Identity Characterization via STRING Degree**
- Hypothesis: SV1-high genes (top 20% at L0) are enriched for low-degree (few known interactions) genes in the STRING network, consistent with them being poorly annotated or peripheral genes, not housekeeping genes.
- Test: Fetch STRING degree (number of interaction partners) for all 2039 gene names. Split by SV1-high/low. Compare STRING degree distribution (MW test). Also check fraction of genes with STRING degree = 0.
- Expected signal: SV1-high genes have significantly lower STRING degree (fewer interactions).
- Null: No STRING degree difference.
- Value: high | Cost: medium (STRING API fetch required).

**D3 — TRRUST Pair Asymmetry: Does TF→Target Direction Matter for SV5-7 Distance?**
- Hypothesis: In TRRUST pairs (TF, target), the direction matters — the TF's SV5-7 position at L0 is closer to its direct targets than an equivalent TF is to random non-target genes of the same expression level.
- Test: For each TF in TRRUST, compute mean SV5-7 distance to its annotated targets vs mean distance to expression-matched non-targets. Paired Wilcoxon test.
- Expected signal: TF-to-target distance significantly lower than TF-to-non-target.
- Null: No difference after expression matching.
- Value: high | Cost: medium.

**D4 — Cell-Type Identity Gradient Along SV5-7 at L0**
- Hypothesis: If scGPT's early layers preserve cell-type–specific regulatory programs, SV5-7 loadings at L0 should correlate with immune cell-type identity markers (e.g., CD3, CD19, CD14) ranked by their TSP14 expression specificity.
- Test: Compute per-gene SV5-7 L0 loading norm. Rank genes by immune specificity score (e.g., fraction of CD3+ cells expressing gene). Spearman correlation of SV5-7 norm vs specificity score.
- Expected signal: Positive Spearman r — cell-type specific genes have higher SV5-7 norm.
- Null: r ≈ 0.
- Value: medium | Cost: medium.

---

## Top 3 for Immediate Execution

### Candidate 1 — High-Probability Discovery: A1 + A2 (SV5-7 Bootstrap + Spectral Decay)
Combine into a single script: (a) 100-resample bootstrap of H03 at L0-L2 for SV5-7; (b) extend residualization to SV8-10 and SV11-20.
Rationale: Validates the most important positive result and characterizes whether the spectral signal is unique to SV5-7 or broader. Extremely low computational cost; uses identical pipeline. Will either confirm a major finding or expose it as noise.

### Candidate 2 — High-Risk/High-Reward: D3 (TF→Target Direction Asymmetry in SV5-7)
Tests whether the regulatory encoding in SV5-7 is directional — a genuinely novel mechanistic claim about how scGPT represents gene regulation. If positive, this is publishable as a standalone result (directed regulatory geometry). If negative, it constrains the interpretation of H03 to symmetric co-regulation rather than directed control.

### Candidate 3 — Cheap Broad Screen: A4 + B2 (Free Re-Analysis of Existing Data)
A4: Re-examine the H03 CSV to test whether rbc decays monotonically from L0→L2, probing embedding-layer vs attention-layer origin — zero new computation.
B2: Compute (rbc_signal / variance_explained) ratio across SV bands to quantify functional information density — all numbers already available.
Together these cost nothing and answer interpretive questions that frame the paper narrative.

---

## Direction Status Table

| Family | Direction | Status |
|--------|-----------|--------|
| persistent_homology | H1 Betti loops for regulatory circuits | retired |
| module_structure | SV2-4 co-expression-confounded regulatory proximity (L0-L7) | retired |
| module_structure | SV5-7 regulatory signal at L0-L2 (H03) | **ACTIVE — top priority** |
| module_structure | SV2-4 regulatory signal at L8 | active (secondary) |
| intrinsic_dimensionality | SV1 anti-TF axis | active (needs better annotation) |
| intrinsic_dimensionality | Housekeeping enrichment (Eisenberg) | rescue_once (switch reference to STRING/GO) |
| manifold_geometry | Geodesic kNN on full embeddings | new — candidate |
| persistent_homology | H0 persistence of TRRUST subgraph | new — candidate |
| biological_anchoring | STRING degree for SV1-high characterization | new — candidate |
| biological_anchoring | TF→Target directional asymmetry in SV5-7 | new — **high reward candidate** |
