# Brainstormer Hypothesis Roadmap: iter_0032 → iter_0033+

**Date:** 2026-02-23

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| STRING PPI distance gradient (generic) | AUROC=0.494 on 195-gene OOV-corrected set (iter_0031) | **RETIRED** |
| Two-trajectory / diverging gene cluster | OOV artifact (14 genes with zero L0 norm), not biology | **RETIRED** |
| Repression anti-pole on SV2 | mean z=-1.41, 0/12 layers (iter_0013) | **RETIRED** |
| TF hub centrality on PC1 | AUROC=0.492 (iter_0027), no geometric centrality of TFs | **RETIRED** |
| PC1 = T-cell vs APC axis | Both marker sets cluster at same pole at L11 (iter_0028) | **RETIRED** |
| Dorothea stable AUROC across all layers | OOV-inflated artifact; revised to early-layer only (iter_0032) | **RETIRED (revised)** |
| Precision@k on SV2 distance | Only 1.20× enrichment — co-pole test is better framing | Deprioritize |
| H1 persistence layer profiling (basic) | Covered; monotonic decline rho=-0.916 confirmed | Deprioritize |

---

## New Hypothesis Portfolio

### CLUSTER A: Explain the 2-Community Split (actionable, close to existing data)

**A1 — Community membership predicts secreted vs membrane-bound localization**
- Hypothesis: The L11 2-community split encodes protein secretion vs surface/intracellular biology — community 1 (B-cell enriched) holds more secreted/paracrine molecules, community 0 holds more intracellular/surface receptors.
- Test: Annotate 195 in-vocab genes with UniProt subcellular location (secreted/extracellular vs membrane vs nuclear/cytoplasmic). Fisher exact test for each location class vs community assignment. Null: 100 random community permutations of same sizes.
- Expected signal: Secreted/ECM genes enriched in one community, OR ≥ 2.0, p < 0.05.
- Control: random partition of same sizes.
- Value: **high** | Cost: **low** (annotation lookup + Fisher test)

**A2 — Community membership predicts TF vs non-TF (functional class split)**
- Hypothesis: The 2-community split encodes transcription factor vs effector/cytokine functional class — TFs dominate one community, cytokines/ligands dominate the other.
- Test: Annotate 195 genes as TF (from TRRUST TF list + GO molecular function TF terms) vs non-TF. Fisher test per community. Also test: cytokines/chemokines vs receptors vs signaling enzymes vs TFs as 4-way contingency.
- Expected signal: OR ≥ 2 for TF vs non-TF partitioning into separate communities.
- Control: random gene label shuffling.
- Value: **high** | Cost: **low**

**A3 — Community boundary is defined by PC1 sign**
- Hypothesis: The 2 communities correspond to positive vs negative PC1 projections at L11 (PC1 explains 26% variance; communities might be PC1+ and PC1− groups).
- Test: Project 195 in-vocab genes on PC1 at L11. Compute AUROC for community assignment from PC1 sign. Kendall tau between community label and ranked PC1 score.
- Expected signal: AUROC ≥ 0.80 for PC1 sign predicting community assignment; if true, all biology of PC1 is community biology.
- Control: PC2, PC3 sign as null comparisons.
- Value: **high** | Cost: **low** (already have both pieces)

### CLUSTER B: Early-Layer Regulatory Geometry (new positive signal to extend)

**B1 — Activation pairs decay faster than repression pairs in Dorothea**
- Hypothesis: The AUROC decay (rho=-0.853 for high-conf pairs) is driven specifically by activation-type regulation, consistent with TRRUST asymmetry (iter_0023). Repression pairs in Dorothea should show flat AUROC or slight opposite trend.
- Test: Separate 205 high-conf Dorothea pairs into activation (stimulates), repression (inhibits), and unknown (flags A or B but unspecified). Compute AUROC per layer for each subset. Spearman rho(layer, AUROC) per subset.
- Expected signal: Activation subset: rho < -0.7; Repression: rho ≈ 0 or positive; confirming that early-layer regulatory geometry encodes activating co-regulation.
- Control: random pairs matched in size.
- Value: **high** | Cost: **low** (Dorothea file already loaded with direction field)

**B2 — Early-layer regulatory geometry clusters TFs by common target overlap**
- Hypothesis: At L0–L2 (peak regulatory proximity), pairs of TFs that share many targets are geometrically closer than TF pairs with no target overlap (co-regulon proximity).
- Test: For TFs in the 195-gene set, compute Jaccard overlap of their Dorothea high-conf target sets (among all 195 genes). Spearman(Jaccard_TF_targets, L2_distance_TF_TF) at L0, L1, L2, L11 as control.
- Expected signal: Significant negative correlation at L0–L2, not significant at L11.
- Control: permuted Jaccard null (shuffle TF-gene assignments).
- Value: **medium** | Cost: **low**

**B3 — GO CC proximity OOV-corrected replication**
- Hypothesis: The GO CC signal (rho=0.106, iter_0024) was computed on the full 209-gene set. On 195 in-vocab genes it should survive — GO CC annotates subcellular localization which is a property of the protein itself independent of expression, so OOV removal should not eliminate it.
- Test: Re-run GO CC + GO BP Jaccard vs L2 distance Spearman at all 12 layers for 195 in-vocab genes. Report rho(layer) profiles and compare to prior 209-gene results.
- Expected signal: rho ≈ 0.08–0.12 at multiple layers, similar to iter_0024.
- Control: chromosomal proximity (already negative at AUROC=0.515, iter_0025).
- Value: **high** | Cost: **low** (straightforward replication)

### CLUSTER C: PC1 Axis Biological Identity (needs more power)

**C1 — Expanded B-cell program on PC1 using all B-cell-related genes in vocab**
- Hypothesis: The PC1 negative pole is a B-cell/plasma cell program axis. Only 3 canonical B-cell markers were in the current 195-gene set (CD19, MS4A1, CD79A). Additional B-cell related genes (JCHAIN, MZB1, IGLL1, CD79B, BLNK, CD27, PAX5 if in vocab, IGHM, IGHG) may be present in the broader 8181-gene vocabulary file.
- Test: Scan the scGPT vocabulary for all B-cell program genes (using a curated list of ~30 B-cell/plasma cell markers). Project in-vocab B-cell markers on PC1 at L5, L11. Mann-Whitney AUROC vs random background genes. Also test T-cell vs myeloid programs.
- Expected signal: B-cell program AUROC < 0.3 (negative pole), T-cell > 0.6 (positive pole), myeloid intermediate.
- Control: random gene sets of same sizes.
- Value: **high** | Cost: **medium** (need to mine vocabulary file)

**C2 — PC2 separates myeloid from lymphoid markers**
- Hypothesis: While PC1 separates B-cell vs TF/T-cell, PC2 at L11 captures myeloid (monocyte/macrophage/dendritic cell) biology as an orthogonal axis.
- Test: Project 195 in-vocab genes on PC2 at L0, L5, L11. Test myeloid markers (e.g., CD14, ITGAM, FCGR3A, LYZ, S100A8 — check which are in vocab) vs lymphoid markers AUROC on PC2. Compare PC1 vs PC2 AUROC profiles for each cell type.
- Expected signal: Myeloid markers at one PC2 pole, lymphoid at other, AUROC ≥ 0.65.
- Control: PC3 as negative (should show no cell-type separation).
- Value: **medium** | Cost: **low**

### CLUSTER D: Cross-Model Validation

**D1 — Geneformer static embedding shows early-layer regulatory proximity analog**
- Hypothesis: The Geneformer static token embedding (iter_0019) encodes regulatory proximity in the same way as scGPT L0–L2. Since Geneformer's static embedding is its "layer 0" equivalent, Dorothea high-conf pairs should show elevated proximity (AUROC > 0.52) in Geneformer static space, matching the scGPT early-layer finding.
- Test: Load Geneformer static token embedding for 195 in-vocab named genes (use the 207 found in iter_0019, intersect with current 195). Compute Dorothea high-conf pair AUROC vs 5000 random pairs. Compare to scGPT L0 result (AUROC=0.564).
- Expected signal: Geneformer AUROC ≥ 0.53, confirming cross-model regulatory geometry in early representation.
- Control: Geneformer AUROC for Dorothea low-conf pairs (should be ~0.50).
- Value: **high** | Cost: **low**

**D2 — Geneformer community structure: does the 2-community split replicate?**
- Hypothesis: The 2-community structure at L11 (z=34) in scGPT reflects something about gene biology, not just scGPT training, and should appear in Geneformer's static embedding as well.
- Test: Apply greedy modularity community detection to the k=10 kNN graph of 195 in-vocab genes in Geneformer static embedding. Test modularity vs 100 random partition nulls. Check if community assignments align with scGPT communities (adjusted rand index).
- Expected signal: Modularity ≥ 0.3 above null; ARI ≥ 0.4 between scGPT L11 communities and Geneformer communities.
- Control: random graph with same degree sequence (configuration model).
- Value: **high** | Cost: **low**

### CLUSTER E: Manifold Geometry Deep Dives

**E1 — Local intrinsic dimension profile (TwoNN estimator) vs layer**
- Hypothesis: The PR collapse measures global effective dimensionality. Local intrinsic dimension (TwoNN estimator) measures the dimensionality of individual gene neighborhoods and may reveal that different gene subsets have different ID profiles — e.g., TFs may live in higher-ID neighborhoods than non-TFs at L0, but collapse to the same low ID by L11.
- Test: Compute TwoNN intrinsic dimension for each gene using its k=20 nearest neighbors at L0, L5, L11. Report mean ID per layer and variance. Mann-Whitney ID for TFs vs non-TFs, community 0 vs community 1 genes.
- Expected signal: Community 1 (B-cell enriched) has lower ID than community 0 at L11; both communities collapse from L0 to L11.
- Control: random gene partitions of same sizes.
- Value: **medium** | Cost: **medium**

**E2 — Persistent homology H0 (connected components) across layers on 195 in-vocab**
- Hypothesis: The kNN spectral gap collapse and 2-community structure both suggest the gene manifold becomes more disconnected with depth. H0 persistent homology captures the birth/death of connected components as a filtration parameter varies — it should show more long-lived H0 components at late layers than early layers.
- Test: Compute H0 persistence (ripser) on pairwise L2 distances for 195 in-vocab genes at each of 12 layers. Count H0 bars with lifetime > 20th percentile. Spearman(layer, n_components). 20-shuffle feature null.
- Expected signal: rho > 0 for n_components vs layer (more disconnection at late layers), consistent with spectral gap.
- Control: feature-shuffled embeddings.
- Value: **medium** | Cost: **medium**

**E3 — Geodesic vs Euclidean distance divergence on kNN graph across layers**
- Hypothesis: As the spectral gap decreases (manifold becomes more modular), geodesic distances on the kNN graph should increasingly diverge from Euclidean distances for cross-community pairs. This measures the "tunneling cost" across community boundaries.
- Test: At each layer, compute kNN graph (k=10) and shortest-path graph distance between all gene pairs. Spearman(Euclidean, Geodesic) per layer. Also compute ratio Geodesic/Euclidean for within-community vs cross-community pairs.
- Expected signal: Within-community Geodesic/Euclidean ≈ 1.0 at all layers; cross-community ratio increases with layer (tunneling becomes more costly at late layers).
- Control: random graph with same adjacency.
- Value: **medium** | Cost: **medium**

### CLUSTER F: Mechanistic / Algorithmic Signatures

**F1 — Attention head specialization: do specific heads encode regulatory vs PPI proximity?**
- Hypothesis: The mechanistic dissociation (attention = TF regulatory, SVD = PPI proximity, iter_0018) may be head-specific. Individual attention heads may be more specialized than the aggregated attention matrix.
- Test: Load per-head attention matrices from scGPT (if available in attention_scores.npy — verify shape). For each head separately, compute TRRUST activation AUROC and STRING PPI AUROC. Identify which heads preferentially encode regulation vs PPI.
- Expected signal: 2–4 heads with TRRUST AUROC ≥ 0.65 and STRING AUROC ≈ 0.50; other heads with opposite profile.
- Control: random attention matrix of same sparsity.
- Value: **high** | Cost: **medium** (depends on artifact format)

**F2 — Cell-type marker clustering OOV-corrected replication (priority replication)**
- Hypothesis: Cell-type marker clustering (AUROC=0.851, iter_0022) was computed on the 209-gene set including OOV genes. On 195 in-vocab genes, the signal should survive as it uses canonical immune markers that are overwhelmingly in-vocabulary.
- Test: Re-run cell-type marker separation test (T-cell, B-cell, fibroblast, epithelial panels) on 195 in-vocab genes per layer. Keep same markers (drop any OOV ones). Mann-Whitney AUROC within-type vs cross-type pairs.
- Expected signal: AUROC ≥ 0.80 at most layers; signal preserved.
- Control: random gene partitions of same sizes.
- Value: **high** | Cost: **low** (priority replication of core claim)

---

## Top 3 for Immediate Execution

### Slot 1: High-Probability Discovery Candidate
**F2 — Cell-type marker clustering OOV-corrected replication + A3 (PC1 vs community alignment)**
Combined execution: (1) Re-run cell-type marker AUROC on 195 in-vocab genes — high probability of confirming core claim. (2) In the same script, test whether PC1 sign predicts community membership (AUROC for PC1 sign vs community label). This resolves two open threads at once.
- Why now: The cell-type clustering result (AUROC=0.851) is the paper's strongest claim and has not been validated on the OOV-corrected gene set. Without this, the claim carries an asterisk. Also, A3 costs zero additional compute given existing embeddings.
- Expected time: Low (single script, existing embeddings).

### Slot 2: High-Risk / High-Reward Candidate
**B1 — Activation vs repression decay rate asymmetry in Dorothea**
The TRRUST activation/repression asymmetry (iter_0023) was one of the most biologically interesting findings. Testing whether this asymmetry survives in the Dorothea early-layer decay could confirm a mechanistic principle: scGPT early-layer geometry specifically encodes *activating* regulatory co-proximity, not repression. If true, this is a publishable mechanistic claim about what the model learns. If false, the decay is direction-agnostic and harder to interpret.
- Why now: The Dorothea data is already loaded with direction annotations; the test is low cost relative to the potential finding.

### Slot 3: Cheap Broad-Screen Candidate
**B3 + A1 — GO CC OOV-corrected replication + secreted/membrane community annotation**
Two cheap tests packaged together: (1) re-run GO CC Jaccard vs distance on 195 in-vocab genes (confirm or retire third biological anchor), (2) annotate the 107/88 community members with UniProt subcellular location and run Fisher test. Total cost: ~1 hour of compute; covers two critical open threads.
- Why now: GO CC is cited as the third independent biological anchor in the paper. If it fails on the OOV-corrected set, the multi-predictor joint model result also needs revision. This must be resolved before the paper can be finalized.
