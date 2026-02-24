# Brainstormer Hypothesis Roadmap — iter_0018 → iter_0019+

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| TRRUST signed direction on SV3/SV4 | 0/12 significant (iter_0017 H03). No rescue path visible. | `retire_now` |
| SV2/SV3 GO enrichment sweeps | Completed. Biologically characterized. Paper-ready. | `retire_now` |
| Axis independence (SV2/SV3/SV4) | Completed. Paper-ready. | `retire_now` |
| SV2 distance as smooth PPI ranking metric | H02 shows 1.2x only; confirmed partition framing is correct. Retire distance framing. | `retire_now` |
| Basic kNN clustering coefficient (iter_0002) | Confirmed positive. Already in paper. | `retire_now` |
| H1 persistence raw count tests | Done in iter_0003, paper-ready. No new signal expected. | `retire_now` |
| Cross-model feature alignment via centered_cosine summaries | Iter_0003 p=0.41 — wrong level of representation. Rescue only via raw Geneformer residuals. | `rescue_once_with_major_change` (use raw embeddings, not summaries) |

---

## New Hypothesis Portfolio

### Family A: Attention Geometry (depth follow-up on H03)

**A1 — TRRUST Repression vs Activation Attention Specificity**
- Hypothesis: TRRUST repression TF-target pairs have weaker co-attention than activation pairs, establishing regulatory-direction specificity in attention heads.
- Test: Run same MW-U co-attention test separately for TRRUST repression pairs (separate tsv column). Compare mean_att(repression) vs mean_att(activation) vs background.
- Expected signal: activation > repression > background (ordered enrichment).
- Null/control: gene-label shuffle within attention matrix.
- Value: high | Cost: low (existing data, one script change)

**A2 — Attention Encoding: Layer Specificity**
- Hypothesis: The TRRUST activation attention signal is concentrated in specific Transformer layers, not uniformly distributed.
- Test: If per-layer attention tensors exist in `invariant_causal_edges/lung/`, load each layer separately and compute MW-U for TRRUST activation pairs per layer. Plot layer profile.
- Expected signal: peak co-attention in middle/late layers (layers 6-10).
- Null/control: same shuffle null per layer.
- Value: high | Cost: medium (depends on data availability)

**A3 — Attention Precision@k for Regulatory Pairs**
- Hypothesis: Attention co-occurrence scores provide a better ranking of TF-target pairs than SVD distance.
- Test: Use mean symmetric attention as similarity score for all pairs in 209-gene universe. Compute precision@k at k=50,100,200 for TRRUST activation pairs. Compare to SV2 co-pole score for same pairs.
- Expected signal: attention P@k > SV2 P@k for TRRUST; SVD P@k > attention P@k for STRING PPI.
- Null/control: random pair ranking.
- Value: high | Cost: low

**A4 — GO Biological Process Enrichment in Attention Clusters**
- Hypothesis: Genes with high mutual co-attention (top 1% attention pairs) are enriched for same GO biological process terms.
- Test: Extract top 0.5% of symmetric attention pairs within 209-gene universe. Run Fisher's exact test for shared GO BP terms vs background. Compare to SV2 co-pole pair enrichment for same GO terms.
- Expected signal: GO BP enrichment for both, with complementary term sets (attention=regulatory processes, SVD=physical binding).
- Null/control: random gene pairs.
- Value: high | Cost: medium

---

### Family B: SVD Co-Pole Partition Refinement

**B1 — Multi-Axis Co-Pole Composite Score for PPI Ranking**
- Hypothesis: A composite co-pole membership score across 12 layers × 3 axes (SV2+SV3+SV4) provides >3x enrichment for STRING pairs in precision@k, far exceeding the 1.2x from SV2 distance.
- Test: For each gene pair (i,j), compute fraction of (layer, axis) combinations where both genes are in same top-K or bottom-K pole (K=52). This yields a [0,1] co-pole rate. Use as ranking score for all 21,736 pairs. Compute precision@k and enrichment vs 14.23% baseline.
- Expected signal: P@100 enrichment ≥ 3x (0.43+ vs 0.142 baseline).
- Null/control: co-pole rate on random Gaussian embeddings.
- Value: high | Cost: low

**B2 — Out-of-Sample Hold-Out Co-Pole Benchmark**
- Hypothesis: The co-pole composite signal generalizes to held-out STRING edges not used during any axis selection.
- Test: Reserve 20% of STRING pairs stratified by confidence score. Train/select nothing on test set. Compute co-pole composite score. Measure precision@k on held-out edges only.
- Expected signal: Similar enrichment as full set (robustness).
- Null/control: random 20% of pairs from background.
- Value: high | Cost: low

**B3 — SV2 Pole Stability Across Biological Contexts (Cell Type Comparison)**
- Hypothesis: Genes assigned to the same SV2 pole in lung scGPT are also co-assigned in immune scGPT (cross-context stability of pole partition).
- Test: Load scGPT immune embeddings (iter_0006 data). Compute SV2 poles for immune layer 10. Compare pole membership overlap with lung layer 10 SV2 poles for matched genes. Measure Jaccard index vs null.
- Expected signal: Jaccard > null, especially for layer 10 where signal peaks.
- Null/control: random pole assignment.
- Value: medium | Cost: medium

---

### Family C: Cross-Model Validation (critical gap)

**C1 — Geneformer Raw Embedding SV2 vs STRING (primary cross-model test)**
- Hypothesis: Geneformer residual-stream embeddings show an analogous SV2 co-pole enrichment for STRING PPI pairs, demonstrating the finding is not scGPT-specific.
- Test: Load Geneformer mean gene embeddings from `cycle12_geneformer_lung_bootstrap/`. If raw residuals not available, use `geneformer_edge_dataset.tsv` per-edge similarities. Run identical SV2 co-pole test (K=52, 300 shuffles, 209 named genes). Report n_sig_layers and mean_z.
- Expected signal: ≥6/N significant layers, mean_z ≥ 2.
- Null/control: random Gaussian embeddings of same shape.
- Value: high | Cost: medium

**C2 — CKA Similarity Between scGPT and Geneformer Residual Manifolds**
- Hypothesis: Centered Kernel Alignment (CKA) between scGPT and Geneformer embedding matrices (matched genes × embedding) is significantly higher than permutation null.
- Test: Load 209 matched gene embeddings from scGPT (layer 10) and Geneformer (final layer). Compute linear CKA. Permutation null: shuffle gene labels on one matrix (1000 permutations).
- Expected signal: CKA > 0.3 and p < 0.01.
- Null/control: gene-label shuffle permutation.
- Value: high | Cost: medium

**C3 — STRING Confidence Gradient in Geneformer: Spearman rho**
- Hypothesis: Geneformer embedding cosine similarity shows a monotonic STRING confidence gradient (Spearman rho ≥ 0.7), matching the scGPT SV2 signal.
- Test: Use Geneformer edge similarities (already computed in cycle12). Group by STRING confidence quintile. Compute Spearman rho between STRING score and Geneformer similarity. Compare to scGPT SV2 rho ≥ 0.90.
- Expected signal: Spearman rho ≥ 0.5 for Geneformer.
- Null/control: permuted STRING scores on fixed gene pairs.
- Value: high | Cost: low (existing Geneformer edge data)

---

### Family D: Manifold Geometry

**D1 — Intrinsic Dimension Profile Across Layers**
- Hypothesis: The intrinsic dimension (ID) of the scGPT residual manifold decreases from early to late layers, reflecting progressive dimensionality collapse onto a task-relevant subspace.
- Test: Compute TwoNN or MLE intrinsic dimension estimator on 209-gene embeddings at each of 12 layers. Plot ID vs layer index.
- Expected signal: Monotonic decrease in ID from L0 to L11, with steeper drop at layers 8-10 (where SV2 signal peaks).
- Null/control: ID of random Gaussian embeddings of same shape.
- Value: medium | Cost: low

**D2 — Geodesic vs Euclidean Distance Ratio for PPI vs Non-PPI Pairs**
- Hypothesis: PPI gene pairs have lower geodesic/Euclidean distance ratio than non-PPI pairs, indicating they lie on the same curved manifold sheet.
- Test: Build k=15 kNN graph on layer-10 embeddings. Compute shortest-path geodesic distances for STRING pairs vs background. Compute ratio = geodesic/Euclidean. Mann-Whitney test.
- Expected signal: Lower ratio for STRING pairs (MW_p < 0.01).
- Null/control: background non-STRING pairs.
- Value: medium | Cost: medium

---

### Family E: Biological Depth

**E1 — Cell Ontology Cell-Type Marker Co-Pole Enrichment**
- Hypothesis: Known lung cell-type marker genes (CellMarker/CellTypist) show stronger co-pole clustering than random, and different cell types occupy distinct SVD poles.
- Test: Load lung cell-type marker list for 5-10 major lung cell types. For each type, compute fraction of marker pairs that are co-pole at SV2 layer 10. Compare to non-marker background.
- Expected signal: Each cell type's markers cluster in same pole; enrichment > 2x.
- Null/control: random gene sets of same size.
- Value: medium | Cost: low

**E2 — KEGG Pathway Co-Pole Enrichment (pathway-level biological structure)**
- Hypothesis: Genes within the same KEGG pathway show higher co-pole rate than genes from different pathways.
- Test: Use KEGG gene sets for 10 well-defined pathways with ≥5 genes in 209-gene universe. Compute intra-pathway vs inter-pathway co-pole rate at SV2 layer 10. Mann-Whitney test.
- Expected signal: MW_p < 0.01, intra-pathway co-pole rate > 2x inter-pathway.
- Null/control: randomized pathway assignments of same sizes.
- Value: medium | Cost: low

---

## Top 3 for Immediate Execution

### 1. High-Probability Discovery Candidate
**C3 — STRING Confidence Gradient in Geneformer** + **C1 — Geneformer SV2 co-pole**
Combine into one script. Use existing `geneformer_edge_dataset.tsv` similarities (already computed). This tests cross-model replication with minimal new computation. If positive, it closes the biggest gap in the paper. If negative (Geneformer doesn't replicate), that itself is a major mechanistic finding (scGPT-specific encoding).

### 2. High-Risk / High-Reward Candidate
**A2 — Attention Layer Specificity**
If per-layer attention tensors exist, this converts the aggregate TRRUST finding into a mechanistic claim: specific layers are responsible for regulatory encoding. This would be a strong interpretability result. Risk: per-layer attention tensors may not be saved; cost scales with data availability check.

### 3. Cheap Broad-Screen Candidate
**B1 — Multi-Axis Co-Pole Composite Score** (run with A1 — TRRUST Repression Attention)
Single script: compute composite co-pole score for all 21,736 pairs → precision@k. Simultaneously run TRRUST repression vs activation comparison on existing attention matrix. Total cost: one 30-minute compute session. Both have high expected signal.

---

## Execution Priority for iter_0019

| Priority | Hypothesis | Family | Expected Time | Paper Impact |
|----------|-----------|--------|---------------|--------------|
| 1 | C3+C1: Geneformer cross-model | cross_model | medium | critical |
| 2 | B1: Co-pole composite P@k | SVD partition | low | high |
| 3 | A1: TRRUST repression attention | attention | low | high |
| 4 | A3: Attention P@k for regulatory | attention | low | high |
| 5 | D1: Intrinsic dimension profile | manifold | low | medium |
