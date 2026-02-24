# Brainstormer Hypothesis Roadmap — iter_0062

**Date**: 2026-02-23

---

## Retire / Deprioritize

| Direction | Reason | Verdict |
|-----------|--------|---------|
| FFL geometric ordering | Permutation null at ~0.50 by symmetry; no rescue path | retire_now |
| L8 CKA boundary hypothesis | Consecutive CKA uniform; no signal | retire_now |
| Dual-role gene ambiguity | All genes target-like; null result | retire_now |
| Topology stability via principal angles | Negative in iter_0056/H03 and iter_0057/H02; no major change possible | retire_now |
| Cross-cycle Procrustes with cosine AUROC | Design-flawed; rescued only by full redesign with LR classifier | rescue_once_with_major_change |
| TF hub centrality (geometric centrality of activation TFs) | iter_0024 H02: 0/12 layers; no rescue path | retire_now |

---

## New Hypothesis Portfolio

### H-A: SV2-4 Repulsion Mechanistic Probe
**Hypothesis**: TF-target pairs have significantly lower cosine similarity in SV2-4 than random pairs (anti-clustering) because SV2-4 encodes subcellular compartment/PPI axes, and TF-target pairs are preferentially drawn from different compartments.
**Test**: At L0 (peak AUROC for SV5-7): compute cosine similarity distribution for TF-target pairs vs random pairs in SV2-4. Compare GO CC Jaccard for TF-target pairs vs random pairs. Regress SV2-4 AUROC deficit on GO CC dissimilarity.
**Expected signal**: TF-target pairs have significantly lower GO CC Jaccard than random pairs, explaining why they anti-cluster in the compartment axis.
**Null**: TF-target GO CC Jaccard = random GO CC Jaccard.
**Value**: high | **Cost**: low (reuse existing embeddings + existing GO CC data from iter_0024)

### H-B: Edge AUROC Bonferroni Correction and Subspace Scan (SV1-10 systematic)
**Hypothesis**: Among SV subspaces SV1–SV10, SV5-7 is the optimal window for TF-target co-embedding, and the signal survives Bonferroni correction across L0–L8.
**Test**: Run edge AUROC for each overlapping triplet SV(k,k+2) for k=1..8 at L0. Apply Bonferroni for 8 × 9 = 72 tests. Report which subspaces survive.
**Expected signal**: SV5-7 and possibly SV6-8 survive Bonferroni; SV1-3 and SV2-4 do not.
**Null**: No subspace AUROC > 0.55 after Bonferroni.
**Value**: high | **Cost**: low

### H-C: LR Classifier Cross-Cycle Procrustes (H03 redesign)
**Hypothesis**: A logistic regression classifier trained on cycle1 TF/target labels in SV5-7 space transfers to Procrustes-aligned cycle4 embeddings above chance.
**Test**: Cycle1 L3: SVD→SV5-7, label genes as TF=1/target=0. Train LR. Procrustes-rotate cycle4 SV5-7 → cycle1 frame. Apply LR to aligned cycle4 embeddings. Evaluate AUROC on cycle4 genes.
**Expected signal**: AUROC > 0.60 on Procrustes-aligned cycle4 (above cycle1-only baseline of 0.54–0.56).
**Null**: AUROC ≤ random (0.50).
**Value**: high | **Cost**: medium (requires LR + careful label alignment)

### H-D: kNN Assortativity by TF/Target Label vs Layer
**Hypothesis**: kNN graph assortativity by TF/target label is highest at L0 in SV5-7 space and declines by L9, mirroring the AUROC trajectory.
**Test**: At each layer: build kNN (k=10) on SV5-7 embeddings of all 2039 nonzero genes. Assign TF/target label. Compute networkx assortativity coefficient. Plot vs layer.
**Expected signal**: Assortativity > 0 at L0, declining to ~0 at L9.
**Null**: Assortativity ≈ 0 at all layers (label-shuffle null).
**Value**: medium | **Cost**: low

### H-E: Persistent Homology of TF+Target Gene Neighborhood at L0 vs L9
**Hypothesis**: The loss of edge-level AUROC between L0 and L9 corresponds to a topological change: Betti-1 loops among TF+target gene neighborhoods collapse or become indistinguishable from random at L9.
**Test**: Extract SV5-7 embeddings for TF+target genes (all 589 valid pairs) at L0 and L9. Compute Vietoris-Rips Betti-0/1 via ripser. Compare to gene-label shuffle null (N=100).
**Expected signal**: Betti-1 > null at L0, Betti-1 ≈ null at L9.
**Null**: No difference between L0 and L9 in Betti-1 relative to shuffle.
**Value**: high | **Cost**: medium (ripser computation on ~200 genes)

### H-F: Effective Rank vs AUROC Cross-Correlation
**Hypothesis**: Layer-wise AUROC decline correlates with effective rank increase (embedding geometry becomes higher-dimensional/less structured in later layers).
**Test**: Compute effective rank at each layer (from SVD eigenvalue spectrum: exp(entropy of normalized squared singular values)). Spearman(AUROC_L, eff_rank_L).
**Expected signal**: Negative correlation: AUROC declines as eff_rank increases.
**Null**: No significant Spearman correlation.
**Value**: medium | **Cost**: low (SVD already computed per layer for M1)

### H-G: Signed Regulation Geometry in SV5-7 (Activation vs Repression Edge Direction)
**Hypothesis**: Activation TF-target pairs have higher cosine similarity in SV5-7 than repression TF-target pairs at early layers.
**Test**: Split TRRUST cycle4 edges into activation (N=?) and repression (N=?). Compute SV5-7 cosine similarity AUROC for each subtype at each layer. MWU activation vs repression similarity.
**Expected signal**: Activation > repression in SV5-7 at L0–L4; both decline toward L9.
**Null**: No activation/repression difference in SV5-7 (sign-blind geometry).
**Value**: high | **Cost**: low

### H-H: L2 Distance vs Cosine Similarity as AUROC Metric (Metric Comparison)
**Hypothesis**: Negative L2 distance in SV5-7 is a stronger or weaker TF-target discriminator than cosine similarity, revealing whether magnitude or angle encodes the signal.
**Test**: At each layer: score each TF-target pair by both (a) cosine similarity and (b) negative L2 distance in SV5-7. Compute AUROC for each. Paired comparison.
**Expected signal**: Cosine similarity AUROC ≥ L2-distance AUROC at L0–L3; may diverge at later layers.
**Null**: Both metrics give AUROC ≈ 0.50.
**Value**: medium | **Cost**: low

### H-I: Cycle Comparison — cycle1 vs cycle4 Edge AUROC Profile Without Procrustes
**Hypothesis**: cycle1 TF-target pairs show a different layer-AUROC profile than cycle4 pairs, reflecting different training conditions, without requiring any alignment.
**Test**: Run identical M1 protocol (SV5-7 AUROC, permutation null N=200) on cycle1 embeddings with cycle1 edge dataset (288 pos pairs). Compare layer profile to cycle4 (735 pos pairs). Spearman correlation between cycle1 and cycle4 AUROC vectors.
**Expected signal**: Both show monotonic decline but with different peak AUROC and layer of signal collapse.
**Null**: Identical profiles (same geometry across training cycles).
**Value**: medium | **Cost**: low

### H-J: Geneformer Cross-Model Edge AUROC Replication
**Hypothesis**: TRRUST TF-target pairs show higher SV5-7 cosine similarity than random in Geneformer residual embeddings, replicating the scGPT finding.
**Test**: Load Geneformer gene embeddings (if available from prior artifacts), run M1 protocol on TRRUST cycle4 pairs. Compare AUROC profile to scGPT.
**Expected signal**: AUROC > 0.55 at early Geneformer layers, with similar decline.
**Null**: AUROC ≈ 0.50 (no regulatory geometry in Geneformer SV5-7).
**Value**: high | **Cost**: medium (requires verifying Geneformer artifact availability)

### H-K: Norm Confound Control per Layer
**Hypothesis**: Embedding L2-norm varies systematically by layer and confounds cosine similarity AUROC (higher-norm layers may show spurious similarity structure).
**Test**: Compute mean embedding L2-norm per layer (2039 genes). Spearman(AUROC, mean_norm). Also: re-run AUROC after L2-normalizing all embeddings before SVD.
**Expected signal**: Norm is constant or its correlation with AUROC is low, confirming no norm confound.
**Null**: Strong Spearman(AUROC, norm), meaning norm drives AUROC.
**Value**: medium | **Cost**: low

### H-L: Multi-Database Joint Edge Model (STRING + TRRUST + Dorothea at Edge Level)
**Hypothesis**: STRING PPI pairs show lower SV5-7 cosine AUROC than TRRUST TF-target pairs; Dorothea high-confidence regulatory pairs show intermediate AUROC. The three databases occupy distinct AUROC bands in SV5-7 space.
**Test**: Run M1 AUROC protocol on three edge sets in the same SV5-7 embedding space: (a) TRRUST cycle4 (existing, AUROC~0.60), (b) STRING >= 0.7 pairs from the 2039-gene vocab, (c) Dorothea A/B pairs from the 2039-gene vocab. Compare layer profiles.
**Expected signal**: TRRUST > Dorothea > STRING in SV5-7 AUROC; STRING peaks in SV2-4.
**Null**: All three AUROC ≈ 0.50 in SV5-7.
**Value**: high | **Cost**: medium

---

## Top 3 for Immediate Execution

### Priority 1 — High-probability discovery candidate
**H-A: SV2-4 Repulsion Mechanistic Probe**
- Rationale: Sub-null AUROC in SV2-4 is the most anomalous signal in iter_0062. If TF-target pairs have lower GO CC Jaccard than random (different compartments), this would explain the repulsion and provide a clean mechanistic account connecting the two major geometric axes (PPI compartment vs regulatory co-expression). Very cheap to test — all data in hand.
- Key artifact needed: GO CC Jaccard matrix (iter_0024 artifacts), cycle4 SV2-4 cosine similarities (compute from existing embeddings)
- Expected compute: < 5 min

### Priority 2 — High-risk/high-reward candidate
**H-G: Signed Regulation Geometry in SV5-7**
- Rationale: Prior results (iter_0012, iter_0017) showed activation/repression distinction in SV2, not SV3/SV4. If SV5-7 also encodes sign (activation>repression), this would be a new axis-specific finding with mechanistic implications. The risk is that the TRRUST cycle4 dataset may not have enough repression pairs for power.
- Key artifact needed: TRRUST cycle4 edge dataset with regulation_type column (check if present in cycle1_edge_dataset.tsv)
- Expected compute: < 10 min

### Priority 3 — Cheap broad-screen candidate
**H-B: Subspace Scan SV1–SV10 + Bonferroni**
- Rationale: We have tested SV5-7 and SV2-4 only. A systematic scan of all overlapping triplets SV(k,k+2) for k=1..8 at L0 will definitively map which subspace windows encode TF-target co-embedding. Cheap (single layer, 8 subspace windows), immediately interpretable, and will either confirm SV5-7 as optimal or reveal a better window.
- Key artifact needed: cycle4 embeddings (already loaded), TRRUST edge set
- Expected compute: < 10 min
