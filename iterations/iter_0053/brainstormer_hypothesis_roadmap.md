# Brainstormer Hypothesis Roadmap — iter_0053 → iter_0054+

**Date**: 2026-02-23

---

## Retire / Deprioritize

| Hypothesis Family | Reason | Action |
|-------------------|--------|--------|
| Betti loop regulatory circuits | Retired iter_0052, consistently negative | `retire_now` |
| TwoNN intrinsic dimensionality | Retired iter_0047, no signal found | `retire_now` |
| kNN purity (cell-type markers) | Retired iter_0045, wrong framing for this data | `retire_now` |
| SV8-14 secondary signal at L8 | Very weak rbc (0.065–0.078), not worth dedicated iteration | `deprioritize` |
| Continuous degree → SV1 correlation | Confirmed negative 3 times; binary membership is the right split | `retire_now` |
| SV1 alone as standalone test | Confirmed 3× as circuit identity proxy; no new info to extract without combining | `deprioritize` |

---

## New Hypothesis Portfolio

### H-A: TF vs Target-only Role Classification from SV5-7 (High value, Low cost)
**Hypothesis**: SV5-7 coordinates at L0 are sufficient to classify a gene as TF (out-degree > 0) vs target-only (in-degree > 0, out-degree = 0) with AUROC > 0.60.
**Test**: Assign TF/target-only labels from TRRUST for circuit genes (n=295). Train logistic regression on SV5-7 [3D] at L0. Evaluate 5-fold CV AUROC. Null: permuted labels, compare AUROC distribution.
**Expected signal**: AUROC 0.62–0.70, p < 0.01 vs permutation. Positive class = TF.
**Null/control**: Permuted label AUROC distribution (1000 permutations).
**Value**: high | **Cost**: low

### H-B: SV2-4 at L8 TF→Target Directional Asymmetry (High value, Low cost)
**Hypothesis**: The directional TF→target asymmetry found in SV5-7 at L0 is reproduced — and stronger — in SV2-4 at L8, consistent with L8 being the dominant regulatory encoding layer.
**Test**: Replicate H02 protocol but project to SV2-4 at L8. Compute signed displacement TF→target per SV axis. One-sample t-test per axis + combined test. Compare Cohen's d to L0 SV5-7 result (d≈0.12).
**Expected signal**: Cohen's d > 0.15, p_combined < 0.001 at L8 in SV2-4.
**Null/control**: Same test at SV2-4 at L0 (expected negative).
**Value**: high | **Cost**: low

### H-C: Multi-axis Joint Classifier (SV1 + SV5-7 → TF Identity) (High value, Medium cost)
**Hypothesis**: Combining SV1 loading and SV5-7 coordinates (4D feature vector) classifies TF identity better than either subspace alone, with AUROC > 0.70.
**Test**: For circuit genes: features = [|SV1_L0|, SV5_L0, SV6_L0, SV7_L0]. Logistic regression with 5-fold CV. Compare to SV5-7 alone (H-A) and SV1 alone.
**Expected signal**: AUROC improvement of ≥0.05 over SV5-7 alone. Feature importance: SV1 separates circuit membership, SV5-7 separates direction role.
**Null/control**: Single-feature baseline (SV1 only), permuted labels.
**Value**: high | **Cost**: medium

### H-D: Cross-seed Replication of SV5-7 L0 Signal (Medium value, Low cost)
**Hypothesis**: The SV5-7 L0 rbc_resid signal (rbc=0.147) replicates in independent seeds (seed43, seed44) of the cycle4_immune model.
**Test**: Load seed43/seed44 embeddings. Replicate full SV5-7 L0 pipeline (SVD, projection, residualized rbc). Report rbc_resid and p for each seed vs main (seed42).
**Expected signal**: Both seeds show rbc_resid > 0.08, p < 0.05.
**Null/control**: Within-seed bootstrap already done (CI=[0.084, 0.199]); seeds below 0.084 = failure.
**Value**: medium | **Cost**: low

### H-E: Regulatory Hierarchy Depth Encoded in SV5-7 (High value, Medium cost)
**Hypothesis**: Genes further downstream in the TRRUST regulatory cascade (higher minimum path length from a master TF) have systematically different SV5-7 positions, encoding hierarchy depth.
**Test**: BFS from master TFs (out-degree ≥ 5, in-degree = 0) in the TRRUST network. Assign hierarchy depth (layer 0 = master TF, layer 1 = direct targets, etc.). Spearman correlation between depth and projection onto the dominant directionality axis (from H02 SVD of displacement matrix).
**Expected signal**: Spearman rho > 0.15, p < 0.01, monotonic depth-axis relationship.
**Null/control**: Random depth permutation; shuffled network.
**Value**: high | **Cost**: medium

### H-F: Edge Prediction from SV5-7 Geometry (High value, High cost)
**Hypothesis**: The SV5-7 distance between two genes, combined with their relative position on the directionality axis, predicts novel TRRUST-positive edges with AUROC > 0.65.
**Test**: For all gene pairs not in training set: compute SV5-7 Euclidean distance + directionality axis projection. Logistic regression with these two features. Evaluate AUROC on held-out positive edges (20% split from cycle1_edge_dataset).
**Expected signal**: AUROC 0.65–0.75. Directionality feature contributes ≥0.03 AUROC over distance alone.
**Null/control**: Random forest with co-expression distance only as baseline.
**Value**: high | **Cost**: high

### H-G: SV5-7 vs SV2-4 Encode Different Edge Subsets (Medium value, Medium cost)
**Hypothesis**: SV5-7 (early) and SV2-4 (deep) encode structurally different subsets of the regulatory network — e.g., SV5-7 captures activator TF→target pairs, SV2-4 captures repressor pairs (or short-range vs long-range regulatory distance).
**Test**: For each positive edge, compute rbc_resid signal strength in SV5-7 at L0 and SV2-4 at L8. Fisher exact test: do edges with top-quartile SV5-7 signal overlap with top-quartile SV2-4 signal? If not, characterize the distinguishing biological property (TF type, edge type, regulatory distance from TRRUST annotation).
**Expected signal**: Low overlap (Jaccard < 0.3) between top-quartile edges in each subspace. Biological correlate: activator vs repressor split, or direct vs indirect regulation.
**Null/control**: Random subsets of equal size; shuffled edge labels.
**Value**: medium | **Cost**: medium

### H-H: SV5-7 Local Curvature Near TF Clusters (Medium value, High cost)
**Hypothesis**: The SV5-7 manifold is locally curved near TF gene clusters (higher local intrinsic dimension or geodesic deviation) compared to non-regulatory gene neighborhoods.
**Test**: In SV5-7 3D space at L0, compute local PCA dimension (number of PCs explaining 90% variance) in k=20 nearest neighbors for each gene. Compare circuit genes vs non-circuit genes. Also compute geodesic vs Euclidean distance ratio for TF-TF pairs.
**Expected signal**: Circuit genes show higher local effective dimension or geodesic/Euclidean ratio > 1.05 vs non-circuit.
**Null/control**: Random matched neighborhood samples; shuffled gene labels.
**Value**: medium | **Cost**: high

### H-I: Topological Persistence in SV5-7 Subspace by Gene Class (Medium value, Medium cost)
**Hypothesis**: TF genes form topologically distinct persistent loops or clusters in SV5-7 3D space compared to target-only genes, detectable via 0-dimensional persistent homology (connected component structure).
**Test**: Compute 0-dim PH (clustering persistence) on SV5-7 L0 coordinates for: (a) TF-only gene subset, (b) target-only subset, (c) random background subset. Compare birth-death diagrams; test if TF persistence entropy differs significantly.
**Expected signal**: TF genes show fewer, more persistent clusters (lower entropy) than targets or background.
**Null/control**: Permuted gene class labels; PH on random subsets of equal size.
**Value**: medium | **Cost**: medium

### H-J: Layer-resolved Directionality Trajectory (Medium value, Low cost)
**Hypothesis**: The mean TF→target displacement vector in SV5-7 space follows a systematic trajectory across layers L0→L11, with magnitude increasing monotonically from L0 to L8 then stabilizing or reversing.
**Test**: At all 12 layers: compute mean displacement vector (3D) for all positive TRRUST edges in SV5-7. Track displacement magnitude and angular change between consecutive layers. Test monotonicity of magnitude increase L0→L8 with page trend test.
**Expected signal**: Monotonic increase in displacement magnitude L0→L8 (confirmed at L0/L1/L2/L8); clear directional convergence.
**Null/control**: Same analysis with shuffled gene-layer assignments.
**Value**: medium | **Cost**: low

### H-K: SV5-7 Coordinate Generalization Across Biological Contexts (Medium value, High cost)
**Hypothesis**: The SV5-7 regulatory proximity signal is not immune-cell-specific; it replicates in a different cell-type dataset (e.g., fibroblast or neural embeddings from subproject_38 reference data).
**Test**: If non-immune embeddings are available in the reference root, extract layer-0 embeddings for overlapping genes. Compute SV5-7 rbc_resid for shared TRRUST edges. Compare to immune rbc=0.147.
**Expected signal**: rbc_resid > 0.08 in non-immune data; signal is context-general.
**Null/control**: Permuted edges in each context separately.
**Value**: medium | **Cost**: high (depends on data availability)

### H-L: SV5-7 Displacement Axis Interpretability via GO Enrichment (Medium value, Low cost)
**Hypothesis**: The dominant directional axis of TF→target displacement in SV5-7 (from H02 displacement SVD) is aligned to biological GO terms — genes with high positive projection are enriched for TF-associated GO terms (DNA binding, transcription regulation).
**Test**: Project all 2039 nonzero genes onto the dominant displacement axis from H02. Compare top-quartile (high positive projection) vs bottom-quartile (low/negative projection) for GO biological process enrichment. Use Fisher exact test with Benjamini-Hochberg correction.
**Expected signal**: Top-quartile enriched for GO:0006355 (regulation of transcription), GO:0003700 (TF activity), or similar.
**Null/control**: Random projection axis of same length; same GO test.
**Value**: medium | **Cost**: low

---

## Top 3 for Immediate Execution

### Candidate 1 — High-probability discovery candidate
**H-A: TF vs Target-only Role Classification from SV5-7**
- Rationale: H02 proves the geometric signal exists. H-A converts it to a predictive AUROC claim. If AUROC ≥ 0.65, this is paper-level result that the embedding geometry functionally distinguishes regulatory roles without any sequence or network features. Low implementation cost (logistic regression on 3D features, n=295 circuit genes). Result available in one iteration.

### Candidate 2 — High-risk/high-reward candidate
**H-E: Regulatory Hierarchy Depth Encoded in SV5-7**
- Rationale: If hierarchy depth correlates with position on the directionality axis, the SV5-7 subspace is a low-dimensional map of the regulatory cascade — a compelling narrative for the paper. High risk: the TRRUST network may be too sparse for robust BFS depth assignment. High reward: a depth-axis relationship would be a qualitatively new finding about how LLMs encode biological hierarchy.

### Candidate 3 — Cheap broad-screen candidate
**H-J: Layer-resolved Directionality Trajectory (L0→L11)**
- Rationale: The displacement magnitude data at L0/L1/L2/L8 is already partially collected from H02. Extending to all 12 layers requires running the same displacement computation with no new methodology. This produces a complete layer-trajectory figure at near-zero incremental cost. If monotonicity holds, it is strong supporting evidence for the directed manifold narrative.
