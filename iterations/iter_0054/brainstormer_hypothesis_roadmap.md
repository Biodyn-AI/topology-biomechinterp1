# Brainstormer Hypothesis Roadmap — iter_0054 → iter_0055+

**Date**: 2026-02-23

---

## Retire / Deprioritize

| Hypothesis Family | Iterations Tried | Reason | Action |
|-------------------|-----------------|--------|--------|
| TwoNN intrinsic dimension | iter_0047+ | No signal, fundamental mismatch with data sparsity | `retire_now` |
| Betti loops / PH on small gene sets | iter_0050-0052 | Consistently negative, sparse gene count kills PH | `retire_now` |
| kNN cell-type purity | iter_0045 | Wrong framing — embedding not trained for cell-type separation | `retire_now` |
| Continuous degree → SV1 correlation | 3× confirmed negative | Binary circuit membership is the right split; graded degree adds nothing | `retire_now` |
| SV8-14 secondary signal at L8 | iter_0050+ | rbc 0.065-0.078, below threshold for standalone claim | `deprioritize` |
| SV1 alone as standalone test | confirmed 3× | Circuit identity proxy confirmed; no new extractable signal without combination | `deprioritize` |

---

## New Hypothesis Portfolio

### H-NEW-01: Cross-Seed Replication of H01 + H03 (High value, Low cost)
**Hypothesis**: The TF/target AUROC (≥0.65) and directionality trajectory (peak at L9) replicate in seeds 43 and 44 of cycle4_immune, establishing cross-seed robustness.
**Test**: Load seed43/seed44 embeddings. Replicate H01 (LR AUROC, SV5-7 L0) and H03 (all-layer displacement trajectory, SV5-7) exactly. Report AUROC and trajectory magnitude per seed. Permutation null for AUROC.
**Expected signal**: Both seeds AUROC > 0.65, p < 0.01; both seeds show peak magnitude at L9±1.
**Null/control**: Same protocol on shuffled gene labels (seed-specific permutation).
**Value**: high | **Cost**: low

### H-NEW-02: SV2-4 Directionality Trajectory Control (High value, Low cost)
**Hypothesis**: SV2-4 at all 12 layers shows a complementary directionality trajectory — weak at L0, growing and peaking L8-L11 — opposite to SV5-7.
**Test**: Replicate H03 protocol using SV2-4 projection at all 12 layers. Compare magnitude trajectories: SV5-7 (from H03) vs SV2-4. Also test SV2-4 at L9 (the SV5-7 peak) for comparison.
**Expected signal**: SV2-4 magnitude monotonically smaller than SV5-7 at L0-L4, then comparable or larger at L8-L11. This creates a complete two-subspace story.
**Null/control**: SV5-7 trajectory (already in hand) as explicit comparison.
**Value**: high | **Cost**: low

### H-NEW-03: Classification at Peak Layer (L9) vs L0 (High value, Low cost)
**Hypothesis**: The SV5-7 directionality peaks at L9 (H03); therefore TF vs target-only classification should be most accurate at L9, not L0.
**Test**: Run H01 protocol (LR on SV5-7, 5-fold CV AUROC) at every layer L0-L11. Track AUROC trajectory. Test if L9 AUROC > L0 AUROC (McNemar test or paired t-test on fold scores).
**Expected signal**: AUROC peaks at L9 (≥ 0.72), declining at L11. L0 AUROC (0.694) is a minimum. Creates classification × directionality unified picture.
**Null/control**: Permuted labels at each layer; SV2-4 AUROC trajectory as control.
**Value**: high | **Cost**: low

### H-NEW-04: Axis Rotation Characterization Across Layers (High value, Medium cost)
**Hypothesis**: The rotation of dominant directionality axis across layers (observed in H03) follows a structured pattern — axes between consecutive layers form predictable angular trajectories — not random rotation.
**Test**: For each pair of consecutive layers (L0-L1, L1-L2, ..., L10-L11): compute cosine similarity of the dominant directionality axis (from per-layer SVD of displacement vectors). Test if axis trajectory is smooth vs discontinuous. Compare to null: random unit vectors in 3D.
**Expected signal**: High consecutive-layer cosine similarity (>0.7) for most transitions, with one or two "break" transitions at functionally important layers (e.g., L4→L5 or L8→L9).
**Null/control**: Random unit vectors in 3D; shuffled layer order.
**Value**: high | **Cost**: medium

### H-NEW-05: Partial Spearman for BFS Depth (Controlling for TF/Target Label) (Medium value, Low cost)
**Hypothesis**: The BFS depth correlation with SV5/SV7 (rho=0.167, H02) persists after controlling for TF vs target-only binary status, indicating cascade position is independently encoded.
**Test**: Partial Spearman: depth ~ SV5 | TF_binary and depth ~ SV7 | TF_binary. Compare partial rho to marginal rho. If |partial rho| > 0.10 and p < 0.05, the effect is not merely TF/target identity confound.
**Expected signal**: Partial rho 0.10-0.15 (reduced but surviving). If partial rho ≈ 0, depth signal is just TF/target identity and H02 is fully accounted for by H01.
**Null/control**: Permuted depth assignments; partial correlation against random covariate.
**Value**: medium | **Cost**: low

### H-NEW-06: GO Enrichment of Displacement Axis (Medium value, Low cost)
**Hypothesis**: Genes with high positive projection onto the dominant TF→target displacement axis in SV5-7 at L0 are enriched for transcription factor GO terms (DNA binding, transcriptional regulation).
**Test**: Project all 2039 nonzero genes onto the dominant displacement axis (first PC of displacement vectors at L0). Fisher exact test comparing top/bottom quartile for GO:0003700 (TF activity), GO:0006355 (transcription regulation), GO:0005634 (nucleus). Benjamini-Hochberg correction.
**Expected signal**: Top quartile enriched for TF GO terms (OR > 2, p_adj < 0.05). This biologically interprets what the geometry encodes.
**Null/control**: Random projection axis; GO test on random gene splits of equal size.
**Value**: medium | **Cost**: low

### H-NEW-07: AUROC Using Full SV1-7 Feature Set (Medium value, Low cost)
**Hypothesis**: Combining SV1 loading with SV5-7 coordinates (7D feature vector: SV1-7) improves TF vs target-only AUROC above 0.70.
**Test**: Features = SV1-7 at L0 (7D). LR with 5-fold CV AUROC. Compare to SV5-7 alone (0.694). Also test SV1-7 at L9 (the peak directionality layer). Ablation: which SV axes contribute most (permutation importance)?
**Expected signal**: AUROC 0.70-0.75. If not, SV5-7 is already near-optimal and additional axes are redundant.
**Null/control**: Permuted labels; SV5-7 baseline (0.694).
**Value**: medium | **Cost**: low

### H-NEW-08: Activator vs Repressor Edge Directionality Asymmetry (Medium value, Medium cost)
**Hypothesis**: TRRUST annotates TF→target edges as activating or repressing; activating and repressing edges show different directionality signatures in SV5-7 space (different magnitude or direction of displacement vector).
**Test**: Subset the 589 positive pairs by TRRUST edge type (activating vs repressing). For each subset, compute displacement vector magnitude in SV5-7 at L0 and L9. Mann-Whitney test: do activating edges show higher directionality than repressing?
**Expected signal**: Activating edges show mean displacement magnitude ≥1.2× repressing edges at L9, p < 0.05.
**Null/control**: Random partition of edges into equal-sized groups.
**Value**: medium | **Cost**: medium

### H-NEW-09: Edge Prediction from L9 SV5-7 Geometry (Medium value, High cost)
**Hypothesis**: At the peak directionality layer (L9), SV5-7 distance plus directionality axis projection jointly predict TRRUST-positive edges with AUROC > 0.65.
**Test**: For all gene pairs in cycle1_edge_dataset: features = [SV5-7 Euclidean distance at L9, projection of gene_pair vector onto dominant displacement axis at L9]. Logistic regression with 20% held-out split. Compare to L0 version (same design).
**Expected signal**: L9 AUROC > L0 AUROC (hypothesized L0 ≈ 0.58-0.62 for distance alone). L9 advantage of ≥0.03 AUROC.
**Null/control**: Shuffled edges; L0 version as comparison.
**Value**: medium | **Cost**: high

### H-NEW-10: 0-dim Persistent Homology on SV5-7 by Gene Class (Medium value, Medium cost)
**Hypothesis**: TF genes form topologically more cohesive clusters in SV5-7 L0 space than target-only genes — detectable as higher persistence in 0-dim PH (more persistent connected components survive longer).
**Test**: Compute 0-dim PH (Rips filtration) on SV5-7 L0 coordinates for: (a) TF genes (n=68), (b) target-only genes (n=215), (c) random background (n=68, matched size). Compare persistence entropy and lifetime of top-5 components.
**Expected signal**: TF genes show lower persistence entropy (more cohesive clustering) than target-only or background.
**Null/control**: Label-permuted subsets; random matched-size samples.
**Value**: medium | **Cost**: medium

### H-NEW-11: Cross-seed Ensemble AUROC (Medium value, Low cost)
**Hypothesis**: Combining SV5-7 L0 features from 3 seeds (main, seed43, seed44) as a 9D feature vector improves TF classification AUROC above any single seed.
**Test**: Concatenate SV5-7 L0 coordinates from 3 seeds per gene (9D feature). LR 5-fold CV AUROC on overlapping genes. Compare to single-seed 3D AUROC (0.694).
**Expected signal**: Ensemble AUROC 0.72-0.78. This tests whether independent model seeds capture complementary geometric information.
**Null/control**: Single-seed baseline; feature-shuffled ensemble.
**Value**: medium | **Cost**: low

### H-NEW-12: Geodesic vs Euclidean Distance Ratio for TF-TF vs TF-Target Pairs (Low value, Medium cost)
**Hypothesis**: TF-TF gene pairs have geodesic/Euclidean distance ratio closer to 1.0 (approximately linear manifold) in SV5-7 at L0, while TF-target pairs show higher ratio (curved manifold region), reflecting different local manifold geometry by edge type.
**Test**: Build k=15 nearest neighbor graph on SV5-7 L0. Compute shortest path (geodesic) vs Euclidean distance for all TF-TF pairs vs TF-target pairs in TRRUST. Mann-Whitney on ratio distributions.
**Expected signal**: TF-target pairs show mean ratio > TF-TF pairs (ratio 1.05-1.15 vs 1.00-1.05).
**Null/control**: Random matched-size gene pairs; permuted edge labels.
**Value**: low | **Cost**: medium

### H-NEW-13: Non-circuit Gene Directionality as Negative Control (Medium value, Low cost)
**Hypothesis**: Genes not in TRRUST (non-circuit genes) show no directional asymmetry in SV5-7 — random displacement vectors that do not systematically displace along any axis — validating that the H03 signal is circuit-specific.
**Test**: Select 589 randomly paired non-circuit gene pairs (matching circuit pair count). Compute displacement magnitude in SV5-7 at all 12 layers. Compare to H03 circuit trajectory.
**Expected signal**: Non-circuit pairs show magnitude ≤ 0.006 at all layers (matching L0 baseline), no amplification trend. Cohen's d ≈ 0 at all layers.
**Null/control**: Circuit pairs (H03 results) as positive comparison.
**Value**: medium | **Cost**: low

### H-NEW-14: Layer-Resolved Classification AUROC × Directionality Joint Plot (Medium value, Low cost)
**Hypothesis**: Classification AUROC (H01 at all layers) and directionality magnitude (H03) are correlated across layers — layers with higher directionality signal support better TF/target classification.
**Test**: Compute H01 AUROC at all 12 layers (same LR protocol). Correlate AUROC vector vs directionality magnitude vector across layers (Spearman). If rho > 0.6, the two measures reflect the same underlying geometry.
**Expected signal**: rho > 0.7 between AUROC-by-layer and magnitude-by-layer. Peak classification at L9 matching peak magnitude. This creates a unified geometric narrative.
**Null/control**: Random permutation of layer order for one vector; correlation at SV2-4 as control.
**Value**: medium | **Cost**: low

---

## Top 3 for Immediate Execution

### Candidate 1 — High-probability discovery candidate
**H-NEW-01: Cross-Seed Replication of H01 + H03**
- Rationale: The AUROC=0.694 and L0→L9 trajectory are the two strongest iter_0054 findings but remain single-seed results. Cross-seed validation is the minimum bar for paper credibility. Low cost (exact protocol replication, just swap embeddings). If both seeds replicate, these become the paper's primary claims with statistical robustness.

### Candidate 2 — High-risk/high-reward candidate
**H-NEW-03: Classification AUROC at All Layers (L0→L11)**
- Rationale: H03 showed directionality peaks at L9. H01 tested only L0 and L8. Testing H01 at all 12 layers costs nearly nothing (same code, loop over layers) but could reveal that L9 classification AUROC is substantially higher than L0 (potentially 0.72+). This would unify H01 and H03 into a single coherent story: "the geometry builds regulatory discriminability from L0 to L9, peaking there." If L9 AUROC is not higher than L0, it means classification and directionality are mechanistically different — also an interesting finding.

### Candidate 3 — Cheap broad-screen candidate
**H-NEW-02 + H-NEW-06 combo: SV2-4 Control Trajectory + GO Enrichment of Displacement Axis**
- Rationale: Both are extremely cheap (SV2-4 trajectory reuses all H03 code, GO enrichment is a Fisher test on already-computed projections). Running them together in one iteration fills the two most conspicuous gaps in the paper: the complementary-subspace story (SV2-4) and the biological interpretation of what the geometry encodes (GO). Together they cost roughly one hour of compute and produce two figures.
