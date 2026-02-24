# Brainstormer Hypothesis Roadmap — iter_0047

---

## Retire / Deprioritize

| Direction | Reason | Status |
|---|---|---|
| ID compression phase transition / breakpoint detection | R²=0.986 linear fit; piecewise breakpoint marginal (L5, slopes -1.61 vs -1.22); no scientific story to pursue | `retire_now` |
| SV1 at 93.4% variance (iter_0046 claim) | Confirmed artifact of zero-norm contamination; corrected to 18.6% | `retire_now` (as stated claim) |
| 14× SVD eff_rank collapse (iter_0046 claim) | Same artifact; real value 4.9×; retire the inflated claim | `retire_now` (as stated claim) |
| Cross-model alignment via low-dim feature-effect summaries | Three iterations inconclusive; method is not powerful enough without matched residual tensors | `deprioritize` |
| TRRUST co-target clustering on full 4941 gene set | Zero-norm contamination means 58.7% of genes are noise points; all prior TRRUST results on full set are suspect | `rescue_once_with_major_change` (nonzero-only + larger TF set) |
| Any test on global PCA of full 4941 gene set | Zero-vector contamination dominates PCA; must always use nonzero-only n=2039 | `retire_now` (the method variant, not the question) |

---

## New Hypothesis Portfolio

### FAMILY A: Gene Identity and Biology (Critical Blocker)

**H-A1: Locate gene names file**
- Hypothesis: The gene names for the 4941-gene embedding are stored in an alternate path, companion artifact, or the scGPT preprocessing output directory, and can be recovered without retraining.
- Test: Search reference root for any `.npy`, `.csv`, `.json`, `.pkl` file containing gene name arrays of length 4941 or compatible with the embedding shape. Check scGPT vocab files, AnnData obs/var fields, and cycle4_immune_main preprocessing logs.
- Expected signal: File found; gene name assignment allows all downstream biological annotation to proceed.
- Null: If not found, attempt to reconstruct from scGPT model vocabulary (known public file).
- Value: HIGH | Cost: LOW (filesystem search + one scGPT vocab lookup)

**H-A2: Zero-norm gene biological identity**
- Hypothesis: The 2902 zero-norm genes are biologically distinguishable — they are either low-expression in the immune training corpus, not in the scGPT vocabulary, or belong to specific functional categories (e.g., tissue-specific, non-coding).
- Test: Once gene names recovered: (1) GO enrichment Fisher-exact for zero vs nonzero set; (2) STRING degree comparison; (3) median expression in matched immune scRNA reference; (4) TRRUST in-degree. Compare to random 2902-gene split.
- Expected signal: Zero-norm genes are enriched for low-expression or non-immune functional categories at FDR < 0.05.
- Null: Random split of 4941 genes; annotation-density controlled permutation.
- Value: HIGH | Cost: LOW (once gene names available)

**H-A3: SV1 and SV2 biological axis identity**
- Hypothesis: The leading singular vectors of the nonzero-gene embedding at L11 (SV1=18.6%, SV2=11.7% variance) correspond to interpretable biological axes — likely expression level (SV1) and a functional/lineage axis (SV2).
- Test: Project 2039 nonzero genes onto SV1/SV2 at L11. Correlate SV1 score with: mean expression, gene length, number of GO terms, TRRUST out-degree. Run GO enrichment on top/bottom quartile of SV1 and SV2 separately.
- Expected signal: |Spearman ρ| > 0.35 for SV1 vs expression level; SV2 shows functional GO enrichment with no expression correlation.
- Null: Permutation of gene-score labels; annotation-density matched random sets.
- Value: HIGH | Cost: LOW

---

### FAMILY B: Spectral Structure

**H-B1: Full spectral decay curve across 12 layers**
- Hypothesis: The spectral decay of the nonzero-gene embedding matrix changes its shape across layers — early layers have a slow power-law decay (diffuse), late layers have a steeper or step-function decay (concentrated), consistent with the 4.9× eff_rank collapse.
- Test: At each layer, compute singular values of the 2039×512 nonzero-gene matrix. Normalize by total variance. Plot cumulative explained variance vs rank. Fit power-law exponent to the spectral distribution. Track exponent across layers.
- Expected signal: Power-law exponent steepens from L0 to L11 (more concentrated spectrum); inflection in exponent matches largest per-layer ID drops (L3-L4, L9-L10).
- Null: Shuffled embedding at each layer.
- Value: HIGH | Cost: LOW

**H-B2: Spectral bulk vs tail dynamics**
- Hypothesis: The 4.9× eff_rank collapse is driven by the top few singular values growing (dominant directions strengthen) rather than the bulk shrinking, indicating active structuring rather than passive compression.
- Test: At each layer, track: (1) top-5 singular values; (2) bulk (rank 6–50 cumulative variance); (3) tail (rank 51+ cumulative variance). Classify collapse mechanism.
- Expected signal: Top-5 SV fraction increases from ~18% to ~40%+ across layers; bulk fraction stable or decreasing.
- Null: Isotropic compression would show all fractions stable and eff_rank drop from ambient noise only.
- Value: MEDIUM | Cost: LOW

**H-B3: Per-layer singular value stability (bootstrap)**
- Hypothesis: The top singular values at L11 are stable across random gene subsamples (indicating they represent global structure, not sampling noise).
- Test: Bootstrap 500 subsamples of 500 nonzero genes at L11. Compute top-5 SVs each time. Report CV of each singular value.
- Expected signal: CV < 5% for top-3 SVs, indicating robust global directions.
- Null: Same bootstrap on shuffled embeddings.
- Value: MEDIUM | Cost: LOW

---

### FAMILY C: Topological Reanalysis (Corrected)

**H-C1: H1 persistent homology on nonzero-only genes**
- Hypothesis: The H1 persistence signal from iter_0003 (positive, p<0.05 in 11/12 layers) is stronger and cleaner when restricted to nonzero-only genes (n=2039 vs the contaminated full set).
- Test: Re-run Ripser H1 on nonzero-only gene embeddings at each of 12 layers (PCA-20, n=2039, 20 shuffle replicates). Compare H1 lifetime delta and z-scores to iter_0003 values.
- Expected signal: Higher z-scores at each layer; cleaner layer-pattern; possibly stronger layer dependence.
- Null: Feature-shuffle of nonzero-only embeddings.
- Value: MEDIUM | Cost: MEDIUM

**H-C2: Betti-0 evolution (connected components) across layers**
- Hypothesis: As ID compresses from L0 to L11, the kNN graph develops fewer connected components (Betti-0 decreases), reflecting increasing global connectivity of the gene manifold.
- Test: At each layer, build kNN graph (k=10) on nonzero-only PCA-20 embeddings. Count connected components. Compare to shuffled null.
- Expected signal: Betti-0 decreases monotonically L0→L11; shuffled null has many more components.
- Null: Feature-shuffle at each layer (same n=2039).
- Value: MEDIUM | Cost: LOW

**H-C3: Wasserstein distance between adjacent layer distributions**
- Hypothesis: The Wasserstein-2 distance between consecutive layer embedding distributions peaks at layers with largest per-layer ID drops (L3-L4, L9-L10), providing a transport-geometry characterization of the compression events.
- Test: At each of 11 consecutive layer transitions, compute sliced Wasserstein-2 distance between nonzero-gene embedding distributions. Correlate with per-layer TwoNN ID drop.
- Expected signal: Wasserstein distance and ID drop are correlated (Spearman ρ > 0.5); largest transport at L3-L4 and L9-L10.
- Null: Adjacent layers in shuffled embeddings (should have ~zero Wasserstein distance by construction).
- Value: HIGH | Cost: MEDIUM

---

### FAMILY D: Biological Anchoring (Corrected Substrate)

**H-D1: TRRUST co-regulatory geometry on nonzero-only genes**
- Hypothesis: TRRUST activation pairs are closer in embedding space than random gene pairs, and this signal is detectable on the nonzero-only gene set with proper null.
- Test: On nonzero-only genes (n=2039), extract all TRRUST TF–target pairs where both genes have nonzero norm. Compute mean pairwise distance for activation pairs vs random size-matched pairs. Test at each layer. Repeat for repression pairs.
- Expected signal: Activation pairs closer than null in late layers (L8–L11); repression pairs not closer or farther.
- Null: Gene-label permutation (500 replicates) within nonzero-only gene set.
- Value: HIGH | Cost: MEDIUM

**H-D2: STRING co-expression clustering in embedding space (nonzero-only)**
- Hypothesis: Genes with high STRING combined-score co-expression are closer in embedding space than random pairs, and this proximity signal increases from L0 to L11 as the model converges to functional geometry.
- Test: Extract STRING co-expression edges at score ≥ 700 for nonzero-only gene pairs. Compute mean embedding distance for high-score pairs vs random pairs at each layer. Test monotone trend.
- Expected signal: High-STRING pairs are closer in embedding space (z < -2) in ≥ 6/12 layers; effect grows across depth.
- Null: Permutation of gene-STRING labels within nonzero-only set.
- Value: HIGH | Cost: MEDIUM

**H-D3: Cell-type marker gene clustering in late-layer embeddings**
- Hypothesis: Known immune cell-type marker genes (CD3D/T-cell, CD19/B-cell, CD14/monocyte) cluster together in L11 embedding space more than randomly selected genes.
- Test: Use canonical immune marker lists (5–15 genes per major lineage). Compute within-lineage vs between-lineage mean embedding distance at L11. Permutation null by random gene-set assignment.
- Expected signal: Within-lineage distance significantly lower than between-lineage (z < -2 for ≥ 3/5 lineages).
- Null: Random gene sets of same size; gene-label permutation.
- Value: HIGH | Cost: LOW (small curated gene sets, no gene-names file needed if markers chosen to be in nonzero set)

---

### FAMILY E: Manifold Geometry

**H-E1: Local linear dimensionality (PCA neighborhood)**
- Hypothesis: The intrinsic dimension estimated locally (k-nearest neighborhood PCA) differs from global TwoNN and from the ambient dimension, revealing that the manifold is locally low-dimensional but globally curved.
- Test: At L0 and L11, for each nonzero gene, compute k-NN PCA dimensionality in k=20 neighborhood. Report distribution of local ID estimates. Compare mean local ID to TwoNN global ID.
- Expected signal: Local ID at L11 < TwoNN ID < ambient dim; distribution narrows at L11 (more homogeneous local structure).
- Null: Shuffled embeddings at same layers.
- Value: MEDIUM | Cost: MEDIUM

**H-E2: Anisotropy of compression across layers**
- Hypothesis: Compression from L0 to L11 is anisotropic — certain directions in embedding space compress faster than others, consistent with selective structuring rather than global scaling.
- Test: Compute ratio of top-N SV to total variance at each layer (N=1,5,10,20). Track how variance concentrates into top directions. Quantify anisotropy as Gini coefficient of singular value spectrum at each layer.
- Expected signal: Gini coefficient increases monotonically L0→L11; top-5 SV capture progressively more variance.
- Null: Isotropic compression would show constant Gini coefficient.
- Value: MEDIUM | Cost: LOW

**H-E3: Neighborhood preservation across layers (kNN stability)**
- Hypothesis: Despite ID compression, the k-nearest neighbor graph for nonzero genes is highly preserved across layers — compression is structure-preserving, not random.
- Test: At each consecutive layer pair, compute Jaccard similarity of kNN graphs (k=10) for nonzero-only genes. Compare to shuffled baseline.
- Expected signal: Mean Jaccard > 0.7 for adjacent layers; gradual decay for non-adjacent layers; near-zero for shuffled.
- Null: kNN graph built on shuffled embeddings.
- Value: MEDIUM | Cost: LOW

---

## Top 3 for Immediate Execution

### Candidate 1: High-probability discovery — H-A1 + H-A2 + H-A3 (Gene Names Recovery + Biological Annotation)

**Justification**: The single largest bottleneck is the missing gene names file. Without it, ~70% of the high-value hypotheses (zero-norm biology, SV1/SV2 identity, TRRUST validation) are blocked. Finding it or reconstructing it from the scGPT vocabulary is a LOW-cost task with HIGH unlock value. Execute as a combined bundle:
1. Search all reference root paths for gene name arrays of length 4941.
2. Fall back to scGPT model vocabulary if not found.
3. Once gene names are resolved, immediately run zero-norm vs nonzero GO enrichment (H-A2) and SV1/SV2 correlation with expression level (H-A3).

**Expected runtime**: LOW. File search + annotation is cheap.

---

### Candidate 2: High-risk/high-reward — H-C3 (Wasserstein distance between layer distributions)

**Justification**: Wasserstein-2 distance provides a genuinely novel geometric characterization not yet explored in this project. If it correlates with per-layer ID drop, it provides a continuous, well-grounded metric of representation change that complements TwoNN and CKA. The transport-geometry framing is novel for this type of model analysis and could be the most citable new result. Sliced Wasserstein is numerically tractable at n=2039.

**Expected runtime**: MEDIUM (sliced Wasserstein at n=2039, 11 layer pairs).

---

### Candidate 3: Cheap broad screen — H-B1 + H-E2 (Spectral decay shape + Anisotropy)

**Justification**: Computing the full singular value spectrum at all 12 layers costs almost nothing (SVD on 2039×512 matrix). Tracking spectral shape (power-law exponent, Gini coefficient, top-N fraction) gives a rich low-cost characterization of HOW the 4.9× collapse happens. This directly addresses the question of "dominant directions growing vs bulk shrinking" and produces publication-quality spectral evolution figures. Can be combined into one script with H-B2 and H-E2 at no extra cost.

**Expected runtime**: LOW.
