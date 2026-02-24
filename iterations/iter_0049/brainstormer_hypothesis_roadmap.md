# Brainstormer Hypothesis Roadmap: iter_0049 → iter_0050

---

## Retire / Deprioritize

| Direction | Reason | Decision |
|-----------|--------|----------|
| SWD / optimal transport as primary metric | Weak effect (17% excess), inconsistent across transitions, limited interpretability | `deprioritize` |
| Generic PH on full 512-dim clouds without biological anchoring | Repeated inconclusive results in prior iterations | `deprioritize` |
| Feature-shuffle null comparisons as standalone analysis | Covered by H02; confirmatory role exhausted | `retire_now` |

---

## New Hypothesis Portfolio

### [NH01] SV1 Identity: Expression-Level Axis
**Hypothesis**: SV1 captures mean expression level across cell types, not regulatory function.
**Test**: Correlate per-gene SV1 loading at each layer with that gene's mean expression across all cells in the adata object. Spearman correlation. Also: rank genes by absolute SV1 loading and test enrichment for housekeeping genes (HK gene lists, e.g., Eisenberg & Levanon 2013).
**Expected signal**: r > 0.5 between SV1 loading and mean expression; housekeeping genes over-represented in top SV1 loadings.
**Null**: Randomly shuffled gene-expression ranks; permutation of SV1 loadings.
**Value**: HIGH | **Cost**: LOW (data already loaded, no new database needed)

### [NH02] SV1 Rotation Drivers: Which Genes Cause the 110° Drift
**Hypothesis**: A small subset of genes (cell-type-specific or immune-activated) drives the SV1 drift from L0 to L11, not the whole population.
**Test**: Rank genes by |SV1_L0_loading - SV1_L11_loading|. Test: (a) progressive removal of top movers — does cumulative cosine similarity recover? (b) enrichment of top 50/100 movers for immune cell markers, TFs, or specific GO terms.
**Expected signal**: ~50-100 genes responsible for >80% of rotation; enrichment for immune-activation or cell-type marker genes.
**Null**: Random gene subsets of same size; random permutation of movers.
**Value**: HIGH | **Cost**: LOW (uses h01_sv1_vectors.npy already computed)

### [NH03] SV2-SV4 STRING Co-expression vs TRRUST: Specificity Test
**Hypothesis**: SV2-SV4 proximity is specific to transcriptional regulation (TF→target), not generic co-expression.
**Test**: Run same H03 Mann-Whitney protocol but with STRING functional association pairs (score > 700) as positive pairs vs. TRRUST TF→target pairs. Compare per-layer effect sizes between TRRUST and STRING.
**Expected signal**: If TF-specific: TRRUST effect > STRING effect. If generic co-expression: STRING effect ≥ TRRUST effect.
**Null**: Random gene pairs (same procedure as H03 negatives).
**Value**: HIGH | **Cost**: LOW (H03 infrastructure reusable; STRING flat file download ~10 MB)

### [NH04] SV5-SV10 Biological Axis Mapping
**Hypothesis**: Higher spectral components (SV5-10) encode additional functional categories (e.g., cell-cycle, metabolism, stress) that are distinct from the TF-regulatory axis in SV2-SV4.
**Test**: For SV5-SV8 and SV9-SV12 subspaces separately, run TRRUST proximity test (H03 clone). Additionally: compute GO term intra-category vs. inter-category distance ratios in each SV-subspace.
**Expected signal**: Lower p-values in specific GO term categories per SV-slice; different biological functions in different spectral bands.
**Null**: Same H03 negative controls + gene-position shuffle.
**Value**: HIGH | **Cost**: MEDIUM (3 separate proximity tests + GO annotation needed)

### [NH05] Layer-Mechanism: SV2-SV4 Subspace Alignment Trajectory
**Hypothesis**: The SV2-SV4 subspace undergoes monotonic convergence to a "regulatory-encoding configuration" peaking at L8, then partially reverting at L9-L11.
**Test**: At each layer compute the principal angles between the 3D SV2-SV4 subspace (as a linear subspace of R^512) and L8's SV2-SV4 subspace. Plot principal angle trajectory across layers.
**Expected signal**: Principal angles decrease L0→L8, increase L8→L11 (V-shape centered at L8).
**Null**: Random 3D subspaces of R^512 for angle baseline.
**Value**: HIGH | **Cost**: LOW (uses SVD already computed; subspace angle via scipy.linalg.subspace_angles)

### [NH06] L10 as Spectral Transition Zone: Simultaneous SV1 + SV2-4 Disruption
**Hypothesis**: L10 is a mechanistic transition layer where both SV1 (rotates 37° at L10→L11) and SV2-SV4 (loses regulatory signal at L10) are simultaneously disrupted, suggesting a non-trivial representational reorganization.
**Test**: At L9, L10, L11: compute (a) SV1 loading correlations with L8, (b) SV2-4 TRRUST proximity effect sizes, (c) effective rank (nuclear norm / Frobenius norm ratio) of the embedding matrix. Compare L10 profile to adjacent layers.
**Expected signal**: L10 shows lower effective rank AND lower proximity signal AND highest SV1 rotation rate — all metrics simultaneously abnormal.
**Null**: Compare to random shuffles of layer order.
**Value**: HIGH | **Cost**: LOW (all data in hand)

### [NH07] SV2-SV4 with DoRothEA/CollecTRI: Extended Regulatory Database
**Hypothesis**: The H03 result replicates with a larger, independent TF-target database (CollecTRI or DoRothEA) using the full 2039-gene vocabulary.
**Test**: Download CollecTRI TF-target interactions (freely available). Match against full 2039-gene nonzero set. Re-run H03 proximity test.
**Expected signal**: Replication of L5-L8 significance peak with larger n and independent evidence source.
**Null**: Same random pair control.
**Value**: HIGH | **Cost**: MEDIUM (new database matching; larger n increases computation slightly)

### [NH08] GO Biological Process Spatial Clustering in SV2-SV4
**Hypothesis**: Genes sharing the same GO Biological Process term are spatially clustered in SV2-SV4 subspace more than random gene sets of the same size.
**Test**: For each GO BP term with 10-200 genes present in the 2039-gene set: compute intra-term mean pairwise distance in SV2-SV4 and compare to 1000 random gene sets of same size (z-score). Identify which GO terms are most strongly clustered.
**Expected signal**: Immune regulation, T-cell activation, cytokine signaling terms showing z-score < -2 (clustered).
**Null**: Random gene set permutation.
**Value**: HIGH | **Cost**: MEDIUM (GO annotation parsing + many distance computations; manageable with scipy)

### [NH09] SV2-SV4 Regulatory Signal in Cell-Type Marker Genes
**Hypothesis**: Cell-type marker gene pairs (from CellMarker or PanglaoDB) are NOT closer in SV2-SV4 than random, confirming the SV2-SV4 signal is specific to regulatory proximity not cell-type similarity.
**Test**: Use CellMarker immune cell marker genes for the gene subset present in n=2039. Construct "within-cell-type" positive pairs and "cross-cell-type" negative pairs. Run H03-protocol on SV2-SV4.
**Expected signal**: No significant proximity (p > 0.1 at most layers), unlike TRRUST which shows p<0.001.
**Null**: Same random pair control.
**Value**: MEDIUM | **Cost**: LOW (specificity control, uses H03 infrastructure)

### [NH10] SV2-SV4 TRRUST Effect Size: Graph-Distance Weighted Analysis
**Hypothesis**: Euclidean distance in SV2-SV4 correlates with regulatory network distance (shortest path in TRRUST directed graph), not just binary TF→target.
**Test**: Construct directed TRRUST graph on the 295 named genes. Compute shortest-path distances between all gene pairs (BFS). Compute SV2-SV4 Euclidean distances between same pairs. Spearman correlation at each layer.
**Expected signal**: Negative Spearman r (closer in graph = closer in SV2-SV4) at L5-L8; peak at L8.
**Null**: Random permutation of graph distances.
**Value**: HIGH | **Cost**: LOW (networkx BFS + correlation; n=295 genes is small)

### [NH11] Intrinsic Dimension of SV2-SV4 Cloud Across Layers
**Hypothesis**: Local intrinsic dimensionality of the gene cloud in SV2-SV4 space decreases at L8 (the regulatory encoding peak), reflecting geometric compression of the regulatory module.
**Test**: Use TWO-NN estimator or PCA variance explained in SV2-SV4 subspace on the 295 named genes at each layer. Plot intrinsic dimension vs. layer alongside H03 effect sizes.
**Expected signal**: Dimension drops at L5-L8 (same window as H03 significance); increases at L10.
**Null**: Shuffled gene positions.
**Value**: MEDIUM | **Cost**: LOW

### [NH12] SV1 Anti-Correlation at L0-L1: Regulatory Repulsion Mechanism
**Hypothesis**: At L0-L1, SV1 positively loads on genes that are TF-regulated (SV1 scores are larger for TF-target genes), creating spatial SEPARATION in SV1 while SV2-SV4 has not yet acquired regulatory structure.
**Test**: At L0 and L1, compare mean SV1 loading of TRRUST-positive genes vs. TRRUST-negative genes (Mann-Whitney). Determine sign and effect size. Cross-reference with H03 result (SV1 at L0: p<0.0001, effect=+0.273 means TF pairs have HIGHER distances in SV1, implying TRRUST genes have LOWER SV1 loading magnitude).
**Expected signal**: TRRUST target genes have systematically lower SV1 loading at L0; this inverts or disappears by L5.
**Null**: Random gene pair assignment.
**Value**: HIGH | **Cost**: LOW (uses h01_sv1_vectors.npy + H03 gene pairs)

---

## Top 3 for Immediate Execution

### Candidate A — High-Probability Discovery: NH03 (SV2-SV4 STRING Specificity)
**Why**: The most critical gap in H03 is the confound between TF-specific regulation vs. generic co-expression. This one test determines whether the headline result is specific or generic. Infrastructure is fully reusable. Result will be binary and interpretable within one iteration.
**Setup needed**: Download STRING human protein-protein functional association file, filter score>700, match to 2039-gene set, run H03 protocol.

### Candidate B — High-Risk/High-Reward: NH04 (SV5-SV10 Multi-Axis Mapping)
**Why**: If higher spectral components encode distinct biological functions, the result is a landmark: scGPT's residual stream organizes biological knowledge in a spectral hierarchy (SV1=expression, SV2-4=regulation, SV5-8=?). This would be a major structural interpretability finding.
**Setup needed**: Extend H03 code to loop over SV-slices [5:8], [9:12], [13:16]; run TRRUST + at least one GO-category proximity test per slice.

### Candidate C — Cheap Broad Screen: NH02 (SV1 Rotation Drivers)
**Why**: We have h01_sv1_vectors.npy already computed. Ranking genes by |SV1_L0 - SV1_L11| and running enrichment is a 30-minute analysis. Will immediately answer "what does SV1 encode" by revealing the biology of the genes that rotate most. No new data or databases needed.
**Setup needed**: Load h01_sv1_vectors.npy, compute per-gene delta, rank, match to gene names, run GO/KEGG enrichment on top 50/100 movers.

---

## Spectral Axis Map (Current State)

| Spectral slice | Biological content | Evidence | Status |
|---|---|---|---|
| SV1 | UNKNOWN (NOT regulation; likely expression level) | H01, H03 | Hypothesis NH01/NH02 |
| SV2-SV4 | TRRUST TF→target regulatory proximity | H03 (8/12 layers, effect=0.268) | Confirmed |
| SV5+ | UNKNOWN | — | NH04 |
