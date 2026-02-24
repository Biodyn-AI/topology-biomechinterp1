# Brainstormer Hypothesis Roadmap — iter_0023

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| Proxy persistence entropy (histogram entropy of distances) | **retire_now** | Confirmed negative, invalid proxy for PH. |
| Quintile-binned STRING effects | **retire_now** | Superseded by continuous Spearman; no additional insight. |
| Standalone SVD direction counting | **retire_now** | Subsumed by co-polarity tests. |
| Macrophage/NK/myeloid cell-type expansion | **deprioritize** | Blocked by named-gene availability; not actionable without curated panel. |
| Generic GO BP at flat effect size | **deprioritize** | rho=-0.077 is confirmed but small. Only pursue via deep stratification, not re-replication. |

---

## New Hypothesis Portfolio

### H-A: Dorothea TF-target activation asymmetry — replication with higher N
**Hypothesis**: The TRRUST activation/repression proximity asymmetry (AUROC=0.640 vs 0.459) replicates in Dorothea dataset, with higher-confidence tier pairs (A/B) showing stronger activation proximity than lower-confidence (D/E).
**Test**: Load dorothea_human.tsv. Extract TF-target pairs with direction=activation / direction=repression. Filter for named scGPT genes. Split by confidence tier. Compute AUROC vs non-STRING background at each layer.
**Expected signal**: Activation AUROC >0.60 (tiers A-B), repression AUROC ≈ chance; confidence tier correlates with AUROC magnitude.
**Null**: Uniform AUROC across direction groups and confidence tiers.
**Value**: high | **Cost**: low (CSV already on disk, same infrastructure as H03)

---

### H-B: TF activation hub centrality — are high-fan-out activators geometrically central?
**Hypothesis**: TFs with many activation targets (degree ≥ 5 in TRRUST activation network) occupy more central positions in the embedding manifold (shorter mean distance to all named genes) than low-fan-out TFs or repressors.
**Test**: Build TRRUST activation degree distribution for TF genes present in named gene set. Compute each TF's mean L2 distance to all 209 named genes at each layer. Rank-correlate TF activation degree vs mean distance (negative rho = more central). Compare activation TFs (N≈20-30) vs repression TFs (same N).
**Expected signal**: Spearman(degree, mean_dist) < -0.3 in later layers; activation TFs systematically closer to center than repression TFs.
**Null**: No correlation between fan-out and centrality.
**Value**: high | **Cost**: low (uses same embeddings, TRRUST already loaded)

---

### H-C: GO ontology comparison — BP vs MF vs CC proximity prediction
**Hypothesis**: GO Biological Process Jaccard predicts embedding proximity more strongly than GO Molecular Function or GO Cellular Component, because scGPT is trained on expression context which is process-driven.
**Test**: Fetch GO MF and GO CC annotations for all 209 named genes (same mygene API). Compute Jaccard for 20K pairs for each ontology. Spearman(Jaccard, L2_distance) per layer for each ontology. Compare rho profiles across three ontologies and across layers.
**Expected signal**: GO BP rho > GO MF rho > GO CC rho; or GO MF wins for enzymatic gene clusters.
**Null**: All three ontologies yield equivalent rho.
**Value**: medium | **Cost**: low (same infrastructure as H02, one API call)

---

### H-D: Activation cascade distance — direct vs indirect TF targets
**Hypothesis**: Direct TF→target pairs (distance 1 in TRRUST activation graph) cluster more tightly than 2nd-order pairs (TF→X→target, distance 2), indicating that geometric proximity encodes regulatory immediacy.
**Test**: Build directed activation graph from TRRUST. Find TF→target 1st-order and TF→X→target 2nd-order pairs within named genes. Compare L2 distance distributions (MW test, AUROC) for 1st-order vs 2nd-order at each layer.
**Expected signal**: 1st-order AUROC > 2nd-order AUROC vs non-STRING background; gap widens with depth.
**Null**: Equal proximity for both cascade levels.
**Value**: high | **Cost**: medium (need to build activation graph, 2nd-order pairs may be sparse)

---

### H-E: Repressor TF anti-clustering — are repression pairs pushed apart?
**Hypothesis**: Repressor TF-target pairs (TRRUST direction=Repression) are not merely at chance distance — they are actively anti-clustered (mean distance > background), meaning scGPT geometry repels repressor pairs.
**Test**: One-sided Mann-Whitney: is repression pair L2 distance *greater* than non-STRING background? Compute effect size and AUROC (expected >0.5 for anti-clustering). Test at each of 12 layers.
**Expected signal**: AUROC > 0.55 (repressors farther than background), sig in ≥ 6/12 layers, growing with depth.
**Null**: Repression pairs at chance distance (AUROC ≈ 0.50).
**Value**: high | **Cost**: low (H03 data already available; just flip the MW direction)

---

### H-F: Cross-model cell-type geometry transfer — scGPT vs Geneformer
**Hypothesis**: Cell-type marker gene clustering (AUROC=0.851 in scGPT) replicates in Geneformer word embeddings, and the two models' pairwise distance matrices are correlated (Mantel test or CKA), indicating architecture-agnostic cell-type geometry.
**Test**: Load Geneformer word embeddings (iter_0019 artifact). Extract embeddings for 17 cell-type marker genes. Compute within-type vs cross-type AUROC. Mantel test between scGPT pairwise L2 matrix and Geneformer pairwise L2 matrix for the 17 genes. CKA between full 209×D matrices.
**Expected signal**: Geneformer cell-type AUROC > 0.75; Mantel r > 0.5, p < 0.01; CKA > 0.5.
**Null**: Geneformer shows no cell-type structure or structure is uncorrelated with scGPT.
**Value**: high | **Cost**: medium (Geneformer embeddings need to be located/recomputed from iter_0019)

---

### H-G: GO BP pathway-level geometric clustering
**Hypothesis**: Specific biological pathways (immune signaling: GO:0002682 / extracellular matrix: GO:0030198 / metabolic: GO:0008152) form distinct geometric clusters in scGPT embedding space, identifiable by inter-cluster distance exceeding intra-cluster distance.
**Test**: Assign named genes to top-3 GO BP pathways by member count. Compute within-pathway vs cross-pathway L2 distance distributions (MW + AUROC) at each layer. Compare AUROC to cell-type AUROC (0.851) — is process-level clustering comparable to cell-type clustering?
**Expected signal**: Within-pathway AUROC > 0.65, substantially below cell-type AUROC, indicating cell-type is a stronger organizer than process.
**Null**: No pathway-level clustering (AUROC ≈ 0.5).
**Value**: medium | **Cost**: low

---

### H-H: Co-activation partner sharing — genes sharing ≥2 common activation TFs cluster
**Hypothesis**: Gene pairs that share 2+ common activation TFs (co-regulated by same activators) are geometrically closer than gene pairs sharing 0 activation TFs, independent of STRING PPI.
**Test**: Build co-activation TF sharing matrix from TRRUST (activation direction only). Define: shared_TF_count per gene pair. Filter pairs not in STRING. Spearman(shared_TF_count, L2_dist) per layer. AUROC for pairs sharing ≥2 TFs vs pairs sharing 0 TFs.
**Expected signal**: Spearman rho < -0.15 in mid/late layers; AUROC > 0.60 for shared-TF pairs.
**Null**: No correlation between shared TF count and proximity.
**Value**: high | **Cost**: low (TRRUST + existing embeddings)

---

### H-I: STRING physical interaction sub-network (CORUM complexes) vs functional STRING
**Hypothesis**: Genes that are members of the same protein complex (CORUM database) are embedded closer than genes with only functional STRING links, because scGPT training captures co-expression of complex subunits more directly than general PPI.
**Test**: Download CORUM human complexes (already in reference root or downloadable). Identify complex-member gene pairs within named genes. Compute AUROC: complex-member pairs vs STRING-only (non-complex) pairs at each layer.
**Expected signal**: CORUM complex AUROC > STRING-all AUROC (currently 0.614); gap widens at later layers.
**Null**: CORUM shows equal or lower proximity vs generic STRING functional.
**Value**: medium | **Cost**: medium (CORUM download + gene name matching)

---

### H-J: Layer-wise local intrinsic dimension profile
**Hypothesis**: Local intrinsic dimension (estimated by two-NN method or PCA on k-nearest neighbors) decreases with layer depth, reflecting progressive manifold compression — and cell-type marker neighborhoods have lower intrinsic dimension than random gene neighborhoods.
**Test**: For each of 12 layers, compute: (a) global intrinsic dimension of all 209 embeddings via two-NN estimator; (b) local intrinsic dimension in 5-NN neighborhood of each cell-type marker gene; (c) compare cell-type neighborhood ID vs random gene neighborhood ID.
**Expected signal**: Global ID decreases from L0 to L11; cell-type marker neighborhoods have ID 2-3 lower than background.
**Null**: Flat ID profile, no cell-type neighborhood compression.
**Value**: medium | **Cost**: medium

---

### H-K: Gene essentiality (DepMap) vs embedding centrality
**Hypothesis**: Essential genes (DepMap gene effect score < -1, pan-essential) occupy more central positions in the embedding manifold (lower mean distance to all other genes) because scGPT training data reflects the fact that essential gene expression is tightly co-regulated with many partners.
**Test**: Load DepMap CRISPR gene effect data. Classify named genes as essential / non-essential. Compute each gene's mean L2 distance to all 208 other named genes at each layer. Mann-Whitney: essential vs non-essential centrality at each layer.
**Expected signal**: Essential genes are significantly more central (lower mean distance) in later layers (L8-L11), AUROC > 0.65.
**Null**: No centrality difference between essential and non-essential genes.
**Value**: medium | **Cost**: medium (DepMap download, gene ID matching)

---

### H-L: Chromosomal proximity vs embedding proximity
**Hypothesis**: Genes co-localized on the same chromosome band (cytogenetic band) are NOT more proximal in scGPT embeddings, confirming that geometric structure reflects functional/regulatory biology not genomic co-localization (negative control test).
**Test**: Fetch chromosome band annotations for 209 named genes via mygene. Define same-band pairs vs different-band pairs. AUROC: same-band proximity vs different-band at each layer.
**Expected signal**: AUROC ≈ 0.50 (no effect) — this is a designed negative control that strengthens the biological anchor story.
**Null** (here the null *is* the expected result): Significant same-band clustering would imply genomic, not biological, structure.
**Value**: medium | **Cost**: low (mygene API already configured)

---

### H-M: GO BP term specificity (IC) modulates proximity prediction
**Hypothesis**: Shared GO BP terms with high information content (rare, specific terms) predict embedding proximity more strongly than shared broad terms, indicating scGPT captures functional specificity not generic pathway membership.
**Test**: Compute information content (IC) for each GO BP term as -log(freq). For each gene pair, compute IC-weighted Jaccard (sum IC of shared terms / sum IC of union terms). Spearman(IC-Jaccard, L2_dist) per layer. Compare to unweighted Jaccard rho (-0.077).
**Expected signal**: IC-weighted Jaccard Spearman rho < -0.10 (stronger than unweighted), peak at L7.
**Null**: IC-weighting does not improve rho vs unweighted Jaccard.
**Value**: medium | **Cost**: low (GO annotations already loaded from H02)

---

### H-N: STRING confidence tier stratification — does physical evidence score predict proximity better than co-expression?
**Hypothesis**: STRING pairs where the evidence is dominated by physical interaction (experiments score > 0.5) cluster more tightly than pairs dominated by co-expression evidence (coexpression score > 0.5), because physical interaction is more direct and expression co-regulation is more diffuse.
**Test**: Download STRING v12 detailed scores for named gene pairs. Split pairs into: experiment-dominant (exp > coexp) vs coexpress-dominant (coexp > exp). AUROC at each layer. Compare the two groups.
**Expected signal**: Experiment-dominant AUROC > coexpress-dominant AUROC in later layers (L8-L11).
**Null**: Equal AUROC for both evidence types.
**Value**: medium | **Cost**: medium (STRING detailed scores download)

---

## Top 3 for Immediate Execution

### Candidate 1 — HIGH PROBABILITY DISCOVERY
**H-A: Dorothea TF-target activation asymmetry replication**
- Consolidates the iter_0023 breakthrough with N ~3-5x larger dataset and confidence tier stratification.
- Low cost: Dorothea CSV + same H03 infrastructure.
- Clear prediction: activation AUROC > 0.60 in tiers A-B; repression ≈ chance.
- If confirmed, this becomes the central mechanistic finding of the paper.

### Candidate 2 — HIGH RISK / HIGH REWARD
**H-B: TF activation hub centrality**
- If high-fan-out activation TFs are geometrically central, this implies scGPT has learned an activation-centric topology where TFs radiate outward to their targets. Highly novel structural claim.
- Risk: N is small (~20-30 activation TFs in named genes); result may not be statistically robust.
- If positive: directly testable claim about model geometry mirroring TF regulatory architecture.

### Candidate 3 — CHEAP BROAD SCREEN
**H-C: GO ontology comparison (BP vs MF vs CC)**
- One additional mygene API call, same Spearman infrastructure as H02.
- Will either confirm GO BP as the best predictor (supporting the "process-driven" interpretation) or reveal that GO MF is stronger (pivoting toward molecular function encoding).
- Either result is informative; no wasted compute if negative.
