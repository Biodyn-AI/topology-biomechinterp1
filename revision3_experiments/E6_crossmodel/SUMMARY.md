# E6 — Cross-model comparison (verified)

Copy of the verifier-rewritten SUMMARY.txt (verify/ holds the checks).

```
E6_crossmodel: cross-model comparison (sensitivity analysis; no pass/fail). Answers R2.4 ("why only 186 genes / selection bias"), R1.7, R2.8, and supports R1.2.
CONVENTION3 sha256 3e3249d1d0a8ce1bc3cf64639fe233b016535456f82432a8bbfcacdd64a05c01 (matches FROZEN.txt; recorded in results.json).
This file was rewritten by the independent verifier. The implementer's text is kept in verify/SUMMARY_implementer_before_verification.txt.
Corrections are marked [verifier correction] or [verifier addition]; the VERIFICATION section at the end lists them and says what was reproduced.

PROVENANCE AND STATUS
Three implementer agents and one verifier did this work. Agent 1 wrote e6lib.py, e6_build.py, e6_run.py, e6_selection.py, e6_legacy.py,
legacy_submitted_repro.py and e6_assemble.py. Agent 2 made three corrections:
  (a) e6_patch1d.py. Estimator (iii) had scored its twelve 1-D bands through stats3.l2norm_rows (|x|+1e-12), which breaks exact sign ties by
      magnitude. The bands were re-scored with exact sign normalisation on the same null draws, and max-T was recomputed.
  (b) e6_tierb_static_maxt.py completed estimator (iii) for the static models in Tier B. Estimators (i)/(ii) came out bit-identical
      (660 rows, 0 differences; tbmaxt_check.log).
  (c) The ESM2 file was relabelled, and delete-d intervals were added in e6_assemble.py.
Agent 3, after the 22:47 reboot, re-ran e6_fixreport.py and e6_assemble.py. It also added e6_supp.py, which gives:
  - a gene-bootstrap CI (B=1,000) and a within-decile permutation p (N=2,000) for the detection-stratified Lambert AUROC (node_detstrat_results.csv);
  - paired deltas vs scGPT for every task (pairwise_deltas_vs_scgpt.csv).
  It regenerates the draws with e6_run's own seeds and rng streams, so they are paired across embeddings, and it reproduces the stored values and
  CIs to 1e-9.
PPI intervals are given twice: the pre-registered sqrt(m/n) interval (ci_lo/ci_hi) and the delete-d sqrt(m/(n-m)) interval (ci_*_deleted;
CONVENTION3_ERRATA item 3). The manuscript should use the delete-d interval. Every PPI interval quoted below is delete-d.
Every number below comes from attrition.csv, tierA_results.csv, tierB_results.csv, node_detstrat_results.csv, pairwise_deltas.csv,
pairwise_deltas_vs_scgpt.csv, legacy_reproduction.csv, selection_bias.csv, projection_control.csv, p1d_fix_effect.csv, results.json, or the
verifier files verify/verify_results.json, verify/reproduced_numbers.csv, verify/centring_sensitivity.csv and verify/centring_delta.csv.

MODELS (static gene vectors; file provenance in results.json "models")
- scGPT static token embedding (whole-human checkpoint, 512-d).
- Geneformer V2-316M (1152-d). Identified from the subproject_54 PLAN: HF snapshot 05fcbeb8, root model.safetensors, 18 layers; word embeddings with
  the embedding LayerNorm applied.
- Geneformer V1-10M (256-d) and Geneformer V2-104M (768-d).
- AIDO.Cell-100M (640-d) and scFoundation (512-d). Neither is a token table: each vector is the output of one forward pass on an HPA mean-expression
  "neutral" cell.
- tGPT (1024-d).
- ESM2-15B (UCE per-protein; 5,120-d). The submitted manuscript called this "ESM2-3B", but 5,120 is the ESM2-15B width; ESM2-3B is 2,560-d.
- ESM2-650M (1280-d, own mean-pooled embedding) and 32 hand-built sequence descriptors.
- HPA co-expression (PCA-128).
- scGPT contextual residual stream, L0-L11, in lung, lung_external, kidney_B, immune and kidney_A.
Estimator (i), the primary one, uses the full dimension and centres the vectors over the gene set being analysed (Tier A, a Tier B set or a legacy
set). value_vocab_centred holds the vocabulary-centred value. Two deviations from CONVENTION3 follow from this; both are allowed by the E6 spec,
which centres over the gene set:
- scGPT contextual vectors are centred over the Tier B set, not over V (CONVENTION3 §2).
- The TRRUST task is scored on labelled genes, but the vectors are centred over every gene in the set. See the [verifier addition] in section 4.
Effect of the first deviation on lung, detection-stratified Lambert AUROC: L1 0.668 in Tier B vs 0.664 on all of V centred over V; L11 0.638 vs 0.625
(verify/verify_results.json).

1. WHY THE SUBMITTED COMPARISON HAD 186 GENES (attrition.csv) -- direct evidence
- Vocabularies are large: 19,264-60,694 genes. The three submitted models (Geneformer V1-10M, V2-104M, ESM2-15B) share 18,907 genes. All 11 static
  vocabularies share 17,693 (Tier A). The vocabulary intersection was therefore never the bottleneck.
- The 186 came from the gene list, step by step:
  - The submitted code started from the 209 genes named in the round-1 kidney_A edge dataset, cycle1_edge_dataset.tsv. That file has 1,152 TRRUST
    edge and non-edge rows, so the list includes genes that appear only in sampled non-edges.
  - 187 of the 209 have a non-zero kidney_A embedding (22 were never observed).
  - Intersecting with the three vocabularies removes only NEAT1, a lncRNA, which leaves 186.
  - Intersecting with all 11 vocabularies also removes 8 HLA genes that tGPT lacks, which leaves 178 ("legacy178").
- [verifier addition] The choice of context was itself limiting. kidney_A observes only 1,972 of its 4,803 genes, and only 180 of the observed genes
  take part in a TRRUST edge whose two ends are both observed (172 within Tier A). A kidney_A-only analysis could therefore not have been much larger
  for the TRRUST tasks. Larger sets need other contexts or the static-only Tier A.
- The same tasks fit on much larger sets:
  - Tier A: 17,693 genes, 1,587 Lambert TFs, 775 TRRUST TFs vs 1,978 target-only genes, 73,430 STRING-physical>=700 pairs, and 7,593 TRRUST edges
    (430 TFs) vs 70,882 matched candidate non-edges.
  - Tier B (each context's V intersected with Tier A): lung 5,226; lung_external 5,154; kidney_B 4,441; immune 1,338; kidney_A 1,298.

2. IS THE SMALL SET SELECTION-BIASED? (selection_bias.csv; 10,000 matched random subsets) -- direct evidence
- The null draws random 186-gene subsets of the 2,753 TRRUST-labelled Tier A genes. Each subset matches the legacy set on size, on the TF:target-only
  ratio (64:122) and on lung mean-expression decile (13 legacy genes outside lung V form an "NA" cell). The statistic is the TRRUST TF/target-only
  AUROC under estimator (i), centred over the subset, which is the same convention as for the legacy value.
- In all 11 static models the legacy value lies at the 98.8th-100th percentile of the matched null:
  - legacy186: percentile >= 0.995, two-sided p 0.0002-0.010, excess over the null mean +0.088 to +0.250;
  - legacy178: percentile >= 0.988, excess +0.075 to +0.248.
  - Examples (legacy vs null mean): scGPT static 0.833 vs 0.739; Geneformer V2-104M 0.874 vs 0.745; ESM2-15B 0.954 vs 0.862.
- The excess survives an additional match on TRRUST-degree quartiles (legacy186_degmatched: percentile 0.979-1.000, excess +0.060 to +0.249). The
  legacy genes are better studied than the pool: mean TRRUST degree for TFs is 16.6 vs 10.8, and for targets 7.0 vs 3.2. Degree matching alone does
  not remove the excess.
- The model ORDERING is largely preserved. Kendall tau, legacy vs Tier A:
  - 0.867 [bootstrap 0.689, 0.956] over 10 models (legacy186);
  - 0.855 [0.673, 0.964] over 11 models (legacy178).
  [verifier addition] When Tier A is centred, like the legacy set, over its TRRUST-labelled genes, legacy178 vs Tier A gives tau = 0.964 (point
  estimate only; verify/centring_sensitivity.csv).
- Interpretation (inference): the submitted 186-gene set is an "easy" set, roughly +0.06 to +0.25 above matched random sets of the same composition.
  It changes the ranking only slightly. Absolute legacy AUROCs should not be quoted as model properties. [verifier correction] The gap between legacy
  and Tier A values is not all selection: for expression-trained models, part of it comes from the centring set (section 4).
- Surprise: scFoundation's subset-null mean is 0.30, far below 0.5. Its Tier A value is 0.487 centred over Tier A and 0.496 centred over the labelled
  genes, so the centring set does not explain the gap. The cause was not investigated.

3. REPRODUCTION OF THE SUBMITTED TABLE 1 (legacy_reproduction.csv; legacy_submitted_repro.log) -- direct evidence
- The submitted code (revision2_experiments/crossmodel.py: SVD fitted on the 186-gene subset, SV2-SV7 band) reproduces the submitted values exactly
  (|diff| < 1e-12). The scGPT rows come from the kidney_A L11 matrix; all 186 genes are observed there, so zero rows play no part on this set.
  - TF/target AUROC: 0.650 / 0.713 / 0.837 / 0.924 for scGPT L11 / Geneformer V1-10M / Geneformer V2-104M / ESM2-15B (the submitted "ESM2-3B").
  - Co-pole z: +2.01 / -0.39 / +4.94 / +1.91, on 809 STRING functional pairs.
- Other estimators on the same 186 genes:

  model             | subset-SVD (submitted) | SVD over full vocab / V, SV2-SV7 | estimator (i) [95% CI]
  scGPT L11         | 0.650                  | 0.611                            | 0.684 [0.587, 0.769]
  Geneformer V1-10M | 0.713                  | 0.592                            | 0.802 [0.676, 0.909]
  Geneformer V2-104M| 0.837                  | 0.460                            | 0.874 [0.781, 0.946]
  ESM2-15B          | 0.924                  | 0.942                            | 0.954 [0.911, 0.987]

  The submitted Geneformer band values therefore depended on fitting the SVD to the 186-gene subset.
- Estimator (i), other tasks, on legacy186 (CIs are wide on 186 genes):
  - Lambert AUROC: 0.698 / 0.816 / 0.855 / 0.968.
  - STRING-physical PPI AUROC, 163 pairs: scGPT L11 0.624 [0.496, 0.754]; Geneformer V1-10M 0.644 [0.538, 0.766]; Geneformer V2-104M 0.688
    [0.582, 0.803]; ESM2-15B 0.772 [0.690, 0.870].
  - Within-TF edge AUROC, 213 edges, 29 TFs: scGPT L11 0.599 [0.496, 0.720]; Geneformer V1-10M 0.567 [0.464, 0.694]; Geneformer V2-104M 0.618
    [0.511, 0.740]; ESM2-15B 0.552 [0.446, 0.657].
- The submitted table compared scGPT at L11, its weakest layer for this task. On legacy186, TRRUST AUROC is 0.79-0.82 at L0-L7 (L2 0.817).

4. BROAD UNIVERSE, TIER A (17,693 genes; tierA_results.csv, node_detstrat_results.csv, pairwise_deltas_vs_scgpt.csv) -- direct evidence
Values are estimator (i), with 95% CIs in brackets. "Delta" means the model minus scGPT static on the same genes (paired bootstrap).
- Lambert TF, detection-stratified (the H2 statistic). Every within-decile permutation p = 1/2001.
  - 32 seq descriptors 0.916 [0.908, 0.924]; ESM2-15B 0.902 [0.893, 0.911]; Geneformer V2-104M 0.864; ESM2-650M 0.861; Geneformer V2-316M 0.855.
  - Geneformer V1-10M 0.776; scGPT static 0.747 [0.734, 0.760]; tGPT 0.690; AIDO 0.623; HPA co-expression 0.617; scFoundation 0.554.
  - Deltas: ESM2-15B +0.155 [+0.141, +0.168]; seq descriptors +0.169 [+0.155, +0.182]; Geneformer V2-104M +0.117 [+0.104, +0.130];
    Geneformer V1-10M +0.029 [+0.015, +0.042]; scFoundation -0.193 [-0.214, -0.175].
- Lambert, unstratified: the same order. ESM2-15B 0.889 [0.882, 0.897]; scGPT static 0.704 [0.694, 0.713].
- TRRUST TF vs target-only, centred over all of Tier A:
  - ESM2-15B 0.870 [0.856, 0.884]; seq descriptors 0.869; ESM2-650M 0.826; Geneformer V2-104M 0.774; Geneformer V2-316M 0.743; Geneformer V1-10M 0.688.
  - scGPT static 0.684 [0.660, 0.708]; HPA 0.615; AIDO 0.579; tGPT 0.520; scFoundation 0.487 (p 0.62).
  - Delta: ESM2-15B +0.186 [+0.159, +0.212].
  [verifier addition: centring sensitivity] For expression-trained models this number depends heavily on the centring set. Centring over the 2,753
  TRRUST-labelled genes only (the convention every legacy and subset value uses) instead of over all 17,693 Tier A genes changes:
    scGPT static 0.684 -> 0.789; AIDO 0.579 -> 0.690; tGPT 0.520 -> 0.639; Geneformer V1-10M 0.688 -> 0.755; V2-316M 0.743 -> 0.802;
    V2-104M 0.774 -> 0.821.
  Protein models, the seq descriptors, HPA and scFoundation change by <= 0.010 (ESM2-15B 0.870 -> 0.877).
  ESM2-15B minus scGPT static becomes +0.088 [+0.066, +0.111] instead of +0.186. Geneformer V2-104M minus scGPT static becomes +0.032
  [+0.014, +0.051] instead of +0.089. The ordering survives but the gaps roughly halve (verify/centring_delta.csv).
  In Tier B lung the convention hardly matters: ESM2-15B minus scGPT static is +0.078 under both.
  Inference: the about 15,000 Tier A genes that carry no TRRUST label shift the mean of the expression-model tables in a way that depresses their
  TRRUST contrast. The Tier A TRRUST gap between protein and expression models is therefore not convention-free. Quote it with both conventions,
  or use Tier B.
- PPI pair-cosine AUROC (delete-d CI):
  - seq descriptors 0.671 [0.660, 0.680]; scGPT static 0.669 [0.658, 0.679]; Geneformer V2-104M 0.662; Geneformer V2-316M 0.659;
    ESM2-15B 0.647 [0.638, 0.657].
  - ESM2-650M 0.610; AIDO 0.603; HPA 0.598; Geneformer V1-10M 0.588; tGPT 0.571; scFoundation 0.559.
  - Deltas: ESM2-15B -0.022 [-0.032, -0.011]; seq descriptors +0.002 [-0.009, +0.012].
- Within-TF edge AUROC (two-way bootstrap CI):
  - Geneformer V2-104M 0.580 [0.549, 0.610]; Geneformer V2-316M 0.578; scGPT static 0.575 [0.534, 0.610]; AIDO 0.552 [0.499, 0.591]; tGPT 0.544.
  - Geneformer V1-10M 0.529; HPA 0.518; scFoundation 0.500; seq descriptors 0.497; ESM2-650M 0.492; ESM2-15B 0.491 [0.466, 0.523] (p 0.97).
  - Delta: ESM2-15B -0.084 [-0.124, -0.042].
- Caveat: detection deciles in Tier A come from lung. 12,467 of the 17,693 genes are not in lung V and are all placed in the median decile, so
  detection stratification is largely nominal in Tier A. Tier B is the detection-controlled comparison.

5. BROAD UNIVERSE, TIER B (scGPT contextual vs static; tierB_results.csv, node_detstrat_results.csv, pairwise_deltas_vs_scgpt.csv) -- direct evidence
scGPT contextual, detection-stratified Lambert AUROC across L0-L11 (best layer in brackets):
  - lung 0.638-0.668 (L1 0.668 [0.652, 0.687]); lung_external 0.631-0.675 (L1); kidney_B 0.616-0.657 (L1);
  - immune 0.616-0.746 (L2 0.746 [0.677, 0.819]); kidney_A 0.637-0.763 (L3).
  In every context the AUROC declines towards L11.
Node identity, deltas on the same genes:
  - ESM2-15B minus scGPT contextual, detection-stratified: lung +0.212..+0.241; lung_external +0.206..+0.250; kidney_B +0.204..+0.245;
    immune +0.187..+0.317; kidney_A +0.151..+0.277. The CI excludes 0 at 12/12 layers in all five contexts.
  - Same comparison on TRRUST TF/target-only: +0.112..+0.300, 12/12 layers in every context.
  - [verifier correction] 32 seq descriptors: detection-stratified Lambert +0.137..+0.296; TRRUST +0.121..+0.301; 12/12 layers everywhere.
  - scGPT's own static token embedding beats its residual stream on Lambert (detection-stratified): +0.083..+0.113 (lung), +0.056..+0.100
    (lung_external), +0.080..+0.120 (kidney_B), +0.065..+0.195 (immune), each at 12/12 layers; +0.016..+0.142 in kidney_A (9/12).
  - Geneformer V2-104M also beats it: +0.157..+0.202 in lung, lung_external and kidney_B (12/12); +0.069..+0.199 in immune (11/12);
    +0.042..+0.167 in kidney_A (7/12).
PPI:
  - Contextual best-layer AUROC (delete-d CI): lung 0.674 [0.658, 0.690]; lung_external 0.678 [0.660, 0.694]; kidney_B 0.671 [0.649, 0.688];
    immune 0.742 [0.701, 0.789]; kidney_A 0.786 [0.741, 0.831].
  - ESM2-15B minus contextual ranges from -0.055 to +0.022. Its delete-d CI lies below 0 at 5/12 (lung), 6/12 (lung_external), 6/12 (kidney_B),
    0/12 (immune) and 2/12 (kidney_A) layers, and above 0 only at 1 lung layer. The per-row flag is ci_excludes_0_deleted [verifier addition].
  - HPA co-expression and scFoundation fall below contextual scGPT at nearly every layer.
Edge:
  - Contextual within-TF AUROC: lung 0.515-0.542; lung_external 0.520-0.539; kidney_B 0.522-0.553; immune 0.576-0.632; kidney_A 0.618-0.664.
  - ESM2-15B is at chance in every Tier B set (0.483-0.529). ESM2-15B minus contextual is negative at every layer (-0.017..-0.149). The CI excludes 0
    at 3/12 (lung), 3/12 (lung_external), 3/12 (kidney_B), 8/12 (immune) and 2/12 (kidney_A) layers.
  - [verifier correction] scGPT static and Geneformer V2-104M lie slightly above contextual on edges in lung, lung_external and kidney_B (+0.015..+0.051;
    CI > 0 at 3-4/12 layers). In immune and kidney_A no layer's CI excludes 0.
Cross-layer max-T over 936 (layer, band) cells for scGPT contextual (estimator bestband_alllayers): adjusted p <= 0.002 for every task and context.
This is a selection-adjusted detection, not an effect size. [verifier addition] For the node tasks, estimator (iii) and this max-T use the
UNSTRATIFIED AUROC, so they do not control detection.

6. ESM2 PROJECTION CONTROL (Tier A; projection_control.csv) -- direct evidence
- ESM2-15B projected 5,120 -> 1,024 by random Gaussian projection (5 seeds), mean (sd) vs full: Lambert 0.887 (0.007) vs 0.889; TRRUST 0.870 (0.007)
  vs 0.870; PPI 0.648 vs 0.647; edge 0.490 vs 0.491.
- scGPT static projected 512 -> 102 (the same 5:1 ratio): Lambert 0.678 vs 0.704 (-0.027); TRRUST 0.635 vs 0.684 (-0.049); PPI 0.634 vs 0.669
  (-0.036); edge 0.556 vs 0.575 (-0.019).
- Interpretation (inference): ESM2's node-identity advantage does not come from its larger dimension. At 1,024-d it still exceeds full 512-d scGPT
  static by about +0.18 on unstratified Lambert. The 32-d hand-built sequence descriptors reach 0.90, so a very low-dimensional protein-sequence
  representation suffices.

7. ESTIMATORS (ii) AND (iii) (Tier A; tierA_results.csv) -- direct evidence
- The fixed bands SV2-SV7 / SV1-SV6, taken from an SVD over each model's whole vocabulary, often recover much less than the full dimension:
  - Geneformer V2-104M TRRUST 0.448 / 0.549 vs full 0.774;
  - Geneformer V2-316M Lambert 0.530 / 0.492 vs 0.817;
  - scGPT static TRRUST 0.513 / 0.493 vs 0.684.
- [verifier correction] The implementer wrote that the best contiguous band "reaches roughly the full-dimension value". That is true only for some
  models. Estimator (iii) is the best of 78 bands within SV1-SV12, with max-T over bands, N=1,000, unstratified AUROC.
  - Close to or above full: the protein/sequence models (ESM2-15B Lambert SV6-SV12 0.913, adjusted p 0.001, vs full 0.889; ESM2-650M 0.916 vs
    0.859; seq descriptors 0.897 vs 0.902) and HPA (0.602 vs 0.576).
  - Well below full: the Geneformers (V2-104M Lambert SV2-SV12 0.630 vs 0.836; V2-316M 0.633 vs 0.817; V1-10M 0.571 vs 0.742), scGPT static
    Lambert (0.675 vs 0.704) and tGPT (0.589 vs 0.639).
  - Inference: in the Geneformer tables, TF identity is spread beyond their first 12 directions.
- Static Tier B rows with max-T p < 0.05: 49/55 (Lambert), 51/55 (TRRUST), 55/55 (PPI), 40/55 (edge). The non-significant rows are mostly
  Geneformer V1-10M and scFoundation, plus ESM2 and the seq descriptors on edges.
- The exact-sign 1-D fix (p1d_fix_effect.csv, 684 rows) changed significance at 0.05 for the best band in 9 rows. Examples: legacy178 Geneformer
  V1-10M PPI 0.001 -> 0.45; legacy186 Geneformer V1-10M edge 0.036 -> 0.64; tierB_kidney_A Geneformer V1-10M edge 0.001 -> 0.20; tierB_immune
  ESM2-650M edge 0.47 -> 0.008. All reported outputs use the fixed values; the pre-fix values are kept in the *_pre1dfix columns.

WHAT THIS MEANS (inference, kept separate from the evidence)
- The submitted claim about ordering on TF identity (protein LM > Geneformer > scGPT) holds on 17,693 genes and in every context's
  detection-controlled Tier B set. The absolute numbers in the submitted Table 1 are inflated by an easy, curated gene set, and two of its four model
  values depended on a subset-fitted SVD band. [verifier correction] The size of the Tier A TRRUST gap depends on the centring convention (+0.186 vs
  +0.088 for ESM2-15B minus scGPT static). The detection-stratified Lambert and Tier B gaps are the more robust numbers to quote.
- TF identity is carried most strongly by protein-sequence representations, including the 32 hand-built descriptors. That suggests "is a TF" is largely
  a protein-family / sequence property (DNA-binding domains, composition), and that recovering it from an scFM residual stream says little about
  learned regulatory biology. This is a hypothesis; E6 does not test mechanism.
- The ranking reverses for regulatory edges. Once TF identity is controlled within TF, protein LMs and sequence descriptors are at chance, while
  expression-trained models (scGPT contextual and static, Geneformer V2, AIDO) are modestly above it (0.52-0.66). Where edge-level signal exists,
  it is therefore expression-derived, not sequence-derived.
- PPI neighbourhood is recovered at broadly similar levels by most models (full-dimension AUROC 0.54-0.79 across Tier A and Tier B sets, scFoundation lowest), and no model family dominates.

WHAT THIS DOES NOT SHOW / LIMITATIONS
- E6 is a sensitivity analysis. Its p-values are not part of the pre-registered H1-H4 families and carry no confirmatory weight.
- Edge permutation p-values (within-TF label permutation) ignore target-level dependence: 49 of the 104 full-estimator edge rows with p < 0.05 have a
  two-way bootstrap CI that includes 0.5. Use the two-way bootstrap interval.
- PPI non-pairs come from ONE Maslov-Sneppen rewiring, and the vertex-subsampling CI does not include variation between rewiring draws.
  [verifier addition] An independent rewiring by the verifier gave Tier A PPI AUROCs 0.003 higher for both ESM2-15B (0.650 vs 0.647) and scGPT static
  (0.672 vs 0.669), within the delete-d CIs; the paired delta was unchanged (-0.022).
- Detection stratification in Tier A is mostly nominal (70% of genes imputed to the median lung decile). Static models' detection deciles are a
  documented proxy: lung n_obs where available, else the median.
- Centring convention: estimator (i) centres over the analysed gene set. For the TRRUST task in Tier A this matters by up to +0.12 for
  expression-trained models (section 4). Vocabulary-centred values (value_vocab_centred) differ again, e.g. scGPT static Lambert 0.660 vs 0.704.
- AIDO.Cell and scFoundation vectors are forward-pass outputs on one neutral cell, not learned token tables, so they are not like-for-like with the
  other static models.
- The scGPT contextual embeddings come from one extraction per context (192-256 cells), so deltas against contextual layers are conditional on that
  extraction.
- Best layers and best bands are selected post hoc, which is why ranges over all 12 layers are reported.
- The selection-bias null matches on size, TF ratio and lung expression decile, plus degree in one variant. It cannot rule out other properties of the
  curated set, such as annotation depth. It shows that the set is unrepresentative, not why. (work/selection_meta.json's "composition" field for the
  degree-matched variant collapses the expression stratum in its keys. This is a display issue only; the draws used the full cells.)
- Identity labels (Lambert, TRRUST) and STRING carry known study bias, and better-studied genes may be easier for every model.
- pairwise_deltas.csv holds node-task deltas for all static model pairs and against every contextual layer. pairwise_deltas_vs_scgpt.csv holds all
  tasks against scGPT. Contextual-vs-contextual layer deltas were not computed.

VERIFICATION (independent adversarial verifier)
Scripts and outputs: verify/verify_e6.py (part1-3), verify/verify_centring.py, verify/verify_centring_delta.py, verify/add_verification.py ->
verify/verify_results.json, verify/reproduced_numbers.csv, verify/centring_sensitivity.csv, verify/centring_delta.csv; logs verify/*.log.
The verifier's own code does not import e6lib, e6_run, stats3 or prep_pairs. It reads the raw model npz files, the Lambert census, the TRRUST and
DoRothEA tsvs, the STRING v12 links and info files, the cache obs_genes/obs_idx/ncells/expr, and the local residual-stream copies. AUROCs use
sklearn, and LOO cosine, stratified pooling, SVD, double-edge-swap rewiring, vertex subsampling and matched subset draws are re-implemented.
Reproduced (verify/reproduced_numbers.csv; |diff| < 1e-6 unless stated):
  - Tier A size 17,693 and the 18,907 three-model intersection. The legacy list rebuilt from the 209 named genes is identical to the stored 186
    (64 TFs).
  - Tier A Lambert AUROC: ESM2-15B 0.8894, scGPT static 0.7044. Detection-stratified 0.9023 / 0.7470 (own rank-decile convention 0.9019 / 0.7463).
    TRRUST 0.8702 / 0.6840 with E6's centring.
  - Tier B lung, detection-stratified: contextual L1 0.6681, L11 0.6384; ESM2-15B 0.8799; scGPT static 0.7512.
  - Legacy186: estimator (i) TRRUST 0.6838 / 0.8744 / 0.9539 (scGPT L11 / Geneformer V2-104M / ESM2-15B); the submitted subset-SVD SV2-SV7 values
    0.650 / 0.713 / 0.837 / 0.924 exactly.
  - Tier A STRING-physical pair set identical (73,430 pairs). PPI AUROC on E6's non-pairs 0.6468 / 0.6692; on an independent rewiring 0.6497 / 0.6718.
    Own delete-d half-width (B=200) 0.0104 vs stored 0.0096-0.011 for scGPT static.
  - Tier A within-TF edge AUROC 0.4907 / 0.5746 from the stored pairs. Construction rules checked: every positive is in TRRUST; no candidate non-edge
    is in TRRUST or any DoRothEA level (either direction); every non-edge target is a target of another TF; >= 3 positives per TF; no duplicate pairs.
  - Selection bias (own 2,000 matched draws, different seed): scGPT static legacy 0.8334, null mean 0.737 (stored 0.739), percentile 0.995;
    ESM2-15B null mean 0.862, percentile 1.000.
  - Paired ESM2-15B minus scGPT static TRRUST delta (own bootstrap): +0.186 [0.160, 0.212] vs stored [0.159, 0.212].
Audit checklist:
  (1) E6 has no H1-H4 primary analysis. Its detection-stratified Lambert statistic is the H2 statistic, but on V intersected with Tier A and centred
      over that set rather than over V (lung L1 0.668 vs 0.664 on V). This is justified by the E6 spec and documented above.
  (2) Zero rows: static all-zero rows are dropped (seq_features 181); contextual vectors use V only (asserted non-zero). The static SVD is
      mean-centred over the model's zero-free vocabulary; the contextual SVD is stats3.svd_observed on V.
  (3) Alignment of obs_genes / obs_idx / pair indices: verified by exact reproduction of the lung contextual values and of the legacy set.
  (4) SV indexing: band SVa-SVb = columns a-1..b-1 of U*S; SV2-SV7 = columns 1..6 (matches the submitted slice(1, 7)).
  (5) Leakage: there are no supervised probes in E6. Best band and best layer are max-T adjusted, but their values are post-hoc maxima.
  (6) Nulls: node label permutation keeps the TF count; within-decile permutation keeps per-decile counts; within-TF permutation keeps per-TF counts;
      the PPI band null uses further degree-preserving rewirings; the selection null keeps size, label ratio and expression cell. All hold.
  (7) Bootstrap units: genes (node), TF x target two-way (edge), vertex subsampling with delete-d also reported (PPI). The PPI delta flag now also
      exists for delete-d (ci_excludes_0_deleted).
  (8) Resolution: node/edge N=2,000 (floor 0.0005), max-T N=1,000 (floor 0.001), selection N=10,000. No claim is finer than its floor.
  (9) SUMMARY numbers were checked against the CSVs. Corrections: estimator (iii) wording (section 7); seq-descriptor delta range (section 5);
      edge static-vs-contextual range and the omitted lung_external (section 5); Lambert order of ESM2-650M vs Geneformer V2-104M (section 4);
      legacy178 percentile 0.988 and excess ranges stated per variant (section 2); PPI range "0.60-0.79" -> 0.54-0.79 (inference). Additions: centring sensitivity of the Tier A TRRUST task
      (sections 2, 4 and inference); the kidney_A TRRUST ceiling (section 1); the unstratified nature of estimator (iii) (section 5); rewiring-draw
      variation (limitations).
  (10) Labels: ESM2-15B (UCE per-protein) everywhere except legacy_reproduction.csv's "model" column, which keeps the submitted name "ESM2-3B" next to
      label_corrected. Geneformer V2-316M / V1-10M / V2-104M, scGPT whole-human, AIDO.Cell-100M, scFoundation and tGPT agree with the subproject PLANs.
Output changes made by the verifier: e6_supp.py aggregate adds the column ci_excludes_0_deleted to pairwise_deltas_vs_scgpt.csv. e6_supp.py aggregate
and e6_assemble.py were re-run; every other CSV is byte-identical to its pre-verification copy (verify/_pre_verification_outputs/).
results.json gained a "verification" block (verify/add_verification.py). No bug changed any stored estimate.
```
