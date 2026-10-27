# Round-3 analysis convention and pre-registered claims (v2, frozen)

Frozen before any primary (H1–H4) analysis was run. `FROZEN.txt` records this file's SHA-256; every
`results.json` written by a primary experiment records the hash it was run under. Anything not listed in
§4 as H1–H4 is a **sensitivity analysis**: reported as an effect size with an interval, with no pass/fail.

Code: `contexts.py` (data), `stats3.py` (estimators and nulls; `python stats3.py` runs both self-test
suites), `prep_cache.py` → `cache/<context>/`. Python:
`/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python` (verified to reproduce round-2
values bit-for-bit). `TMPDIR=/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2`; BLAS threads ≤ 2 per job.
Local embedding copies: `~/biomechinterp-local/subproject41_r3/emb/<context>.npy` (identical bytes).

## 1. Contexts and their roles

Each context is one scGPT residual-stream extraction: a sample of cells, per-layer mean-pooled gene
embeddings, 12 layers × G genes × 512.

| tier | context | source dataset | donors | screen usage |
|---|---|---|---|---|
| **primary** | `lung` | Tabula Sapiens lung | 4 (shared with immune) | not used by the 183; explored after iteration 63 |
| **replication** | `lung_external` | Krasnow lung atlas | 3, fully independent study | not used by the 183; explored after iteration 63 |
| second tissue | `kidney_B` | Tabula Sapiens kidney | 1 (TSP2) | barely used |
| continuity / discovery | `immune` | Tabula Sapiens immune | multi-donor | **discovery**: where the screen's SV2–4 / SV5–7 / SV2–7 bands were found |
| continuity / discovery | `kidney_A` | Tabula Sapiens kidney, same donor as kidney_B, different HVG axis | 1 (TSP2) | **discovery**: the matrix of the submitted version and of screen iterations ~1–38 |
| technical robustness | `immune_seed43/44`, `immune_TSP14/21/25`, `immune_Bcell/CD4T/CD8T` | TS immune subsets | — | partly used |
| model robustness | `kidney_ckpt_{whole_human,kidney,brain}` | TS kidney, same cells, three scGPT checkpoints | 1 | 1 iteration |

The contexts are separate extractions of four datasets (TS immune, TS kidney, TS lung, Krasnow lung) —
not independent samples: kidney_A/kidney_B share one donor; TS immune and TS lung share donors. **No context
was untouched by the screen**; lung and lung_external were not used to generate the 183 audited hypotheses.

**Replication criterion** (pre-declared): the primary statistic in `lung_external` has the same sign as in
`lung` and its interval excludes the null value. Counting "significant contexts" is not used.

## 2. Genes, geometry, labels

- **Observed universe V**: genes with a non-zero embedding (≥ 1 sampled cell's input contained the gene).
  Never-observed genes are all-zero rows and are excluded everywhere (the zero-row SVD artefact: with them
  included, SV1 is the observed/unobserved indicator, r = −1.00). Sensitivity: V restricted to n_obs ≥ 3, 5, 10.
- **Embedding**: each layer's V × 512 matrix, mean-centred over V. Full-dimensional statistics use cosine
  similarity of these centred vectors.
- **SVD**: thin SVD of the centred V × 512 matrix (`stats3.svd_observed`). "SVk" always means the k-th
  direction of this observed-gene decomposition, of that context and layer. Directions are compared across
  layers or contexts only by principal angles / |cos| on matched genes, never by index. Top-k truncations use
  U_k S_k. Adjacent singular values with gap ratio S_k/S_{k+1} whose bootstrap CI contains 1 are a 2-D subspace.
- **Detection** n_obs = number of sampled cells whose input contained the gene (`cache/<ctx>/ncells.npy`,
  verified against the reproduced extraction by A2).
- **Node labels (primary)**: Lambert et al. (2018) TF census (`humantfs_v1.01.txt`, "Is TF?" = Yes) vs every
  other gene in V. **Secondary / legacy**: TRRUST TF vs TRRUST target-only.

## 3. Reference networks

- **TRRUST v2**: distinct (TF, target) pairs (rows deduplicated; TRRUST lists one row per PMID).
- **STRING v12 physical subnetwork** (`9606.protein.physical.links.v12.0`, verified from the fetch script;
  1,477,610 lines = 738,805 pairs × 2 directions), combined score ≥ 700, symbol-mapped, undirected. Chosen as
  primary because the functional network's co-expression and text-mining channels are circular for an
  expression-trained model. The submitted version used functional-network pairs; a kidney_A bridge row
  quantifies the switch. Sensitivity: ComplexPortal and hu.MAP 2.0 co-complex pairs.
- **DoRothEA**: its curated tier integrates TRRUST, so it is never used as independent positives. Independent
  positives = ChIP-seq-file level-D pairs (binding evidence without curated support) not in TRRUST and not in
  the main file's A–C. Every DoRothEA level from both files is excluded from every candidate non-edge set.

## 4. Pre-registered confirmatory claims

| | **H1 — network neighbourhood** | **H2 — node identity** | **H3 — edge relation** | **H4 — rank vs full** |
|---|---|---|---|---|
| question | do spectral coordinates place physical interaction partners together beyond degree? | is TF identity linearly recoverable beyond detection frequency? | is a specific TF→target relation recoverable once TF identity and target identity are controlled? | does truncating to the top-k spectral directions improve recovery over the full representation? |
| primary context | lung (replication: lung_external) | lung (replication: lung_external) | lung (replication: lung_external) | k* chosen on `immune_seed43`; Δ evaluated in lung, lung_external, kidney_B, immune |
| universe | V; STRING-physical ≥ 700 edges within V | V; Lambert TF vs non-TF | TRRUST edges within V; TFs with ≥ 3 positives and ≥ 3 candidate non-edges | identical genes and pairs as H1–H3 |
| primary statistic | scalar assortativity r of the SVk coordinate over edges (K-free) | detection-stratified AUROC (strata = deciles of log n_obs) of the leave-one-out mean cosine to the TF set, full dimension | within-TF stratified AUROC of signed cosine (full dimension) | Δ = metric(k*) − metric(512) per task; metrics: STRING-pair cosine AUROC vs degree-preserving rewired non-pairs; H2 statistic; H3 statistic |
| null | Maslov–Sneppen degree-preserving rewiring of the V-induced STRING subgraph (burn-in 10m, thinning m accepted swaps), N = 10,000 joint draws in lung / lung_external, 2,000 elsewhere | label permutation within detection strata, joint across layers, N = 10,000 | within-TF label permutation, joint across layers, N = 10,000 | paired bootstrap of Δ (gene-level for H2, two-way TF × target for H3, vertex subsampling for H1), B = 2,000 |
| family / correction | SV1–SV10 × 12 layers = 120; Westfall–Young max-T | 12 layers; max-T | 12 layers; max-T for detection; intersection–union equivalence for absence | 3 tasks × 4 contexts; reported with CIs (no test) |
| decision rule | supported if max-T p < 0.05 in lung and in lung_external | supported if max-T p < 0.05 in lung and in lung_external | **detected** if max-T p < 0.05; **AUROC ≥ 0.60 excluded** if every layer's one-sided 95% upper bound (two-way cluster bootstrap) < 0.60; **not above co-expression** if the upper bound of Δ vs \|co-expression\| < 0.05 at every layer; otherwise **inconclusive** | one sentence: "rank reduction improves / does not improve X by Δ [CI] in n of 4 contexts" |
| smallest effect of interest | — | — | AUROC 0.60; Δ vs \|co-expression\| 0.05 | — |

Candidate non-edges for H3 ("TRRUST-negative candidate pairs", never called negatives without qualification):
(TF, g) where g is a TRRUST target of at least one *other* TF in V; (TF, g) is absent from TRRUST, DoRothEA
(any level, both files); (g, TF) is absent from every regulatory set; up to 10 per positive, each matched to its
positive on target TRRUST in-degree decile and target detection decile. Within-TF stratification removes TF
identity; target matching removes "is a commonly documented target". Symmetric scorers (cosine,
\|co-expression\|) cannot carry direction; direction is tested separately with directed scorers.

## 5. Sensitivity analyses (effect sizes with intervals; no pass/fail)

E1 compression on V (all contexts; zero-row values only as a documented artefact) · E2 extras: label-shuffle,
degree-stratified (within-universe and genome-wide degree) and detection-stratified nulls; co-expression-matched
non-pairs; co-pole at K/n ∈ {0.10, 0.20, 0.27, 0.33, 0.50} as observed/null ratios; hub exclusion; ComplexPortal /
hu.MAP 2.0 · E3 extras: 78 contiguous bands × 12 layers (max-T adjusted, N = 2,000), single directions
(sign-invariant), cross-fitted gene-half band estimate, legacy bands mapped by `E9_screen/legacy_band_mapping.json`,
TRRUST TF vs target-only, baselines (detection, mean expression, scGPT static token embedding, ESM2) ·
E4 extras: pooled AUROC with all-V non-edges (legacy only), macro per-TF AUROC, |cos|, rank-k and bands,
GENIE3 / GRNBoost2, DoRothEA-ChIP independent positives, strata (PMID count, activation / repression, TF
detection ≥ 10 %), label-noise attenuation (π = 0.05, 0.1, 0.2), a priori MDE, direction test, cell-type
contexts, supervised probes (C1 random pairs, C2 held-out TFs, C3 pair-disjoint; rewired-label null) ·
E5 conditional logistic confound adjustment · E6 cross-model · E7 full curves · E8 Word2Vec ·
E10 cross-context variance decomposition and diagnosis · E11 compartment re-test · checkpoints.

**Retrospective family**: the screen's 183 hypotheses (iterations 1–63, the set reported in the preprint),
Bonferroni 0.05/183 and BH; sensitivity over all records.

## 6. Inference rules

- Empirical p = (b+1)/(N+1). `stats3.assert_resolution` is called for every test that carries a significance
  claim. Stouffer combination is disabled (layers share genes; contexts share donors).
- Intervals for edge statistics: two-way (TF × target) pigeonhole cluster bootstrap, B = 2,000 (TF-only as a
  check); for network statistics: vertex subsampling without replacement (m = n/2, √(m/n) rescaling); for node
  statistics: gene bootstrap.
- Comparisons between contexts or models are tested as paired differences on matched genes and labels, never as
  "significant here, not there".

## 7. Outputs

Each experiment writes `revision3_experiments/<ID>_<name>/`: script(s), log, CSVs, `results.json`
(including the CONVENTION3 hash) and `SUMMARY.md` (what was run, numbers with file provenance,
interpretation, caveats). Null distributions are stored as summaries (b, N, mean, sd, quantiles).
