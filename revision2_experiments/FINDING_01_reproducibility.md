# Finding 1 — the headline numbers do not reproduce from any deposited artefact

Established by `convention_matrix.py` (output: `convention_matrix_results.txt`), which sweeps the
16 combinations of the conventions actually present in the deposited material:

* embedding matrix: `cycle1_main` (12, 4803, 512) vs `cycle4_immune_main` (12, 4941, 512)
* zero-norm ("out of vocabulary") filter applied at L0: on / off
* SVD computed on the curated SUBmatrix (loop convention, `iter_0010`) vs on the FULL
  gene matrix and then restricted (revision convention, `e3_full_vocab/run.py`)
* TRRUST pairs deduplicated / not (TRRUST lists one row per supporting PMID, so the same
  TF-target pair appears several times)

## What the manuscript claims vs what any convention yields

| statistic | manuscript | best reproduction | convention that gives it |
|---|---|---|---|
| curated gene set | **195** | 209 (unfiltered) or 187 (L0 nonzero) | none gives 195 |
| TRRUST pairs | **326** | **326 exactly** | cycle1_main, no nz filter, **not deduplicated** |
| STRING >= 0.7 pairs | **1,022** | 968 | none gives 1,022 |
| TF count | **51** | 69 (209 genes) / 65 (187 genes) | none gives 51 |
| SV2 TRRUST co-pole | **0.206** | **0.206 exactly** | cycle1_main, no nz, subset SVD, not deduplicated |
| SV1 extracellular OR | **6.37** | 3.39 | none comes close |
| TF/target AUROC (12-layer mean) | **0.744** | 0.717 (a *different* convention) | none reproduces it with the co-pole convention |

## The three problems this exposes

1. **The headline numbers come from mutually inconsistent conventions.** The row that reproduces
   the co-pole rate exactly (0.206) simultaneously gives OR = 3.39, not 6.37, and AUROC = 0.583, not
   0.744. No single row of the sweep produces more than one headline number.

2. **326 is a count of TRRUST *rows*, not distinct pairs.** Deduplicating gives 288. Using rows means
   a TF-target pair supported by five PMIDs contributes five times to the co-pole rate, so the
   statistic is implicitly weighted by citation count — a literature-bias confound that is not
   disclosed and not controlled.

3. **Methods 4.4 misdescribes the null.** It states a "gene-label shuffle null with N = 200-1,000
   permutations". `iter_0010/H02`, the source of 0.206 vs 0.122, states: *"Null: N=1000 random
   same-size gene pairs sampled from named genes"* — a pair-resampling null, which tests a different
   hypothesis (is this pair set special?) than a gene-label shuffle (is this axis special?).

## Additional estimator defect, independent of the above

`tf_target_auroc_on_subset` scores every gene by its mean cosine similarity to the TF set, and each
TF's own self-similarity (= 1) is included in its mean. This adds 1/|TF| to every TF score and
nothing to any target score — a systematic upward bias. Measured here:

* 209-gene curated set: **+0.021 AUROC**
* 4,941-gene matrix, 295-gene set: **+0.016 AUROC**
* 187-gene L0-filtered set: **+0.029 AUROC**

The leave-one-out estimator (`common.tf_target_scores(exclude_self=True)`) removes it.

## Consequence for the plan

The blocking problem is no longer "re-run the permutations at higher N". It is that **every headline
statistic has to be recomputed from the embeddings under one declared convention**, and the resulting
effect sizes are materially smaller than those currently in the manuscript. Adding controls (GENIE3
baselines, matched ESM2, all-layer confounds) to numbers that do not reproduce would be wasted work.
