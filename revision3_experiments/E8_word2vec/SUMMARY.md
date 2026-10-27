# E8 — Word2Vec co-occurrence instrument (verified)

Copy of the verifier-corrected SUMMARY.txt (verify/ holds the checks).

```
E8 - Word2Vec co-occurrence-only instrument on H1/H2/H3 (sensitivity; CONVENTION3 §5 "E8 Word2Vec"; no pass/fail)
CONVENTION3 sha256 3e3249d1d0a8ce1bc3cf64639fe233b016535456f82432a8bbfcacdd64a05c01 (matches FROZEN.txt; results.json hash_match = true)
Answers R2.1 (the word-embedding model from the Introduction is now used) and feeds the instrument comparison.
Status: independently verified (verify/; see VERIFICATION at the end). No bug that changes a result was found; 12 wording
corrections were applied to the implementer's text (verify/v3_wording_corrections.csv). The implementer's original text is kept in
verify/SUMMARY_implementer_before_verification.txt.

QUESTION
Does an instrument that sees only which genes share a cell (the scGPT token sets, order shuffled, no expression values, no attention)
already recover the three signals the scGPT residual stream carries? If it does, those signals in scGPT need no transformer to explain
them. If it does not, the transformer (or its pretraining) adds something. scGPT = whole-human checkpoint (the extraction script's
default), pretrained on ~33M cells.

DESIGN (what was run)
- Corpora (e8_corpus.py): every cell of each context's extraction h5ad (lung 15,000; lung_external 9,000; kidney_B 8,000; immune 20,000;
  kidney_A 11,376 cells), tokenised with A2's rule (non-zero genes, top 512 by expression). Check: per-gene token counts over the
  A2-reproduced extraction cells equal layer_gene_embedding_count.npy exactly in all 5 contexts (corpus/<ctx>_meta.json).
- gensim 4.4.0 Word2Vec, sg=1, vector_size=256, min_count=1, negative=10, workers=2, 10 epochs, 5 seeds; token order and cell order
  re-shuffled on every pass; other parameters gensim defaults (recorded in models/*.json; sample = 1e-3 frequent-token downsampling
  affects only immune: 20 V genes, 0.3% of V tokens, minimum keep probability 0.88; verify/v1_inputs.json). "full" = every cell;
  "extr" = the 192/256 extraction cells only (lower bound). Every V gene is in every model's vocabulary (0 missing).
- Evaluation (e8_eval.py), on V (n = 1,972-5,648 genes), every representation mean-centred over V: H1a assortativity of SV1-10
  (STRING-physical >= 700 edges within V), Maslov-Sneppen null N = 1,000, max-T over 10 SVs (over 120 SV x layer tests for scGPT);
  H1b basis-free pair-cosine AUROC vs one MS-rewired non-pair set, vertex subsampling m = n/2, B = 2,000; H2 detection-stratified LOO TF
  AUROC (Lambert), within-decile permutation N = 2,000, gene bootstrap B = 2,000 with scores held fixed; H3 within-TF stratified AUROC
  on h3_pairs, two-way cluster bootstrap B = 2,000. Paired deltas vs scGPT layers 0-11 on the same genes, pairs and resampling draws.
- Every vertex-subsampling interval (H1b) is given both ways: the pre-registered sqrt(m/n) interval (ci_lo/ci_hi) and the delete-d
  interval (*_deleted, half-widths x sqrt 2; CONVENTION3_ERRATA item 3). Below, H1b intervals are delete-d.
- H2 note: the LOO estimator is not centred at 0.5 under the within-stratum permutation null (Word2Vec null means 0.27-0.50; scGPT
  0.475-0.496). Reason (inference, consistent with the data): a TF's LOO score omits one same-stratum TF, so when the geometry
  is dominated by detection (as for under-trained Word2Vec, where rarely seen genes stay near their shared initial direction), TFs
  score systematically below same-stratum non-TFs. Raw AUROC < 0.5 therefore does not mean anti-signal; significance is judged only
  against the permutation null. W2V-vs-scGPT H2 deltas are quoted null-centred (delta_nullcentred*: each representation's own null mean
  subtracted); raw deltas are in the CSV too.
- This is a sensitivity analysis. The scGPT H1a/H2 p-values here (N = 1,000 / 2,000) are not the primary H1/H2 tests (N = 10,000)
  and must not be quoted as such.

DEVIATIONS / REDUCTIONS (stated exactly)
1. Window. The spec says window = 512. Measured cost: 1,912 s for 10 epochs on 256 lung cells (models/lung__extr__w512__e10__s0.json),
   i.e. ~31 h per seed on the full lung corpus (x 5 seeds x 5 contexts). Not feasible on this machine. Production models use window = 5
   (gensim default) with per-cell, per-pass shuffling, so each pass sees a random subset of each token's within-cell partners. This
   weights each within-cell pair by ~1/(L-1) (L = tokens in the cell) instead of uniformly; with L = 512 in most lung/kidney_B/
   lung_external cells this matters little, more in immune (median 228 tokens). The spec's alternative ("token order shuffled per cell
   per epoch") is met. The spec's window is met only in a calibration (extraction cells, seed 0, lung and immune): w512 x 10 epochs vs
   w5 x 10 epochs vs w5 with the pair-update count matched to w512 (571 / 365 epochs).
2. Training budget. Full-corpus w5 x 10-epoch models get about as many pair updates as the w512 calibration gets on 256 cells (lung
   4.5e8 vs 4.5e8). That is ~1/58 of what the specified w512 full-corpus model would get. Saturation check (e8_train_sat.py): seed 0,
   40 epochs, lung and immune (saturation.csv).
3. Not reduced: N_MS = 1,000 and N_perm = 2,000 as the E8 spec says (the primary H1/H2 use more). assert_resolution passes for every
   significance claim (floor 1/1001 < 0.05/10; 1/2001 < 0.05).
4. gensim workers = 2 is not bit-reproducible; the seeds fix initialisation and the corpus stream only. gensim/Cython 3 prints
   "Exception ignored in our_dot_float" when a float dot product is exactly -1.0 (the sentinel BLAS sdot uses for errors); that one dot
   product then counts as 0. 23-84 occurrences per training-process log (each log covers 1-30 models), against roughly 1e9 dot products
   per model. Benign.
5. Fixes to the previous agent's code (by the implementer): e8_merge.py read macOS "._*.json" AppleDouble files on the exFAT drive and
   crashed (now excluded); delete-d H1b columns and the saturation table added. No other code changed. Adding the e40 models and
   re-running eval for lung and immune reproduced every earlier row exactly (max abs diff 0.0; reproduction_check.csv; earlier outputs
   in eval_v1/).

DIRECT EVIDENCE (5-seed means of the full-corpus models unless stated; seed SD <= 0.008 for every H1b/H2/H3 statistic; <= 0.016 for
the extraction-only models)
H1b STRING pair-cosine AUROC (paired_vs_scgpt.csv, w2v_results.csv)
  W2V: lung 0.605 [0.589, 0.621], lung_external 0.597 [0.581, 0.614], kidney_B 0.605 [0.586, 0.622], immune 0.704 [0.656, 0.751],
  kidney_A 0.739 [0.690, 0.794].
  scGPT layer range: 0.619-0.680, 0.630-0.677, 0.621-0.675, 0.668-0.743, 0.749-0.814 (best at L1/L2).
  Delta W2V - scGPT(best layer): lung -0.075 [-0.090, -0.059], lung_external -0.080 [-0.097, -0.062], kidney_B -0.071 [-0.089, -0.052],
  immune -0.040 [-0.082, +0.005], kidney_A -0.074 [-0.107, -0.042].
  scGPT is higher (delta CI below 0) at 12/12 layers in lung, lung_external and kidney_B, at 8/12 in kidney_A, and at 0/12 in immune
  (4/12 with the pre-registered narrower intervals). The gap shrinks with depth (lung: -0.074 at L0 to -0.014 at L11). Delta at L11:
  lung -0.014 [-0.027, -0.00003] (upper bound touches 0), lung_external -0.032 [-0.049, -0.015], kidney_B -0.016 [-0.029, -0.002],
  kidney_A -0.015 [-0.037, +0.008], immune +0.036 [-0.005, +0.078]. In immune, W2V > scGPT at L10-L11 only under the pre-registered
  sqrt(m/n) intervals, not under delete-d.
H2 detection-stratified TF AUROC
  W2V: lung 0.632, lung_external 0.602, kidney_B 0.601 (each permutation p = 0.0005, the floor; null mean ~0.49). immune 0.512
  (p = 0.22) and kidney_A 0.469 (p = 0.28; null mean 0.444): no detectable TF signal.
  scGPT range: lung 0.625-0.665, lung_external 0.627-0.662, kidney_B 0.615-0.654, immune 0.566-0.720, kidney_A 0.678-0.783.
  Null-centred delta vs best layer: lung -0.033 [-0.060, -0.007], lung_external -0.061 [-0.082, -0.039], kidney_B -0.051
  [-0.075, -0.025], immune -0.204 [-0.281, -0.120], kidney_A -0.284 [-0.354, -0.208].
  scGPT is higher at 7, 12, 8, 10 and 12 of 12 layers (null-centred; lung, lung_external, kidney_B, immune, kidney_A). W2V is never
  higher at 10 epochs. The scGPT advantage is concentrated in early/middle layers: the null-centred CI covers 0 at L7-L11 in lung,
  L8-L11 in kidney_B and L10-L11 in immune; never in lung_external or kidney_A.
H3 within-TF stratified AUROC
  W2V: lung 0.512 [0.481, 0.549], lung_external 0.534 [0.498, 0.570], kidney_B 0.528 [0.476, 0.582], immune 0.603 [0.533, 0.696],
  kidney_A 0.577 [0.470, 0.695].
  scGPT range: 0.516-0.542, 0.518-0.538, 0.524-0.555, 0.577-0.638, 0.601-0.652.
  The paired-delta CI covers 0 at all 60 context x layer comparisons (best-layer deltas -0.005 to -0.076).
H1a assortativity of SV1-10
  W2V full: in lung, lung_external, kidney_B and immune, all 10 SVs have max-T p < 0.05 in every seed; in kidney_A, 9/10. Max mean r:
  0.274, 0.295, 0.265, 0.620, 0.693. scGPT: 120/120 SV x layer tests at max-T p < 0.05 in four contexts, 118/120 in immune (L7 SV1
  and L8 SV1 not); max r 0.435, 0.407, 0.404, 0.566, 0.702 (h1_assortativity_scgpt.csv). The two instruments' SV indices do not
  correspond, so the H1a comparison is descriptive only.
Extraction-cells-only lower bound (5 seeds, w5 x 10 epochs): lower than the full-corpus model on every task in every context (e.g.
  H1b lung 0.549; H2 lung 0.487, p = 0.66). Its H2 permutation p is nonetheless < 0.05 in lung_external (0.539, p = 0.0005) and in
  kidney_A (0.361, p = 0.026, null mean 0.274: below 0.5 yet above its own null). The calibration shows these models are under-trained:
  CKA with the w512 model is 0.26 (lung) and 0.40 (immune); for the matched-update w5 model it is 0.80 and 0.83. Matched-update / w512
  extraction-only models (seed 0): lung H1b 0.571 / 0.571, H2 0.611 / 0.594 (p = 0.0005 both), H3 0.516 / 0.520; immune H1b
  0.689 / 0.656, H2 0.526 / 0.492 (p = 0.12 / 0.43), H3 0.592 / 0.569 (calibration.csv).
Saturation (full corpus, seed 0, 10 -> 40 epochs; saturation.csv; CKA e10 vs e40 = 0.944 lung, 0.945 immune)
  lung: H1b 0.606 -> 0.612, H2 0.630 -> 0.664, H3 0.512 -> 0.511. immune: H1b 0.705 -> 0.702, H2 0.520 -> 0.515, H3 0.608 -> 0.614.
  At 40 epochs, lung H2 cannot be distinguished from scGPT at L0-L10 (null-centred delta CI covers 0; vs best layer L5: +0.000
  [-0.024, +0.022]) and W2V is higher than L11 (+0.039 [+0.014, +0.063]). Single seed.
Representation similarity (cka.csv): W2V-full vs scGPT linear CKA is 0.17-0.58 (context means 0.22-0.42). Seed pairs: >= 0.968 for
  full-corpus models, >= 0.963 for extraction-only models.

INFERENCE
- The H1-type signals (STRING neighbours lie close; SVs are assortative on STRING) do not need a transformer. A model that sees only
  within-cell co-occurrence of the same tokens yields STRING pair AUROCs of 0.60-0.74 and assortative SVs; even on the extraction cells
  alone, with enough training (matched-update / w512), it gives 0.57 (lung) and 0.66-0.69 (immune). scGPT is better by 0.071-0.080
  AUROC at its best (early) layers in 4/5 contexts. At its last layer scGPT exceeds co-occurrence by 0.014-0.032 in those four contexts;
  that gap excludes 0 in lung_external and kidney_B, only touches it in lung, and does not exclude 0 in kidney_A. In immune, the two
  cannot be told apart at any layer.
- H2 (TF identity) splits by context. In the three larger-V contexts (lung, lung_external, kidney_B), co-occurrence alone gives a clear
  TF signal (0.60-0.63), and in lung, with 4x the training (seed 0, 40 epochs), it cannot be distinguished from scGPT at L0-L10 and
  exceeds L11. In immune and kidney_A, scGPT's strong TF signal (up to 0.72 / 0.78) has no counterpart in this co-occurrence instrument
  (W2V at null at 10 epochs, at 40 epochs in immune, and in the immune w512 / matched-update calibration models). That part of H2 is
  what this co-occurrence instrument does not reproduce.
- H3 (edge relation): both instruments are weak (<= 0.65), and they cannot be told apart anywhere. This experiment gives no evidence
  that scGPT encodes TF-target relations beyond what co-occurrence gives. Equivalence is not claimed: the delta CIs are +/- 0.04-0.10.

HYPOTHESIS (not tested here)
- The immune / kidney_A H2 gap might reflect that those observed universes are small (V = 2,039 / 1,972 genes from 4,941 / 4,803-gene
  HVG axes) and that few detected genes are lost to truncation (A2: 99.6% / 96.5% of detected entries tokenised), so co-occurrence
  there is dominated by cell-type programmes rather than TF expression patterns. Or scGPT's TF signal there might come from pretraining
  across ~33M cells, not from these cells. Neither is tested.
- The late-layer convergence of scGPT to co-occurrence-level H1b might mean that late layers keep mostly context-co-occurrence
  information. This is untested. CKA (0.2-0.6) shows the geometries are far from identical anyway.

WHAT THIS DOES NOT SHOW / LIMITATIONS
- It is not the specified window-512 full-corpus model. The full-corpus models are window 5 and get ~1/58 of the specified pair
  updates (lung). The saturation check covers 2 of 5 contexts and 1 seed. H1b, H3 and immune H2 had plateaued, but lung H2 was still
  rising (+0.034 from 10 to 40 epochs). So the 10-epoch H2 deltas in lung, lung_external and kidney_B may overstate scGPT's advantage.
  Only lung was checked.
- The comparison is not data-matched. W2V-full sees 8,000-20,000 cells of the context. scGPT's embeddings come from 192/256 cells, but
  its weights come from ~33M pretraining cells. The extraction-cells-only W2V is the other bound, and it depends strongly on training
  budget (calibration).
- Single-seed calibration and saturation. Seed SD is small (<= 0.008 for full-corpus models, <= 0.016 for extraction-only).
- All intervals condition on the representation; seed variation is reported separately (sd_seeds), not propagated into the paired CIs.
  The H2 gene bootstrap holds the LOO scores fixed, so variability of the TF centroid itself is not propagated either.
- The H2 null-centring is a descriptive adjustment (a constant shift by each representation's permutation-null mean); it is not part of
  the pre-registered statistic.
- H1a uses N = 1,000 MS draws, so the max-T p floor is 0.001. The SV comparison is not index-matched.
- kidney_A and kidney_B share one donor's cells (errata item 2). They are not independent replications.
- Token sets at the 512 cut-off: many cells have ties at the 512th value (lung ~85% of checked cells). The corpus uses the same numpy
  argsort tie-break as A2's reproduction of the extractor; exact count equality supports, but does not strictly prove, identical
  per-cell token sets (A2 caveat).
- Surprises: (i) the extraction-only w5 x 10-epoch models are badly under-trained, so the "lower bound" must be read from the
  matched-update models; (ii) in immune, co-occurrence H1b cannot be distinguished from scGPT at any layer (delete-d); (iii) the size of
  the LOO-null shift for Word2Vec in kidney_A (null mean 0.27 for the extraction models).

VERIFICATION (independent adversarial verifier; scripts and outputs in verify/)
- Code audit: e8_corpus.py, e8_train.py, e8_train_sat.py, e8_eval.py, e8_merge.py, e8_driver.sh read in full. Statistics use stats3
  (svd_observed on V, V-centred; assortativity; maslov_sneppen_chain burn 10m / thin m; vertex_subsample_ci m = n/2; tf_scores_loo +
  stratified_auroc; permute_within_strata; StratifiedPairs two-way bootstrap; max_t_adjust; assert_resolution). The two speed-only
  re-implementations (IndexedStratifiedPairs.auroc, batched H2 null) are asserted equal to stats3 inside e8_eval.py. Paired deltas use
  identical genes, pairs and resampling draws for every representation. SV k = column k-1 of US (sign-invariant statistic). No
  never-observed rows. Gene alignment: model row r = V gene obs_idx[r] = h5ad var index = embedding row (A2). No supervised probes,
  so no fold leakage is possible. No bug that changes a result was found.
- Inputs (verify/v1_inputs.py, own anndata tokenisation, not a2_lib): token sets equal the corpus in all 2,716 checked cells (all
  extraction cells + 300 random others per context); extraction counts equal layer_gene_embedding_count.npy in all 5 contexts; local
  embedding copies equal the source files on V rows (layers 0 and 11, SHA-256).
- Recomputed from raw inputs with own code (verify/v2_stats.py; verify/v2_stats.csv, mine vs stored): exactly reproduced (|diff|
  < 1e-5): H1b W2V seed-mean AUROC lung 0.605 and immune 0.704, scGPT lung L1 0.680 / immune L2 0.743 / immune L11 0.668 and the
  three deltas; H2 W2V seed-mean lung 0.632, immune 0.511, kidney_A 0.469 and scGPT best layers lung L5 0.665, immune L2 0.720,
  kidney_A L3 0.783; H3 W2V seed-mean kidney_A 0.577 and immune 0.603, scGPT kidney_A L1 0.652, immune L0 0.638; H1a r for SV1-10 of
  lung and immune W2V seed 0 and scGPT L2 SV2 in both; calibration CKA lung w512 vs w5_e10 0.262. Monte-Carlo quantities agree within
  resampling error using my own draws and seeds: H1b delete-d delta CIs (B = 400) lung L1 [-0.091, -0.061] vs stored [-0.090, -0.059],
  immune L2 [-0.085, +0.006] vs [-0.082, +0.005], immune L11 [-0.001, +0.075] vs [-0.005, +0.078]; H2 permutation null means within
  0.004 and p-values (immune s0 0.149 vs 0.156; kidney_A s0 0.294 vs 0.327; lung at floor); H3 two-way bootstrap delta CIs (B = 300)
  kidney_A L1 [-0.170, +0.015] vs [-0.174, +0.021], immune L0 [-0.099, +0.035] vs [-0.100, +0.033]; an independent networkx
  double-edge-swap MS null (300 draws) for immune W2V s0 gives null means -0.004 to -0.006 and SDs 0.027-0.030 vs stored -0.001 to
  -0.005 and 0.028-0.029, with observed z = 7.6-21.2 against a 95% max-null z of 2.4. The stored rewired non-pair sets contain no
  STRING edge and no self-pair.
- Number check: every number in this summary was checked against w2v_results.csv, paired_vs_scgpt.csv, h1_assortativity_scgpt.csv,
  calibration.csv, saturation.csv, cka.csv, models/*.json and the training logs. Twelve wording corrections were made
  (verify/v3_wording_corrections.csv); none changes a conclusion's direction, but (6) and (7) weaken two statements: the lung L11 H1b
  gap only touches 0, and 40-epoch lung H2 is "not distinguishable from" rather than "equal to" scGPT.
- Labels: no ESM2 model is used in E8. scGPT embeddings are from the whole-human checkpoint.

FILES
w2v_results.csv (per seed and seed-mean, every task); paired_vs_scgpt.csv (per layer deltas, both H1b scalings, null-centred H2);
h1_assortativity_scgpt.csv; calibration.csv; saturation.csv; cka.csv; reproduction_check.csv; results.json (hash, spec, deviations,
per-model metadata, headline, null summaries via eval/<ctx>.json); eval/<ctx>.json (null summaries for every representation);
corpus/*_meta.json; models/*.json. Scripts: e8_corpus.py, e8_train.py, e8_train_sat.py, e8_eval.py, e8_merge.py, e8_driver.sh.
Verification: verify/v1_inputs.py (+ .json), verify/v2_stats.py (+ v2_stats.csv, v2_extra.json, v2_stats.log),
verify/v3_wording_corrections.csv, verify/SUMMARY_implementer_before_verification.txt.
```
