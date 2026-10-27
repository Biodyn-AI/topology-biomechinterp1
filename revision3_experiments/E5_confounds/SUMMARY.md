# E5 — Confound adjustment of the edge signal (verified)

Copy of the verifier-corrected SUMMARY.txt (verify/ holds the checks).

```
E5_confounds: confound adjustment of the edge-level (H3) signal. Sensitivity analysis (CONVENTION3 §5), answers R1.6 and R1.3.
CONVENTION3 sha256 3e3249d1d0a8ce1bc3cf64639fe233b016535456f82432a8bbfcacdd64a05c01 (matches FROZEN.txt; recorded in results.json).
Every number below comes from clogit_<ctx>.csv, residual_auroc.csv, continuity_kidneyA.csv, continuity_hashseed_spread.csv,
inference_scale.csv, stratum_weights.csv or results.json in this directory, from verify/verify_summary.json (the independent
verifier's checks), or, where stated, from E4a_edges/h3_primary_<ctx>.csv and E4a_edges/sens_macro.csv.
This version incorporates the independent verifier's corrections (see VERIFICATION at the end; the implementer's original text
is in verify/SUMMARY_pre_verify.txt). No result changed. Wording was corrected in eight places, listed at the end.

QUESTION
Does the within-TF scGPT edge signal survive adjustment for pair-level and node-level nuisance covariates? The signal is the
signed cosine of V-centred full-dimensional embeddings on the H3 pairs (TRRUST positives plus matched TRRUST-negative
candidate pairs). The covariates are log|co-expression|, log mean expression, log n_obs, genome-wide STRING-physical (>=700)
degree and GO annotation count, each for TF and target.

METHOD (as specified; deviations under LIMITATIONS)
- Conditional logistic regression, stratum = TF, for each context x layer: y ~ z(score) + 5 target/pair covariates (z-scored
  over all candidate pairs; the OR is per pooled SD of the cosine over all pairs). The 4 TF-level covariates are constant
  within a stratum, so they cancel from the conditional likelihood. The stratification adjusts for them, and they cannot be
  estimated. TRRUST target in-degree is not used (label leakage; candidate non-edges are already matched on it).
- statsmodels ConditionalLogit (the specified engine) overflows and exceeds the recursion depth in the large strata. lung,
  lung_external and kidney_B have TFs with up to 367 positives among about 1,600 candidates. The same exact conditional
  likelihood is evaluated with a binomially normalised recursion with an exact gradient and Hessian (e5_lib.py).
  It was validated against statsmodels on 14 fits (validation_statsmodels.csv). The verifier also reproduced it independently
  (log-space DP in torch) to <= 2e-14 in beta (VERIFICATION).
- Two-way (TF x target) pigeonhole cluster bootstrap, B = 500, refit each time, the same draws for every layer. LRT of the
  score term (chi2, 1 df). Unadjusted conditional OR for comparison. Joint within-TF label permutation, N = 500, one
  permutation for all 12 layers, statistic = mean over layers of the adjusted score coefficient (assert_resolution passed,
  floor 1/501 = 0.002, also at alpha/5). The permutation keeps each TF's positive count and breaks the label's association
  with the score and the covariates together. It therefore tests the global within-TF null with the adjusted coefficient as
  statistic. It is not a conditional (Freedman-Lane-type) test of the score given the covariates; with |r(score, covariate)|
  <= 0.35 the difference is expected to be small.
- AUROC scale: the score is residualised on the covariates by OLS after within-TF centring. The within-TF stratified AUROC
  of the residual is computed with the same two-way bootstrap draws; r_rb = 2*AUROC - 1.
- Sensitivity: degree and GO count entered as log1p (point estimate + LRT).
- Continuity (kidney_A): round-2 run_confounds_sv57.py re-run verbatim (variant A, zero-row SVD), under PYTHONHASHSEED=0.
  Then only the SVD was changed to the observed-gene convention (B with identity bands, and B with bands mapped by
  E9_screen/legacy_band_mapping.json). Then the paper benchmark end to end (C: H3 pairs, E5 covariates, within-TF residual
  AUROC, two-way bootstrap).

DIRECT EVIDENCE
1. Adjusted OR per SD of the score, two-way bootstrap 95% CI, range over layers (largest at L0-L2 and smallest at L10-L11 in
   every context):
   lung          1.12 [1.06,1.20] (L11) to 1.24 [1.15,1.36] (L1); unadjusted 1.13-1.25
   lung_external 1.13 [1.05,1.23] (L10) to 1.22 [1.12,1.34] (L1); unadjusted 1.14-1.23
   kidney_B      1.15 [1.07,1.25] (L11) to 1.26 [1.16,1.40] (L0); unadjusted 1.15-1.28
   immune        1.37 [1.08,1.74] (L10) to 1.65 [1.32,2.12] (L2); unadjusted 1.39-1.68
   kidney_A      1.37 [1.00,2.01] (L11) to 1.78 [1.41,2.62] (L1); unadjusted 1.43-1.86
   The bootstrap CI excludes 1 in all 60 context x layer cells. The lowest lower bound is 1.005 (kidney_A L11); the largest
   fraction of bootstrap draws with beta <= 0 is 0.024. Covariates attenuate the log-OR by 2.7-6.7% in lung, 4.0-5.0% in
   lung_external, 2.2-5.3% in kidney_B, 3.0-5.8% in immune and 6.6-11.2% in kidney_A. log1p coding of degree and GO count
   changes the OR by < 0.01.
2. Mean over layers (joint), adjusted OR [bootstrap CI] vs unadjusted:
   lung 1.189 [1.114,1.281] vs 1.198; lung_external 1.184 [1.098,1.294] vs 1.193; kidney_B 1.212 [1.121,1.332] vs 1.223;
   immune 1.539 [1.230,1.978] vs 1.569; kidney_A 1.602 [1.259,2.354] vs 1.667.
   Joint within-TF permutation: p = 0.002 in all five contexts. This is the resolution floor: 0 of 500 null draws reached the
   observed value. The observed statistic is 6.0 (kidney_A) to 12.3 (lung) null SDs above the null mean. Per-cell LRT p:
   5e-47 to 3e-14 in lung, 2e-39 to 1e-15 in lung_external, 2e-39 to 1e-14 in kidney_B, 3e-21 to 9e-8 in immune and 7e-13 to
   1.6e-4 in kidney_A.
3. Residualised within-TF AUROC (r_rb = 2*AUROC - 1), range over layers:
   lung 0.512-0.537 (r_rb 0.024-0.075); lung_external 0.515-0.533 (0.030-0.066); kidney_B 0.521-0.548 (0.042-0.096);
   immune 0.561-0.628 (0.123-0.255); kidney_A 0.579-0.634 (0.159-0.268).
   Residualisation lowers the AUROC by 0.002-0.007 in the large contexts and by 0.006-0.022 in immune and kidney_A.
   The one-sided 95% bootstrap upper bound of the residual AUROC is < 0.60 at every layer in lung (max 0.577), lung_external
   (0.579) and kidney_B (0.591), but not in immune (max 0.695) or kidney_A (0.716). The two-sided residual CI excludes 0.5 at
   0/12 layers in lung and lung_external and 2/12 in kidney_B (L6, lower bound 0.503, and L5, lower bound 0.500007, i.e. on
   the boundary). It excludes 0.5 at 10/12 layers in immune and 8/12 in kidney_A.
4. Covariate effects (reduced covariate-only model, per SD). Target GO count is the most consistent effect: positive in all
   five contexts (+0.26 lung, +0.29 lung_external, +0.32 kidney_B, +0.20 immune, +0.35 kidney_A) and the largest in four of
   them. In immune, target log mean expression (-0.73) and target log n_obs (+0.72) are larger, with opposite signs (the two
   are strongly collinear). STRING degree of the target is +0.10 to +0.17, except kidney_A (-0.07). log|co-expression| is
   +0.04 to +0.30. Score-covariate correlations are small at every layer: |r| <= 0.17 in the large contexts, <= 0.35 in immune
   and <= 0.25 in kidney_A. This is why adjustment barely moves the score coefficient.
5. Continuity, kidney_A (continuity_kidneyA.csv):
   A. Submitted zero-row convention, re-run under PYTHONHASHSEED=0: SV2-SV7 r_rb 0.094-0.163, p < 0.05 at 12/12 layers;
      SV5-SV7 -0.083 to 0.066, 1/12. The submitted CSV has 0.054-0.130 with 8/12, and 0/12. The joint permutation
      (replacing Stouffer) gives SV2-SV7 mean r_rb 0.125, p = 0.0016, and SV5-SV7 -0.015, p = 0.64.
      The difference from the submitted values is a hash-seed effect, not a change of method. The round-2 script draws its
      negatives from a list built by iterating a Python set, so its output depends on PYTHONHASHSEED. Over seeds 0-20, the
      layer-mean SV2-SV7 r_rb spans 0.077-0.150; seed 0 gives 0.125 and the submitted run 0.093. SV5-SV7 spans -0.039 to
      0.016. The submitted per-layer values lie inside the 21-seed range at every layer, but no seed in 0-20 reproduces them
      exactly (closest: max |diff| 0.024). Per-layer p-values exist only for seed 0, so "12/12 vs 8/12" must not be read as a
      change.
   B. Observed-gene SVD with the same band indices: SV2-SV7 0.007-0.188 (9/12, joint p = 0.0023); SV5-SV7 0.006-0.106
      (2/12, joint p = 0.066). With mapped bands: SV1-SV6 0.092-0.158 (12/12, joint p = 0.0015); SV4-SV6 -0.075 to 0.043
      (0/12, p = 0.70).
   C. Paper convention (within-TF, H3 pairs, two-way bootstrap): the full-512 residual r_rb is 0.159-0.268, with the CI
      excluding 0 at 8/12 layers (L0-L7). No SV band (SV1-SV6, SV2-SV7, SV4-SV6, SV5-SV7) has a CI excluding 0 at any layer.
      SV2-SV7 is borderline at L0-L1: r_rb 0.177 and 0.189, lower bounds -0.0013 and -0.0014. These point estimates are as
      large as variant A's (0.163, 0.155). From L5 to L11 the C point estimates fall to 0.03-0.08 (A: 0.09-0.13). Within C,
      residualisation lowers SV2-SV7 r_rb by 0.02-0.05 (e.g. L0 raw 0.210 to residual 0.177).

INFERENCE
- The covariates named by the reviewers (co-expression, expression level, detection, PPI degree, annotation depth) do not
  explain the within-TF scGPT edge signal. Adjustment removes 2-11% of the log-OR and <= 0.022 of AUROC. In the three large
  contexts the edge signal is small (residual AUROC 0.51-0.55), but it is not a by-product of these covariates.
- In lung, lung_external and kidney_B, the residual AUROC's upper bound stays below the 0.60 smallest effect of interest at
  every layer. Adjustment therefore does not change the "small" reading of H3 there. This is descriptive: B = 500 here,
  while the pre-registered decision uses B = 2,000 in E4a, and kidney_B's maximum (0.591) is close to the margin. immune and
  kidney_A are too small (36 and 29 TFs) to bound the effect below 0.60.
- Continuity. The covariates alone do not remove the band signal: in A and B the round-2 covariates are already
  residualised and SV2-SV7 still has p < 0.05 at most layers. Under C, which is pre-registered pairs and negatives,
  within-TF stratification and a two-way cluster CI, no SV band's CI excludes 0, while the full representation's does at
  8/12 layers. At L0-L1 the SV2-SV7 effect size under C is no smaller than in the submitted analysis. What is lost is
  precision: the cluster-level CI half-width is about 0.18 in r_rb. At later layers the estimate also shrinks. C changes the
  pairs, the negatives, the covariates, the stratification and the interval method all at once, so these data do not
  isolate which of them is responsible. The submitted SV5-SV7 null remains null: no joint test has p < 0.05 and no C
  interval excludes 0. B with identity bands has 2/12 nominal layers.

SURPRISES
- Inference scale (inference_scale.csv). The two-way bootstrap SD of the adjusted coefficient is 1.9-3.1x the model-based SE
  and 1.8-3.0x the SD of the within-TF permutation null. On the bootstrap scale the per-layer z = beta/SD is 3.5-5.4 (lung),
  3.4-4.4 (lung_external), 3.4-4.6 (kidney_B), 2.5-4.0 (immune) and 1.8-3.7 (kidney_A); model-based z is 3.7-14.6.
  The three quantities answer different questions. The model SE assumes pairs are independent given the TF. The permutation
  null is the distribution under the sharp null of within-TF exchangeability, for these TFs and pairs. The bootstrap
  approximates sampling variability over TFs and targets, including heterogeneity of the effect between TFs. In the large
  contexts most of the extra width is TF-level, not target recurrence. For the residual AUROC, the TF-only bootstrap CI is
  84-99% as wide as the two-way CI in lung, 79-92% in lung_external and 80-88% in kidney_B; target resampling adds more in
  immune (60-76%) and kidney_A (47-74%) (verify/verify_summary.json). The permutation and LRT p-values are therefore valid
  answers to the sharp-null question but overstate precision for generalisation to other TFs and targets. The bootstrap CIs
  are the intervals to report. The pigeonhole bootstrap may itself be somewhat conservative. The pre-registered H3 within-TF
  permutation null has the same scale relation: in E4a lung the null SD is 0.0058 and the two-way bootstrap SD is
  0.018-0.024 (E4a_edges/h3_primary_lung.csv). Its detection p-values answer the sharp-null question; its bootstrap
  intervals carry the generalisation uncertainty. This was not re-analysed here.
- OR versus AUROC precision. In lung the OR CI excludes 1 at every layer, while the raw AUROC CI includes 0.5 at 11/12
  layers (bootstrap z of the AUROC 0.8-1.8). Bootstrap draws of the two statistics correlate at 0.69-0.82, so they track the
  same signal; the difference is weighting (stratum_weights.csv). The pooled within-TF AUROC weights a TF by n_pos x n_neg:
  the top 5 TFs carry 51% of the weight in lung (effective number of TFs 15). The conditional logit weights roughly by
  k(n-k)/n: the top 5 carry 19% (effective number 76). The pooled AUROC is therefore dominated by a few heavily documented
  TFs, where the signal is weaker. In E4a lung the macro per-TF AUROC is 0.554-0.580 against a pooled 0.516-0.542 (L0: 0.580
  vs 0.541; E4a_edges/sens_macro.csv, h3_primary_lung.csv).

WHAT THIS DOES NOT SHOW
- It does not show that the score encodes regulation. Only the listed covariates are adjusted for. Not adjusted: shared
  annotation between TF and target (e.g. GO-term or pathway overlap, as opposed to counts), genomic proximity, literature
  co-mention, and any confounder not in the specification.
- TF-level confounding is removed by design (stratification). The TF-level covariate effects are not estimated.
- No direction claim: cosine and |co-expression| are symmetric.
- The joint permutation p = 0.002 is the resolution floor of N = 500. The permutation p could be smaller; it is not evidence
  of a large effect. Like the LRT, it addresses the sharp within-TF null and not generalisation over TFs and targets (see
  SURPRISES). Effect sizes and the bootstrap CIs carry the inferential weight.
- immune and kidney_A are small (36 and 29 TFs, 547 and 213 positives), and kidney_B shares cells with kidney_A (errata item
  2), so the five contexts are not five independent replications.
- The continuity comparison does not say which change between A/B and C (pairs, negatives, covariates, stratification,
  interval method) causes the loss of significance.

LIMITATIONS AND IMPLEMENTATION NOTES
- Model engine: our own implementation of the exact conditional likelihood (statsmodels cannot fit the large strata). This
  is a justified deviation from the specified engine, validated against statsmodels and independently reproduced by the
  verifier. Bootstrap and permutation refits used full Newton for immune, kidney_A, lung L0-L1, lung_external L0-L1 and
  kidney_B L0-L3, and chord Newton (fixed Hessian) for all other layers. Chord Newton reaches the same maximiser
  (|beta diff| <= 2.5e-8, check_fast.json; the verifier's lung L11 refits differ by <= 5e-10). Point estimates and SEs
  always use full Newton. All 60 x (1,000 bootstrap + 500 permutation) refits converged, with no NaN.
- Co-expression = Pearson r on the extraction cells (expr_extraction.npy), entered as log(|r| + 1e-3); the floor is arbitrary.
- The OR is per pooled SD of the cosine over all candidate pairs, not per within-TF SD.
- Bootstrap target multiplicity is implemented by replicating rows, and TF multiplicity by stratum weights (the verifier's
  physical stratum duplication gives identical refits). The percentile pigeonhole bootstrap tends to overestimate
  variance, so the CIs may be somewhat conservative.
- The residualisation removes covariate-explained variance from the score only, not from the label (a semi-partial
  measure). The conditional-logit OR is the fully adjusted estimate.
- B = 500 and N = 500, as the E5 specification states. This is lower than the CONVENTION3 confirmatory B = 2,000 and
  N = 10,000; nothing was reduced relative to the E5 specification. Errata item 3 (delete-d scaling) does not apply: E5 uses
  no vertex subsampling, only the two-way TF x target bootstrap.
- Variant A could not be reproduced exactly (the round-2 PYTHONHASHSEED was not recorded). The 21-seed spread is reported
  instead.
- Resumption: the first agent's background jobs completed all 60 context x layer checkpoints and the continuity runs. The
  merge (e5_merge.py) and diagnostics (e5_summarize.py, no new fits) were re-run afterwards.

VERIFICATION (independent adversarial verifier; scripts and outputs in verify/; collected in verify/verify_summary.json)
Reproduced from raw inputs with the verifier's own code (not e5_lib):
- Covariates: STRING-physical degree re-parsed from the gzipped links file (1,477,610 data lines, i.e. the physical file;
  86,519 symbol pairs >= 700), GO counts, log mean expression, log n_obs and Pearson r. All agree with E5's to <= 1.1e-13.
  obs_idx equals the non-zero embedding rows in all 5 contexts.
- Raw and residualised within-TF AUROC at all 60 context x layer cells (own V-centring and cosine; residualisation by OLS
  with explicit TF dummies; per-TF Mann-Whitney). Maximum difference vs residual_auroc.csv 0.0, and vs E4a's H3 primary
  statistic 0.0.
- Conditional logit with an independent engine (log-space elementary-symmetric DP in torch, autograd gradient, BFGS + Newton):
  adjusted beta vector, unadjusted beta, log-likelihood, LRT and SE at lung L0, lung L11, immune L2 and kidney_A L11. Max
  |beta diff| <= 1.9e-14, |SE diff| <= 2e-12, LRT p agrees to <= 2e-11 relative.
- Bootstrap refits (TF multiplicity by physical stratum duplication, draws regenerated from the documented RNG recipe):
  kidney_A L11 all 500 draws (max |diff| 5e-14). The OR CI [1.0048, 2.0098] is reproduced exactly; this is the borderline
  "excludes 1 in all 60 cells" case. immune L2 20 draws (1e-14); lung L11 1 draw (5e-10, chord-Newton layer).
- Permutation refits: kidney_A L11 20 draws, immune L2 10 draws and lung L11 1 draw, max |diff| 1.3e-11.
- Residual-AUROC two-way bootstrap (own weighted OLS + weighted Mann-Whitney, 500 draws): kidney_A L0 and kidney_B L0 (the
  closest to the 0.60 margin, upper bound 0.5910), every draw identical (max |diff| 0.0).
- Joint permutation p, z and joint bootstrap CIs re-derived from the stored draws for all 5 contexts: identical.
- Continuity C point estimates for all 4 SV bands + full_512 x 12 layers (own SVD of the V-centred observed matrix): max
  |diff| 1e-16. The C full_512 rows equal residual_auroc.csv kidney_A. The main-run variant A/B values equal the hash-seed-0
  spread file exactly (so the main run used PYTHONHASHSEED=0). Stratum-weight table reproduced.
Audit: no zero rows (observed-gene SVD, V-centred; variant A uses the zero-row SVD only as the documented legacy artefact);
SV indexing (0-based slices; mapped bands SV4-SV6 / SV1-SV6) correct; gene/pair index alignment correct; no supervised
probes (no fold leakage possible); the TRRUST in-degree covariate is correctly excluded; the within-TF permutation keeps
per-TF counts; the bootstrap unit is TF x target; permutation resolution supports the stated p. No model labels other than
scGPT appear (no ESM2 is used in E5). Not independently re-implemented: the round-2 variant A/B pair construction (checked
only for internal consistency and against the submitted values' hash-seed range).
Corrections to the implementer's text (no number in any CSV changed):
 (1) Boot SD / model SE ratio is 1.9-3.1x (not 1.8-3.1x).
 (2) GO count is the most consistent covariate, but not the largest in immune (mean expression -0.73, n_obs +0.72).
 (3) "LRT p ~ 1e-14" replaced by the per-context ranges.
 (4) "p = 0.002 ... anti-conservative" is wrong as stated: 0.002 is the resolution floor. The caveat is about the inferential
     target (sharp null vs generalisation), not the floor.
 (5) Surprise 1's mechanism: the extra bootstrap width is mainly TF-level (between-TF heterogeneity), not target recurrence,
     in the large contexts. The permutation and LRT answer the sharp-null question rather than being invalid.
 (6) The variant A re-run is hash-seed 0. Its 12/12 vs the submitted 8/12 is a seed effect.
 (7) Continuity C: SV2-SV7 is borderline at L0-L1, with effect sizes as large as the submitted ones. The loss is precision
     plus later-layer shrinkage. Attributing it to "the evaluation design" is not isolated by the data.
 (8) kidney_B's second residual-CI layer (L5) sits on the boundary (0.500007). The E4a macro vs pooled comparison is
     specified as L0 0.580 vs 0.541 (range 0.554-0.580 vs 0.516-0.542), and the E4a bootstrap SD is 0.018-0.024.
```
