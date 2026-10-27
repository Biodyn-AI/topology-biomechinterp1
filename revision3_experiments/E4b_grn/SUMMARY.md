# E4b — GENIE3 / GRNBoost2 baselines and direction test (verified)

Copy of the verifier-rewritten SUMMARY.txt (verify/ holds the checks).

```
E4b_grn -- GENIE3 / GRNBoost2 baselines on the H3 benchmark, plus a direction test
(sensitivity analysis; answers R1.1, R1.3, R1.2-direction). CONVENTION3 sha256 3e3249d1d0a8ce1bc3cf64639fe233b016535456f82432a8bbfcacdd64a05c01
(= FROZEN.txt; re-verified against CONVENTION3.md in results.json and verify/verification.json). No pass/fail rule attaches.
This SUMMARY was rewritten by the independent verifier (2026-10-27). The implementer's text is kept in
verify/SUMMARY_implementer_before_verification.txt. No bug changed any result. The corrections are rounding and wording,
plus added caveats. They are listed in section V.

STATUS: complete, full design. All 70 planned fits are done (grn_runs.csv, done=True for all).
- Tiers: T1 = the exact extraction cells (256; 192 in kidney_B), seeds 0-4. T2 = cache expr.npy (3,000 cells), seeds 0-1.
- Coverage: 2 methods x 5 contexts. The task's fallback (1 T2 seed in lung/lung_external/kidney_B) was not needed.

WHAT WAS DONE
- Regulators = TFs in cache/<ctx>/h3_pairs.npz plus Lambert TFs in V (input features). Targets = all h3 targets plus the
  direction-test endpoints. Sizes (regulators / targets fitted / h3 targets / direction pairs):
    lung           676 / 1781 / 1700 / 855
    lung_external  716 / 1785 / 1701 / 886
    kidney_B       590 / 1406 / 1340 / 669
    immune         118 /  261 /  250 /  62
    kidney_A        94 /  152 /  145 /  21
- Methods. Both are sklearn re-implementations; arboreto/dask was not used.
  - "GENIE3" = ExtraTreesRegressor, 300 trees, max_features='sqrt'. This is the extra-trees variant of GENIE3;
    arboreto's default is a random forest with 1,000 trees.
  - "GRNBoost2" = GradientBoostingRegressor: lr 0.01, up to 1,000 stages, depth 3, subsample 0.9, max_features 0.1,
    early stopping after 25 stages without improvement on a 10% validation split. Arboreto instead stops on
    out-of-bag improvement.
  - Importances are normalised per target. A regulator is never its own feature.
- Zero-importance rows:
  - 20 target fits had zero-variance expression and get importance 0: 4 lung T2 targets x 4 runs and 1 lung_external T2
    target x 4 runs.
  - [verifier] A further 70 GRNBoost2 fits of non-constant targets had all-zero importances (no informative split):
    0-9 per fit. GENIE3 has none. These rows hold <= 1.8% of h3 pairs in any fit (verify/v4_zero_rows.csv).
- Headline score per (method, tier) = seed ensemble (mean of per-seed normalised importances). Per-seed AUROCs give the
  run-to-run SD. Pair score = importance(TF -> target). Exact zeros are kept and ties count 1/2 (mid-rank).
- Metric = the H3 statistic: within-TF stratified AUROC (stats3), with:
  - a two-way TF x target pigeonhole bootstrap, B = 2,000 (95% CI and one-sided 95% upper bound);
  - a TF-only bootstrap as a check;
  - a within-TF permutation p (N = 2,000; floor 1/2001 = 5e-4).
- Paired deltas use identical bootstrap weight draws for every score in a context.
- Comparators on identical pairs:
  - |Pearson co-expression| on the same cells as the GRN tier;
  - the scGPT primary scorer: signed cosine of V-centred embeddings, all 12 layers. Its point AUROCs match E4a's
    primary H3 table exactly (max |diff| = 0 in all 5 contexts).
- Errata item 3 (delete-d rescaling) does not apply: E4b has no vertex-subsampling interval. No ESM2 embedding is used.

DIRECT EVIDENCE (grn_summary.csv, grn_deltas.csv, direction_test.csv, results.json; verifier additions in verify/)
1) Within-TF AUROC per context. Each GRN cell is the seed-ensemble AUROC [two-way 95% CI], followed by the SD over seeds.
   |coexpr| is given as T1 / T2. scGPT is the cosine at L0, with the L0-L11 range in parentheses.
   lung
     GENIE3     T1 .516 [.488,.550] sd .003    T2 .513 [.480,.550] sd .003
     GRNBoost2  T1 .518 [.494,.546] sd .006    T2 .512 [.483,.544] sd .002
     |coexpr| .508 / .515    scGPT .541 (.516-.542)
   lung_external
     GENIE3     T1 .535 [.508,.565] sd .003    T2 .538 [.501,.585] sd .0003
     GRNBoost2  T1 .521 [.493,.551] sd .007    T2 .535 [.504,.570] sd .001
     |coexpr| .524 / .551    scGPT .538 (.518-.538)
   kidney_B
     GENIE3     T1 .529 [.496,.564] sd .002    T2 .528 [.498,.563] sd .0003
     GRNBoost2  T1 .519 [.492,.548] sd .003    T2 .528 [.501,.557] sd .0003
     |coexpr| .515 / .522    scGPT .554 (.524-.554)
   immune
     GENIE3     T1 .578 [.508,.664] sd .005    T2 .610 [.532,.693] sd .002
     GRNBoost2  T1 .574 [.502,.667] sd .007    T2 .605 [.529,.685] sd .009
     |coexpr| .559 / .562    scGPT .638 (.577-.638)
   kidney_A
     GENIE3     T1 .557 [.460,.677] sd .004    T2 .596 [.504,.701] sd .003
     GRNBoost2  T1 .554 [.449,.670] sd .013    T2 .595 [.502,.700] sd .011
     |coexpr| .559 / .577    scGPT .647 (.601-.652)
   - T1 has 5 seeds and T2 has 2.
   - One-sided 95% upper bounds of the GRN ensembles: <= .576 in the three large contexts; .647-.683 in immune and kidney_A.
   - Run-to-run SD is 0.0003-0.013. The GRN bootstrap half-widths are 0.026-0.110, at least 4x the SD in every cell
     (lowest ratio 4.1, lung_external GRNBoost2 T1). Pair/TF/target sampling, not seed noise, limits these estimates.
2) Ties.
   - GENIE3 scores exactly 0 for 0.06-4.6% of pairs.
   - GRNBoost2 scores exactly 0 for 4.5-33% of pairs: lung_external 32% / 33% (T1 / T2), lung 22% / 25%, kidney_B 26% / 27%.
   - All ties are handled by mid-rank.
   - [verifier] GRNBoost2 stopped at stage 26 (no validation improvement after the first stage) for 20-31% of targets on
     T1, where the validation split is about 19-26 cells, and for 5-16% on T2 (verify/v5_early_stop.csv). Depth-3 trees
     with max_features 0.1 and few stages leave most regulator columns unused. These two effects were not separated.
3) Paired delta, GRN minus |co-expression| on the same cells.
   - Point deltas range from -0.017 to +0.048, and every 95% CI includes 0.
   - The largest is immune GENIE3 T2, +0.048 [-0.031, +0.121].
   - [verifier] The one-sided 95% upper bound of the delta is < 0.05 for every GRN ensemble in lung, lung_external and
     kidney_B (max .044, kidney_B GENIE3 T2). In immune and kidney_A it is .070-.121 (verify/v5_h3_margins.csv).
4) Paired delta, GRN minus scGPT signed cosine.
   - Versus L0 the point deltas are mostly negative: lung -0.023 to -0.029; lung_external -0.016 to +0.001;
     kidney_B -0.026 to -0.035; immune -0.028 to -0.064; kidney_A -0.051 to -0.094.
   - Versus late layers (L8-L11) they are near 0: lung -0.011 to +0.002; lung_external -0.004 to +0.021;
     kidney_B -0.014 to +0.005; immune -0.012 to +0.032.
   - kidney_A is negative at every layer (max -0.005).
   - No two-sided 95% CI excludes 0 in any of the 240 context x method x tier x layer comparisons.
   - [verifier] One one-sided bound does exclude 0: kidney_A GENIE3 T1 vs L1, upper95 = -0.0008. It is unadjusted and
     1 of 240.
5) T2 minus T1 (same method; GENIE3 / GRNBoost2): lung -0.003 / -0.006; lung_external +0.003 / +0.013;
   kidney_B -0.001 / +0.008; immune +0.032 / +0.031; kidney_A +0.040 / +0.041. All CIs include 0.
6) scGPT cosine minus |co-expression| (a by-product of the same draws).
   - At L0 the point delta is positive in 9 of 10 context x tier cells (+0.014 to +0.089). It mostly shrinks with depth.
   - In lung_external T2 it is negative at every layer (-0.013 to -0.034).
   - No two-sided CI excludes 0.
   - [verifier] In immune, at L0-L3 and in both tiers, the one-sided 95% lower bound is > 0 (+0.0025 to +0.0069). This is
     unadjusted (8 of 120 comparisons, verify/v5_onesided.csv).
7) Direction test. Pairs: TRRUST i->j with both genes in R and j->i in no regulatory set. Statistic: seed-ensemble
   fraction of untied pairs with importance(i->j) > importance(j->i).
   - lung: GENIE3 .511 / .509 (T1 / T2); GRNBoost2 .528 / .523. n untied = 834 / 845 / 762 / 723. All sign-test p >= .13.
   - lung_external: GENIE3 .525 / .505; GRNBoost2 .501 / .487. All p >= .15.
   - kidney_B: GENIE3 .523 (p .24) / .551 (p .0095, Clopper-Pearson 95% [.512, .589], source-regulator cluster CI
     [.492, .609]). GRNBoost2 .502 / .505, p >= .84.
   - immune (about 60 pairs): .467-.544, p >= .59. kidney_A (21 pairs): .429-.619, p >= .38.
   - Per-seed fractions span .38-.62 in the small contexts and .46-.55 in the large ones.
   - |co-expression| and cosine are symmetric: n_untied = 0 for every pair by construction, so they cannot carry direction.
8) [verifier] Within-TF permutation p for the 20 GRN ensembles: all <= .019, and 14 are at the 5e-4 floor. This holds even
   where the two-way CI includes .5 (e.g. lung GRNBoost2 T1, .518 [.494, .546], p = 5e-4). Example for |co-expression|:
   immune T1, .559 [.479, .652], p = 5e-4. Verifier recomputation: permutation-null SD 0.015 vs bootstrap SD 0.048.

INFERENCE
- The expression-based GRN baselines are as weak as every other scorer on this benchmark.
  - Within-TF AUROC is .51-.54 in the three large contexts and .55-.61 in immune and kidney_A.
  - They show no detectable gain over |co-expression| and no detectable difference from scGPT cosine. The scGPT point
    estimates at early layers are higher in 4 of 5 contexts (lung_external is a tie), but no interval separates them.
  - So scGPT's H3 result is not "below what an expression GRN method recovers". Equally, there is no evidence that scGPT
    encodes more edge information than co-expression or GRN methods.
- [verifier] Applied descriptively (no decision rule attaches to E4b), the H3 margins give the following. In lung,
  lung_external and kidney_B, every GRN ensemble has:
  - permutation p < .05 (a single score, so no max-T);
  - one-sided upper bound < .60 ("AUROC >= .60 excluded");
  - delta vs |co-expression| upper bound < .05 ("not above co-expression").
  For scGPT, E4a reports "detected" and "AUROC >= .60 excluded" in these contexts, but "not above co-expression" only
  in lung_external with the 3k comparator. In E4b's draws, the maximum over layers of the scGPT-minus-|coexpr| upper
  bound is .035-.088 (grn_deltas.csv); for GRN-minus-|coexpr| the upper bounds are .015-.044. The scGPT bounds are
  higher because its point deltas are larger at early layers. This is not a tested difference between the scorers.
  In immune and kidney_A neither margin is excluded, because the small designs give wide intervals.
- The direction test is null. 19 of 20 ensemble tests have p > .13. The exception, kidney_B GENIE3 T2 (p = .0095), should
  not be read as signal:
  - it fails Bonferroni over the 20 ensemble tests (.0025);
  - its source-regulator cluster CI includes .5;
  - the same method on T1 (5 seeds) gives .523 (p .24);
  - kidney_A, the only other kidney context, consists of overlapping cells of the same donor (not a replication) and has
    only 21 pairs.
  Tree-importance asymmetry does not recover TRRUST direction here.
- [verifier, reworded; consistent with E4a section C] The within-TF permutation and the two-way bootstrap answer different
  questions.
  - The permutation tests the sharp null that, in this fixed benchmark, documented and candidate targets are exchangeable
    within each TF.
  - The two-way bootstrap describes uncertainty over a population of TFs and targets.
  - A floor-level permutation p therefore shows a within-benchmark association. It is not evidence that the
    population-average edge AUROC exceeds .5. The paper should quote the two-way intervals for that question.

HYPOTHESES (not tested here)
- The ceiling may be set by the labels. TRRUST positives, compared against candidate non-edges matched on in-degree, may
  include many context-inactive or indirect edges that no expression-derived scorer can separate in 192-3,000 cells.
  E4a additionally found that other-TF in-degree is not neutralised by the matching.
- The higher AUROCs of every scorer in immune and kidney_A may reflect fewer, better-studied TFs rather than easier biology.

WHAT THIS DOES NOT SHOW / LIMITATIONS
- This is not a test of arboreto's GENIE3/GRNBoost2. The methods are sklearn re-implementations with the specified
  hyperparameters: GENIE3 uses extra-trees with 300 trees; GRNBoost2 stops early on a validation split, which ends 20-31%
  of T1 fits at 26 stages. Label them "GENIE3 (extra-trees, sklearn)" and "GRNBoost2-style boosting (sklearn)" in the paper.
- Only 2 T2 seeds, so each T2 seed SD rests on n = 2. Seed variance is in any case small relative to sampling variance.
- immune (36 TFs, 547 positives) and kidney_A (29 TFs, 213 positives) have wide CIs. Their direction tests (62 and 21
  pairs) are uninformative.
- kidney_A and kidney_B are overlapping cells of one donor, not independent replicates.
- The direction test uses only TF-TF pairs documented in one direction. A reverse edge missing from TRRUST/DoRothEA is not
  evidence that the edge is absent.
- The Clopper-Pearson interval treats pairs as independent. The cluster bootstrap resamples source regulators only, but
  pairs also share target genes.
- No conclusion about causal or mechanistic encoding in scGPT follows, because all scorers are compared on one label set.
- Some importance rows are all zero: 20 constant-target T2 fits and 70 GRNBoost2 fits with no informative split. They
  become ties and are counted 1/2.

SURPRISES
- 12-16x more cells (T2) did not help the GRN methods in lung, lung_external or kidney_B.
- |co-expression| gained most in lung_external (T1 .524 -> T2 .551), where it is the highest-scoring method point-wise.
- GRNBoost2 leaves up to a third of pairs at exactly 0; GENIE3 leaves almost none.

CORRECTIONS MADE BY THE IMPLEMENTER ON RESUME (kept)
- e4b_eval.py: a 'layer' column was added to the scGPT-minus-|coexpr| delta rows (backup e4b_eval.py.bak_pre_layerfix).
  No number changed.
- make_results.py: the stale compute-plan text was replaced, and seed_design_final was added.
- All five contexts were re-evaluated. grn_runs.csv wall_seconds_this_process is the wall time of the process that
  assembled each fit, not total CPU time.

V. VERIFICATION (independent adversarial verifier; scripts and outputs in verify/)
Checks against the audit list:
(1) Primary analysis. E4b is a sensitivity analysis and runs no primary analysis. Its metric is exactly the H3 statistic
    of CONVENTION3 §4: same pairs, within-TF stratified AUROC, two-way pigeonhole bootstrap with B = 2,000, one-sided
    upper bound, and the within-TF permutation null. N = 2,000 rather than 10,000 is acceptable for a sensitivity
    analysis, and stats3.assert_resolution(2000, .05, 6) passed. The fits follow the task's GENIE3/GRNBoost2 specification.
(2) Geometry. The embedding rows used are exactly the non-zero rows: obs_idx = the non-zero rows at L0 in all 5 contexts
    (v1). Cosine is computed on V-centred vectors; no SVD is used.
(3) Alignment.
    - obs_genes = h5ad var_names[obs_idx] in all 5 contexts.
    - expr_extraction.npy, rebuilt by gene name from the h5ad rows at extraction_cells.npy, matches with max |diff| = 0
      in all 5 contexts.
    - The T2 expr.npy gene axis agrees with T1: per-gene mean Pearson r = .992-.998, against |r| <= .021 for a shuffled axis.
    - The regulator, target and direction-pair sets rebuilt from the raw TRRUST/DoRothEA/Lambert files are identical to
      E4b's in all 5 contexts.
    - h3 pairs: no self pairs; every positive is in TRRUST; 0 candidates are in any regulatory set in either direction.
    - All 35 grn_scores_*.npz pair scores equal an independent lookup into partial/*.npz (max |diff| = 0).
(4) SV indexing: not applicable.
(5) Leakage: none possible. The GRN fits never see labels, and labels are used only in evaluation.
(6) Nulls. The within-TF permutation preserves each TF's positive count, as claimed. Its meaning is restated in INFERENCE.
(7) Bootstrap unit. Two-way TF x target, as CONVENTION3 §6 specifies. There is no vertex subsampling, so errata item 3
    does not apply. The direction-test cluster bootstrap is one-way (added as a limitation).
(8) Permutation resolution. Floor 5e-4 < .05/6. No claim rests on resolution beyond the floor. The direction test is exact.
(9) Numbers. Every SUMMARY number was checked against the CSVs at 4-5 decimals. The corrections are in "Wording/rounding
    corrections" below.
(10) Labels. scGPT contexts as in CONVENTION3 §1. No ESM2 is used. The method labels were clarified (see LIMITATIONS).

Independent recomputations (own code, no stats3/e4b functions for the statistics):
- Refits with the verifier's own code, the spec hyperparameters and the same seeds reproduce E4b's normalised importances
  to max |diff| <= 1.3e-8, which is float32 rounding (v2_refit.json):
  - kidney_A T1 GENIE3 s0, all 152 targets;
  - kidney_A T1 GRNBoost2 s0, 150/152 identical. The other 2 rows have zero total importance: E4b stores 0 and the
    verifier's 0/0 gave NaN. This led to the zero-row finding above;
  - lung T2 GRNBoost2 s1, 34 targets including the 4 constant ones;
  - immune T2 GENIE3 s1, 30 targets.
- Headline point estimates. Per-TF sklearn roc_auc_score pooled with n_pos*n_neg weights reproduces every value exactly
  (v3_headlines.json):
  - lung GENIE3 T1 ensemble .51634;
  - lung GRNBoost2 T1 per-seed SD .00588;
  - lung_external |coexpr| T2 .55126;
  - lung_external GRNBoost2 T2 zero fraction .3292;
  - kidney_B scGPT L0 .55446 and L11 .52398, from the raw embeddings;
  - immune GENIE3 T2 minus |coexpr| T2 delta +.04799;
  - immune |coexpr| T1 .55898, with permutation p = 5.0e-4 under the verifier's own permutation.
- Bootstrap CIs from the verifier's own pigeonhole code (B = 500-1,000, different seeds) agree within Monte-Carlo error:
  - lung GENIE3 T1: [.488, .546] vs [.488, .550];
  - lung_external |coexpr| T2: [.517, .592] vs [.517, .585];
  - immune delta: [-.023, .120] vs [-.031, .121];
  - immune |coexpr| T1: [.481, .667] vs [.479, .652].
- Direction test, recomputed with own lookup and scipy binomtest:
  - kidney_B GENIE3 T2: 368/668 = .5509, p = .00948, CP [.5123, .5891];
  - lung GRNBoost2 T1: 402/762 = .5276, p = .137.
  Both are exact matches.

Wording/rounding corrections to the implementer's SUMMARY:
- kidney_B scGPT L0: .555 -> .554; its range .524-.555 -> .524-.554 (.55446; E4a also reports .554).
- lung GENIE3 T2 CI lower bound: .479 -> .480 (.47950).
- lung_external GRNBoost2 T1: .522 -> .521 (.52149).
- kidney_B GENIE3 T1 CI upper bound: .565 -> .564 (.56446).
- kidney_B GENIE3 T2 CI lower bound: .497 -> .498 (.49750).
- kidney_B GRNBoost2 T1: .520 -> .519 (.51949).
- kidney_A GENIE3 T2 CI lower bound: .503 -> .504 (.50353).
- immune GENIE3 T2 delta vs |coexpr|: CI [-0.032, +0.122] -> [-0.031, +0.121] (-.03147, .12145).
- "bootstrap half-widths 0.026-0.12, SD at least ~3x smaller" -> "GRN half-widths 0.026-0.110, at least 4x the SD".
  The value 0.12 was a |coexpr| half-width.
- GRNBoost2 zero fractions: lung_external T1 (32%) and lung T1 (22%) were added for completeness.
- DE4 omitted lung_external and did not say that kidney_A is negative at every layer. Both were added, as were the
  late-layer ranges and the single one-sided exclusion.
- "12x more cells" -> "12-16x": kidney_B T1 has 192 cells.
- "20 target fits ... get importance 0" was correct, but 70 further GRNBoost2 all-zero rows were unreported. Added.
- The method labels were made precise: extra-trees GENIE3; validation-split early stopping in GRNBoost2, which ends
  20-31% of T1 fits at 26 stages.
- INFERENCE on permutation p: "the paper should not quote these p-values as evidence of edge signal" was reworded. They
  show a within-benchmark association (as H3 "detected" does in E4a), not a population-level AUROC > .5. Added: all 20
  GRN ensembles have p <= .019.
- The kidney_A/kidney_B argument against the kidney_B direction result was reworded. Shared cells do not by themselves
  discount kidney_B; the point is that kidney_A (21 pairs, overlapping cells) cannot replicate it.
- Added a descriptive application of the H3 margins to the GRN baselines (verify/v5_h3_margins.csv).
No E4b output CSV/JSON was changed by the verifier, and no re-run was needed.

FILES: grn_runs.csv, grn_summary.csv, grn_deltas.csv, direction_test.csv, results.json, grn_scores_<ctx>_<tier>_<seed>.npz
(35 files), eval/<ctx>_{summary,deltas,direction}.csv and _meta.json, partial/ (fits and batch checkpoints),
eval_final_*.log. Verifier files are in verify/:
  scripts   v1_inputs.py, v2_refit.py, v3_headlines.py, v4_zero_rows.py, v5_extra.py, v6_collect.py
  outputs   v1_inputs.json, v2_refit.json, v3_headlines.json, v4_zero_rows.csv, v5_early_stop.csv, v5_onesided.csv,
            v5_h3_margins.csv, v5_extra.json, verification.json
  backup    SUMMARY_implementer_before_verification.txt
```
