# E4a_chipfix — ChIP-seq benchmark rebuilt with corrected candidates (verified)

Copy of the verifier-rewritten SUMMARY.txt (verify/ holds the checks).

```
E4a_chipfix - corrected ChIP-seq-only independent-positive benchmark (E4a sensitivity E7), 2026-10-27
This file: an independent adversarial VERIFICATION section (V), then the implementer's report with the verifier's corrections applied
inline and marked [verifier ...]. The unmodified implementer text is in verify/SUMMARY_implementer.txt.

V. VERIFICATION (independent; scripts verify/verify_build.py and verify/verify_auroc.py; outputs verify/verify_build.csv,
   verify/verify_build_vs_audit.csv, verify/verify_auroc.csv; logs verify/*.log)
Verdict: no code bug was found and nothing was re-run by the implementer's scripts. The fixed pairs and every number in
sens_chip_fixed.csv that I checked are reproduced exactly. Four statements in the interpretation are corrected (V5): one is a
"new result" that is neither new nor robust, and one result was omitted.

V1. Rebuild, all 16 contexts (verify_build.py).
- How it was done:
  - The reference sets are parsed directly from the raw TSVs (TRRUST, DoRothEA main, DoRothEA ChIP).
  - The candidate builder is my own re-implementation. No prep_pairs.build or load_refs, and nothing from e4a_chipfix.py.
  - The RNG replay is my own: default_rng(20261026) -> h3 build -> legacy draws -> ChIP draw.
- Result: my builder reproduces ../cache h3_pairs.npz, h3_legacy_pairs.npz and the OLD h3_chip_pairs.npz (old rule anyreg - set(ch_v))
  bit for bit. With the full anyreg it reproduces cache/<ctx>/h3_chip_pairs_fixed.npz bit for bit (tf, tg, y, match_level). This holds in
  all 16 contexts.
V2. Candidates in the fixed files (all 16 contexts), each count checked in both directions (a,b) and (b,a):
  - 0 candidates in anyreg;
  - 0 equal to any ChIP-independent pair of ANY TF within V;
  - 0 equal to a benchmark positive;
  - 0 in TRRUST, 0 in DoRothEA main (any level), 0 in DoRothEA ChIP (any level);
  - 0 self pairs, 0 duplicate (tf, tg) rows under either label.
V3. Positives:
  - The y==1 rows are identical to the old file, same order, in all 16 contexts.
  - They equal exactly the ChIP-independent pairs within V for the kept TFs.
  - 0 TFs were dropped (lung 57,017 / 225; lung_external 54,332 / 218).
V4. Match levels:
- The shares recomputed from the files equal chip_match_levels.csv and build_audit.csv (to 1e-12). Lung goes from 0.511/0.057/0.432 to
  0.427/0.062/0.511.
- Consistency: every level-0 candidate matches the (in-degree decile, detection decile) of at least one positive of its TF, and every
  level-1 candidate matches an in-degree decile. There are 0 violations in all contexts.
- The level-2 assignment is covered by the bit-exact replay.
- The audit counts are reproduced exactly: old own-positive copies (lung 52,908), old anyreg rows (54,427), retained (416,561 = 74.7%) and
  new n_neg.
V5. AUROC and bootstrap (verify_auroc.py; no stats3 / e4a_lib / e4a_main import).
- Method, all written from scratch:
  - Signed cosine of the V-centred 512-d embedding.
  - |Pearson r| rounded to 12 decimals.
  - Within-TF stratified AUROC via scipy mid-ranks.
  - Two-way TF x target bootstrap with the stats3 integer stream, rng [seed_of(ctx), 40], B = 2,000. My weighted evaluator was
    cross-checked against a brute per-TF evaluator.
- Fixed pairs, lung: L0 0.528384, L11 0.536046, |co-expr| extraction 0.541997, 3k 0.558652.
  - These equal sens_chip_fixed.csv to |diff| = 0.
  - CI, UB95, SD and both deltas (estimate, CI lower bound, UB) also match, max |diff| < 1e-16.
  - Rounding the comparator moves the |co-expr| AUROC by 1.1e-6 only.
- Old pairs, lung (code-path check): L0 0.525321, L11 0.532027, |co-expr| 0.536369 / 0.550657. These equal E4a_edges/sens_chip.csv,
  including CI, UB and deltas (< 1e-16).
- lung_external, fixed pairs (point estimates): L0 0.533305, L11 0.534353, |co-expr| 0.515166 / 0.527280, all equal to sens_chip_fixed.csv.
V6. Claims checked against the CSVs:
- Confirmed: the lung, lung_external, kidney, immune and range numbers in section 3; the delta ranges and UBs; "UB < 0 at 11/12" for the
  lung 3k delta; the lung_external CI excluding 0 only at L3; the match-level-0 numbers; the shifts +0.003..+0.012 (cosine) and
  +0.002..+0.013 (comparator); and that the deltas change by <= 0.007 (0.0072).
- Corrected (applied inline below):
  (a) "redraw adds <= 0.002": true in lung (<= 0.0011) but <= 0.003 over all contexts (immune 0.0030).
  (b) "Every scorer falls when fallback candidates are removed": not in kidney_ckpt_brain (cosine minimum 0.501 -> 0.508).
  (c) The "new descriptive result" (lung: 3k |co-expr| beats scGPT, delta UB < 0 at 11/12) is NOT new: the old file had UB < 0 at 10/12.
      It is NOT robust either: it reverses in the match-level-0 block (delta +0.0004..+0.0097). It reflects unmatched fallback candidates
      and should not be reported as co-expression outperforming the embedding.
  (d) Omitted by the implementer: in the match-level-0 block of lung_external, scGPT exceeds extraction-cell |co-expr| by +0.019..+0.024,
      with the CI lower bound > 0 at 8/12 layers (L0-L7; unadjusted, UB <= 0.040). In lung the level-0 delta is positive as well
      (+0.005..+0.014, CI contains 0). So the sign of the scGPT-vs-co-expression delta in the lung pair depends on candidate matching.
- The E4a conclusions still stand. The CI excludes 0.5 at 12/12 layers in both lung contexts in both blocks; AUROC >= 0.60 is excluded
  in the lung pair (UB <= 0.561); and no >= 0.05 scGPT advantage over |co-expr| is seen in the lung pair (max delta UB 0.047 all-candidate,
  0.040 level 0).
V7. Not verified independently:
- The bootstrap CIs outside lung. They rest on the same code path as lung, which matches exactly, and on the implementer's 16-context
  repro.
- The scGPT embeddings themselves (used as given from ~/biomechinterp-local/subproject41_r3/emb).
V8. Still open (not in scope here): E4a_edges/SUMMARY.md section E7, the manuscript ChIP numbers and tracking/prompt.md still show the
    old values. The replacements are in section 4 below, with the verifier corrections.

==================================================================================================== implementer report (corrected)
CONVENTION3.md sha256 3e3249d1d0a8ce1bc3cf64639fe233b016535456f82432a8bbfcacdd64a05c01 (= FROZEN.txt; in results.json)
Script: e4a_chipfix.py (build | sens | repro | aggregate); runner run_all.sh; logs build.log, w1-w5*.log.
Nothing outside E4a_chipfix/ was written. prep_pairs.py, stats3.py, contexts.py, ../cache/ and ../E4a_edges/ were only read.

0. THE DEFECT
- prep_pairs.main() built cache/<ctx>/h3_chip_pairs.npz with build(ch_v, ..., anyreg - set(ch_v), ...).
- So a TF's own ChIP-only positives, and their reverses, were not excluded from its candidate non-edges.
- In lung, 52,908 of 557,447 candidate rows (9.5%) were copies of a positive row of the same TF. Another 1,769 had the reverse pair as a ChIP
  positive; 54,427 rows were in anyreg in either direction.
- The contaminated share is 9.5-13.9% in every context (build_audit.csv: old_neg_is_own_chip_positive / old_n_neg).
- CONVENTION3 sections 3-4 require every pair in any regulatory resource to be excluded. anyreg = TRRUST + DoRothEA main + DoRothEA ChIP, any
  level, both directions.

1. REBUILD (build_audit.csv / .json; cache/<ctx>/h3_chip_pairs_fixed.npz, same layout as the old file: tf, tg, y, match_level)
- RNG replay (direct evidence). I replayed prep_pairs.main() per context: default_rng(20261026); then the h3_pairs build and the legacy draws,
  which consume the generator first; then the ChIP draw.
  - The replay reproduces ../cache h3_pairs.npz, h3_legacy_pairs.npz and h3_chip_pairs.npz bit for bit in all 16 contexts.
  - The fixed draw starts from a deep copy of the identical generator state. It calls the unmodified prep_pairs.build with the full anyreg.
    That is the only change.
- Positives are IDENTICAL to the old file: same rows in the same order, in all 16 contexts. No TF and no positive was dropped (lung 57,017 / 225 TFs;
  lung_external 54,332 / 218).
- Contamination after the fix: 0 candidate rows are in anyreg in either direction, in every context.
- How the candidates changed:
  - The candidate pools shrink, so rng.choice draws from different sets. The random stream therefore diverges at the first affected TF (AHR in
    the lung/kidney contexts, BACH2 in immune/kidney_A). From there on, the remaining draws are re-randomised under the same procedure.
  - Lung: 555,997 new vs 557,447 old candidates. 416,561 (74.7%) of the old candidate pairs are retained, 140,886 are dropped (54,427 of them
    contaminated) and 139,436 are new.
  - lung_external: 74.7% of the old candidates retained.
  - kidney contexts: 77-78% retained. Immune contexts and kidney_A: 84-86% retained.
  - The new candidate counts are smaller because some pools are exhausted: immune 21,456 -> 20,180; kidney_A 16,851 -> 14,841.
- Share of candidates at each match level (chip_match_levels.csv; level 0 = in-degree AND detection decile, 1 = in-degree only, 2 = unmatched
  fallback):
    lung           old 0.511 / 0.057 / 0.432  ->  new 0.427 / 0.062 / 0.511
    lung_external  old 0.503 / 0.071 / 0.426  ->  new 0.424 / 0.072 / 0.504
    kidney_B       old 0.493 / 0.064 / 0.444  ->  new 0.410 / 0.068 / 0.522   (the three ckpt contexts: 0.414 / 0.069 / 0.517)
    immune group   old 0.43-0.48 level 0      ->  new 0.35-0.41 level 0; level 2 rises to 0.49-0.55
    kidney_A       old 0.455 / 0.099 / 0.447  ->  new 0.354 / 0.111 / 0.535
  The contaminated rows were disproportionately fully matched. A factor's own ChIP targets share the in-degree and detection deciles of its
  positives, which is what makes them the best matches. So the corrected benchmark is LESS well matched: about half the candidates are now
  unmatched fallbacks.

2. CODE-PATH CHECK (repro_old_sens_chip.csv)
- I ran the same functions on the OLD pairs with E4a's settings:
  - B = 2,000, or 1,000 for kidney_B and the ckpt contexts;
  - seed [seed_of(ctx), 40];
  - raw float comparator outside the lung pair, tie-corrected comparator in lung and lung_external.
- This reproduces every value in E4a_edges/sens_chip.csv: AUROC, CI, UB, SD and both deltas, in all 16 contexts, max |diff| < 1e-16.
- The new numbers therefore differ from the old ones only through the pairs, the uniform B = 2,000 and the uniform tie-corrected comparator.

3. RESULTS ON THE FIXED PAIRS (sens_chip_fixed.csv: same columns as E4a_edges/sens_chip.csv; old vs new per row in chip_old_vs_new.csv)
Settings:
- Scorer: signed cosine of the V-centred full 512-d embedding, per layer.
- Statistic: within-TF stratified AUROC.
- Intervals: two-way TF x target pigeonhole bootstrap, B = 2,000 in all contexts, seed [seed_of(ctx), 40].
- Comparator: |r| rounded to 12 decimals (the verifier's tie fix), on the extraction cells and on the 3,000-cell tier. Deltas use the same
  resamples.
Results:
- lung:
  - AUROC 0.525-0.534 -> 0.528-0.538. CI lower bound 0.506-0.514, > 0.5 at 12/12 layers. Max UB 0.552 -> 0.559.
  - |co-expr| 0.536 / 0.551 -> 0.542 / 0.559 (extraction / 3k).
  - Delta vs extraction cells: -0.014 to -0.004, max UB 0.018 (old 0.016).
  - Delta vs 3k: -0.030 to -0.021, UB <= 0.0001 at every layer, and < 0 at 11/12 layers.
- lung_external:
  - AUROC 0.530-0.537 -> 0.533-0.541. CI > 0.5 at 12/12 layers. Max UB 0.555 -> 0.561.
  - |co-expr| 0.513 / 0.524 -> 0.515 / 0.527.
  - Delta vs extraction cells: +0.018 to +0.026, max UB 0.047 (old 0.042). The two-sided CI excludes 0 only at L3 (lower bound +0.0004;
    unadjusted, Monte-Carlo-borderline).
  - Delta vs 3k: +0.006 to +0.014, UB <= 0.034.
- kidney_B 0.534-0.546 (CI > 0.5 at 12/12). ckpt_whole_human 0.531-0.541 (12/12). ckpt_kidney 0.523-0.529 (0/12). ckpt_brain 0.501-0.514 (0/12).
  Max UB 0.539-0.568. Delta UB <= 0.037 (extraction) and <= 0.014 (3k).
- Immune group and kidney_A: 0.526-0.594. Max UB 0.585-0.642: the immune UBs now exceed 0.60 (they did before as well, 0.601-0.624).
  - Delta vs extraction-cell |co-expr|: the delta CI contains 0 at every layer in every context.
  - Delta UB: 0.023 (kidney_A) to 0.109 (immune_Bcell).
- Every new scGPT AUROC is above the old one: +0.003 to +0.012. The only exception is kidney_ckpt_brain, which changes by -0.0001
  to +0.0013. The comparator moves the same way (+0.002 to +0.013).
- Range over all contexts and layers: 0.501-0.594 (old 0.50-0.58).
Decomposition (chip_old_vs_new.csv: old_minus_contaminated_auroc = old pairs with the contaminated candidate rows simply deleted):
- Deleting the contaminated rows reproduces most of the shift: lung 0.528-0.537 vs fixed 0.528-0.538. The redraw adds <= 0.0011 in lung and
  <= 0.003 across all contexts (immune 0.0019-0.0030). [verifier correction: was '<= 0.002']
- Mechanism (inference, consistent with the numbers): a contaminated candidate is an exact copy of a positive of the same TF. Any symmetric
  scorer ties it with that positive, which shrinks the AUROC toward 0.5 by about f x (A - 0.5), where f is the contaminated share.
  - Lung: 0.095 x 0.035 = 0.003, which is the observed shift.
  - This affects scGPT and |co-expression| alike, so the deltas barely move (<= 0.007).
Match-level-0 block (sens_chip_fixed_matchlevel0.csv: fully matched candidates + all positives; seed [seed_of(ctx), 41]; B = 2,000):
- lung: 0.521-0.531, CI > 0.5 at 12/12, UB <= 0.544. |co-expr| 0.516 / 0.521. Delta UB 0.030 / 0.027.
- lung_external: 0.526-0.531, CI > 0.5 at 12/12, UB <= 0.546. Delta UB 0.040 / 0.030.
- kidney_B / ckpt_whole_human: 12/12 > 0.5. Immune: CI > 0.5 at 0-7 layers, UB <= 0.595.
- Almost every scorer falls when the fallback candidates are removed, |co-expr| most of all (lung 0.542 -> 0.516). The exception is
  kidney_ckpt_brain, where the cosine minimum rises 0.501 -> 0.508. Part of the all-candidate ChIP AUROC is therefore carried by
  unmatched target properties (detection / in-degree), for both scorers.
- [verifier addition] Because |co-expr| falls more than scGPT, the comparator ordering depends on matching:
  - lung: the delta vs extraction cells goes from -0.014..-0.004 (all candidates) to +0.005..+0.014 (level 0; CI contains 0, UB <= 0.030),
    and the delta vs the 3k tier goes from -0.030..-0.021 to +0.0004..+0.0097 (UB <= 0.027).
  - lung_external: in the level-0 block the delta vs extraction cells is +0.019..+0.024, with the CI lower bound > 0 at 8/12 layers
    (L0-L7; unadjusted) and UB <= 0.040.

4. WHAT CHANGES IN THE E4a TEXT (E7 sensitivity; no pass/fail)
- The qualitative conclusion is unchanged. With ChIP positives the effect is small, and its CI excludes 0.5 in the lung pair (12/12 layers in
  both). AUROC >= 0.60 is excluded in the lung pair (UB <= 0.561), and scGPT shows no >= 0.05 advantage over |co-expression| in the lung pair
  (delta UB <= 0.018 lung, <= 0.047 lung_external, extraction cells).
- The numbers to replace:
  - lung 0.525-0.534 -> 0.528-0.538; UB 0.552 -> 0.559; |co-expr| 0.536 / 0.551 -> 0.542 / 0.559; delta -0.011..-0.003 -> -0.014..-0.004.
  - lung_external 0.530-0.537 -> 0.533-0.541; UB 0.555 -> 0.561; delta +0.017..+0.024 -> +0.018..+0.026; delta UB 0.042 -> 0.047.
  - all contexts 0.50-0.58 -> 0.50-0.59.
- The old E4a_edges/SUMMARY 'CI lower bound 0.505-0.512' for lung becomes 0.506-0.514.
- The lung_external delta UB moves closer to the 0.05 margin (0.047). The match-level-0 block puts it at 0.040.
- [verifier correction] The implementer listed this as new: in lung the scGPT delta vs the 3,000-cell |co-expression| is negative
  (UB < 0 at 11/12 layers, max 0.0001). It is neither new nor robust.
  - Not new: the OLD sens_chip.csv already had UB < 0 at 10/12 layers (max 0.0008).
  - Not robust: it reverses in the match-level-0 block (+0.0004..+0.0097, UB <= 0.027), so it is carried by the unmatched fallback
    candidates.
  - If mentioned at all, say it is matching-dependent. Do not write that co-expression outperforms the embedding on ChIP positives.

5. LIMITATIONS
- Fixing the exclusion forced a re-randomisation of about 25% of candidates (lung) once the RNG stream diverged. Old and new candidate sets are
  therefore not a pure deletion. The 'old minus contaminated' points bound this: the redraw contributes <= 0.0011 AUROC in lung and <= 0.003
  in any context. [verifier correction]
- The corrected benchmark is less well matched (level-0 share 0.51 -> 0.43 in lung), so part of the signal reflects target properties. The
  match-level-0 block is reported for this reason; its positives keep TFs with >= 1 level-0 candidate (lung loses 3 positives).
- B = 2,000 in all contexts, where E4a used 1,000 in kidney_B / ckpt. The deltas share resamples within a context only.
- The old sens_chip.csv rows outside the lung pair used the float-tie-broken comparator; the new rows are all tie-corrected. Per the verifier's
  audit this shifts the comparator by <= 0.0011.
- Candidate non-edges remain 'ChIP-negative candidate pairs', not confirmed negatives. Level-D ChIP pairs are noisy binding evidence. The contexts
  are not independent: the three ckpt contexts share pairs and cells, and the immune subsets share cells.
- work/ (825 MB of per-pair score checkpoints) can be regenerated with `python e4a_chipfix.py sens <ctx>` / `repro <ctx>`.
- tracking/prompt.md was not updated (this task was restricted to E4a_chipfix/).
```
