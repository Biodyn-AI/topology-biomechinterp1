# E2 — H1 network neighbourhood (confirmatory) and interaction sensitivity analyses

Written by the orchestrator from the implementer's and verifier's returned findings, with verifier Corrections A–D
applied. CONVENTION3 hash verified in `results.json`; stats3.assert_resolution called for every family test.

## Primary test (CONVENTION3 §4, H1) — SUPPORTED
Scalar assortativity of SVk (k = 1–10) × 12 layers over STRING-physical ≥ 700 edges within V (lung 14,474 edges,
lung_external 14,132, kidney_B 12,434, immune 1,251, kidney_A 1,390), Maslov–Sneppen null (burn-in 10m, thinning m,
joint draws), max-T, N = 10,000 (lung, lung_external) / 2,000 (others).
- Family min p_adj = 1e-4 (floor) in lung and lung_external; 120/120 cells p_adj < 0.05 and r > 0 in both
  (`h1_summary.csv`). Null-max |z| 95th pct 3.41 (lung) < smallest observed |z| 4.25.
- Top cells: lung L3 SV1 r = 0.435 (null mean −0.001, sd 0.008, z 53.1); lung_external L8 SV2 r = 0.406 (z 49.2);
  median r over 120 cells 0.181 / 0.170. Secondary: kidney_B L9 SV3 0.404; immune L4 SV6 0.566 (118/120);
  kidney_A L5 SV2 0.702 (all at their N = 2,000 floors).
- Intervals: pre-registered √(m/n) vertex-subsampling CIs are ≈ 0.71× too narrow (CONVENTION3_ERRATA item 3).
  Delete-d intervals (verifier `verify/verify_c.json`): lung [0.388, 0.487], lung_external [0.343, 0.462],
  kidney_B [0.351, 0.454], immune [0.465, 0.688], kidney_A [0.602, 0.809].
- Replication: lung top-3 cells each matched in lung_external by right-singular-vector |cos| (random 0.035): L3 SV1 ↔
  L3 SV1 (0.979; gene-space r 0.985 on 4,535 genes; ext r 0.407); L6 SV2 ↔ L6 SV2 (0.946; 0.358); L8 SV3 ↔ L8 SV2
  (0.912; 0.406). Same sign, intervals exclude null under either scaling (`h1_direction_match.csv`).
- Chain diagnostics: acceptance 0.96 (lung, lung_external, kidney_B), 0.93 immune, 0.90 kidney_A; lag-1
  autocorrelation ≈ 0.13 (effective N ≈ 0.77N); co-pole null mean vs configuration expectation within 0.27–0.39%
  (large graphs), 1.9% (kidney_A).

## Sensitivity analyses
- Other nulls (`sens_nulls_summary.csv`), lung top cell: label −0.001 (z 53.1); within-V degree-stratified 0.005
  (z 50.4); genome-wide STRING-degree-stratified 0.021 (z 41.0). Family min p_adj at floor in all 5 contexts.
- Hub exclusion (`sens_hub_exclusion.csv`; corrected scope): removing the top 5% / 10% by degree RAISES r in the three
  large graphs (lung 0.435 → 0.515 / 0.511; lung_external 0.406 → 0.518 / 0.487; kidney_B 0.404 → 0.499 / 0.488) and
  LOWERS it in the two small graphs (immune 0.566 → 0.537 / 0.528; kidney_A 0.702 → 0.569 / 0.512), all z 6.7–11.9.
  Leave-one-gene-out: max influence 0.0045 vs jackknife SE 0.026 (lung). Top-influence lung genes are mostly hubs
  (8/10 top-5%: FN1, CCL11, CXCR4, DCN, CCR2, CXCL12, CCR1, CCR5; not hubs: CXCR3, CXCL1).
- Detection-decile null: 117/120 (lung), 120/120 (lung_external) cells significant; the 3 lung failures are SV1 at
  L9–L11. kidney_A top cell detection-null mean 0.434 vs observed 0.702 (excess 0.268, z 8.7).
- Co-expression-matched non-edges (no physical link at any score; deciles of |co-expression| on extraction cells;
  `sens_coexpr_matched.csv`): lung cosine (L3, full-dim) AUROC 0.725 [0.711, 0.738], top-cell −|Δx| 0.631, residual
  |co-expression| 0.560; lung_external 0.679, 0.617, 0.534. 50-bin version (post hoc, `sens_coexpr_fine.csv`):
  residual co-expression 0.53–0.56, cosine 0.69–0.74; 3,000-cell tier within 0.02.
- Co-pole at K/n ∈ {0.10, 0.20, 0.27, 0.33, 0.50} (`sens_copole_K.csv`), observed / MS-null ratio, lung L3 SV1: 2.27,
  1.89, 1.70, 1.60, 1.32; lung_external 2.10, 1.85, 1.70, 1.56, 1.30. **Correction A:** label-null co-pole null means
  are close to MS in the large graphs but not in kidney_A: L5 SV2 MS null mean 0.284 vs label 0.146 (ratio 3.93 →
  2.02); L2 SV2 4.23 → 2.09; L11 SV1 2.78 → 1.58; immune L4 SV6 2.37 → 1.88; kidney_B L11 SV1 label 0.79 (z −3.9) vs MS
  1.40 (z 15.8). Degree concentration inflated the label-null co-pole effect size ≈ 2× in kidney_A; significance
  survives (MS z ≈ 29). Label-null co-pole SD ≈ 3× MS SD.
- Co-complex pairs, MS N = 2,000 (`sens_complexes_summary.csv`), lung L3 SV1: ComplexPortal 1,934 pairs (87% already
  STRING) r 0.588 (z 26.1); ComplexPortal minus STRING (246) r 0.423 (z 7.0); hu.MAP 2.0 minus STRING (5,159) r 0.196
  (z 14.1; 95/120 cells). lung_external L8 SV2: 0.524, 0.278, 0.172.
- Robustness contexts (N = 1,000; `sens_robustness_summary.csv`): all 11 at floor; immune draws/donors/cell types top
  r 0.50–0.59 (100–117/120 cells); checkpoints whole-human 0.383, kidney 0.418, brain 0.302 (120/120 each).
- kidney_A bridge (`bridge_kidneyA.csv`): S0 screen value reproduced exactly (L5 "SV3" co-pole 0.333 vs 0.141,
  z 7.5); S1 observed SVD: identical value at observed L5 SV2 (zero-row SVD shifts low indices by one; verified for
  k = 1–5, directions mix at SV6–SV8 for L5–L7); S2 universe V 0.656 (z 15.7; label null, partly degree-inflated, see
  Correction A); S3 physical STRING 0.614 (z 25.3); S4 assortativity r 0.702 (36/36); S5 MS z 28.2 (36/36 at floor).
  Off-path: physical STRING on the 187 legacy genes (163 pairs) best z 3.9, max-T over 36 p = 0.017.

## Interpretation
- Direct: physical partners have more similar spectral coordinates than degree-preserving rewiring allows, in every
  extraction and checkpoint; the effect is not produced by degree, detection (outside kidney_A and late lung SV1) or
  pairwise co-expression.
- Inference: the screen's co-pole statistic under the label shuffle was exposed to degree (≈ 2× effect inflation in
  kidney_A) but its conclusion survives a degree-preserving null.
- Hypothesis (untested): partners share coordinates because they share cellular programmes (the lung L3 SV1 axis runs
  from secreted/signalling genes to ubiquitous nuclear machinery — also the compartment axis of E11); H1 does not show
  that the model encodes binding.

## Not shown / limitations
SV indices not comparable across contexts; no paired cross-context test; V-restriction (n_obs ≥ 3/5/10) not run;
co-expression matching pairwise only; p at floor everywhere (use z, CIs); chain mildly autocorrelated, starts from the
observed graph; vertex subsampling assumes exchangeable nodes; compartment co-membership not conditioned on (follow-up).
