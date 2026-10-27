# Provenance of the SV2–4, SV5–7 and joint SV2–7 bands

Sources: `iterations/iter_XXXX/executor_hypothesis_screen.json` (method and result strings), the
`executor_iteration_report.md`, `executor_next_steps.md` and `brainstormer_*.md` files of the same
iteration, the iteration scripts where they exist (iterations 45–63 left no script), and this folder's
`S2_hypotheses.csv`, `crosswalk.csv` and `legacy_band_mapping.json`. Quotes are verbatim from those files.

**Bottom line.** None of the three bands was prespecified. Each was chosen after the screen had seen
results:
- **SV2–4** was put together from single-axis PPI tests on the kidney matrix. One axis was added at a
  time, and the cut-off after SV4 rested on a knife-edge. The brainstormer's stated criterion was
  "z > 3 in ≥6/12 layers" (iter_0014). SV4 met it at exactly 6/12 layers, its layer-1 z being
  3.00003. SV5 failed it at 2/12. By the looser emp_p < 0.05 count, both pass: SV4 at 7/12, SV5 at
  6/12. The executor also evaluated the criterion on 3,092 pairs at STRING ≥ 0.4, whereas the brief had
  specified 1,022 pairs at ≥ 0.7. (Corrected by the verifier: an earlier version of this file said SV5
  met the same criterion as SV4.)
- **SV5–7** is simply "the next three" directions, tried because SV2–4 had a signal.
- **SV2–7** is the union of the two, formed after their layer profiles were seen to complement each
  other. Its success threshold was then relaxed to fit the result.

When the screen itself scanned every three-direction window with a Bonferroni correction
(iter_0063), no window survived.

## Matrices and conventions involved

| phase | iterations | matrix (true tissue) | gene rows | SVD convention |
|---|---|---|---|---|
| A | 8–26 | `cycle1_main` = **kidney_A** (TS kidney; called "immune" in the preprint) | 209 edge-dataset genes at their true rows; **22 of them are never-observed all-zero rows** | fitted on the 209-row subset, mean-centred over it, zero rows included |
| A′ | 29–40 | kidney_A | labels on **rows 0–208** (gene_list.txt read with blank lines dropped); "195 in-vocab" = those rows minus 14 zero rows | fitted on the misaligned subset |
| B | 47–58, 61–63 | `cycle4_immune_main` = **immune** (TS immune) | 2,039 observed (non-zero) rows | observed rows only, mean-centred over them; identical to `stats3.svd_observed` (crosswalk \|cos\| = 1.000 for SV1–SV12 at all layers) |
| B′ | 59–60 | kidney_A, **summed** embeddings (`layer_gene_embedding_sum.npy`) | 1,972 observed rows | observed rows of the sum matrix; poorly aligned with the mean-embedding basis (layer-averaged mean cos² of SV5–7 with the best observed 3-band = 0.21) |

## Timeline

| step | iteration | matrix / gene rows / convention | what happened | result that motivated the next step |
|---|---|---|---|---|
| 1 | iter_0008 H01 | kidney_A, 209 true rows, subset-fitted | first use of SV2 and SV3: "SVD of layer-11 mean-centered 209-gene embedding. Project onto SV2 and SV3." | SV2 IL-4 term emp_p=0.000; SV3 emp_p=0.088 |
| 2 | iter_0010 H02 | same | TRRUST co-pole on SV1 and SV2 at L11 (326 TRRUST rows, i.e. 288 distinct pairs) | "SV1: obs=0.163 vs null=0.123, emp_p=0.016; SV2: obs=0.206 vs null=0.122, emp_p=0.000" |
| 3 | iter_0013 H01 | same | STRING (score ≥ 0.7, 1,022 API pairs) co-pole on SV1, SV2 and SV3 | "SV1 z=0.39 (2/12 sig), SV2 z=10.04 (12/12 sig), SV3 z=7.18 (12/12 sig)" |
| 4 | iter_0014 brief | — | brainstormer proposes extending to the next axes: "H01: SV4/SV5 PPI co-pole test", with the success criterion "SV4 or SV5 z > 3 in ≥6/12 layers" | — |
| 5 | iter_0015 H01 | same, but 3,092 pairs at STRING ≥ 0.4, not the brief's 1,022 at ≥ 0.7 | SV4 and SV5 tested | "SV4 mean_z=2.21, sig=7/12; SV5 mean_z=1.65, sig=6/12" ("sig" = emp_p < 0.05). Against the stated criterion "z > 3 in ≥6/12" (per-layer z in `iter_0015/h01_sv45_ppi_copole.json`), SV4 passes at 6/12 with L1 z = 3.00003, and SV5 fails at 2/12 |
| 6 | iter_0015 brief | — | SV5 dropped: "SV5 (deprioritized; rescue only if SV6+ shows new pattern)" | SV4 kept by a z-margin of 3×10⁻⁵ at one layer; the group becomes SV2, SV3, SV4 |
| 7 | iter_0016 H02, iter_0019 H02, iter_0020 | same | "SV2/SV3/SV4 axis independence"; "multi-axis co-polarity count (SV2+SV3+SV4)" | the SV2–SV4 triple becomes a unit, defined by PPI co-pole on the kidney 209-gene subset |
| 8 | iter_0047 H01/H02 | immune, 2,039 observed | the screen finds the zero rows: "Found 2902 genes have zero norm at ALL layers", and that the earlier "93.4%" SV1 figure "was dominated by zero vectors" | switch to observed rows only |
| 9 | iter_0048 brief | immune, observed | the kidney band is carried over unchanged: "H02 (HIGH VALUE): SV2-SV4 PPI Geometry in 2039-Gene Active Set", comparing with "prior rho values (SV2: 1.0, SV3: 0.9, SV4: 0.7)" from the kidney matrix | — |
| 10 | iter_0049 H03 | immune, observed; coordinates = left singular vectors | first band-level test of SV2–SV4 on TRRUST pairs | "Best: L8 p<0.0001 effect=0.268; significant (p<0.05) at 8/12 layers" |
| 11 | iter_0049 next steps | — | **origin of SV5–7**: "Q: Is the effect specific to SV2-SV4, or does it extend further into the spectrum? Test SV5-SV7, SV8-SV10 subspaces" | fixed-width windows of three, chosen after the SV2–4 result |
| 12 | iter_0050 H03 | immune, observed | first recorded SV5–7 (with SV1, SV2–4, SV8–10) | "SV2-4: rbc=+0.229 (10/12 layers p<0.05); SV5-7: rbc=+0.152 (7/12 layers); SV8-10: rbc=+0.123 (5/12 layers)" |
| 13 | iter_0050 H01, iter_0051 H01 | immune, observed | SV2–4 TRRUST proximity tested against "co-expression" (in practice the cosine of the full 512-d embeddings, not expression correlation) | "0/12 layers specific"; residualised: "1/12 layers significant (L8: rbc=0.083, p=0.0016)" (recomputed exactly here: claims_recompute.json, `residualised_immune_iter51_52`) |
| 14 | iter_0051 brainstormer | — | SV5–7 promoted as the alternative: "SV5-7 regulatory signal (independent of co-expression) \| UNTESTED \| High potential" | — |
| 15 | iter_0052 H03 | immune, observed | SV5–7 residualised | "rbc_resid=0.148 (L0), 0.119 (L1), 0.083 (L2); all p<0.002. L3+ lose significance". Recomputed exactly (0.148 / 0.119 / 0.083); layers 5–11 have *negative* rbc |
| 16 | iter_0053 H01 | immune, observed | other bands scanned only after SV5–7 had been chosen: "SV2-4, SV5-7, SV8-10, SV11-14, SV15-20, SV21-30" at L0/L1/L2/L8 | "SV5-7 dominant at L0-L2; SV2-4 dominant at L8" |
| 17 | iter_0055 H02 | immune, observed | layer profiles of SV5–7 and SV2–4 classifiers | "Complementary encoding: SV5-7 high early+late, SV2-4 high mid-depth" |
| 18 | iter_0055 brainstormer | — | **origin of SV2–7**: "can we dissect the joint SV2-7 space to produce a layer-stable classifier?"; decision threshold "joint AUROC ≥ 0.72 at 10+ of 12 layers" | — |
| 19 | iter_0056 H01 | immune, observed | joint 6-D classifier | "joint >= 0.72 at 9/12 layers". The threshold was then relaxed: "Lower threshold for criterion since 9/12 vs 10/12 is minor; note claim should state 'most layers'" |
| 20 | iter_0063 H63-B | immune, observed | the screen's own systematic scan of SV(k, k+2), k = 0…9, at L0 with Bonferroni over 10 windows | "Best raw AUROC: SV8-10=0.619 (perm_p=0.67); SV5-7=0.611 (perm_p=0.28); no window survives Bonferroni" |

## What this means for the round-3 analysis

- **The SV2–4 group was defined on another tissue, with a defective gene set, on another task.** It
  was defined on kidney_A, on a 209-gene subset that included 22 all-zero rows, for PPI co-location.
  It was then reused unchanged on immune for TF/target tasks. Measured against the observed-gene basis
  of kidney_A, the kidney 209-subset SV2–4 matches observed **SV1–3** best: layer-averaged mean cos² is
  0.49, against 0.45 for observed SV2–4. Observed SV1–3 wins at 8 of 12 layers
  (`legacy_band_mapping.json`). Even its best match is weakly aligned.
- **The submitted manuscript's convention moves every band down one index.** That convention is the
  full 4,803-row matrix with its zero rows (D2). Under it, legacy SVk is observed SV(k−1):
  - for k = 2…12 the best match is observed SV(k−1) in 94.7% of (layer, k) cases, with layer-averaged
    best-match |cos| between 0.945 and 0.998 (`crosswalk.csv`, `results.json`);
  - legacy SV2–4 maps to observed **SV1–3** (mean cos² 0.99, all 12 layers);
  - SV5–7 maps to **SV4–6** (0.94, all 12 layers);
  - SV2–7 maps to **SV1–6** (0.97, all 12 layers; 0.82 for the identity band SV2–7).
- **On immune the screen's bands are the observed-basis bands.** The screen's convention there
  ("centered nonzero embeddings [2039,512]") equals `stats3.svd_observed`, so SV2–4, SV5–7 and SV2–7
  map to themselves (mean cos² = 1.000).
- **The bands need a selection correction.** They were chosen after the results were seen, so any
  round-3 test of them should be reported against a scan of all bands and layers with a max-statistic
  correction (E3), not as a confirmatory test of a prespecified hypothesis. The mapping in
  `legacy_band_mapping.json` was computed from geometry alone; no outcome data were read. It is the
  fixed translation to use.
