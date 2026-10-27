# E1 — Compression of the observed-gene representation (sensitivity analysis; CONVENTION3 §5)

Written by the orchestrator from the implementer's and verifier's returned findings (subagents cannot
write .md files). CONVENTION3 hash 3e3249d1…a05c01 recorded in `results.json` and verified.

## What was run
All 16 contexts, L0–L11, on the mean-centred observed-gene matrix V × 512 (`stats3.svd_observed`).
Effective rank (ER, exp of the entropy of normalised S²), participation ratio (PR), SV1 variance share,
cumulative SV1–10. Gap-ratio bootstrap B = 200 (genes with replacement); feature shuffle (each column
permuted independently across genes), 20 draws; detection-matched subsampling n = 1,300 on common pooled
deciles × 200 draws; zero-row (full-matrix) reproduction; detection tracking; n_obs threshold sensitivity.
Seed 20261026. Scripts: `e1_compression.py`, `e1_finalize.py`.

## Direct evidence (file provenance)
Main-text table (`core_table.csv`), order lung | lung_external | kidney_B | immune | kidney_A:

| | lung | lung_external | kidney_B | immune | kidney_A |
|---|---|---|---|---|---|
| ER L0 | 457.4 | 459.4 | 456.2 | 442.4 | 437.4 |
| ER L11 | 324.1 | 340.8 | 324.1 | 309.6 | 297.1 |
| ER fold | 1.41 | 1.35 | 1.41 | 1.43 | 1.47 |
| PR L0 | 86.4 | 87.8 | 88.5 | 90.1 | 64.8 |
| PR L11 | 12.3 | 17.8 | 14.8 | 15.2 | 10.7 |
| SV1 var L0 | 7.3% | 7.3% | 7.0% | 5.6% | 8.1% |
| SV1 var L11 | 24.0% | 17.9% | 17.6% | 18.6% | 23.1% |

- ER decreases strictly at every layer step in all 16 contexts (ρ(ER, layer) = −1.000). PR strictly decreasing
  in lung and kidney_B; one small rise in lung_external (L10→L11), immune and kidney_A (L0→L1, ≤ 1.1)
  (`depth_trend.csv`). Non-core exception kidney_ckpt_brain (ρ = −0.909).
- Feature shuffle: shuffled ER at L11 467.5–479.8 (1.40–1.57× observed); at L0 1.09–1.12× (`shuffle.csv`).
- Detection-matched (n = 1,300): ER L0→L11 lung 432.0→296.3; lung_external 435.7→318.3; kidney_B 433.1→302.8;
  immune 431.4→299.5; kidney_A 427.4→290.0 (fold 1.37–1.47; sd 0.3–1.1) (`matched.csv`).
- Zero-row artefact: full kidney_A matrix (2,831 zero rows) gives ER 383.4→133.0 (2.9-fold), SV1 57.2%→93.2%,
  shuffled ER 360.9; |r(SV1, observed indicator)| 0.998/0.9998. Immune full matrix: 93.2% at L11
  (`zero_row_artefact.csv`).
- Gap ratios: 429/2,112 adjacent pairs (≈20%) near-degenerate by the basic-CI rule (4% at k = 1, 3% at k = 2,
  12–13% at k = 3–4, 30–39% at k = 8–11). Core flags among SV1–SV4 (all 9 confirmed at B = 1,000 by the
  verifier, `verify/gap_B1000.csv`): immune SV1/SV2 at L0/L1/L2/L10 and SV3/SV4 at L5; kidney_A SV2/SV3 at L11;
  kidney_B SV2/SV3 at L7; lung SV3/SV4 at L11; lung_external SV3/SV4 at L0.
- Detection tracking (`detection_corr.csv`, verifier `verify/detection_by_layer.csv`): median |ρ(SV1, log n_obs)|
  0.51–0.62 at L0–L10 and SV1 is the most detection-correlated direction in 13/16 contexts at L0; **at L11
  the median for SV1 falls to 0.23 and SV2 (0.56) becomes the top detection direction in 10/16 contexts**
  (lung SV1 0.38, immune 0.20, kidney_B 0.09; kidney_A SV1 0.75 is the exception).
- n_obs ≥ 3 restriction (lung): L11 SV1 share 24.0%→13.7%, PR 12.3→20.5, ER 324→335; other core contexts move
  ≤ 4 points / ≤ 10 ER (`detection_threshold.csv`).
- Checkpoints on the same kidney cells: kidney-fine-tuned 474.9→190.9 (2.49-fold; PR 13-fold); whole-human
  1.40-fold; brain 1.64-fold.
- Exactly one singular value < 1e-3·S1 in all 192 context-layers; aligned with the uniform vector (|cos| 0.999)
  where checked (5 core contexts, L0 and L11). ER ceiling 511.

## Verification
Independent recomputation (no stats3 import): 115/115 checks pass to < 1e-12; stochastic outputs reproduced
under a different seed; data integrity (non-zero rows identical across layers, equal to obs_idx) in all 16.

## Interpretation
- Direct: ER falls ≈1.4-fold, PR ≈6-fold, SV1 share roughly triples with depth; feature-marginal shuffle does not
  reproduce it.
- Inference: the "2.9-fold" ER fall and "93.2% on one direction" are properties of the zero-row matrix, not of the
  representation. Depth concentrates variance into the leading ~10 directions; it is not a collapse to low rank.
- Hypothesis (untested): low-detection genes load on the leading directions because their mean-pooled vectors
  average fewer cells; the kidney checkpoint's narrow fine-tuning makes late layers more anisotropic; the null
  direction reflects a LayerNorm-type constraint.

## Caveats
Gene bootstrap is biased for ER (duplicates), so no ER bootstrap CI; ER sampling sd from matched draws.
Near-degeneracy rule is a heuristic (sorted bootstrap ratios ≥ 1), and non-core flags are Monte-Carlo-sensitive
at B = 200 (≈20% count is approximate). Signs of SV-indicator correlation are fixed by stats3's sign rule; |r|
reported. Contexts share donors and datasets; absolute ER comparable across contexts only in `matched.csv`.
