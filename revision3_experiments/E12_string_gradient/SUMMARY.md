# E12 — STRING confidence gradient (current convention)

Sensitivity analysis added during write-up; no pass/fail. CONVENTION3 hash unchanged (3e3249d1…a05c01).
Script: `e12_gradient.py` (one process per context); output: `gradient.csv` (per context x layer), `gradient_<ctx>.csv`, `e12_<ctx>.log`.

## Question
The screen reported that co-location of STRING partners rises with STRING confidence (iter_0015/H02: co-pole z over five
score quintiles on one direction of the kidney_A 209-gene subset, Spearman rho = 1, reported p = 1.4e-24, exact p = 0.0167).
Does a confidence gradient exist on the observed genes of each extraction, for the full-dimensional cosine, beyond endpoint
degree, detection and co-expression?

## Design
- Links: STRING v12 physical subnetwork, combined score >= 150, both ends in V (lung 132,077; lung_external 127,663;
  kidney_B 112,823; immune 7,803; kidney_A 7,938).
- Per layer: centred (over V), L2-normalised 512-d vectors; cosine of each link.
- rho = Spearman(score, cosine). rho_adj = partial Spearman controlling for rank(log deg_i + log deg_j) (degree in the
  >= 150 network within V), rank(log1p n_obs_i + log1p n_obs_j), rank(|co-expression|) (Pearson on <= 3,000 cells of
  the population, `cache/<ctx>/expr.npy`).
- Bands 150-399 / 400-699 / 700-899 / 900-999: AUROC of link cosine vs one degree-matched non-linked pair per link
  (endpoint degree deciles matched; no physical link at any score), drawn once (seed 20261026).
- Intervals: vertex subsampling, m = n/2, delete-d scaling (= 1 at m = n/2; ERRATA item 3), B = 500; subsample
  statistics use the full-sample ranks restricted to retained links (Pearson of ranks; partial by residualisation on the
  restricted covariate ranks). Check: the first two lung layers computed with full re-ranking in every subsample gave
  identical point estimates and intervals within 0.001 (rho 0.098 [0.081, 0.116] at L0 both ways).
- Descriptive: Spearman(score, -|x_i - x_j|) for SV1-SV3 coordinates.

## Direct evidence
| context | rho (12 layers) | rho_adj | min lower bound of rho_adj | AUROC by band, L0 | AUROC by band, L11 | layers monotone over 4 bands |
|---|---|---|---|---|---|---|
| lung | 0.050-0.102 | 0.051-0.096 | 0.036 | 0.569 / 0.597 / 0.642 / 0.692 | 0.545 / 0.558 / 0.584 / 0.630 | 12/12 |
| lung_external | 0.051-0.095 | 0.060-0.099 | 0.044 | 0.568 / 0.590 / 0.643 / 0.683 | 0.550 / 0.561 / 0.596 / 0.629 | 12/12 |
| kidney_B | 0.043-0.099 | 0.046-0.098 | 0.028 | 0.556 / 0.584 / 0.639 / 0.694 | 0.544 / 0.551 / 0.587 / 0.633 | 12/12 |
| immune | 0.092-0.153 | 0.055-0.115 | 0.006 | 0.624 / 0.637 / 0.736 / 0.767 | 0.590 / 0.585 / 0.636 / 0.691 | 8/12 |
| kidney_A | 0.085-0.188 | 0.115-0.190 | 0.050 | 0.661 / 0.659 / 0.749 / 0.820 | 0.647 / 0.613 / 0.699 / 0.761 | 0/12 |

Band sizes in lung: 91,165 / 26,438 / 6,833 / 7,641. Score correlates with degree product (rho 0.09-0.19) and only weakly with
|co-expression| (0.03-0.10). Lung band AUROC ranges over layers: 0.545-0.569, 0.558-0.600, 0.584-0.643, 0.630-0.698.
Mean cosine at lung L3: 0.056 / 0.074 / 0.105 / 0.151. Single-direction rho (SV1-SV3) <= 0.06 in the large networks.

## Inference
- A confidence gradient exists on the observed genes in every extraction and at every layer: higher-confidence physical
  links join more similar vectors, beyond endpoint degree, detection and co-expression. It is monotone over the four
  STRING bands at all 12 layers in the three large networks. In the small immune and kidney_A networks the two lowest
  bands are not separated (non-monotone), but the two highest bands are clearly above them.
- The pair-level association is weak (rho ~ 0.1); the screen's rho = 1 is a five-point summary and its p = 1.4e-24 an artefact.
- It is strongest at early layers and declines with depth, like the other neighbourhood statistics.

## Not shown / limitations
- STRING confidence partly reflects study intensity; the degree adjustment controls this only in part.
- One degree-matched negative draw per link (no resampling of negatives); the band AUROCs carry no intervals.
- Subsample statistics use restricted full-sample ranks (checked against full re-ranking on two layers).
