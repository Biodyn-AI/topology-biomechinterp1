# Round-2 rebuild — final results

Convention: `CONVENTION.md`. Scripts: `rebuild_headline.py`, `sensitivity.py`, `run_controls.py`,
`run_trees.py`, `regen_audit.py`. Logs and CSVs alongside.

## Verdict on every headline claim

| # | claim | manuscript | round 2 | verdict |
|---|---|---|---|---|
| 1 | compression "effective rank" | 23.6 -> 1.6 (14.4x) | effective rank 383 -> 133 (**2.9x**) | **mislabelled**: 14.4x is a participation ratio on a different gene set |
| 2 | participation ratio | 58 -> 9.5 (6.1x) | 59.5 -> 10.3 (5.8x) | reproduces |
| 3 | SV1 variance fraction | 53.7% -> 93.4% | 57.2% -> 93.2% | reproduces |
| 4 | compression is learned | shuffle raises rank 17.6x | 133 -> 361 under feature shuffle | reproduces |
| 5 | SV1 subcellular localisation | OR 6.37, p = 2.6e-4 @L11 | **OR 0.46, p = 0.13** (reversed); best layer L8 OR 2.41, p = 0.031 | **FAILS** |
| 6 | SV2 STRING co-pole | 0.226 vs 0.124, 12/12 layers | 0.246 vs 0.141; 8/12 at p<0.05, **3/12 Bonferroni** | partially survives |
| 7 | SV3 STRING co-pole | 0.198, 12/12 | 0.297 best; 11/12 at p<0.05, **9/12 Bonferroni** | **SURVIVES — strongest axis** |
| 8 | SV4 STRING co-pole | 0.158, 7/12 | 3/12 at p<0.05, 2/12 Bonferroni | weak |
| 9 | SV2 TRRUST co-pole | 0.206 vs 0.122, survives Bonferroni | 0.220 vs 0.141, best p = 1.2e-2, **0/12 Bonferroni** | **FAILS Bonferroni** |
| 10 | TF/target AUROC SV2-SV7 | **0.744** (max 0.789 @L3) | **0.612** (max 0.640 @L11), null 0.496 | **halves, but survives** (Stouffer Z = 8.69, p = 1.8e-18) |
| 11 | SV5-SV7 regulatory residual | r_rb 0.061, p = 0.080 (L0 only) | all 12 layers: r_rb +0.038 @L0 (p = 0.359) -> **-0.073**; **0/12 significant**, Stouffer Z = -1.01, p = 0.84 | **FAILS — null across depth** |
| 12 | SV2-SV7 regulatory residual | not reported | r_rb +0.054 to +0.130, 8/12 layers significant, Stouffer Z = 6.80, p = 5.3e-12 | **new positive result** |
| 13 | GRN utility, full-512 | AUROC 0.860 vs co-expression 0.793 | **0.591 [0.499, 0.684]** vs **0.556 [0.468, 0.641]**; paired delta +0.034 [-0.073, +0.144] | **FAILS — does not separate** |
| 14 | ESM2 TF/target | 0.889, different gene set, projected | **0.924** un-projected on the matched set | **survives and strengthens** |
| 15 | ESM2 lacks SV1/SV2: artefact? | "untested rather than settled" | matched control run: projection does not destroy structure (scGPT 512d z=+2.00 -> 102d z=+2.72) | **SETTLED — biological, not artefact** |
| 16 | kidney regulatory transfer | confounded, "cannot be separated" | matched TF set (n=59 both sides): AUROC 0.49-0.55 | **SETTLED — genuinely at chance** |

## The GRN benchmark, with the baselines round 1 never ran

Held-out 20% TRRUST split: 49 positives / 147 negatives. B = 2,000 bootstrap.

| method | AUROC | 95% CI | vs \|co-expression\| |
|---|---|---|---|
| scGPT full-512 (L0) | 0.591 | [0.499, 0.684] | +0.034 [-0.073, +0.144] |
| GENIE3 | 0.590 | [0.500, 0.680] | +0.034 [-0.065, +0.136] |
| scGPT spectral SV2-SV7 | 0.589 | [0.495, 0.669] | +0.032 [-0.078, +0.142] |
| \|co-expression\| | 0.556 | [0.468, 0.641] | -- |
| GRNBoost2 | 0.555 | [0.465, 0.643] | -0.001 [-0.103, +0.105] |
| scGPT spectral SV5-SV7 | 0.541 | [0.451, 0.630] | -0.014 [-0.130, +0.100] |
| random | 0.513 | [0.414, 0.605] | -0.040 [-0.170, +0.106] |

**No method separates from any other, and none separates from co-expression.** With 49 positives the
benchmark cannot resolve differences of this size; the honest statement is that it is underpowered,
not that scGPT wins. GENIE3 and GRNBoost2 — the field-standard comparators — do no better.

## Multiplicity audit, rebuilt

| | round 1 | round 2 |
|---|---|---|
| hypotheses | 183 | 183 |
| extractable p | 112 | 112 |
| p recorded as exactly 0.0 | 18 | 0 (replaced by 1/(N+1)) |
| BH q<0.05 passes | 70 | 70 |
| **Bonferroni passes** | **43** | **25** |
| **headline survivors** | **8 of 11** | **5 of 11** |

Flipped to failure: SV1 secretory-pathway layering, SV2 TRRUST co-pole, TF/target cross-seed
validation (a 50-permutation test whose p-floor is 2.0e-2).

## What the paper now says

Three claims die: **subcellular localisation** (5), the **SV5-SV7 regulatory residual** (11), and
**functional utility on regulatory edges** (13). Two strengthen: **sequence-derivability of TF
identity** (14, 15) and the **failure of regulatory geometry to transfer across tissue** (16). The
protein-interaction geometry (7) is the robust positive result, and it is SV3, not SV2, that carries
it.

The thesis is intact and cleaner than before. What is fixed below the level of transcription —
whether a gene's product is a transcription factor — is read best by a model that sees only protein
sequence (ESM2 0.924 vs scGPT 0.650 on identical genes). What is relational — which factor acts on
which target — is recovered by **nothing**: not the residual stream, not its spectral summary, not
co-expression, and not GENIE3 or GRNBoost2. That is a stronger and more useful negative result than
the one the current draft reports, because it now rests on the standard baselines rather than on
their absence.


## Addendum: completed confound analysis (supersedes the partial run)

The first controls run was killed at layer 7, leaving SV5-SV7 with 8 of 12 layers, and
`make_figures.py` silently fell back to a stale CSV from an aborted run -- so an earlier draft of
Figure 2b plotted SV2-SV7 twice under two different labels. `run_confounds_sv57.py` re-ran both
subspaces to completion against the **same** negative-pair sample, which is what makes the contrast
between them interpretable:

| subspace | r_rb range | layers p<0.05 | Stouffer |
|---|---|---|---|
| SV2-SV7 | +0.054 to +0.130 | **8 / 12** | Z = 6.80, p = 5.3e-12 |
| SV5-SV7 | +0.038 to -0.073 | **0 / 12** | Z = -1.01, p = 0.84 |

The SV5-SV7 retraction is therefore firmer than the partial run suggested: not merely
non-significant, but null with a slightly negative central tendency from layer 1 onward.
