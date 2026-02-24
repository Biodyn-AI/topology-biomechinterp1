# Brainstormer Structured Feedback — iter_0017

## Gate Status
`passed_min_research_gate: true` — 3 hypotheses tested, 2 positive, 1 negative.

## What iter_0017 achieved

The iteration completes the multi-axis confidence gradient story:
- SV4 quintile gradient: rho=0.900, Q5 z=6.07. Matches SV3 (rho=0.90) and SV2 (rho=1.00).
- Biological orthogonality: GO Jaccard ≤ 0.081 across all axis pairs. Axis-dominant sets encode distinct biological programs: SV2=innate immune/exosome, SV3=membrane/cell-surface, SV4=receptor-binding/apoptosis.
- TRRUST signed regulation on SV3/SV4: cleanly negative. Confirms SV2 is the unique regulatory-polarity axis.

The executor's call—"project is ready for paper consolidation"—is warranted for the current narrative. However, the paper as written has a major gap: **zero cross-model validation**. Every experiment has been on scGPT lung, seed42 (or seed42–44 for topology). The biological claims are strong within-model, but the paper cannot credibly claim this is a general property of biological language models without at least one replication on Geneformer.

## Completed narrative (no further testing needed)

| Claim | Evidence | Status |
|---|---|---|
| kNN clustering coefficient elevated in residual stream | iter_0002, all 12 layers, all seeds | Done |
| H1 persistent homology elevated | iter_0003, 11/12 layers | Done |
| Cross-layer CKA ≈ 1.0 (stable residual stream) | iter_0004 | Done |
| Low intrinsic dimension vs shuffle | iter_0004 | Done |
| SV2/SV3 encode compartment biology | iter_0008–0010 | Done |
| SV2 signed regulatory geometry (activation 12/12) | iter_0011–0012 | Done |
| STRING PPI co-pole SV2+SV3 12/12 layers | iter_0012–0013 | Done |
| Hub-degree confound ruled out | iter_0014 | Done |
| GO co-annotation confound ruled out | iter_0015 | Done |
| Confidence gradient SV2/SV3/SV4 rho ≥ 0.90 | iter_0015–0017 | Done |
| Axis independence (near-zero inter-axis r) | iter_0016 | Done |
| Biological orthogonality (GO Jaccard ≤ 0.081) | iter_0017 | Done |
| SV2 regulatory polarity is axis-specific | iter_0017 | Done |

## Critical gaps (paper vulnerability)

1. **No cross-model replication.** Geneformer SV2/SV3 PPI co-pole has been listed as "highest-priority next experiment" in iters 0010–0016 and never executed. This is the paper's biggest weakness.
2. **No quantitative predictive claim.** Does geometric distance predict STRING score as a continuous variable? A regression R² number would substantially strengthen the story.
3. **No out-of-sample PPI prediction test.** Does the geometric structure generalize to interactions not seen during training? This would be the highest-impact claim.
4. **Single dataset.** All experiments are on scGPT lung. No test on immune-only or a different tissue.
5. **No random initialization control.** No test that a randomly initialized scGPT fails to show PPI geometry.

## Retired / deprioritized directions

- TRRUST signed regulation on SV3/SV4: retired (iter_0017).
- GO BP enrichment in SV2 poles: retired (iter_0011).
- Repression anti-pole hypothesis: retired (iter_0013).
- SV5 extension: weak/inconsistent signal (iter_0015), deprioritized.
- Bootstrap CI on act/rep differential: underpowered, superseded by per-axis empirical p approach.

## Scoring of iter_0017

| Dimension | Score |
|---|---|
| Hypotheses tested | 3/3 planned |
| Positives | 2 |
| Negatives (informative) | 1 |
| Paper update | Yes |
| New open question generated | Yes (cross-model gap) |

Overall: strong consolidation iteration. The negative result (H-E) is clean and informative.
