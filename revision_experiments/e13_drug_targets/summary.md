# E13 Summary — Drug-target candidate ranking

**Status:** Done.
**Code:** `run.py`.
**Outputs:** `results.csv`, `per_target.csv`, `results.json`.

## Headline finding

For drug-target candidate ranking on 8 in-vocabulary drug targets, the spectral PPI subspace ($\SV{2}$--$\SV{4}$) **does not** outperform co-expression. In fact, co-expression Pearson correlation is the best baseline, and scGPT full-512 cosine matches it.

| Method | mean P@5 | mean P@10 | mean P@25 |
|---|---|---|---|
| **Co-expression Pearson** | **0.025** | **0.0375** | 0.035 |
| scGPT full-512 cosine | 0.025 | 0.025 | 0.025 |
| **scGPT $\SV{2}$--$\SV{4}$** | 0.000 | 0.0125 | 0.015 |
| Random | 0.000 | 0.000 | 0.000 |

(8 drug targets evaluated; median 12 known STRING interactors per target.)

## Interpretation

Three points:

1. **The spectral PPI subspace is real but practically modest.** The original paper's claim that "geometric proximity in $\SV{2}$--$\SV{4}$ can serve as an ordinal predictor of protein interaction likelihood" is supported in the *enrichment* sense (Section 2.3) but **fails to translate into a competitive drug-target ranking method** in this specific benchmark. The PPI subspace's P@5 is 0; the full 512-dim cosine gets 0.025; co-expression gets 0.025-0.0375.
2. **The full 512-dim scGPT cosine matches co-expression** but does not exceed it. This is consistent with prior work showing that single-cell foundation model representations primarily encode information that is also recoverable from raw expression statistics.
3. **The benchmark is small** (8 drug targets in vocab + STRING). The denominator of in-vocab + STRING-covered drug targets is small because the immune-relevant HVG vocabulary doesn't include many canonical pharmacology targets. A larger benchmark with broader vocabulary (E3 work) would give a more reliable estimate.

## Manuscript implication

The §"Practical Applications" subsection §"Drug target prioritization" in the original draft claims that "geometric proximity in $\SV{2}$--$\SV{4}$ can serve as an ordinal predictor of protein interaction likelihood" and that "this enables prioritizing candidate interacting partners of a target protein by ranking genes along the relevant spectral axes."

The honest reframe: the PPI subspace shows enrichment for STRING co-pole co-occurrence (Section 2.3), but **at the level of ranking known interactors of specific drug targets, the spectral approach is neither better nor worse than co-expression**. The application is therefore better characterised as "interpretable prioritisation" (one knows which axis is doing the work) than as "high-precision prioritisation" (it is not state-of-the-art).

Suggested manuscript text replacement:

> "We evaluated the practical utility of the PPI subspace ranking for drug-target candidate prioritisation. For 8 in-vocabulary drug targets, ranking genes by cosine similarity in $\SV{2}$--$\SV{4}$ at L11 yields P@5 = 0 (vs. 0.025 for co-expression Pearson and 0.025 for scGPT full-512 cosine; Table~\ref{tab:drug_targets}). The spectral PPI ranking is therefore **not** state-of-the-art for drug-target prioritisation in this small benchmark. We interpret the spectral structure as a real and *interpretable* organisational principle of scGPT's geometry rather than as a competitive ranker; downstream applications would need to combine spectral proximity with co-expression or other signals to be useful in practice."

## Caveats

- 8 drug targets is a small denominator. Replicating with broader gene vocabulary (E3) and external drug-target databases (e.g., ChEMBL, DrugBank) would strengthen the conclusion.
- Median 12 interactors per target → max possible P@5 ≈ 12/968 ≈ 0.012 if hits were uniform random; the observed 0.025 is ~2x random. This is a real-but-modest signal.
- The benchmark uses STRING ≥ 0.7 as the ground-truth interactor set, which is restrictive; loosening to STRING ≥ 0.4 would change the absolute numbers.
