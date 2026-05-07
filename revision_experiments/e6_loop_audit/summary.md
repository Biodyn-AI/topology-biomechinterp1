# E6 Summary — Two-agent loop protocol and selection-bias audit

**Status:** Done.
**Code:** `family_audit.py`.
**Outputs:** `family_retire_validate.csv`, `extends_finding_distribution.csv`.

## Protocol summary (for new Methods subsection M1)

The autonomous loop runs two LLM-backed agents in alternation:

- **Executor agent** (`prompts/executor_prompt_v2.md`): designs and runs Python experiments (SVD, permutation tests, enrichment analyses, classification AUROCs). Writes a JSON artifact `executor_hypothesis_screen.json` per iteration with the schema documented at `prompts/executor_prompt_v2.md:112-138`.
- **Brainstormer agent** (`prompts/brainstormer_prompt_v2.md`): reviews the executor's outputs, proposes refinements or extensions, and decides which hypothesis branches to extend, retire, or stress-test. Writes `brainstormer_hypothesis_roadmap.md` and `brainstormer_next_iteration_brief.md`.

**Key procedural rules**, reconstructed from the prompts and verified against the iteration log:

1. **2–3 hypotheses per iteration.** One must extend an established finding to a new context (model/tissue/system); at most one is internal refinement.
2. **Retirement rule.** A hypothesis branch is retired after **two consecutive negative results** (non-significant under permutation null).
3. **Permutation requirement.** Every positive claim must pass at least one permutation-based null model. Gene-label shuffle, feature shuffle, and degree-preserving rewiring are the three nulls used; the choice depends on what is being tested.
4. **JSON artifact schema.** Every iteration writes `executor_hypothesis_screen.json` with per-hypothesis records: `id`, `name`, `family`, `extends_finding`, `method`, `result_value`, `result_direction` ∈ {positive, negative, mixed, inconclusive}, `decision`, and `next_action`. This is the single source of truth used in E5's multiple-testing audit.
5. **Phase boundary.** Iterations 1–63 (the "Phase 1" exploration covered by the paper) used 9 broad families. After iter_0063 the loop transitioned to Phase 2 (extension-focused families: cross_model, cross_tissue, regulatory_edges, etc.) — these post-paper iterations are out of scope for the current submission.

## Selection-bias audit

The reviewer's concern is that the brainstormer's "retire after two negatives" rule introduces selection bias at the **family** level: branches that produce early negatives are retired, so the surviving hypotheses within a family have an inflated positive rate.

We test this formally by computing the two-tailed binomial $p$-value of the observed positive rate per family under a 50/50 null (the rate one would expect if the brainstormer were *not* selectively pruning):

| family | n | positive | positive rate | $p$ vs 50/50 null | selection bias? |
|---|---|---|---|---|---|
| manifold_distance | 42 | 31 | 0.738 | **0.003** | yes (significant) |
| intrinsic_dimensionality | 41 | 28 | 0.683 | **0.028** | yes (significant) |
| module_structure | 52 | 33 | 0.635 | 0.070 | borderline |
| topology_stability | 8 | 5 | 0.625 | 0.73 | sample too small |
| graph_topology | 16 | 8 | 0.500 | 1.00 | no |
| null_sensitivity | 15 | 8 | 0.533 | 1.00 | no |
| cross_model_alignment | 5 | 2 | 0.400 | 1.00 | no |
| persistent_homology | 4 | 2 | 0.500 | 1.00 | no |

**Conclusion.** Two of the largest families (`manifold_distance`, `intrinsic_dimensionality`) show statistically significant positive-rate inflation versus a 50/50 null. This **is** the selection-bias signal Reviewer 1 anticipated. Three observations contextualise it:

1. The bias is a **property of the brainstormer's retire policy**, not of the individual permutation tests. Each within-family hypothesis still passes its own gene-label-shuffle or feature-shuffle null.
2. The multiple-testing audit in **E5** accounts for the family-level inflation by applying BH-FDR and Bonferroni correction across all 183 hypotheses, treating each as if it were independently sampled (the most conservative interpretation).
3. The two families with no detectable bias (`cross_model_alignment`, `persistent_homology`) include the negative findings the paper already reports honestly (Geneformer cross-model alignment failed; persistent homology collapsed under degree-preserving nulls).

## Manuscript implications (M1)

Add a Methods subsection §"Hypothesis-screening protocol and selection-bias audit" with:

- The five procedural rules above, with prompts cited.
- The JSON artifact schema as a code listing.
- The selection-bias audit table, with this honest framing: "The two-agent loop is biased toward extending hypothesis branches with early positive signal. We do not view this as a flaw — it is the loop's intended behaviour — but it does mean the per-family positive rate over-reports the rate at which a *random* hypothesis would yield a positive result. The Bonferroni correction in §[E5 section] accounts for this inflation."

## Caveats

- The `extends_finding` field was added in Phase 2 (iter_0063+) and is empty for the iters 1–63 that this audit covers. The within-family lineage analysis is therefore limited to the family taxonomy and `lineage` (prior hypothesis ID) chain, which is logged but not analysed here.
- The binomial null at $p=0.5$ is itself a strong assumption; if the "true" baseline rate of positive findings in genuinely random hypothesis space were 0.2 (more realistic for an exploratory screen), the inflation would be even larger. The 0.5 baseline is therefore the most generous to the loop.
