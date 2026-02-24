# Brainstormer Structured Feedback — iter_0012

## Gate Status
`passed_min_research_gate: true`. Full brainstorming pass proceeds.

---

## Per-Hypothesis Assessment

### H01 — Bootstrap CI on act-rep co-pole differential
**Verdict: RETIRE as standalone claim.**
The qualitative direction (act > rep) is consistent across 10/12 layers but the bootstrap CI is wide because n_rep=64 and binary pole-membership is a coarse metric. The per-edge-type empirical p-values from iter_0011 (activation 12/12 significant, repression 1/12) already contain the actual signal. A CI test on the differential is not the right tool at this sample size. Do not re-run this test in this form. Instead, redirect toward signed/directional geometry tests (see roadmap H-NEW-04 and H-NEW-05).

### H02 — SV2 spatial concentration (mean pairwise distance)
**Verdict: CONSOLIDATE into confirmed evidence.**
Activation pairs more concentrated than null at 10/12 layers; repression at 6/12. This is a continuous-geometry complement to the discrete co-pole result and uses the same null. The finding that repression is also concentrated (though less so) weakens the "activation-only" signed geometry claim from iter_0011. Important next step: test whether activation pairs are significantly more concentrated than repression pairs (direct paired comparison, not just each vs. null). Without this the act > rep claim is not fully established on continuous geometry.

### H03 — STRING PPI co-pole enrichment
**Verdict: STRONGLY POSITIVE — primary finding this iteration.**
12/12 layers significant, z=3.3–6.5, ~2× enrichment over null. This is the highest z-score of any biological network test across all iterations. However, two confounds require urgent attention before this is a clean result:

1. **Hub degree confounder**: High-degree PPI hubs have many within-set partners. If hubs cluster in one pole (which they might, as hubs are often functionally specialized), the co-pole enrichment could be driven almost entirely by a few high-degree nodes. This must be tested.

2. **STRING channel decomposition**: STRING combines co-expression, experimental, database, text-mining. The co-expression channel is partially circular with scGPT's training signal (both derived from expression data). The experimental channel (binding, co-complex) is the cleanest. If the co-pole enrichment is primarily driven by co-expression links, interpretability is weakened.

Both of these are cheap analyses on the already-cached STRING data and should be done before this result is cited as primary evidence.

---

## Cumulative State Assessment

**Established (fully validated):**
- SV1 = secretory/extracellular axis (12/12 layers, gene-label-shuffle, emp_p < 0.004)
- SV2 = cytoskeletal vs. EV export axis (12/12 layers compartment, OR up to 20x)
- TRRUST activation co-pole in SV2: 12/12 layers (signed geometry)
- STRING PPI co-pole in SV2: 12/12 layers, z=3.3–6.5 (convergent validation)
- SV2 spatial concentration (continuous): activation at 10/12 layers

**Retired:**
- GO BP enrichment in SV2 poles
- Bootstrap CI on act-rep differential as a standalone claim
- Feature-column shuffle as a null (degenerate for L2 norms)
- Drift TF enrichment (fails gene-label-shuffle specificity test)

**Open confounds requiring resolution:**
- STRING hub degree confound
- STRING channel decomposition (co-expression vs. experimental channel)
- Act > rep signed distance (continuous form not yet established)

**Gaps worth exploiting:**
- SV1/SV3 PPI co-pole specificity (is STRING enrichment SV2-specific?)
- Cross-model validation (Geneformer never meaningfully tested)
- Repression anti-pole geometry (does rep encode opposite-pole?)
- Layer trajectory of STRING z-scores: peak at L1 and L11, trough at L9-10 — unexplained
- Geodesic vs Euclidean distance in regulatory subspace
- Gene complex membership of SV2 poles (CORUM database)
