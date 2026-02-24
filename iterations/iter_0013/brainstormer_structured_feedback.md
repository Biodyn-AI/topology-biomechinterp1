# Brainstormer Structured Feedback — iter_0013

## Gate Status
`passed_min_research_gate: true` — Execution was valid. 3 hypotheses tested with code artifacts, quantitative results, and paper updated.

## Per-Hypothesis Assessment

### H01: SV Axis Specificity (Mixed/Positive)
Strong result. SV2 (z=10.04) and SV3 (z=7.18) both significant across all 12 layers; SV1 is noise (z=0.39). This upgrades the working model: scGPT encodes PPI proximity across a 2D subspace (SV2+SV3), not a single axis.

Next lever: What is SV3 biologically? Is SV2 encoding something like "metabolic/cytoskeletal" and SV3 "signaling/nuclear"? GO annotation of SV3 poles is the obvious next step. Also worth asking: do SV2 and SV3 together create a 2D geometry (e.g., ring or arc) for PPI neighborhoods?

### H02: Repression Anti-Pole (Definitive Negative)
Clean negative. rep_xpole z=−1.41 uniformly across all layers. Retired correctly. The encoding model is: both activators and repressors of common target genes co-localize (mild positive z for both), but regulatory sign is not encoded as spatial opposition. This is actually an interesting constraint on the representation — it suggests scGPT captures "functional neighborhood" but not "regulatory direction."

### H03: Hub-Degree Control (Inconclusive/Underpowered)
The critical confound question (hub degree driving PPI co-pole signal) cannot be answered with N=35 non-hub pairs. This is not a dead end — it is a data-expansion problem. Lowering STRING threshold to score≥0.4 is the right fix. Also worth trying: use a degree-matched null instead of random null (sample null pairs with matching degree distribution to the test pairs).

## Cumulative Picture After iter_0013

Confirmed positives:
- SV2 (primary) and SV3 (secondary) jointly encode PPI neighborhood geometry — strong, 12/12 layers
- TRRUST activation pairs cluster in SV2 (continuous distance test, iter_0012 H02)
- Activation > repression co-pole qualitative signal (iter_0011)

Confirmed negatives:
- Regulatory sign (act vs rep) NOT encoded as opposite poles in any SV axis
- Bootstrap CI underpowered for direct act/rep differential (iter_0012 H01)

Open questions:
- Is hub degree a confound for PPI co-pole signal (need more non-hub pairs)
- What does SV3 encode biologically?
- Does the SV2+SV3 2D embedding form interpretable geometry (clusters, arcs, gradients)?
- Are these effects scGPT-specific or generic to any transformer trained on gene expression?

## Strategic Direction for iter_0014

Primary opportunity: characterize the SV2+SV3 joint geometry. A 2D projection scatterplot with GO/pathway annotation would tell us whether we have discrete clusters, continuous gradients, or manifold structure. This is publishable characterization, not just a binary test.

Secondary opportunity: resolve the hub confound definitively — lower STRING threshold and use a degree-matched null.

Tertiary opportunity: begin cross-model or cross-layer structure tests to establish whether this is scGPT-specific or a general feature of gene language models.
