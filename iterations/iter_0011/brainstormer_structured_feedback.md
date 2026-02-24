# Brainstormer Structured Feedback — iter_0011

## Research Gate
`passed_min_research_gate: true` — all systems go.

---

## Signal Assessment

### H01 — TRRUST co-pole × 12 layers + Signed stratification
**Verdict: Major positive.** This is the strongest finding in the project so far.

- Co-pole enrichment at ALL 12 layers (emp_p<0.05) confirms this is a structural property of the model, not a layer-11 artifact.
- Signed stratification (activation 12/12 sig; repression 1/12 sig) is the first evidence that the transformer residual stream geometry reflects regulatory sign — a genuinely novel result.
- Layer 4 peak (emp_p=0.000, obs=0.234 vs null=0.130) suggests a specific architectural amplification point worth characterizing.
- **Immediate follow-up needed**: The act−rep difference is visible but not yet formally tested. A bootstrap test of Δ = obs_act − obs_rep > 0 across layers would convert this observation to a publishable claim with error bars.

### H02 — Mito pole-flip at L3→L4
**Verdict: Productive negative.** Hypothesis falsified but with useful byproducts.

- The L3→L4 mito compartment OR shift (iter_0010) was a rank-based enrichment effect, not a sign flip. This is now understood. No further L3→L4 flip work needed.
- 12/14 mito genes stably bottom-pole from L1 onward reinforces the compartment structure finding.
- HIF3A (3-flip outlier at L0→1, L8→9, L10→11) is anomalous. HIF3A is a hypoxia-inducible factor — its dynamic behavior may reflect context-sensitive regulatory state and is worth one focused characterization iteration.

### H03 — GO BP enrichment in SV2 poles
**Verdict: Clean negative. Retire.** 0/591 terms at 3 layers. SV2 encodes cellular compartment (CC) structure, not biological process programs. This is itself an interpretable result: the axis is localization-based, not function-based.

---

## Synthesis of Running Picture

After 11 iterations the core positive findings are:

1. **SV2 = cellular localization axis**: CC compartments (cytoskeleton OR=20, mito OR=6.4, EVs OR=4-5) are significantly enriched in SV2 poles across all 12 transformer layers.
2. **Regulatory network geometry**: TRRUST TF-target pairs co-localize in SV2 space above null at all 12 layers. Peak L4.
3. **Signed geometry**: Activation edges are geometrically segregated in SV2 space (12/12 layers); repression edges are not (1/12 layers). This is a new structural asymmetry finding.

The current main gap: **we don't know what SV1 encodes, what SV3+ encode, or whether the regulatory geometry is specific to SV2 or present in a multi-SV subspace.**

---

## Priority Signal for iter_0012

The activation/repression differential needs a formal statistical test immediately. This is the highest-value low-cost action. After that, STRING PPI co-pole provides independent cross-graph validation.
