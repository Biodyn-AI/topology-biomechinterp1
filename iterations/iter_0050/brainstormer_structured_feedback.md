# Brainstormer Structured Feedback: iter_0050

## Gate Status
PASSED. 2/3 hypotheses positive. Paper updated with iter_0050 marker. Artifacts present.

---

## Per-Hypothesis Assessment

### H01: SV2-4 Co-expression Specificity Control
**Decision: CRITICAL NEGATIVE — claim revision required**

The specificity test failed decisively: co-expression pairs (top-300 by embedding cosine, mean cos=0.885) are CLOSER in SV2-4 than TRRUST pairs at all 12 layers. TRRUST < co-expression: 0/12. This collapses the iter_0049 claim that SV2-4 encodes "regulatory proximity" specifically — it encodes co-expression proximity, within which TRRUST pairs happen to sit.

Implication: The regulatory signal in SV2-4 is mediated by co-expression. TFs tend to be co-expression network hubs; their targets co-express with them. The model may be encoding co-expression structure (which is biologically real), not causal regulatory topology per se.

This is scientifically informative — the finding is not worthless, it's refined. But the paper's framing must be updated.

**Next action options**:
1. Test whether residual TF enrichment survives after regressing out co-expression similarity (cosine) from SV2-4 distance.
2. Test whether TF hub degree (in the inferred co-expression network from embeddings) explains SV2-4 loading.
3. Reframe the claim: "SV2-4 encodes co-expression structure enriched for transcription factor regulatory modules."

### H02: SV1 Identity — Correlation Sign Flip
**Decision: POSITIVE — mechanistic insight, follow-up needed**

Clear sign flip: rho(SV1, norm) goes -0.33 (L0) → ~0 (L3-L8) → +0.27 (L11). Aligns with geometric rotation at L9-L10. The SV1 axis switches from "low-norm" to "high-norm" gene selection across the network depth.

This quantifies the rotation mechanistically: it's not just a geometric rotation in abstract space — it's a semantic inversion in what type of gene anchors the primary variance axis.

**Gap**: We have correlations but not gene identities. We don't yet know which genes are high-SV1 at L0 (low expression) vs. L11 (high expression). The biological identity of these two populations is the key open question.

### H03: SV5-10 Spectral Screen
**Decision: POSITIVE — extends spectral hierarchy**

Spectral hierarchy confirmed: SV2-4 (rbc=0.229) > SV5-7 (rbc=0.152) > SV8-10 (rbc=0.123) > SV1 (rbc=-0.11). SV5-7 is a new positive (7/12 layers, p<0.05 for TFs). The regulatory signal doesn't collapse sharply after SV4 — it decays gradually across SV5-10.

**Open question**: Is SV5-7 TF enrichment redundant with SV2-4 (same TFs, same co-expression structure), or does it capture an independent regulatory subspace? This matters for whether the spectral hierarchy is scientifically meaningful or just reflects that TF-related co-expression explains variance across many SVs.

---

## Narrative Revision

The current paper anchor — "SV2-4 encodes regulatory proximity" — must be revised to:

> "SV2-4 encodes co-expression structure that is preferentially enriched for transcription factor regulatory modules. The model learns co-expression, not causal regulation, as its primary organizational principle."

This is a more defensible and arguably more interesting claim: scGPT's intermediate representations reflect co-expression network topology, which is the ground truth the model sees during training (single-cell expression matrices). TF-target enrichment within SV2-4 follows directly from TFs being co-expression hubs.

The SV1 story (semantic inversion from L0 to L11) is now the strongest mechanistic result with no competing confound identified.

---

## Key Open Questions Entering iter_0051

1. **Gene identity**: Which specific genes have highest SV1 loading at L0 (anti-norm) vs L11 (pro-norm)?
2. **Residual TF signal**: After controlling for co-expression proximity, does any TF-specific signal survive in SV2-4?
3. **SV5-7 independence**: Is SV5-7 TF enrichment orthogonal to SV2-4, or fully explained by co-expression variance spilling into higher SVs?
4. **Topological structure**: Do TRRUST TF-target pairs occupy distinct topological positions in SV2-4 space (0-dim PH clusters) beyond mere proximity?
5. **Biological identity of SV5-7**: What GO categories are enriched in top-50 genes by SV5-7 loading magnitude at peak layer?
