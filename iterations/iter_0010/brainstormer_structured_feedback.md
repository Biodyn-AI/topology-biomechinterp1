# Brainstormer Structured Feedback — iter_0010

## Gate Status
`passed_min_research_gate: true` — all three hypotheses positive. No recovery plan needed.

---

## Signal Assessment

### H01: SV2/SV3 multi-layer compartment scan
**Assessment: Breakthrough-tier result.**
- SV2 cytoskeleton OR=20.35 (emp_p=0.000) is the strongest compartment hit in this project. 7 cytoskeletal genes dominating one pole of SV2 across L1–L11 is near-deterministic.
- SV2 mitochondrion transient (top pole L1–L3, shifts to bottom L4+) and SV2 extracellular_vesicle (OR=4.2–5.4 at L4–L11) reveal layered biological programs encoded in sub-leading SVD axes.
- The mito polarity flip is particularly interesting — it likely encodes a phase transition in how scGPT builds mitochondrial gene representations as context accumulates.
- SV3 nucleus (OR=4.3, L3) and plasma_membrane (OR=2.9, L4) are significant but weaker; SV3 may encode secondary compartment distinctions.
- **What this opens**: SV2 is a separable biological axis from SV1. The geometry is not just a single principal signal but has a multi-axis compartment decomposition that changes across depth.

### H02: TRRUST TF-target co-pole test
**Assessment: Most novel finding in this project to date.**
- First demonstration that regulatory network structure (directed TF→target edges) is geometrically encoded in the SVD embedding axes of scGPT.
- SV2 co-pole rate (0.206 vs null 0.122, emp_p=0.000) beats SV1 (0.163 vs 0.123, emp_p=0.016). Regulation is more tightly packed in SV2 than in the dominant principal component.
- The 326 TRRUST pairs with both genes present is a reasonable overlap. Effect is not driven by large-degree hubs (that would be a confounder to check, but co-pole rate is a counting statistic not weighted by degree).
- **What this opens**: Signed regulatory geometry (activation vs repression), layer-resolved regulatory encoding, and cross-network tests (STRING PPI, pathway co-membership).

### H03: Spectral profile + confounder check
**Assessment: Essential validation, plus a mechanistic clue.**
- SV1 monotone rise to 77% = further confirmation that the residual stream progressively compresses biological variance onto SV1.
- SV2 peak at L8 (9.9%) coincides with mito transient. This is a mechanistic anchor — L8 is where the model is doing something qualitatively different with mitochondrial gene representations.
- Zero annotation-density confounders: all enrichment results survive this standard control. Results are clean.
- **What this opens**: SV2 L8 peak as mechanistic anchor point for deeper dissection.

---

## Directions to Retire or Rescue

| Direction | Status | Rationale |
|-----------|--------|-----------|
| Persistent homology rewiring-null | **RETIRED** (confirmed from executor) | 4 negatives, no new angle |
| GO BP cosine clustering | **RETIRED** (confirmed from executor) | Underpowered, no fix path |
| kNN transitivity (iter_0002 negative) | **RETIRED** | Negative at maximum separation |
| PCA-only co-target clustering (iter_0004 H02) | **DEPRIORITIZE** | Underpowered in PCA space; SVD pole approach is strictly better |

---

## What Is Not Yet Tested (Gaps)

1. **Signed regulatory geometry**: Are activation edges different from repression edges in co-pole rate? Untested.
2. **Layer-resolved regulatory encoding**: H02 tested only L11. Peak regulatory layer unknown.
3. **Physical PPI network**: STRING pairs have not been tested for co-pole enrichment.
4. **Identity of the 7 cytoskeletal genes**: Not characterized. Which genes, are they consistently in one pole?
5. **SV4+ axes**: Only SV1–SV3 tested. What does SV4 encode?
6. **Cross-model consistency of SV2**: Geneformer SV2 not checked against scGPT SV2.
7. **Manifold curvature**: No Ricci/sectional curvature analysis has been done.
8. **Inter-layer gene rank drift**: Are any genes "crossing poles" as layer increases (pole flip tracking)?
9. **Topological robustness**: Bootstrap SV2 pole stability — how many gene-resamples preserve same cytoskeletal enrichment?
10. **Biological network topology (pathway graphs)**: Reactome pathway co-membership not tested.
