# Next Steps — iter_0011 → iter_0012

## Priority Queue

### P1 — Activation vs Repression differential significance (refinement, family=module_structure)
**Why**: H01 shows activation significant at 12/12 layers vs repression at 1/12. Need direct statistical test of the difference.
- Bootstrap test: is act_copole_rate − rep_copole_rate > 0 across all layers?
- Compute per-layer effect size: (obs_act − null_mean_act) / null_std_act vs (obs_rep − null_mean_rep) / null_std_rep.
- Also stratify by "Unknown" mode edges to see if they track with activation or repression.
**Target**: Quantified evidence that activation geometry > repression geometry across the transformer stack.

### P2 — STRING PPI co-pole test (new_family, family=graph_topology)
**Why**: TRRUST regulatory co-pole is significant. Physical interaction (STRING PPI) co-pole would establish convergent evidence from a completely different graph type.
- Download or locate STRING human PPI (high-confidence >=700). Check /Volumes/Crucial X6/.../external/networks/.
- Map to 209-gene set, compute co-pole rate for high-confidence PPI pairs in SV1 and SV2.
- Null: N=1000 gene-label shuffles.
**Target**: If PPI and regulatory edges both show co-pole enrichment, this is model-independent evidence of functional-network geometry.

### P3 — HIF3A multi-flip characterization (refinement, family=intrinsic_dimensionality)
**Why**: HIF3A is the only mito gene that flips SV2 pole at layers 0, 8, and 10. This is anomalous.
- Extract HIF3A's layer-wise embedding trajectory (L2 norms, SV1 vs SV2 projections).
- Check HIF3A's TRRUST regulatory context: which TFs regulate HIF3A in our set? Are they co-located?
- Compare HIF3A trajectory to CCR7 and PMAIP1 (the two stable-top mito genes) to understand what makes them different.

### P4 — Cross-model SV2 alignment: check Geneformer availability (new_family, family=cross_model_alignment)
**Why**: Given strong SV2 results in scGPT, testing whether Geneformer encodes a similar SV2 axis would be a major finding.
- Check subproject_38 outputs for Geneformer embeddings.
- If available: compute Procrustes alignment / CKA of scGPT vs Geneformer SV2 axes.

## Retired Directions
- GO BP → SV2 enrichment (0/591 significant at 3 layers, NEGATIVE).
- Persistent homology rewiring-null (retired since iter_0010).
- GO BP cosine clustering (retired since iter_0010).
- kNN graph transitivity (retired iter_0002).

## Confirmed Strong Findings for Paper
1. SV2 compartment structure (CC: cytoskeleton OR=20.35, mito OR=6.36, extracellular vesicle OR=4-5) — stable across 11+ layers.
2. TRRUST TF-target co-pole enrichment — consistent across ALL 12 layers (not a layer-11 artifact). Peak emp_p=0.000 at L4.
3. Signed regulatory geometry: activation edges (12/12 layers significant) vs repression edges (1/12 layers) — first evidence of sign-specific regulatory geometry in transformer residual stream.
