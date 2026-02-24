# Brainstormer Structured Feedback — iter_0058

**Date**: 2026-02-23
**Iteration**: iter_0058
**Research gate**: PASSED

---

## Iteration Assessment

iter_0058 produced two strong positives and one weak neutral.

### H01 (TF boundary gene identity + family enrichment): STRONG POSITIVE
Top anchors BCL11A, NFKB1, FOXO3 are biologically coherent for an immune dataset. Family stratification (Forkhead/STAT high-margin, bZIP/C2H2-ZF negative-margin) is a clean finding. The cross-layer consistency (L2 and L3) confirms this is not noise. The biological story is compelling: master immune TFs are geometrically farthest from their target cloud. **Next: cross-seed validation to claim full robustness.**

### H02 (TRRUST Laplacian alignment): WEAK NEUTRAL — PARTIAL RETIRE
Z=7.81 is statistically real but the 3° absolute improvement over random is biologically vacuous. The hypothesis that TRRUST graph structure imprints on SV subspace geometry is refuted at the whole-graph level. Rescue option: restrict to the 295 circuit-gene subgraph where the adjacency is denser. If the circuit-restricted version also fails, retire fully.

### H03 (Effective rank → AUROC, rho=0.855): STRONG POSITIVE
This is the strongest single-metric predictor to date. The gap between full-rank (rho=0.855) and SV2-7 rank (rho=0.539) is the key mechanistic clue: the full spectrum drives AUROC, not just the TF-separating subspace. Open question is whether this is a causal mechanism (more dimensions = more room for separation) or a confound (early scGPT layers happen to do both: preserve dimensionality and preserve gene identity). The next priority should disentangle these.

---

## Direction Status

| Direction | Status | Rationale |
|-----------|--------|-----------|
| Joint 6D SV2-7 classifier | Active | Mean AUROC 0.751, fully cross-seed validated, central result |
| TF boundary anchor identity | Active | BCL11A/FOXO3/STAT3 geometry, needs cross-seed validation |
| Effective rank → AUROC | Active | rho=0.855, strongest predictor, needs mechanistic decoupling |
| SV2-4 vs SV5-7 encoding regimes | Active | Complementary encoding confirmed; L9 crossover needs seed43 follow-up |
| TRRUST Laplacian alignment | Rescue-once | Restrict to circuit subgraph only; if still weak, retire |
| SV basis permutation | RETIRED | Definitively ruled out (iter_0057) |
| SV rotation angle → AUROC | RETIRED | rho=-0.27, p=0.42 (iter_0056) |
| PH on circuit genes | RETIRED | iter_0052 negative |
| Raw graph rewiring null | RETIRED | iter_0050 retired |

---

## Key Open Questions for Steering

1. **Mechanistic decoupling**: Is the effective rank → AUROC link (rho=0.855) mediated by global capacity, or by something specific about early-layer representations preserving fine-grained gene identity? Test: partial out layer depth from both variables.

2. **Boundary anchor causality**: Are BCL11A/NFKB1/FOXO3 boundary anchors because they are *highly connected* hubs in TRRUST (regulatory breadth → central position in embedding space), or because they have *distinctive expression profiles* that scGPT encodes as geometrically extreme?

3. **bZIP/C2H2-ZF proximity**: Why are these TF families near the target centroid? Two hypotheses: (a) high TRRUST degree (they regulate many of their own targets) or (b) co-expression structure (they share expression patterns with their targets). Testable with TRRUST degree distribution.

4. **Cross-model stability**: Does the effective rank–AUROC relationship hold in Geneformer embeddings? If yes, this is a model-agnostic property of transformer-based gene language models.

5. **Individual pair geometry**: All results so far are class-level (TF centroid vs target centroid). Can we predict individual TF→target edge probability from pairwise distance in the 6D SV2-7 space?
