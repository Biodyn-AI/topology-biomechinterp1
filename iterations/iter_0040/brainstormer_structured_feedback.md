# Brainstormer Structured Feedback — iter_0040

## Gate Status
`passed_min_research_gate: true`. Clean iteration with 2 strong positives and 1 informative negative.

---

## Signal Quality Assessment

### Strong confirmed findings entering iter_0041
1. **GC attractor onset = Layer 3** — confirmed by two independent measures (rank scan rho=-0.951; p@10 onset). This is now a publication-grade claim.
2. **GC subspace → B-cell subspace convergence** — principal angle rho=-0.888, p=0.0001. Directionally specific (plasma = null). Geometrically clean.
3. **Revised GC attractor = {BATF, BACH2, PAX5}** — BCL6 excluded. PAX5 is the earliest (rank 36 at L0).
4. **CD19 criticality** — confirmed for 3rd time across different methods.

### What the BCL6 exclusion means
BCL6 ranks 750-1500 throughout. This is not noise — BCL6 does something different. It may cluster with another cell-type program (plasma, memory B) or exist in a near-orthogonal direction. This is unexplored.

### What the Geneformer negative means
The B-cell/GC-TF proximity is scGPT-specific (not a generic gene co-expression artifact). This strengthens the claim that scGPT's learned layer representations encode lineage-regulatory logic. However: we tested INPUT token embeddings only. Geneformer LAYER-WISE representations remain untested.

---

## What to Retire

- **GC-plasma subspace convergence**: rho=-0.098, null for 2 iterations. Retire.
- **BCL6 as GC attractor member**: corrected. Paper claims being updated.
- **Geneformer INPUT token embedding cross-model test**: negative; not worth repeating. Requires a different approach (layer-wise forward pass).

---

## Primary Open Questions for Next Iterations

1. **Mechanistic: what changes at L3?** — This is the most important unanswered question. Both rank and subspace measures agree onset = L3. We need structural evidence for why.
2. **PAX5 pre-wiring at L0** — rank 36 even before any contextual processing. Is this a training-data co-occurrence artifact, or does scGPT assign PAX5 a structural role in input representations?
3. **BCL6 neighborhood** — If BCL6 is NOT in the B-cell/GC-TF cluster, where does it live? Its rank ~900 near B-cell centroid needs a positive-control neighborhood test.
4. **Other lineage attractors** — We have only tested the B-cell/GC-TF system. Do T-cell or myeloid TF attractors exist? Do they have the same L3 onset?
5. **Topological signature of L3 transition** — Persistent homology on gene embedding point cloud. Do topological features (connected components, loops) change at L3?
