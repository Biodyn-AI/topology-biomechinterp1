# Executor Next Steps — iter_0007

## Confirmed positive claims (publishable with caveats)

1. **SV1 extracellular/cytosol axis** — confirmed with gene-label shuffle null (empirical p=0.004). Layer-11 SV1 top-pole = extracellular space (GO:0005615, OR=6.37), bottom-pole = cytosol (GO:0005829, OR=2.96). The SV1 explains 93.2% of named-gene variance at layer 11.

2. **Layer-wise SVD trajectory** — SV1/SV2 ratio rises 4.07→7.70 across 12 layers. 5x dominance threshold crossed at layer 2. Extracellular axis detected in 8/12 layers. Effective rank 19.54→1.63 (12x compression). Monotonic compression correlated with SV1 axis strengthening.

## Retired/demoted claims

- **Drift GO enrichment for TF activity** — does NOT survive gene-label shuffle null (empirical p=0.124). Retire as primary evidence. Retain as supporting trend only with this caveat.

## Top 3 priorities for iter_0008

### Priority 1: Layer-wise label-shuffle null (extend H01 across all 12 layers)
- Run gene-label shuffle null (N=200) at each of 12 layers
- Test whether extracellular axis specificity is maintained at all layers or only emerges late
- Expected: specificity is strongest at layers 9-11 where SV1/SV2 ratio is highest
- This would confirm the "trajectory" as biologically meaningful throughout, not just at layer 11

### Priority 2: SV2/SV3 biological axes (layer 11 and trajectory)
- Currently only SV1 analyzed. SV2/SV3 may encode additional biological axes.
- GO enrichment for SV2 top/bottom quartile at layer 11
- Label-shuffle null to confirm specificity
- Expectation: SV2 may encode immune vs. metabolic axis or proliferation axis

### Priority 3: STRING network distance vs. SV1 axis proximity
- Test whether gene pairs that are closer in the SV1 axis (smaller |SV1 score difference|) are also more connected in the STRING protein-protein interaction network
- Metric: Spearman correlation between |SV1 score difference| and STRING co-expression/interaction score for all 209-gene pairs
- Null: permute gene labels to break the correspondence
- This would directly connect embedding geometry to protein functional networks

## Novel high-priority explorations (new families)

### Manifold topology — persistent homology at layer 11 (NEW)
- Prior PH attempts retired (rewiring-null failed). Fresh approach: use 209 named genes only (not full vocab), with raw Euclidean distance filtration, and compare Betti-0 curve vs. random label shuffle null.
- Small dataset makes this computationally feasible.

### Cross-model alignment with layer trajectory
- iter_0004 H03 found strong CKA cross-model alignment. Now test whether SV1 alignment trajectory mirrors the layer compression trajectory (does SV1 emergence coincide with cross-model alignment peak?).

## Methodological notes for iter_0008

- Always use gene-label shuffle (not feature shuffle) as the null for L2-norm or SVD-based statistics.
- Feature shuffle is only valid for statistics that genuinely depend on the relative positions of features (e.g., directional projections onto fixed unit vectors).
- N=500 label-shuffle reps is adequate for empirical p-values down to ~0.002.
