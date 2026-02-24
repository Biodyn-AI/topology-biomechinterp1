# Executor Next Steps — iter_0012

## Status Summary (iter_0012)

- H01 NEGATIVE: Bootstrap CI on act-rep differential underpowered (n_rep=64). Retire this specific sub-claim.
- H02 PROMISING: TRRUST activation pairs significantly concentrated in SV2 at 10/12 layers (continuous distance). Repression also concentrated at 6/12 layers.
- H03 STRONGLY POSITIVE: STRING PPI co-pole enrichment at 12/12 layers (z=3.3–6.5). Cross-graph-type convergent validation.

## Priority Queue for iter_0013

### H-A′ (HIGH PRIORITY): STRING PPI co-pole on SV1 and SV3
- H03 established STRING PPI co-pole is strong on SV2 (12/12 layers).
- Test same pipeline on SV1 and SV3 to ask: is PPI clustering SV2-specific or a general manifold property?
- If SV2-specific: supports axis identity theory. If SV1 also strong: SV1 may encode different biological organization.
- Family: module_structure. Novelty: new SV axis.

### H-B′ (HIGH PRIORITY): SV2 vs SV1 co-pole differential for PPI and TRRUST
- At each layer, compare co-pole rate on SV1 vs SV2 for both STRING PPI and TRRUST edges.
- This directly tests the axis-specificity of the biological organization.
- Family: module_structure. Novelty: axis comparison.

### H-C (MEDIUM PRIORITY): SV1 axis identity
- What does SV1 encode? The SV1 axis has been shown to have large variance fraction (iter_0007).
- Run GO CC enrichment (as done for SV2) on SV1 poles to characterize SV1's biological content.
- Family: module_structure. Novelty: SV1 characterization.

### H-D′ (HIGH PRIORITY): Repression edge spatial geometry — detailed characterization
- H02 shows repression pairs are also concentrated (6/12 layers) in SV2. Is this at the same or different SV2 positions as activation pairs?
- Compute mean SV2 projection for activation-pair genes vs repression-pair genes (signed, not just distance).
- Test if activation vs repression gene sets are on opposite sides of SV2 axis (signed geometry).
- Family: manifold_distance. Novelty: signed comparison of pair types.

### H-E (MEDIUM PRIORITY): STRING PPI co-pole by edge subtype (coexpression, experimental, database)
- STRING API returns subscores: nscore (neighborhood), fscore (fusion), escore (experiments), dscore (databases), tscore (textmining).
- Split edges by dominant evidence type and test if co-pole signal is driven by specific interaction types.
- Family: module_structure. Novelty: evidence-type stratification.

### H-F (LOWER PRIORITY): Cross-layer co-pole stability for STRING PPI
- For each pair of layers (i,j), compute how often the same PPI pair is co-pole in both.
- Quantify cross-layer consistency. If PPI co-pole is stable across layers, this further strengthens the claim.
- Family: topology_stability. Novelty: cross-layer co-pole correlation.

## Retired Directions
- **GO BP → SV2 enrichment**: 0/591 terms significant (iter_0011 H03). Retired.
- **Mito pole flip hypothesis**: No sign flips observed (iter_0011 H02). The "transient" was a rank artifact. Retired.
- **Bootstrap CI on act-rep differential** (iter_0012 H01): Underpowered at n=64; retire this specific framing. The per-set empirical p approach (iter_0011) remains valid.

## Paper Update Notes
- Add STRING PPI co-pole result (12/12 layers, z=3.3-6.5) to main results section.
- Add SV2 spatial concentration result (H02, 10/12 layers, act pairs concentrated).
- Note Bootstrap CI limitation (H01) as a methodological note.
- Expand "SV2 axis" section to cover convergent evidence from three graph types.
