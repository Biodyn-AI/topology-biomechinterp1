# Brainstormer Structured Feedback — iter_0031

## Gate Status
`passed_min_research_gate: true` — clean pass, no recovery needed.

## Iteration Quality Assessment

**H01 (spectral gap):** Strong result. Spearman rho = −1.000, real/null gap ratio = 5.6×. The kNN graph grows progressively more modular across layers. This is now a second confirmed graph-topology law (alongside clustering coefficient from iter_0002). Ready to exploit biologically.

**H02 (STRING→embedding distance):** Definitively null. rho ≈ −0.015, AUROC ≈ 0.494 across all layers. Combined with the positive STRING PPI co-pole results (iter_0012, 0013), this creates an important dissociation: STRING partners co-localize in SVD-pole space but do NOT cluster in L2-distance space. This is itself a finding worth articulating — the geometry is axis-aligned (SVD poles), not radially symmetric (L2 distance). Retire the STRING→L2 direction completely.

**H03 (participation ratio):** Strong result. PR collapses 6.1× from L0→L11, PC1 variance rises 8%→26%. Perfectly monotonic (rho = −1.000). Combined with H01, the picture is: as layers deepen, the gene embedding manifold simultaneously (a) collapses onto fewer dimensions and (b) becomes more modular. This dual compression+fragmentation narrative is the key new structural story to develop.

## Convergent Picture After iter_0031

The two new strong results (H01 + H03) integrate with prior findings to reveal a coherent layer-depth narrative:
- L0: High-dimensional (PR~58), weakly modular (spectral gap 0.145), norms uniform
- L11: Low-dimensional (PR~9.5), strongly modular (spectral gap 0.023), norms still uniform

This is consistent with a model that compresses the gene space into a few dominant directions while simultaneously separating gene clusters. The SVD biological axes (SV2/SV3 encoding PPI, compartment, regulatory geometry) were found on the full 209-gene set but the 195-gene corrected set now provides the proper substrate for downstream community/annotation tests.

## Key Open Questions for Next Iterations

1. What biology do the kNN communities at L11 correspond to? (P1 from executor next steps)
2. Does GO Jaccard (functional similarity) predict embedding distance? (P2 — STRING failed, GO may succeed)
3. Does PC1 at L11 align with any known biological partition? (P3)
4. Is the spectral gap / PR collapse consistent across different cell-type contexts (not just lung)?
5. Do the SVD PPI co-pole results (iter_0012–0013) replicate on the corrected 195-gene set?
6. Does the dual compression+fragmentation signal appear in Geneformer? (cross-model validation)

## Directions to Retire

- **STRING→L2 embedding distance**: Definitively null (H02 this iteration + negative in iter_0028 on different gene sets). No rescue potential.
- **GO BP enrichment in SVD poles (H03 iter_0011)**: Already retired — zero significant terms.
- **Repression anti-pole geometry (H02 iter_0013)**: Already retired — no cross-pole signal.
- **Bootstrap CI on act/rep differential**: Underpowered; do not revisit unless n_repression_pairs grows substantially.

## Directions to Rescue

- **TRRUST co-pole on 195 in-vocab genes**: iter_0012 used 209 genes including OOV. The SV2 co-pole tests should be rerun on 195 genes to validate that the biological signal survives the OOV correction.
- **Cross-model validation (Geneformer)**: Multiple positive findings have accumulated without Geneformer replication. This is a major gap. The dual compression+fragmentation finding is a good candidate for cross-model test.
