# Executor Next Steps: iter_0006 → iter_0007

## Priority Actions

### High priority: Strengthen H01 + H02 with null controls

Both H01 and H02 show promising enrichment p-values (0.0078 and 0.0003) but no FDR significance (n=209 genes is borderline). The critical next step is to add explicit null controls to determine whether the enrichment is meaningful.

**H01 null control (drift enrichment):**
- Shuffle gene labels (assign random GO terms to random drift-rank positions) → compute null enrichment distribution → compare to observed.
- Permute drift ranks within the 209 named genes → test whether the observed enrichment pattern holds under label permutation.

**H02 null control (SV1 enrichment):**
- Feature-shuffle the embeddings → recompute SVD → test whether SV1 GO enrichment is preserved (it should not be if the signal is real).
- Column-permute the embedding matrix → SV1 should become random noise.

### High priority: SV1 directionality test

The SV1 axis separates extracellular (positive) from presumably nuclear/intracellular (negative) genes. Verify:
- Explicitly identify the bottom-quartile SV1 genes and confirm they are enriched in nuclear/chromatin/TF terms (complementary to the top-quartile extracellular signal).
- This would confirm the axis is a coherent extracellular vs nuclear gene-function axis.

### Medium priority: Extend SVD analysis across layers

We have SVD data for layer-0 and layer-11. Run SVD for all 12 layers to:
- Track when the extracellular-space axis emerges (does it appear early or late?).
- Correlate SV1 singular value magnitude across layers with ER (does SV1 dominance increase monotonically?).
- Identify at which layer the top-loading gene set stabilizes.

### Medium priority: Cross-domain SVD replication

SVD biology was computed on lung-domain embeddings. Replicate on:
- Immune domain embeddings (if available) to see if the same SV1 axis (extracellular vs nuclear) emerges.
- If the axis is domain-invariant, it's a strong universality claim.

### Low priority: Increase GO analysis power with expanded gene list

The 209-gene TRRUST-derived set limits FDR power. Options:
1. Use all named genes in the full scGPT vocabulary (need the index-to-gene mapping for the full 4803 genes).
2. Search for a full gene_idx_to_symbol file in subproject_38 codebase.

## Retired Directions (do not revisit without rescue rationale)

- Rewiring-null PH survival: negative across 4+ iterations.
- Cross-model alignment (scGPT vs Geneformer): no Geneformer embeddings in accessible workspace.
- Distance permutation null (too adversarial, not biologically informative).

## Suggested iter_0007 Hypothesis Packet

1. **H01-null**: Feature-shuffle null control for SV1 GO enrichment → confirm biological specificity (decisive test).
2. **H02-fullSVD**: SVD across all 12 layers → track SV1 singular value and top-loading gene stability layer by layer (new family: spectral dynamics).
3. **H03-complement**: SV1 bottom-quartile enrichment → confirm nuclear/TF enrichment as complementary pole of extracellular-vs-nuclear axis.
