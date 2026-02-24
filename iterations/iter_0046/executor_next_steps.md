# Next Steps: iter_0047

## Standing retirement list
- kNN lineage purity (underpowered: <3 genes/lineage)
- BCL6 distance-to-centroid variants (repeated negatives)
- STRING/GO proximity (repeated negatives)
- Lineage centroid orthogonality (H03 this iter: z<0.3 at all layers, negative)
- manifold_distance family generally (too many negatives; only rescue with strong rescue rationale)

## Primary validated finding (carry into paper)
- **Monotone TwoNN ID compression**: L0~33 → L11~19, Δ≈-14, ratio≈0.57, confirmed across 3 seeds, non-overlapping bootstrap 95% CIs (L0=[27,37], L11=[16,20]).
- **SVD effective rank collapse**: L0=23.6 → L11=1.64 (14× collapse), top-1 variance 53.7% → 93.4%, null eff_rank=28.86.

## Priority hypotheses for iter_0047

### H1 (high priority): Top singular vector biological identity
**Family**: intrinsic_dimensionality / cross_model_alignment
**Rationale**: SV spectrum result shows 93.4% variance at L11 in top SV. Is this SV biologically meaningful? Test: project each in-vocab gene onto top SV at each layer; check if loadings correlate with known gene properties (expression level, essentiality, lineage marker status).
**Method**: Compute SVD of L11 embedding; for each of 361 in-vocab genes, record its projection onto U[:,0] (top left singular vector). Compare projections for B-cell vs T-cell vs Myeloid vs background. If loadings differ systematically, the dominant direction has biological content.
**Null**: permuted gene labels.

### H2 (high priority): Cross-tissue generalization of ID compression
**Family**: split_robustness
**Rationale**: ID compression replicated across seeds (same tissue). Does it generalize across tissues? Test on cycle6_lung embeddings and cycle7_external_lung.
**Method**: TwoNN ID at each of 12 layers for cycle6_lung_main and cycle6_lung_seed43, compare profile to cycle4_immune_main. Key: do both show monotone compression? Do final-layer IDs agree?
**Novelty**: new biological context (lung vs immune).

### H3 (medium priority): Bootstrap persistence of topological features
**Family**: topology_stability
**Rationale**: Brainstormer suggested bootstrap persistence of topological features. Test: run TwoNN on 30 bootstrap subsamples (n=500) for each layer of cycle4_immune_main; report mean ± std ID per layer. This gives distribution of ID estimates, quantifying how much sampling noise matters.
**Method**: Already partially done (H01 this iter). Extend to more subsamples and report layer-wise variance to assess whether monotone ordering is stable across all bootstrap draws (not just CI endpoints).

### Fallback: Persistent Betti-0 / connected components on kNN graph
**Family**: persistent_homology
**Rationale**: Untested this project. Simple graph: build k=5 kNN graph on in-vocab gene embeddings (n=361), vary k from 2 to 20, count connected components. Compare across layers. If biologically coherent clusters appear at small k, this is evidence for community structure.
**Method**: scipy sparse kNN graph + connected_components. Fast and cheap.

## Notes on SV spectrum finding
The SVD eff_rank result (1.64 at L11) is stronger than expected. Note: this collapses global linear structure, while TwoNN captures local non-linear geometry. The two measures are complementary and both show compression. For the paper, report both and note that the global linear structure collapses faster (14× for SV eff_rank vs 1.8× for TwoNN ID). The dominant direction at L11 does NOT discriminate lineages (H03 negative) — it may reflect a shared "gene embedding" axis learned from the masked gene prediction task.

## Architecture notes (for iter_0047 scripts)
- Embeddings: `[OUTPUTS]/cycle4_immune_main/layer_gene_embeddings.npy` shape [12, 4941, 512]
- In-vocab genes: `[OUTPUTS]/cycle4_immune_main/cycle1_edge_dataset.tsv` (361 unique genes, mapped by source_idx/target_idx)
- Valid gene rows: ~2039 (non-zero norm); 2902 zero-norm rows (padding)
- Conda env: `subproject40-topology`
