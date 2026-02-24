# Next Steps — iter_0032

## Priority Portfolio

### P1 (high): Community detection on kNN graph → biological annotation
- Family: graph_topology / module_structure
- Method: For 195 in-vocab genes at L11 (most modular layer), build k=10 kNN graph. Apply greedy modularity (networkx) or Louvain. Compute NMI between graph communities and:
  - GO BP family clusters (Jaccard similarity threshold)
  - Cell-type marker groups (T_cell/B_cell/fibroblast/epithelial from iter_0023)
  - STRING hub genes (degree >= median)
- Sweep layers L8–L11 to check stability
- This tests if the modular kNN structure (H01) maps to known biology

### P2 (high): GO Jaccard → embedding distance on 195 in-vocab genes
- Family: manifold_distance
- Method: Use GO BP Jaccard similarity from iter_0023/iter_0024 cache. For 195 in-vocab named genes, compute pairwise GO Jaccard. Spearman rho(GO_Jaccard, L2_embedding_dist) per layer.
- This is the clean replacement for STRING→distance (H02 was negative)
- Prior results (iter_0024 H03) were promising but on full 209-gene set including OOV

### P3 (medium): PC1 alignment at L11 vs biological labels
- Family: intrinsic_dimensionality / module_structure
- Method: Project 195 in-vocab genes onto PC1 of L11 embeddings. Test if PC1 separates:
  - Cell-type markers (T_cell vs B_cell vs other)
  - Transcription factors vs targets
  - High-STRING-degree vs low-degree genes
- This follows up H03 finding (PC1 accounts for 26% of variance at L11)

## Retired Directions
- STRING→embedding distance (H02 this iteration): definitively null at all layers
- OOV artifact analysis: resolved

## Data notes
- Use 195 in-vocab genes throughout (exclude 14 OOV)
- GO annotation cache likely in iter_0024 directory
- STRING cache: `iterations/iter_0015/string_ppi_score04_cache.json` (format: d['pairs'] list of {g1,g2,score})
