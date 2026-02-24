# Executor Next Steps — iter_0035

## Iteration outcomes
- H01 (TRRUST/GO proximity): **negative** — retire annotation→proximity family
- H02 (B-cell kNN precision@10): **strong positive** — B-cell cluster 7x enriched at L2, persists to L11
- H03 (B-cell centroid separation): **positive** — normalized rho=+0.972, increasing separation relative to manifold scale

## Priority for iter_0036

### H-A (HIGH PRIORITY): Cross-model validation — Geneformer B-cell cluster
- Test same B-cell kNN precision@10 on Geneformer embeddings (if available)
- If Geneformer data accessible: compare precision@10 per layer to scGPT
- Family: cross_model_alignment | novelty: new_family

### H-B (HIGH PRIORITY): Expanded B-cell marker set with broader annotation
- Retrieve full B-cell marker list from mygene/GO (not just manually curated 5 genes)
- Re-run H02 with n=15-20 B-cell markers for statistical power
- Expected: stronger precision@10 signal with more markers
- Family: manifold_distance | novelty: refinement of H02

### H-C (MEDIUM): Topology stability — B-cell precision@10 vs kNN scale
- Test precision@10 for B-cell cluster at k={5,10,15,20,30}
- Does cluster coherence scale with k or plateau?
- Family: topology_stability | novelty: new_method

### H-D (MEDIUM): Functional annotation of B-cell geometric neighborhood
- At L11: for each B-cell marker, what genes appear most often in top-10 neighbors?
- Annotate recurring non-B-cell neighbors (GO/STRING lookup)
- Tests whether B-cell cluster co-localizes with functionally related genes
- Family: module_structure | novelty: new_method

## Retired directions
- STRING PPI → embedding distance (iter_0031, iter_0015): negative ×2
- TRRUST/GO co-annotation → embedding proximity (iter_0035 H01): negative ×1 new direction
- PC2/PC3 T-cell/Myeloid axes (iter_0034 H02): negative (insufficient marker count)
- Monotonic community emergence (iter_0034 H01): inconclusive → subsumed by kNN precision

## Key biological anchor
B-cell marker geometric coherence (PC1 signal + kNN clustering + centroid separation) is the strongest reproducible signal. Next iteration should:
1. Cross-validate in Geneformer
2. Expand marker set
3. Characterize neighborhood functional content
