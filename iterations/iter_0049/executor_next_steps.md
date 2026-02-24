# Executor Next Steps: iter_0049 → iter_0050

## Headline finding from iter_0049

**H03: SV2-SV4 subspace encodes TRRUST regulatory proximity** (8/12 layers significant, best L8 p<0.0001, effect=0.268). SV1 is explicitly disassociated from TF regulation (p>0.1 at most layers). This is publication-quality and mechanistically novel.

---

## Top priorities for iter_0050

### P1 (Mandatory): Expand H03 to SV5-SV10 range
- Q: Is the effect specific to SV2-SV4, or does it extend further into the spectrum?
- Test SV5-SV7, SV8-SV10 subspaces with same Mann-Whitney proximity test
- Expected: diminishing effect at higher SVs → identifies the "regulatory encoding" spectral window

### P2 (High): Replicate H03 at seed replicates (cycle4_immune_seed43, seed44)
- If seed43/seed44 dirs contain layer_gene_embeddings.npy, run identical SV2-SV4 proximity test
- This would establish reproducibility across seeds → meets evidence standard
- Command: `ls /Volumes/.../cycle4_immune_seed43/` then run same analysis

### P3 (High): SV1 biological identity via top gene movers
- From H01: SV1 gene loading vector changes from L0 to L11 (cosine = 0.359)
- Compute: for each gene, |SV1_score_L11 - SV1_score_L0| → rank by drift magnitude
- Map to known biology using 295 named genes + GO term enrichment of top/bottom genes
- This addresses the open question: what does SV1 actually capture?

### P4 (Medium): SV2-SV4 FDR correction and effect size profile
- Apply Benjamini-Hochberg FDR to p-values across 12 layers (H03 SV2-4 column)
- Plot effect size by layer to show the L5-L8 peak pattern
- Report effect sizes vs. layer as a discovery narrative

### P5 (Medium, new family): Bootstrap persistence of H03 signal
- Subsample 80% of positive pairs 50 times, recompute Mann-Whitney at each layer
- Confirms H03 is not driven by outlier pairs

---

## Retired directions (do not repeat without rescue)
- Graph topology (kNN purity, centroid cosine similarity): retired iter_0045-0046
- Manifold distance (BCL6 attractor, lineage centroid tracking): retired iter_0043
- TwoNN breakpoint in NK cells: retired iter_0047

## Quick wins for iter_0050 (if blocked)
- Check seed replicate embeddings exist: `ls /path/cycle4_immune_seed43/`
- Use 295 named genes to compute full SV1-SV20 range test in one script
- Run persistent homology (Betti curves) on the SV2-SV4 projected gene point cloud — this is UNTESTED and family=persistent_homology, novel

## Key data paths
- Main embedding: `.../cycle4_immune_main/layer_gene_embeddings.npy` [12, 4941, 512]
- Seed43: `.../cycle4_immune_seed43/` (check if exists)
- Seed44: `.../cycle4_immune_seed44/` (check if exists)
- Edge dataset (gene name → index): `.../cycle4_immune_main/cycle1_edge_dataset.tsv`
- SV1 vectors saved: `iterations/iter_0049/h01_sv1_vectors.npy` [12, 2039]
- H03 results: `iterations/iter_0049/h03_sv_ppi_proximity.csv`
