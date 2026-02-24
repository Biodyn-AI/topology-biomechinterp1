# Next Steps — iter_0008

## Confirmed findings requiring follow-up

### SV2 immune-signaling axis (HIGH PRIORITY)
- SV2 top pole: IL-4 positive regulation (GO:0032753, emp_p=0.000)
- SV2 bottom pole: extracellular vesicles (GO:0070062, OR=7.99)
- **Next**: Cross-layer SV2 trajectory (12 layers). Does SV2 emerge gradually or suddenly? What biological term is top at each layer?
- **Next**: SV2 top-pole genes vs TRRUST TF targets (immune TFs: STAT6, GATA3, IRF4).

### Layer-specific compartment transients (HIGH PRIORITY)
- Mitochondrion: dominant at L3 (OR=23); presumably weakens by L11.
- ER lumen: dominant at L7/L8 (OR=10-7); ER lumen was present at L11 too but weaker.
- **Next**: Plot compartment-term enrichment p-value at each of 12 layers for: mitochondrion, ER lumen, extracellular space, nucleus, cytoplasm. This creates a "compartment trajectory" plot — key figure for paper.
- **Next**: Test cytoplasm (GO:0005737) and nucleus (GO:0005634) transients.

### Signal peptide biological anchor (MEDIUM PRIORITY)
- GO-proxy confirmed (emp_p=0.002).
- **Next**: Get UniProt signal peptide annotations (signal_peptide.tsv from UniProt) to test with ground-truth gene set, not proxy.
- **Next**: Test whether signal peptide genes concentrate specifically at SV1 (vs SV2, SV3) — this would prove SV1 is THE secretory axis.

## Hypothesis families still untested
- **persistent_homology**: Betti curves on 209-gene neighborhood at layer-11. Would complement SVD findings.
- **graph_topology**: kNN graph modularity vs GO module co-membership. Retired in iter_0002 (transitivity) but modularity not yet tested.
- **cross_model_alignment**: CKA of scGPT vs Geneformer SV axes. Not yet executed.
- **split_robustness**: All results so far on full 209-gene set. Need target-disjoint or source-disjoint split validation.

## Recommended iter_0009 portfolio
1. **Compartment trajectory across 12 layers** (12 Fisher tests × 5 compartments) — closes main mechanistic narrative
2. **SV2 cross-layer trajectory** (12 SVD + GO enrichment calls) — new finding expansion
3. **Persistent homology (Betti-0/1) on 209-gene embedding at layer-11** — new family, tests topology directly
