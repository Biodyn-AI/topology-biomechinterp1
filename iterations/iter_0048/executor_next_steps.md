# Executor Next Steps: iter_0048 → iter_0049

## High priority (strong upside)

### 1. SV1 axis biological identity — what IS it?
- All immune cell-type markers score near zero on SV1 at L11 (~-0.01).
- The top SV1 genes are unnamed (not in edge annotation). Need to identify them.
- Approach: Use ALL available gene name sources (scGPT model vocab, AnnData, other tsv files) to name the ~2039 nonzero genes. Then characterize SV1 top/bottom genes via GO or STRING.
- Priority: HIGH — SV1 explains 18.6% variance at L11 (largest single axis); understanding what it captures is scientifically critical.

### 2. Spectral decay biological split — by known cell-type groups
- H02 showed global eff_rank 236.9→48.7 (4.87×). Does compression rate differ by gene type?
- Compute spectral decay separately for: (a) B-cell markers, (b) T-cell markers, (c) GC TFs, (d) myeloid, (e) random matched-size controls.
- Compare eff_rank trajectories across groups.
- Priority: HIGH — would be first demonstration of biologically differentiated spectral compression.

### 3. SV1 trajectory across 12 layers for annotated genes
- Track each named gene's SV1 score across layers (not just L11). Do B-cell genes converge on a shared SV1 value? Do myeloid genes?
- Analogous to the manifold distance tracking done in iter_0040-0043 but in spectral space.
- Priority: MEDIUM.

### 4. Wasserstein + null comparison
- H03 SWD is inconclusive without a null. Compute SWD for shuffled embeddings (random permutation of gene assignments within each layer) and compare to real SWD.
- If real SWD < shuffled SWD, the layer-to-layer transport is constrained by embedding structure.
- Priority: MEDIUM — necessary to interpret H03.

## Retired / do not revisit
- manifold_distance (centroid-proximity tests): retired (>=2 negative/inconclusive outcomes), confirmed retired
- graph_topology (kNN purity, inter-lineage cosine): retired (>=2 negative), confirmed retired
- ID phase transition: retired (linear, no story)
- Inflated SVD claims from iter_0046: corrected and retired

## Notes on SV1 interpretation
- SV1 NOT a cell-type axis (all immune markers score ~0)
- Top SV1 genes are in the unnamed majority — they are biological real genes but not in the edge annotation (TRRUST-based immune gene set)
- Possible interpretation: SV1 at L11 = metabolic/housekeeping axis, or gene length/GC-content confound, or training data frequency
- Need gene name recovery to distinguish
