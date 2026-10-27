# A2 — gene-axis verification (gate) — PASSED, independently verified

**Question.** Is each embedding's gene axis the var index of the h5ad it was extracted from?
**Answer.** Yes, for all 16 contexts.

**Evidence (direct).**
- A model-free re-implementation of the extractor's cell sampler
  (`np.sort(default_rng(seed).choice(n_obs, max_cells, replace=False))`) and of the ScGPTDataset tokeniser
  (non-zero genes in the scGPT vocabulary, top `max_genes = 512` by expression) reproduces
  `layer_gene_embedding_count.npy` **exactly at every gene, in every layer, in all 16 contexts**
  (`per_context.csv`: `n_genes_differ = 0`). An independent verifier using `anndata` instead of raw `h5py`
  (`verify/verify_a2.py`) obtained the same exact match.
- Each run's `cycle1_edge_dataset.tsv` name→index map agrees with the context h5ad var index (0 mismatches;
  2,304–65,560 index pairs per context). The non-zero-embedding set equals the count > 0 set everywhere.

**Extraction parameters recovered.** seed 42 / 256 cells: immune, kidney_A, lung, lung_external, three donors,
three cell types; seeds 43 / 44 / 256 cells: immune resamples; **kidney_B: seed 67 / 192 cells**;
**three kidney checkpoints: seed 87 / 192 cells** (recovered by exhaustive seed search, unique match).

**Tissue truth (reproduced cells).** kidney_A: Tabula Sapiens kidney, donor TSP2 (epithelium 213, immune 41,
endothelium 2 of 256). kidney_B and all checkpoints: kidney, TSP2. lung: TS lung, donors TSP1 107, TSP14 68,
TSP25 67, TSP2 14. lung_external: Krasnow lung, donors 1/2/3 = 131/88/37 (Smart-seq2). immune: 17 donors.

**Non-independence.** kidney_B's 8,000-cell file is exactly the first 8,000 rows of kidney_A's file (head
truncation; it therefore keeps all 352 Smart-seq2 cells at rows 0–351). kidney_A and kidney_B are two
extractions of one donor's cells with different HVG axes, not two samples. Immune draw overlaps with seeds,
donors and cell types are at chance levels (largest: seed44, 8 shared vs 3.3 expected, p = 0.018 unadjusted).

**The 512-token truncation matters outside immune.** Share of detected (cell, observed-gene) entries that were
tokenised: immune ≥ 99.6%, kidney_A 96.5%, kidney_B 46.7%, checkpoints 48.2%, lung 32.4%, lung_external 51.5%.
In lung, 1,114 genes detected in the drawn cells were never tokenised.

**Caveats.** For cells with ties at the 512th value, exact count equality makes identical per-cell token sets very
likely but does not strictly prove them (kidney_A 64, kidney_B 188, lung 223 such cells; ≤ 5 per immune context).

**Cache files written:** `cache/<ctx>/extraction_cells.npy`, `nobs_reproduced.npy`, `expr_extraction.npy`;
token masks in `A2_axis_verification/token_masks/`.
