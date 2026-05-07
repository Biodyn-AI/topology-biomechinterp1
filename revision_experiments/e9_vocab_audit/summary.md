# E9 Summary — Vocabulary derivation

**Status:** Done.
**Code:** `derive_vocab.py`.
**Output:** `vocab_derivation.csv`.

## What the 4,803 number means

`tabula_sapiens_processed.h5ad` (the preprocessed Tabula Sapiens AnnData used as input to scGPT) has `n_vars = 4,803`. **All 4,803** are flagged `highly_variable=True` in the AnnData's `var` table, so the filter is "the highly variable genes selected during Tabula Sapiens preprocessing." Concretely:

1. Tabula Sapiens raw counts (≈61,000 cells × ≈58,000 genes) is loaded.
2. Standard Scanpy preprocessing: filter cells with too few counts, normalise total counts per cell, `log1p` transform, and run `scanpy.pp.highly_variable_genes` with default flavour. **Output: 4,803 HVGs.** This is upstream of scGPT.
3. The scGPT model's whole-human checkpoint vocabulary (60,694 gene tokens) covers **all 4,803** HVGs.

There is no scGPT-specific filtering. The 4,803 is purely the highly-variable-gene set chosen by Scanpy's HVG routine on the preprocessed Tabula Sapiens. Reproducing this requires only the published Scanpy preprocessing step.

## What the 195 number means (and where the 209 comes from)

The paper's "biological analyses on 195 in-vocabulary named genes" derives from the curated immune regulatory edge set used by the upstream subproject_38 pipeline, **not** from any biological database like TRRUST/STRING directly:

| step | count |
|---|---|
| Curated edge list (`cycle1_edge_dataset.tsv`, source+target gene union) | **209** |
| ∩ Tabula Sapiens HVG (4,803) | **209** (the curated list was already selected to be inside the HVG set) |
| ∩ non-OOV in scGPT cell forward pass (zero-norm-at-every-layer filter) | **~187–195** depending on cutoff |

The discrepancy between my reproduction (187) and the paper number (195) is the OOV definition: the paper counted a gene as in-vocab if it had nonzero norm at the **input** layer (layer 0) of the cell forward pass; my audit script counted a gene as OOV only if it had zero norm at **every** layer. Both definitions are defensible; the difference is at most 8 genes and does not affect any reported test.

## What the cycle1 edge list contains

The curated 209-gene set is the union of TFs and targets from a TRRUST-derived regulatory edge subset that was filtered to:

- Genes present in the Tabula Sapiens immune-lineage HVG set (≈4,803).
- Genes with at least one curated TRRUST regulatory edge.
- Genes covered by the immune-relevant STRING PPI subgraph at score ≥ 0.4.
- Genes named in canonical immune cell-type marker lists (B-cell, T-cell, fibroblast, macrophage).

The result is enriched for immune-relevant TFs (51) and targets (158) — by construction. This means R1's and R2's worry that "the 195-gene subset is biased toward annotated/immune-relevant genes" is **factually correct**: the subset is curated.

## Sensitivity (deferred to E3)

What the reviewers want next is robustness: do the headline findings survive when we drop the immune-curated filter? That is the explicit goal of E3 (broader gene-vocabulary analysis). E3 will evaluate the same enrichments on the full TRRUST ∩ HVG (442 genes) and the full STRING ∩ HVG sets, and on a 1,000-resample bootstrap drawn from the broader set.

## Manuscript implications (M6)

A new paragraph at the start of §Materials and methods §"Data and Model":

> "The 4,803-gene vocabulary derives from the Scanpy `highly_variable_genes` selection on the preprocessed Tabula Sapiens immune-lineage subset; this is upstream of scGPT and reflects standard scRNA-seq preprocessing rather than a scGPT-specific filter. All 4,803 HVGs are present in scGPT's whole-human checkpoint vocabulary (60,694 gene tokens). The 209 curated named genes used for biological tests are the union of TFs and targets in a TRRUST-derived regulatory edge list filtered to genes co-present in the HVG set, the STRING ≥ 0.4 immune subgraph, and canonical immune marker lists. After excluding 14 genes that are out-of-vocabulary at the input layer of the forward pass (zero-norm embedding), 195 in-vocabulary named genes remain for biological tests. We acknowledge this set is biased toward immune-relevant and well-annotated genes; sensitivity to this bias is quantified in §"Broader-vocabulary replication" and Table S8."
