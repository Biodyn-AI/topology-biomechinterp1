# Methods and Reproducibility Guide

This document provides detailed instructions for reproducing all analyses in the paper.

## 1. Embedding Extraction

All analyses operate on pre-extracted residual-stream embeddings from scGPT. The extraction procedure is:

1. **Model**: scGPT (12 transformer layers, 512-dimensional hidden states) from [scGPT repository](https://github.com/bowang-lab/scGPT).
2. **Input data**: Immune-lineage cells from Tabula Sapiens. For each cell, scGPT processes a rank-ordered list of gene expression values as input tokens.
3. **Forward pass**: Extract the residual-stream hidden state at the output of each transformer layer, obtaining a 512-dimensional embedding vector per gene per layer.
4. **Averaging**: Average gene embeddings across all input cells to obtain a single representative `G x 512` matrix per layer, where `G = 4,803` (number of gene positions in the input vocabulary).
5. **Output format**: NumPy array of shape `[12, 4803, 512]` saved as `layer_gene_embeddings.npy`.

Three independent extractions were performed with different random seeds for cell sampling (`seed42`, `seed43`, `seed44`).

### Gene Vocabulary

Of the 4,803 positions, 209 correspond to named genes present in biological annotation databases (TRRUST, STRING, Gene Ontology). Of these, 14 were discovered to be out-of-vocabulary (OOV) tokens producing identical zero-norm embeddings, leaving **195 in-vocabulary named genes** for biological analyses.

## 2. Path Configuration

The analysis scripts in `iterations/` contain absolute paths from the original execution environment. To run them locally, you need to set up a data directory and update paths.

### Required Data Layout

```
data/
├── embeddings/
│   ├── cycle4_immune_main/
│   │   └── layer_gene_embeddings.npy    # [12, 4803, 512]
│   ├── cycle4_immune_seed43/
│   │   └── layer_gene_embeddings.npy
│   └── cycle4_immune_seed44/
│       └── layer_gene_embeddings.npy
├── networks/
│   ├── trrust_human.tsv                 # TRRUST v2 regulatory edges
│   ├── cycle1_edge_dataset.tsv          # Pre-processed TF-target edges
│   └── string_pairs.tsv                 # STRING protein interactions
├── annotations/
│   ├── gene2go_all.pkl                  # Gene Ontology annotations
│   └── go_terms.json                    # GO term metadata
└── geneformer/
    └── gene_embeddings.npy              # Geneformer static embeddings (optional)
```

### Updating Paths in Scripts

In each `run_*.py` script, look for path definitions like:

```python
SUBP38_OUTPUTS = Path("/Volumes/Crucial X6/.../subproject_38.../outputs")
```

Replace with your local data directory:

```python
SUBP38_OUTPUTS = Path("/path/to/your/data/embeddings")
```

The key path variables used across iterations:
- `SUBP38_OUTPUTS` or `BASE`: Points to embedding directory
- `TRRUST_PATH`: Points to `trrust_human.tsv`
- `EDGE_TSV`: Points to `cycle1_edge_dataset.tsv`
- `GENE2GO_PATH`: Points to `gene2go_all.pkl`

## 3. Core Analyses and Their Locations

### Spectral Collapse (Section 2.1)

- **Effective rank**: Computed in multiple iterations; see `iter_0035/` onward for the 195-gene corrected version
- **TwoNN intrinsic dimensionality**: `iter_0043/` and `iter_0044/`
- **Participation ratio**: `iter_0035/` onward

### Biological Axes (Section 2.2)

- **SV1 compartment enrichment**: `iter_0009/`, `iter_0010/` (initial), `iter_0031/` (OOV-corrected)
- **SV2 co-pole enrichment**: `iter_0010/`, `iter_0015/`
- **SV3 immune signaling**: `iter_0011/`
- **Cell-type clustering**: `iter_0048/`, `iter_0049/`
- **Gene family clustering**: `iter_0049/`

### PPI Encoding (Section 2.3)

- **Multi-axis PPI co-pole**: `iter_0015/`, `iter_0016/`
- **Confidence gradient**: `iter_0017/`
- **PPI-only vs GO-only**: `iter_0018/`
- **Hub-degree confound**: `iter_0019/`

### Regulatory Geometry (Section 2.4)

- **TF-vs-target classification**: `iter_0056/` (joint SV2-7), `iter_0057/` (cross-seed)
- **Co-expression residualization**: `iter_0041/`, `iter_0042/`
- **Dual-regime discovery**: `iter_0042/`

### Edge-Level AUROC (Section 2.5)

- **Edge-level depth decay**: `iter_0062/`
- **Signed regulation**: `iter_0063/`
- **Mechanistic probes**: `iter_0063/`

### B-Cell Attractor (Section 2.6)

- **B-cell precision@10**: `iter_0048/`
- **GC-TF convergence**: `iter_0050/`, `iter_0051/`
- **GC-plasma orthogonality**: `iter_0052/`
- **BCL6 metabolic isolation**: `iter_0051/`
- **Lineage-specific compression**: `iter_0053/`

### Negative Findings (Section 2.7)

- **Persistent homology**: `iter_0003/`, `iter_0004/`, `iter_0005/`
- **Cross-model alignment**: `iter_0020/`, `iter_0061/`
- **ER-AUROC confound**: `iter_0055/`
- **GO BP enrichment**: `iter_0012/`
- **Feed-forward loops**: `iter_0039/`

## 4. Statistical Methods

All analyses use the following null models:

1. **Gene-label shuffle**: Permute gene-to-embedding-row assignments (N = 200-1,000 permutations)
2. **Feature shuffle**: Independently permute values within each of 512 embedding dimensions
3. **Degree-preserving rewiring**: Rewire kNN graph edges while preserving degree distribution
4. **Co-expression residualization**: OLS regression of spectral proximity on pairwise co-expression

Key statistical metrics:
- **AUROC**: Area under ROC curve for binary classification tasks
- **Co-pole rate**: Fraction of gene pairs where both genes fall in the same SVD pole (top-K or bottom-K)
- **z-score**: (observed - null_mean) / null_std from permutation distribution
- **Rank-biserial correlation**: Effect size measure for two-group comparisons on ranks

## 5. Software Environment

All experiments were executed in the `subproject40-topology` conda environment. Key packages:

- Python 3.10
- NumPy 1.26.4
- SciPy 1.12+
- scikit-learn 1.4+
- matplotlib 3.8+
- UMAP-learn 0.5+

See `requirements.txt` for the complete dependency list.
