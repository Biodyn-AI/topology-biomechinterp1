# External Data Sources

This document describes all external datasets required to reproduce the analyses.

## 1. scGPT Model

- **Source**: [scGPT GitHub](https://github.com/bowang-lab/scGPT)
- **Reference**: Cui, H., Wang, C., Maan, H., et al. (2024). scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods*, 21(8):1470-1480.
- **Model**: Pre-trained whole-human model (12 transformer layers, 512-dimensional hidden states)
- **Usage**: Forward pass to extract residual-stream hidden states at each layer

## 2. Tabula Sapiens

- **Source**: [Tabula Sapiens](https://tabula-sapiens-portal.ds.czbiohub.org/)
- **Reference**: Jones, R.C., Karkanias, J., Krasnow, M.A., et al. (2022). The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. *Science*, 376(6594):eabl4896.
- **Usage**: Immune-lineage cells used as input to scGPT for embedding extraction
- **Format**: H5AD (AnnData)
- **Note**: The `feature_name` column has Categorical dtype; use `.astype(str)` before `var_names_make_unique()`

## 3. TRRUST v2

- **Source**: [TRRUST](https://www.grnpedia.org/trrust/)
- **Reference**: Han, H., Cho, J.W., Lee, S., et al. (2018). TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. *Nucleic Acids Research*, 46(D1):D199-D202.
- **File**: `trrust_human.tsv`
- **Contents**: 9,396 signed human TF-target regulatory edges
  - 589 unique pairs involving in-vocabulary genes
  - 270 activation-only, 141 repression-only, 178 mixed/ambiguous
- **Format**: TSV with columns: TF, target, regulation_type, PMID

## 4. STRING v12.0

- **Source**: [STRING database](https://string-db.org/)
- **Reference**: Szklarczyk, D., Kirsch, R., Koutrouli, M., et al. (2023). The STRING database in 2023. *Nucleic Acids Research*, 51(D1):D99-D105.
- **Usage**: Protein-protein interaction pairs
  - Primary analyses: combined score >= 0.7 (1,022 in-vocabulary pairs)
  - Confidence gradient: combined score >= 0.4 (3,092 pairs)
- **Access**: Via STRING API or bulk download

## 5. Gene Ontology

- **Source**: [Gene Ontology](http://geneontology.org/)
- **Reference**: Ashburner, M., Ball, C.A., Blake, J.A., et al. (2000). Gene Ontology: tool for the unification of biology. *Nature Genetics*, 25(1):25-29.
- **Usage**: Cellular Component (CC) and Biological Process (BP) terms
  - CC terms for subcellular compartment enrichment
  - BP terms for negative control (591 terms tested, 0 significant)
  - Jaccard index of shared annotations for pairwise gene similarity

## 6. Geneformer (Optional)

- **Source**: [Geneformer on HuggingFace](https://huggingface.co/ctheodoris/Geneformer)
- **Reference**: Theodoris, C.V., Xiao, L., Chopra, A., et al. (2023). Transfer learning enables predictions in network biology. *Nature*, 618(7965):616-624.
- **Usage**: Cross-model alignment analysis (negative finding #2)
- **Model**: Static gene embeddings from the pre-trained token embedding layer

## Pre-processed Edge Dataset

The file `cycle1_edge_dataset.tsv` is a pre-processed dataset that maps TRRUST TF-target pairs to the scGPT gene vocabulary:
- 589 positive pairs (known TRRUST regulatory edges)
- 2,351 negative pairs (same TFs, non-target genes)
- Columns: TF_gene, target_gene, label (1=regulatory, 0=non-regulatory), regulation_type

This file is derived from TRRUST v2 intersected with the scGPT vocabulary (195 in-vocabulary named genes, 51 TFs, 144 targets).
