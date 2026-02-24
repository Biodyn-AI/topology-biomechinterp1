# Multi-Dimensional Spectral Geometry of Biological Knowledge in Single-Cell Transformer Representations

This repository contains the complete experimental codebase, results, and agent transcripts for the autonomous hypothesis-screening loop described in:

> Kendiukhov, I. (2026). Multi-Dimensional Spectral Geometry of Biological Knowledge in Single-Cell Transformer Representations.

## Overview

We systematically decode the geometric structure of scGPT's internal representations using an automated two-agent loop (executor + brainstormer) that iteratively proposes, tests, and retires geometric hypotheses. Over **63 iterations**, the loop tested **183 hypotheses** across **13 families**, discovering that scGPT organizes genes into a structured biological coordinate system:

- **SV1**: Subcellular localization along the secretory pathway (mitochondria -> ER -> extracellular)
- **SV2-SV4**: Protein-protein interaction networks with quantitative fidelity to STRING confidence
- **SV5-SV7**: Transcriptional regulatory relationships, with co-expression-independent signal at early layers

Key findings include:
- 14.4-fold effective rank collapse across 12 transformer layers
- Perfect monotonic correlation between STRING confidence quintiles and geometric proximity
- TF-vs-target classification AUROC = 0.744 in a 6D spectral subspace
- B-cell GC master regulators converge toward PAX5 across transformer depth
- Attention and residual-stream geometry encode complementary biological relationships

## Repository Structure

```
.
├── README.md                     # This file
├── METHODS.md                    # Detailed methods and reproducibility guide
├── DATA_SOURCES.md               # External data dependencies and access
├── requirements.txt              # Python dependencies
├── .gitignore
│
├── loop/                         # Autoloop orchestration infrastructure
│   ├── run_claude_topology_autoloop.py   # Main two-agent loop runner
│   ├── config_claude.json        # Claude-specific configuration
│   ├── config.json               # Default configuration
│   ├── start_claude_autoloop.sh  # Start the loop
│   ├── stop_claude_autoloop.sh   # Stop the loop
│   └── status_claude_autoloop.sh # Check loop status
│
├── prompts/                      # Agent instruction templates
│   ├── executor_prompt_topology_hypothesis_screening.md
│   └── brainstormer_prompt_template.md
│
├── iterations/                   # All 63 iteration artifacts
│   ├── iter_0001/                # First iteration
│   │   ├── iteration_meta.json   # Execution metadata (timing, return codes)
│   │   ├── executor_hypothesis_screen.json  # Structured hypothesis results
│   │   ├── executor_iteration_report.md     # Methods and results narrative
│   │   ├── run_*.py              # Executable analysis scripts
│   │   ├── *.csv / *.json        # Numerical results
│   │   └── brainstormer_*.md     # Strategic direction documents
│   ├── iter_0002/
│   │   ...
│   └── iter_0063/                # Final iteration
│
├── reports/                      # Aggregated analyses
│   ├── autoloop_master_log.md    # Comprehensive log of all 63 iterations
│   ├── make_claude_visualizations.py  # Figure generation script
│   └── figures/                  # Report-level visualizations
│
└── figures/                      # Publication-quality figures
    ├── fig1_joint_vs_single_auroc.png
    ├── fig2_cross_seed_robustness.png
    ├── fig3_edge_auroc_depth_decay.png
    ├── fig4_signed_regulation_split.png
    ├── fig5_spectral_and_trrust.png
    └── fig6_subspace_seed_angles.png
```

## Quick Start

### Prerequisites

1. **Python environment** (conda recommended):
   ```bash
   conda create -n topology-screen python=3.10
   conda activate topology-screen
   pip install -r requirements.txt
   ```

2. **External data** (see [DATA_SOURCES.md](DATA_SOURCES.md) for details):
   - scGPT pre-trained model weights
   - Tabula Sapiens immune-lineage cells
   - TRRUST v2 regulatory network
   - STRING v12.0 protein interactions
   - Gene Ontology annotations

3. **Pre-computed embeddings**: The analysis scripts in `iterations/` reference pre-extracted scGPT residual-stream embeddings stored as NumPy arrays (`layer_gene_embeddings.npy`, shape `[12, 4803, 512]`). See [METHODS.md](METHODS.md) for the embedding extraction procedure.

### Reproducing Individual Analyses

Each iteration directory contains self-contained analysis scripts. For example, to reproduce the graph topology screen from iteration 1:

```bash
cd iterations/iter_0001
python run_graph_topology_screen.py
```

> **Note**: Scripts contain absolute paths from the original execution environment. You will need to update the `SUBP38_OUTPUTS` and similar path variables to point to your local data directory. See [METHODS.md](METHODS.md) for path configuration instructions.

### Running the Full Autoloop

The complete hypothesis-screening loop requires the Claude CLI tool:

```bash
cd loop
bash start_claude_autoloop.sh
```

The loop alternates between executor (running experiments) and brainstormer (proposing new hypotheses) agents, producing all artifacts in `iterations/`.

## Hypothesis Families

The 13 hypothesis families tested across 63 iterations:

| Family | Outcome | Key Result |
|--------|---------|------------|
| Persistent homology | Partial | Positive under feature-shuffle; negative under rewiring nulls |
| Graph topology | Positive | Clustering coefficient z > 9 at all 36 tests |
| Intrinsic dimensionality | Validated | 44.6% TwoNN reduction; 14.4x ER collapse |
| Cross-model alignment | Partial | Geneformer PPI replication (p = 7.8e-127); B-cell absent |
| SVD biological axes | Validated | Three orthogonal biological axes; null-controlled |
| PPI network encoding | Validated | rho = 1.000 confidence gradient; multi-dimensional |
| Cell-type/family clustering | Validated | AUROC 0.851; HLA-I perfect (1.000) |
| Attention-SVD dissociation | Validated | Attention = TF regulation; SVD = PPI; complementary |
| Regulatory geometry (SV2-4) | Revised | Co-expression confounded; encodes gene class |
| Regulatory geometry (SV5-7) | Validated | Co-expression-independent; bootstrap-confirmed |
| Edge-level geometry | Validated | AUROC 0.602; monotonic depth decay |
| Signed regulation | Positive | Repression > activation; pending replication |
| B-cell attractor dynamics | Validated | GC-TF convergence; GC-plasma orthogonality |

## Key Iteration Map

The 63 iterations can be grouped into research phases:

- **Iterations 1-10**: Broad exploration (graph topology, persistent homology, manifold distances)
- **Iterations 11-20**: SVD axis discovery (SV1 localization, SV2 PPI, SV3 immune signaling)
- **Iterations 21-30**: Confound detection (OOV genes, co-expression, annotation density)
- **Iterations 31-45**: Regulatory geometry (TF classification, edge-level AUROC, signed regulation)
- **Iterations 46-55**: Cross-seed replication and B-cell attractor dynamics
- **Iterations 56-63**: Mechanistic probes (subspace dissociation, attention comparison, final validation)

See `reports/autoloop_master_log.md` for the complete narrative of all iterations.

## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{kendiukhov2026spectral,
  title={Multi-Dimensional Spectral Geometry of Biological Knowledge
         in Single-Cell Transformer Representations},
  author={Kendiukhov, Ihor},
  year={2026}
}
```

## Related Work

This project builds directly on:

> Kendiukhov, I. (2026). Systematic evaluation of single-cell foundation model interpretability reveals attention captures co-expression rather than unique regulatory signal. arXiv:2602.17532.

## License

MIT License. See [LICENSE](LICENSE).
