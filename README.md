# Multi-dimensional spectral geometry of biological knowledge in single-cell transformer representations

This repository contains the complete experimental codebase, results, and agent transcripts for the autonomous hypothesis-screening loop described in:

> [Author redacted for double-blind review] (2026). Multi-dimensional spectral geometry of biological knowledge in single-cell transformer representations.

> **Note (October 2026).** The headline findings listed below are those of the automated screen. They are re-examined under a frozen
> analysis convention in `revision3_experiments/` (see [Re-analysis and confirmatory tests](#re-analysis-and-confirmatory-tests-2026)),
> which supersedes them where they differ: the matrix the screen called "immune" is Tabula Sapiens kidney, never-observed genes entered its
> decompositions as zero vectors, the protein-language-model embedding is ESM2-15B (UCE per-protein file), and several values below
> (for example the GRN-edge AUROC of 0.860 and the kidney non-transfer) do not hold.

## Overview

This codebase systematically decodes the geometric structure of scGPT's internal representations using an automated two-agent loop (executor + brainstormer) that iteratively proposes, tests, and retires geometric hypotheses. Over **63 iterations**, the loop tested **183 hypotheses** across **13 families**.

The compiled manuscript (with all post-revision analyses) is at `paper/research_paper_plos.pdf`. The LaTeX source is at `paper/research_paper_plos.tex`. The reviewer-response letter is at `paper/response_letter.md`.

### Headline findings (post-revision, after multiple-testing correction)

Three orthogonal spectral directions carry distinct biological information:

- **SV1** correlates with subcellular localization along the secretory pathway (mitochondria → ER → extracellular).
- **SV2–SV4** correlate with protein–protein interaction networks; pair-level cosine similarity in SV2–SV7 has a small but significant partial correlation with STRING combined score (r ≈ 0.10 at L0 after multi-feature confound control).
- **SV5–SV7** correlate with transcriptional regulatory relationships at early layers; the residual signal after stricter confound control is small (r_rb ≈ 0.06, borderline at p = 0.080).

Eight of eleven audited headline claims survive Bonferroni correction at α = 0.05/183.

### Cross-model and cross-tissue replication

- **TF-vs-target distinction is convergent across foundation models**: AUROC = 0.651 (scGPT-L11), 0.748 (Geneformer V1), 0.754 (Geneformer V2) on the 1,672-gene shared subset; 0.889 (ESM2-3B, on a 1,877-gene scGPT∩ESM2 subset).
- **SV1 subcellular and SV2 PPI structures scale with model size**: full strength in scGPT, partial in Geneformer V2, absent in Geneformer V1.
- **Cross-tissue (kidney) replication**: SV1 (OR = 2.39) and SV2 PPI (z = 11.49) replicate; TF-vs-target AUROC drops to chance, indicating the regulatory distinction is immune-specific given the immune-curated TRRUST set.

### Functional utility

- **GRN-edge inference benchmark**: scGPT's full residual-stream cosine at L0 achieves AUROC = 0.860 on a held-out 20% TRRUST split, beating |co-expression Pearson| (0.793). The spectral SV5–SV7 projection alone underperforms co-expression (AUROC = 0.602).
- **Drug-target ranking**: spectral PPI subspace ranking does not outperform co-expression at the small (n=8) drug-target benchmark; the spectral structure is interpretable but not a competitive standalone ranker.

## Re-analysis and confirmatory tests (2026)

These folders re-derive the screen's claims and test four confirmatory hypotheses fixed in advance.

| folder | contents |
| --- | --- |
| `revision3_experiments/` | Frozen analysis convention `CONVENTION3.md` (SHA-256 in `FROZEN.txt`: `3e3249d1d0a8ce1bc3cf64639fe233b016535456f82432a8bbfcacdd64a05c01`, frozen before any confirmatory test) and `CONVENTION3_ERRATA.md`; shared code (`contexts.py`, `stats3.py`, `prep_cache.py`, `prep_pairs.py`); one folder per analysis, each with its scripts, CSV/JSON outputs, a `SUMMARY.md` and, where applicable, an independent re-derivation in `verify/` |
| `revision2_experiments/` | The earlier round of re-analysis scripts (reproduced and compared in `revision3_experiments/`) |
| `sagmb_revision/` | Figure scripts and figures, the supplementary-table builders (`supplement/build_S*.py`), the supplementary tables as CSV (`supplementary_tables/`) and `supplementary_material.pdf` |

Analyses in `revision3_experiments/`:

| folder | question |
| --- | --- |
| `A2_axis_verification` | gene-axis and tissue identity of every extraction; exact reproduction of each stored per-gene cell count |
| `E1_compression` | spectral concentration with depth (three effective-rank definitions, controls) |
| `E11_compartment` | the early-layer compartment axis; linear probes; sequence-model comparison |
| `E2_ppi`, `F1_h1_compartment`, `E12_string_gradient` | H1: physical-interaction partners vs degree-preserving rewiring; degree, detection, co-expression, compartment and confidence controls |
| `E3_node` | H2: transcription-factor identity (detection-stratified); single directions, band scan, cross-fit |
| `E4a_edges`, `E4a_chipfix`, `E5_confounds`, `E4b_grn`, `E4c_probes` | H3: within-factor edge benchmark, power and equivalence; ChIP-seq standard (`E4a_chipfix` rebuilds its candidates so that they exclude every documented pair); conditional logistic regression; GENIE3/GRNBoost2; supervised probes |
| `E7_rank` | H4: rank-reduced vs full-dimensional representations |
| `E6_crossmodel`, `E8_word2vec`, `E10_crosscontext` | static gene embeddings of other models; a co-occurrence-only Word2Vec model; variance components across extractions |
| `E9_screen` | audit of the 183 screen records (tissue, gene labels, decomposition convention, band provenance, multiple testing) |

**Not included (regenerable).** Per-pair score checkpoints (`work/`), per-context intermediates (`partial/`, `_ckpt/`), Word2Vec corpora and
models, run logs, the cached expression matrices (`cache/*/expr*.npy`) and other `.npy` files, and four ChIP-seq pair files over 10 MB. Each is
produced by the scripts in the same folder from the public data sources listed in `DATA_SOURCES.md` and in the paper's data-availability
statement (Tabula Sapiens; the Travaglini et al. 2020 lung atlas, Smart-seq2 portion, via CZ CELLxGENE; scGPT, Geneformer, tGPT, scFoundation,
AIDO.Cell and ESM2 checkpoints; TRRUST v2; DoRothEA; the Lambert et al. 2018 TF census; STRING v12.0 physical links; ComplexPortal; hu.MAP 2.0;
Gene Ontology; Human Protein Atlas). Scripts resolve paths relative to a working copy laid out as in the original project; adjust the
paths in `revision3_experiments/contexts.py` to point at local copies of these sources.

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
├── figures/                      # Publication-quality figures
│   ├── fig1_joint_vs_single_auroc.png
│   ├── fig2_cross_seed_robustness.png
│   ├── fig3_edge_auroc_depth_decay.png
│   ├── fig4_signed_regulation_split.png
│   ├── fig5_spectral_and_trrust.png
│   ├── fig6_subspace_seed_angles.png
│   ├── fig7_cross_tissue.png            # Revision: cross-tissue (immune vs kidney) bar chart
│   ├── fig8_cross_model.png             # Revision: scGPT vs Geneformer V1/V2 vs ESM2
│   └── fig9_grn_benchmark.png           # Revision: GRN AUROC benchmark with bootstrap CIs
│
├── paper/                        # Manuscript and revision deliverables
│   ├── research_paper_plos.tex   # LaTeX source (anonymized)
│   ├── research_paper_plos.pdf   # Compiled PDF (anonymized)
│   ├── response_letter.md        # Response to peer reviewers
│   └── figures_medium/           # Manuscript-embedded figures
│
└── revision_experiments/         # Post-review additional analyses
    ├── e1_cross_tissue/          # Kidney replication
    ├── e2_cross_model/           # Geneformer V1/V2 cross-model comparison
    ├── e3_full_vocab/            # Broader-vocabulary replication
    ├── e4_string_continuous/     # Continuous PPI regression (replaces n=5 quintile)
    ├── e5_multiple_testing/      # BH-FDR + Bonferroni audit across 183 hypotheses
    ├── e6_loop_audit/            # Two-agent loop selection-bias audit
    ├── e7_confounds/             # Multi-feature confound residualisation
    ├── e8_dynamics/              # Pseudotime dynamics test
    ├── e9_vocab_audit/           # 4,803 / 209 / 195 vocabulary derivation
    ├── e10_plm_validation/       # ESM2 protein-LM cross-validation
    ├── e11_baselines/            # Word2Vec and other simpler baselines
    ├── e12_grn_benchmark/        # Held-out TRRUST-edge AUROC benchmark
    └── e13_drug_targets/         # Drug-target ranking case study
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
@article{spectral2026anonymous,
  title={Multi-Dimensional Spectral Geometry of Biological Knowledge
         in Single-Cell Transformer Representations},
  author={[Author redacted for double-blind review]},
  year={2026}
}
```

## Related Work

This project builds directly on:

> [Author redacted for double-blind review] (2026). Systematic evaluation of single-cell foundation model interpretability reveals attention captures co-expression rather than unique regulatory signal. arXiv:2602.17532.

## License

MIT License. See [LICENSE](LICENSE).
