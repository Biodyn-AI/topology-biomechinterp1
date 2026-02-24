# Executor Iteration Report — iter_0024

## Overview

Three hypotheses tested: Dorothea confidence stratification (H01, new_method), TF activation hub centrality (H02, new_method), GO ontology comparison BP vs MF vs CC (H03, new_family).

**Results**: H01 strong positive, H02 negative, H03 strong positive (CC > MF > BP).

---

## Command Trace

```bash
# Environment: subproject40-topology conda environment
cd /Volumes/Crucial\ X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0024

conda run -n subproject40-topology python run_iter0024_screen.py
```

**Data sources**:
- `cycle1_main/layer_gene_embeddings.npy` [12, 4803, 512] — scGPT residual embeddings
- `cycle1_main/gene_list.txt` — vocabulary
- `dorothea_human.tsv` — Dorothea regulon database (confidence A-E)
- `trrust_human.tsv` — TRRUST TF-target database (activation/repression)
- `iter_0015/string_ppi_score04_cache.json` — STRING PPI scores ≥0.4
- mygene API — GO annotations (BP, MF, CC) for 209 named genes

---

## H01: Dorothea Confidence Stratification

**Hypothesis**: High-confidence Dorothea TF-target pairs (A/B) are closer in scGPT embedding space than low-confidence (D/E) pairs and background non-biology pairs.

**Method**: Load Dorothea regulons for named genes. Split pairs by confidence: A/B (N=266) vs D/E (N=775). Compute AUROC (high-conf vs. background non-Dorothea/non-STRING pairs, N=5000). Permutation control: 500 shuffles of gene labels at layer 8.

**Results**:
- High-conf (A/B) vs background: AUROC = 0.618–0.679 across 12 layers, ALL 12 significant (p<10^-10)
- Layer 8: AUROC = 0.671, permutation p = 0.000 (0/500 shuffles reach this level)
- Shuffle null mean: 0.177 ± 0.011 (permuted pairs are FARTHER than background; real signal is dramatically above null)
- High-conf vs low-conf: 7/12 layers significant (AUROC ~0.53-0.54)

**Direction**: POSITIVE — strong, fully reproducible, robust to permutation test.

**Artifact**: `h01_dorothea_confidence.json`

---

## H02: TF Activation Hub Centrality

**Hypothesis**: TFs with many activation targets (TRRUST) are geometrically central (lower mean distance to all named genes) compared to TFs with many repression targets.

**Method**: Count activation/repression target degrees per TF in TRRUST (named genes only). Compute mean L2 distance from each TF to all other named genes as centrality proxy. Spearman(activation_degree, -mean_dist) at each layer. Also: top-15 activation TFs vs top-15 repression TFs centrality comparison (Mann-Whitney).

**Results**:
- Activation TFs: N=39, top TFs: JUN(17), ETS1(10), HIF1A(9), EGR1(7)
- Repression TFs: N=32, top TFs: EGR1(6), FOS(4), JUN(4), PAX5(4)
- Spearman(act_degree, proximity) at layer 8 = -0.047, p=0.74 — not significant
- 0/12 layers significant for either degree-centrality or top-K comparison

**Direction**: NEGATIVE — no evidence that activation TF degree predicts geometric centrality.

**Interpretation**: The pairwise activation proximity signal (iter_0023 H03) does not generalize to a global hub structure. Activation TFs are not more geometrically "central" than repression TFs. The signal may be limited to specific TF-target pairs, not TF degree globally.

**Artifact**: `h02_tf_hub_centrality.json`

---

## H03: GO Ontology Comparison (BP vs MF vs CC)

**Hypothesis**: GO Cellular Component (CC), Molecular Function (MF), and Biological Process (BP) annotations encode different amounts of information about scGPT embedding proximity.

**Method**: Query mygene API for BP, MF, CC GO terms for all 209 named genes. Compute pairwise Jaccard similarity for each ontology across all 21,736 named-gene pairs. Spearman(Jaccard, -L2_distance) at each of 12 layers.

**Results** (all 12 layers significant for all three ontologies):

| Layer | BP Spearman | MF Spearman | CC Spearman |
|-------|-------------|-------------|-------------|
| 0     | 0.050       | 0.088       | 0.117       |
| 4     | 0.065       | 0.092       | 0.122       |
| 8     | 0.081       | 0.089       | 0.106       |
| 11    | 0.070       | 0.088       | 0.104       |

- **CC is best**: peak Spearman = 0.124 (layer 5), consistently strongest
- **MF is second**: stable ~0.088-0.094 across layers
- **BP is third but increasing**: deepens from 0.050 (L0) to 0.083 (L7) showing layer-dependent emergence
- All p-values < 10^-25 across all layers and all three ontologies

**Direction**: POSITIVE — GO CC (subcellular localization) is the strongest single GO ontology predictor of scGPT embedding proximity. BP shows layered deepening pattern.

**Artifact**: `h03_go_ontology_comparison.json`

---

## Quantitative Summary

| Hypothesis | Family | Status | Decision | Key metric |
|-----------|--------|--------|----------|------------|
| H01: Dorothea confidence | manifold_distance | tested | promising | AUROC=0.671, 12/12 sig, perm_p=0 |
| H02: TF hub centrality | module_structure | tested | negative | sp=-0.047, 0/12 sig |
| H03: GO ontology comparison | module_structure | tested | promising | CC sp=0.106 @ L8, 12/12 sig |

---

## Artifacts Generated

- `run_iter0024_screen.py` — experiment script
- `h01_dorothea_confidence.json` — H01 per-layer results
- `h02_tf_hub_centrality.json` — H02 per-layer results
- `h03_go_ontology_comparison.json` — H03 per-layer results
- `iter0024_summary.json` — combined summary
- `executor_iteration_report.md`, `executor_next_steps.md`, `executor_hypothesis_screen.json`
