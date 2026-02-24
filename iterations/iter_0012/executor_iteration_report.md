# Executor Iteration Report — iter_0012

## Summary

Three hypotheses tested. H01 negative (bootstrap CIs on act vs rep co-pole differential do not exclude zero at any layer — sample size n_rep=64 is too small for CI-level resolution). H02 strongly positive (TRRUST activation pairs are significantly more spatially concentrated in SV2 than random at 10/12 layers; continuous distance test validates discrete pole-membership result). H03 strongly positive (STRING PPI co-pole enrichment: 12/12 layers significant, z-scores 3.3–6.5 — cross-graph-type convergent validation of SV2 geometric organization by physical protein-protein interactions).

---

## Command Trace

```bash
# H01, H02 (bootstrap + pairwise distance):
conda run -n subproject40-topology python \
  /Volumes/Crucial\ X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0012/run_iter0012_screen.py

# H03 (STRING PPI, corrected score threshold — score in 0-1 not 0-1000):
conda run -n subproject40-topology python3 -c "
  [inline script: load STRING cache, apply score>=0.7, run co-pole enrichment N=500 shuffles × 12 layers]
"
```

Data inputs:
- `layer_gene_embeddings.npy` [12, 4803, 512] (scGPT residual stream)
- `cycle1_edge_dataset.tsv` (209 named genes with embedding indices)
- `trrust_human.tsv` (9396 TF-target edges)
- STRING v12.0 API response (9905 interactions for 209 genes, cached to `string_ppi_named_genes.json`)

---

## H01: Bootstrap CI on Activation vs Repression Co-pole Differential

**Hypothesis**: The activation > repression co-pole enrichment from iter_0011 (qualitative) becomes a CI-backed quantitative claim when bootstrap-resampling both edge sets independently.

**Method**: Bootstrap resample TRRUST activation (n=116) and repression (n=64) pairs separately at each of 12 layers (N=2000 bootstrap draws). Compute 95% CI on obs_act − obs_rep co-pole rate differential. Report CI coverage of zero and one-sided p (fraction boot diffs ≤ 0).

**Results**:

| Layer | obs_act | obs_rep | obs_diff | 95% CI diff | p(act>rep) | CI excl. 0 |
|-------|---------|---------|----------|-------------|------------|------------|
| 0  | 0.250 | 0.219 | 0.031 | [-0.097, 0.160] | 0.325 | No |
| 1  | 0.250 | 0.188 | 0.062 | [-0.063, 0.184] | 0.170 | No |
| 2  | 0.224 | 0.172 | 0.052 | [-0.064, 0.167] | 0.193 | No |
| 3  | 0.216 | 0.172 | 0.044 | [-0.073, 0.158] | 0.225 | No |
| 4  | 0.241 | 0.234 | 0.007 | [-0.126, 0.139] | 0.447 | No |
| 5  | 0.233 | 0.188 | 0.045 | [-0.081, 0.168] | 0.249 | No |
| 6  | 0.241 | 0.203 | 0.038 | [-0.095, 0.161] | 0.288 | No |
| 7  | 0.233 | 0.203 | 0.030 | [-0.100, 0.154] | 0.307 | No |
| 8  | 0.224 | 0.219 | 0.005 | [-0.121, 0.130] | 0.466 | No |
| 9  | 0.241 | 0.188 | 0.054 | [-0.074, 0.174] | 0.211 | No |
| 10 | 0.224 | 0.172 | 0.052 | [-0.071, 0.166] | 0.187 | No |
| 11 | 0.224 | 0.188 | 0.037 | [-0.083, 0.161] | 0.266 | No |

**Key findings**:
- 0/12 layers have 95% CI excluding zero.
- 0/12 layers show p < 0.05 for the differential.
- The observed differences are consistently positive (act > rep at 10/12 layers) but with wide CIs due to small repression set (n=64) and binary pole membership.
- **Root cause**: Binary co-pole rate is a coarse signal; n_rep=64 gives insufficient power to resolve the act-rep differential with CIs. The per-set empirical p-values from iter_0011 (act enriched 12/12, rep enriched 1/12) are more informative than the differential CI.

**Decision**: NEGATIVE for this specific claim. The signed-geometry qualitative observation survives but bootstrap CI method is underpowered for n=64.

---

## H02: Mean Pairwise SV2 Distance — Spatial Concentration Test

**Hypothesis**: TRRUST activation pairs are geometrically concentrated (shorter SV2 distance) compared to a random-pair null, providing continuous-distance evidence complementary to discrete pole membership.

**Method**: For each of 12 layers: compute SV2 projections for 209 named genes. Compute mean |SV2_i − SV2_j| for TRRUST activation (116) and repression (64) pairs. Null: N=500 gene-label shuffles. Empirical p = fraction null ≤ obs (low p = more concentrated than null).

**Results**:

| Layer | act_dist | null_mean | z_act | p_act | rep_dist | null_mean | z_rep | p_rep |
|-------|----------|-----------|-------|-------|----------|-----------|-------|-------|
| 0  | 2.851 | 3.350 | -1.77 | 0.024 | 3.018 | 3.335 | -0.84 | 0.220 |
| 1  | 2.466 | 2.980 | -1.87 | 0.024 | 2.617 | 2.950 | -1.09 | 0.150 |
| 2  | 2.376 | 2.906 | -2.10 | 0.016 | 2.476 | 2.885 | -1.36 | 0.086 |
| 3  | 2.529 | 3.074 | -1.95 | 0.016 | 2.561 | 3.069 | -1.56 | 0.058 |
| 4  | 2.353 | 2.979 | -2.46 | 0.004 | 2.297 | 2.946 | -2.02 | 0.016 |
| 5  | 2.322 | 2.957 | -2.51 | **0.000** | 2.180 | 2.939 | -2.43 | 0.006 |
| 6  | 2.318 | 2.994 | -2.68 | **0.000** | 2.109 | 2.989 | -2.76 | **0.000** |
| 7  | 2.670 | 3.243 | -1.98 | 0.020 | 2.514 | 3.238 | -2.19 | 0.010 |
| 8  | 3.128 | 3.466 | -1.10 | 0.130 | 2.946 | 3.445 | -1.38 | 0.080 |
| 9  | 2.711 | 3.094 | -1.39 | 0.086 | 2.569 | 3.071 | -1.50 | 0.054 |
| 10 | 2.185 | 2.652 | -1.99 | 0.016 | 2.052 | 2.634 | -2.03 | 0.014 |
| 11 | 1.698 | 2.135 | -2.33 | 0.006 | 1.632 | 2.116 | -2.15 | 0.012 |

**Key findings**:
- Activation pairs significantly more concentrated than null at **10/12 layers** (emp_p < 0.05).
- Repression pairs concentrated at **6/12 layers** — both TRRUST edge types show spatial proximity in SV2.
- Activation pairs have smaller mean pairwise distance than repression (act < rep) at only 4/12 layers, so activation is not uniquely more concentrated than repression.
- **Peak concentration**: layer 5–6 (both act and rep) and layer 11 (final layer) show strongest effect.
- **Interpretation**: Regulatory pairs (both activation and repression) are geometrically proximate in SV2 space. This complements the co-pole result: it is not just a discrete pole co-assignment effect — the continuous SV2 geometry shows regulatory co-localization.

**Decision**: PROMISING (positive, null-controlled, consistent across most layers, continuous complement to pole test).

---

## H03: STRING PPI Co-pole Enrichment

**Hypothesis**: STRING physical protein-protein interaction partners are co-localized in SV2 poles, providing cross-graph-type validation independent of transcriptional regulation.

**Method**: Download STRING v12.0 via API for 209 named genes (9905 total edges returned; score is in 0-1 range). Filter to score ≥ 0.7 (high-confidence): 1022 unique pairs. Run co-pole enrichment test (top-K=52 / bottom-K=52 SV2 poles, N=500 gene-label shuffles) at each of 12 layers. Note: initial script used threshold 700 (wrong units); corrected to 0.7.

**Results**:

| Layer | obs_copole | null_mean | z | emp_p | sig |
|-------|-----------|-----------|---|-------|-----|
| 0  | 0.252 | 0.121 | 5.83 | 0.000 | Yes |
| 1  | 0.266 | 0.121 | 6.52 | 0.000 | Yes |
| 2  | 0.226 | 0.121 | 4.68 | 0.000 | Yes |
| 3  | 0.229 | 0.122 | 4.87 | 0.000 | Yes |
| 4  | 0.220 | 0.123 | 4.35 | 0.000 | Yes |
| 5  | 0.230 | 0.122 | 4.98 | 0.000 | Yes |
| 6  | 0.233 | 0.121 | 5.17 | 0.000 | Yes |
| 7  | 0.213 | 0.122 | 3.98 | 0.000 | Yes |
| 8  | 0.210 | 0.122 | 3.98 | 0.000 | Yes |
| 9  | 0.198 | 0.121 | 3.42 | 0.002 | Yes |
| 10 | 0.198 | 0.123 | 3.26 | 0.002 | Yes |
| 11 | 0.239 | 0.122 | 5.45 | 0.000 | Yes |

**Key findings**:
- **12/12 layers significant** (emp_p ≤ 0.002).
- Z-scores range 3.26–6.52; obs co-pole rate 0.198–0.266 vs null 0.121–0.123 (2× enrichment).
- Peak signal at early layers (L1: z=6.52) and final layer (L11: z=5.45).
- This is **independent of TRRUST** — physical PPI interactions show the same SV2 geometric co-localization as transcriptional regulatory relationships.
- **Converging evidence**: TRRUST TF-target co-pole, mito/ER compartment enrichment, and now STRING PPI co-pole all point to SV2 as a geometrically organized axis of biological function.

**Decision**: PROMISING (strongly positive, large effect size, 12/12 layers, cross-graph-type convergent validation).

---

## Artifacts Generated

| File | Type | Contents |
|------|------|----------|
| `h01_bootstrap_act_rep_diff.json` | JSON | Bootstrap CI × 12 layers, act/rep/diff |
| `h02_sv2_pairwise_dist.json` | JSON | Mean pairwise SV2 distance × 12 layers × act/rep |
| `h03_ppi_copole.json` | JSON | STRING PPI co-pole enrichment × 12 layers |
| `string_ppi_named_genes.json` | JSON | Raw STRING API response (9905 edges) |
| `iter0012_results.json` | JSON | Master results summary |
| `run_iter0012_screen.py` | Python | Experiment script (H01, H02, H03) |

---

## Overall Assessment

- **H01 NEGATIVE**: Bootstrap CI on act-rep differential is underpowered (n_rep=64). The discrete per-set empirical p-value approach from iter_0011 is retained.
- **H02 PROMISING**: Continuous SV2 distance test shows activation pairs are spatially concentrated at 10/12 layers. Repression also concentrated at 6/12 layers. Complementary evidence to co-pole.
- **H03 STRONGLY POSITIVE**: STRING PPI co-pole is significant at all 12 layers (z=3.3–6.5). This is the most important new finding: SV2 geometric organization extends to physical protein interactions, not just TF regulatory edges.

**Top cumulative signal**: SV2 axis organizes at least three types of gene relationships simultaneously: (1) TRRUST transcriptional co-regulation (12/12 layers), (2) subcellular compartment membership (mito/ER/secretory), and (3) STRING physical PPI (12/12 layers, z up to 6.5). This is convergent multi-graph evidence for meaningful SV2 geometry.
