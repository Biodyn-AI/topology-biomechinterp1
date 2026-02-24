# Executor Iteration Report — iter_0013

## Overview

Three hypotheses tested, all using cached data (STRING PPI, TRRUST, scGPT layer embeddings).
Total runtime: ~3 minutes.

## Command Trace

```bash
conda run -n subproject40-topology python \
  iterations/iter_0013/run_iter0013_screen.py
```

Working directory: `/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop`

Input data:
- `subproject_38/.../layer_gene_embeddings.npy` — [12, 4803, 512]
- `subproject_38/.../cycle1_edge_dataset.tsv` — 209 named genes with embedding indices
- `iter_0012/string_ppi_named_genes.json` — STRING v12.0 PPI (9905 edges raw)
- `single_cell_mechinterp/external/networks/trrust_human.tsv` — TRRUST TF-target annotations

Processed inputs:
- 209 named genes × 12 layers × 512 dims
- STRING pairs (score ≥ 0.7, both named): 1022
- TRRUST activation pairs: 113; repression pairs: 62
- SVD at each layer: top/bottom K=52 pole classification

## H01: SV Axis Specificity — STRING PPI Co-pole on SV1/SV2/SV3

**Method**: For each layer (0–11), compute SVD of 209-gene mean-centered embeddings.
Project all genes onto SV1, SV2, SV3 separately. For STRING PPI pairs (score≥0.7),
compute co-pole rate (both genes in top-52 or both in bottom-52).
Compare to null (N=500 random pairs of same size).

**Results**:
| SV Axis | Mean z | Significant layers (p<0.05) |
|---------|--------|------------------------------|
| SV1     | 0.39   | 2/12                         |
| SV2     | 10.04  | 12/12                        |
| SV3     | 7.18   | 12/12                        |

**Interpretation**: SV2 is the strongest axis (mean z=10.04), but SV3 is also consistently
significant (mean z=7.18). SV1 shows no consistent signal. This means the PPI geometry
effect spans at least two principal axes (SV2 and SV3), with SV2 dominant.
The SV2 signal from iter_0012 (z≈3–6) was likely underestimated due to using smaller N_shuffle.
With N=500, SV2 z values are 7.3–14.7.

**Decision**: **Mixed/positive** — SV2 is the primary axis but not uniquely so. SV3 carries
substantial information. This refines the claim: scGPT SV2+SV3 jointly encode PPI proximity.

**Artifact**: `h01_sv_axis_specificity.json`

## H02: Repression Anti-pole (Cross-pole Rate Test)

**Method**: For each layer, using SV2 projections, compute:
- Activation co-pole rate (both TF and target in same pole)
- Repression co-pole rate
- Repression cross-pole rate (TF in top pole, target in bottom or vice versa)
- Activation cross-pole rate (control)
Compare each to null (N=500 random pairs).

**Results** (mean z across 12 layers):
| Metric            | Mean z | N_layers z>1.5 |
|-------------------|--------|-----------------|
| ACT co-pole       | 3.18   | 12/12           |
| REP co-pole       | 1.26   | 5/12            |
| REP cross-pole    | -1.41  | 0/12            |
| ACT cross-pole    | -1.11  | 0/12            |

**Interpretation**: The repression anti-pole hypothesis is **negative**. Repression pairs do NOT
preferentially occupy opposite poles. The directional regulatory sign (activation vs repression)
is not encoded as opposite pole placement in SV2. Both ACT and REP show mild same-pole
tendencies (both positive z), with activation stronger. Cross-pole rates are uniformly suppressed
(negative z), suggesting both regulatory types are geometrically co-localized rather than opposed.

**Decision**: **Negative** — retire repression anti-pole hypothesis. Regulatory sign is not
encoded as SV2 pole opposition.

**Artifact**: `h02_repression_antipole.json`

## H03: Hub-Degree Control for STRING PPI Co-pole

**Method**: Compute STRING PPI degree for each named gene. Split STRING edges into:
- Hub edges: at least one gene has degree > median (7.0); N=987
- Non-hub edges: both genes degree ≤ median; N=35

Run SV2 co-pole test separately on hub vs non-hub edges at each layer.

**Results** (mean z across 12 layers):
| Edge subset | Mean z | Sig layers (z>1.96) | N pairs |
|-------------|--------|----------------------|---------|
| ALL         | 9.89   | 12/12                | 1022    |
| HUB         | 10.01  | 12/12                | 987     |
| NON-HUB     | 1.43   | 3/12                 | 35      |

**Interpretation**: The PPI co-pole signal is heavily concentrated in hub edges. Non-hub edges
show only marginal signal (mean z=1.43, only 3/12 layers borderline). This is a confound concern:
the high-degree hubs (highly connected proteins) drive the enrichment. However, 987/1022 edges
are hub edges, leaving only 35 non-hub edges — insufficient statistical power to distinguish
true signal from noise. The test is underpowered for non-hub edges.

**Conclusion**: Signal is not definitively ruled out for non-hub edges (low N). But hub bias
is a legitimate confound. The STRING co-pole result (iter_0012 H03) should be qualified:
"PPI co-pole enrichment holds for high-degree genes; low-degree gene pairs insufficient to test."

**Decision**: **Inconclusive** (underpowered for non-hub) — confound concern identified but
cannot be resolved with 35 non-hub pairs. Lower STRING threshold would expand non-hub set.

**Artifact**: `h03_hub_degree_control.json`

## Quantitative Summary

| Hypothesis | Decision | Key Metric |
|-----------|----------|-----------|
| H01: SV axis specificity | Mixed/positive | SV2 z=10.04, SV3 z=7.18, SV1 z=0.39 (12/12 layers each for SV2/SV3) |
| H02: Repression anti-pole | Negative | rep_xpole mean z=-1.41, 0/12 layers significant |
| H03: Hub degree control | Inconclusive | non-hub mean z=1.43 (N=35 pairs only), hub mean z=10.01 |

## Artifacts Generated

- `iterations/iter_0013/h01_sv_axis_specificity.json`
- `iterations/iter_0013/h02_repression_antipole.json`
- `iterations/iter_0013/h03_hub_degree_control.json`
- `iterations/iter_0013/iter0013_results.json`
- `iterations/iter_0013/run_iter0013_screen.py`
