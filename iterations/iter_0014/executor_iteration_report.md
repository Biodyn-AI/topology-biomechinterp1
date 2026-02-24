# Executor Iteration Report — iter_0014

## Summary

Three hypotheses tested, all producing clean quantitative artifacts:
- **H01 (GO enrichment of SV3 poles)**: Positive across all 12 layers. SV3 encodes a kinase/signaling vs immune axis that shifts biological content with depth.
- **H02 (Layer-depth Spearman on SV2/SV3 z-scores)**: Mixed. SV2 shows a moderate declining trend (rho=-0.56, p=0.059). SV3 is non-monotonic (peak at L3). Neither definitively monotonic.
- **H03 (Hub-degree control, STRING score>=0.4)**: Strong positive. Non-hub edges (N=138) show SV2 co-pole enrichment in 12/12 layers (mean z=3.11), exceeding hub signal (mean z=2.60, 9/12 sig). Hub-degree confound definitively ruled out.

---

## H01: GO Enrichment of SV3 Poles

### Motivation
iter_0013 established SV3 is a significant PPI geometry axis (mean z=7.18, 12/12 layers). The biological identity of SV3 was unknown.

### Command Trace
```bash
conda run -n subproject40-topology python iterations/iter_0014/run_iter0014_screen.py
```

### Method
- SVD of 209-gene mean-centered embedding at each of 12 layers
- Top-K=52 and bottom-K=52 genes by SV3 projection
- Fisher exact test against 594 GO BP+CC terms (size 3-60)
- Gene-label shuffle null N=300 per best term per layer

### Results
All 12 layers significant for both poles (emp_p <= 0.017):

| Layer | TOP pole best term | BOT pole best term |
|-------|-------------------|-------------------|
| L0 | Integrated stress response (emp_p=0.000) | Neuron apoptotic process (emp_p=0.000) |
| L1 | Peptidyl-serine phosphorylation (emp_p=0.000) | Integrated stress response (emp_p=0.000) |
| L2 | Peptidyl-serine phosphorylation (emp_p=0.003) | Integrated stress response (emp_p=0.000) |
| L3 | Protein Ser/Thr kinase activity (emp_p=0.000) | Integrated stress response (emp_p=0.000) |
| L4 | Protein Ser/Thr kinase activity (emp_p=0.003) | Inflammatory response (emp_p=0.000) |
| L5 | Inflammatory response (emp_p=0.003) | Protein Ser/Thr kinase activity (emp_p=0.007) |
| L6 | Cell surface (emp_p=0.000) | Protein tyrosine kinase activity (emp_p=0.000) |
| L7 | Cellular response to calcium ion (emp_p=0.003) | MHC class II antigen assembly (emp_p=0.000) |
| L8 | MHC class II antigen assembly (emp_p=0.000) | Cell surface (emp_p=0.017) |
| L9 | Chromatin binding (emp_p=0.007) | MHC class II antigen assembly (emp_p=0.000) |
| L10 | DNA binding (emp_p=0.000) | MHC class II antigen assembly (emp_p=0.000) |
| L11 | DNA binding (emp_p=0.000) | Negative regulation of T cell proliferation (emp_p=0.003) |

**Interpretation**: SV3 encodes a depth-progressive axis. Early layers (L0-5): kinase/signaling (top) vs stress/inflammatory (bottom). Mid-layers (L6-8): cell surface vs immune signaling. Late layers (L9-11): DNA binding/chromatin vs MHC-II/immune regulation. The axis transitions from kinase biology to nuclear/immune regulatory biology with increasing depth.

### Artifact
`iterations/iter_0014/h01_sv3_go_enrichment.json`

---

## H02: Layer-depth Spearman trend on SV2/SV3 z-scores

### Motivation
Characterize whether PPI geometry strengthens or weakens monotonically across transformer depth. Uses already-computed data from iter_0013.

### Command Trace
```bash
conda run -n subproject40-topology python iterations/iter_0014/run_iter0014_screen.py
# H02 uses iter_0013/h01_sv_axis_specificity.json — no new embedding computation needed
```

### Method
- Load SV1/SV2/SV3 z-scores per layer from iter_0013 artifact
- Spearman correlation between layer index [0-11] and z-score
- Identify peak layer per axis

### Results

| Axis | Spearman rho | p-value | Peak layer |
|------|-------------|---------|-----------|
| SV1 | +0.168 | 0.602 | L11 |
| SV2 | -0.559 | 0.059 | L1 |
| SV3 | -0.259 | 0.417 | L3 |

**SV2**: Moderate declining trend (rho=-0.56, p=0.059, marginal). Peaks at L1, weakens toward output. Not definitively monotonic.
**SV3**: Non-monotonic — bell shape, peak at L3, falls off in both directions. Not monotonic.
**SV1**: No trend.

### Artifact
`iterations/iter_0014/h02_layer_depth_trend.json`

---

## H03: Expanded Hub-Degree Control (STRING score>=0.4)

### Motivation
iter_0013 H03 found non-hub control underpowered (N=35 pairs at score>=0.7). Expanding to score>=0.4 yields 138 non-hub edges.

### Command Trace
```bash
conda run -n subproject40-topology python /tmp/h03_fix.py
# Loads iterations/iter_0012/string_ppi_named_genes.json (all STRING scores 0.15-0.999)
# Filters to score>=0.4 -> 3364 edges
# Hub: degree > median_degree_04=24 -> N=2347; Non-hub: degree <= 24 -> N=138
```

### Method
- STRING JSON from iter_0012 contains scores 0.15-0.999 for all 209-gene named pairs
- Filter to score>=0.4: 3364 edges (vs 1022 at score>=0.7)
- Degree computed in score>=0.4 graph; median=24
- Hub/non-hub split; run SV2 co-pole Fisher test + gene-label shuffle null N=300 per layer

### Results

| Layer | Hub z (N=2347) | Non-hub z (N=138) |
|-------|--------------|-----------------|
| L0 | 4.21* | 3.36* |
| L1 | 4.54* | 2.30* |
| L2 | 3.02* | 3.36* |
| L3 | 2.81* | 3.50* |
| L4 | 2.31* | 3.02* |
| L5 | 2.92* | 2.48* |
| L6 | 3.15* | 3.12* |
| L7 | 1.84* | 3.26* |
| L8 | 1.67 | 1.95* |
| L9 | 0.91 | 3.26* |
| L10 | 0.98 | 3.24* |
| L11 | 2.78* | 4.44* |
| **Mean** | **2.595 (9/12 sig)** | **3.107 (12/12 sig)** |

**Key finding**: Non-hub edges (low-degree gene pairs) show stronger and more consistent SV2 co-pole enrichment than hub edges. This definitively rules out hub-degree confound. PPI geometry in SV2 reflects the general network structure, not just high-connectivity hubs.

### Artifact
`iterations/iter_0014/h03_hub_degree_control_04.json`

---

## Overall Assessment

| Hypothesis | Direction | Decision |
|-----------|-----------|---------|
| H01: GO enrichment SV3 poles | positive | promising |
| H02: Layer-depth Spearman | mixed | neutral |
| H03: Hub-degree control (score>=0.4) | positive | promising |

**Portfolio**:
- 2 promising (H01, H03), 1 neutral (H02)
- H01 is a new biological characterization of SV3 — novel finding
- H03 resolves the critical hub confound from iter_0013 with adequate power

**Novelty check**:
- H01: new_family (GO enrichment on SV3 is new)
- H02: new_method (Spearman on existing z-scores; zero new data)
- H03: new_method (score threshold change gives adequate power)
- No retired directions explored

---

## Compile information

```bash
latexmk -pdf -interaction=nonstopmode paper/autoloop_research_paper.tex
```
