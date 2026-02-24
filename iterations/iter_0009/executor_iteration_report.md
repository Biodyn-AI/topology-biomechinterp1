# iter_0009 Executor Iteration Report

**Date:** 2026-02-22
**Iteration:** iter_0009
**Status:** PASSED — 3/3 hypotheses tested, all positive

---

## Summary

Three hypotheses tested:

| ID | Name | Decision | Direction |
|----|------|----------|-----------|
| H01 | SV2 bottom-pole extracellular vesicle null (N=1000) | **promising** | positive |
| H02 | 12-layer × 8-compartment systematic scan | **promising** | positive |
| H03 | Cross-layer SV1 rank stability (Spearman) | **promising** | positive |

---

## Command Trace

```bash
# Environment: subproject40-topology
conda run -n subproject40-topology python run_iter0009_screen.py 2>&1 | tee run_stdout.log
```

Script: `iterations/iter_0009/run_iter0009_screen.py`

---

## H01: SV2 Bottom-Pole Extracellular Vesicle Null

**Motivation:** iter_0008 found SV2 bottom-pole enrichment for GO:0070062 (extracellular vesicle, OR=7.99, p=8.7e-9) but did not run a gene-label shuffle null. This was the most critical missing control.

**Method:** Gene-label shuffle null, N=1000 replicates. Top-K=52 genes at SV2 bottom (lowest projection). Fisher exact test vs GO:0070062 (48 annotated genes in universe of 209).

**Results:**
- Observed: OR=7.99, p=8.69e-09, a=28 genes in pole from 48 annotated
- Empirical p (gene-label null, N=1000): **emp_p=0.000** (0/1000 shuffles achieved as small a p-value)
- EV genes in bottom pole (28): AKR1B1, ALDH1A1, APP, ASAH1, ASS1, ATP1B1, CA2, CLU, CTSB, EPCAM, FXYD2, GSTP1, HLA-A, HLA-B, HLA-C, HLA-DRA, HPGD, ITGB8, LDHA, LGALS1, MIF, MMP7, PCBD1, RHOB, RHOC, TNFSF10, VIM, WFDC2

**Interpretation:** SV2 bottom-pole extracellular vesicle enrichment is highly significant and null-controlled. The gene list is biologically coherent: MIF, APP, VIM, LGALS1, CLU are canonical exosome/EV cargo; HLA class I/II are EV surface markers. This confirms SV2 encodes a genuine secretion/EV axis orthogonal to SV1 (secretory pathway).

**Artifacts:** `h01_sv2_bot_null_ps.npy`, `iter0009_results.json`

---

## H02: 12-Layer × 8-Compartment Systematic Scan

**Motivation:** Prior iterations found isolated compartment-layer associations. This hypothesis runs a comprehensive scan to produce a "layer compartment map."

**Method:** For each of 12 layers and 8 GO compartments (mitochondrion GO:0005739, ER_lumen GO:0005788, extracellular_vesicle GO:0070062, cytoskeleton GO:0005856, plasma_membrane GO:0005886, nucleus GO:0005634, secreted GO:0005615, ribosome GO:0005840), compute best-pole Fisher OR on SV1 top/bottom-52 with N=200 gene-label shuffles each.

**Results (96 cells tested; 33 significant at emp_p ≤ 0.05):**

Key significant cells (emp_p=0.000):
| Layer | Compartment | OR | p |
|-------|-------------|-----|---|
| 0 | ER_lumen | 6.95 | 2.0e-3 |
| 0 | secreted | 3.01 | 1.9e-3 |
| 1 | ER_lumen | 4.73 | 1.1e-2 |
| 1 | secreted | 3.87 | 1.5e-4 |
| 2–11 | ER_lumen | 6.95–18.45 | consistently significant |
| 2–11 | secreted | 3.0–5.2 | consistently significant |
| 10 | mitochondrion | 13.77 | 1.8e-5 |
| 11 | ER_lumen | 18.45 | 1.9e-5 |

**Pattern:**
- **ER lumen** and **secreted** compartments: persistent across all 12 layers — the secretory axis is present from layer 0
- **Mitochondrion**: strongest at layers 2–4 and 10 (transient enrichment at L2–4, resurges at L10)
- **Extracellular vesicle, cytoskeleton, plasma membrane, nucleus**: sporadic or absent on SV1 (EV is on SV2 as confirmed by H01)

**Artifacts:** `h02_layer_compartment_scan.json`, `h02_layer_compartment_map.csv`

---

## H03: Cross-Layer SV1 Rank Stability

**Novelty:** New method — tests whether the secretory axis gene ranking is stable across all 12 transformer layers.

**Method:** For each of 12 layers, compute SVD of mean-centered 209-gene embedding, extract SV1 projection. Compute Spearman r between all adjacent layer pairs. Compare against null (row-permutation control, N=200).

**Results:**

Adjacent-layer Spearman r values:
```
L0-L1: 0.951
L1-L2: 0.945
L2-L3: 0.980
L3-L4: 0.977
L4-L5: 0.970
L5-L6: 0.966
L6-L7: 0.952
L7-L8: 0.954
L8-L9: 0.958
L9-L10: 0.899
L10-L11: 0.663  ← largest transition
```

- **Mean adjacent Spearman r = 0.929** (range: 0.663–0.980)
- **Null mean = 0.003**, empirical p = **0.000** (0/200 null reps exceed observed)
- **SV1 variance explained rises monotonically:** 19.1% (L0) → 76.7% (L11)

**Interpretation:** The secretory axis is remarkably stable across layers 0–9 (r > 0.89) with a single transition at L10-L11 (r=0.663) where ER lumen enrichment also strengthens (OR=18.45 at L11 vs ~7 at L0-L9). The monotonically increasing SV1 variance explained suggests the model progressively sharpens this axis. This provides strong evidence that SV1 is not a random artifact but a layer-coherent biological signal.

**Artifacts:** `h03_sv1_projections_12layers.npy`, `h03_sv1_spearman_mat.npy`, `h03_sv1_stability.csv`

---

## Evidence Quality

All three hypotheses:
- ✅ Explicit command trace
- ✅ Quantitative metrics with effect sizes
- ✅ Gene-label null controls
- ✅ Biologically interpretable
- ✅ Machine-readable artifacts generated

---

## Blocker Log

No blockers. Minor JSON serialization bug (numpy bool_) fixed with trivial post-hoc script.
