# Executor Iteration Report — iter_0027

## Summary

Three hypotheses tested: PC1 biological identity at L11 (H01), TRRUST activation/repression
geometric polarity (H02), and Dorothea confidence threshold sensitivity sweep (H03).
H01 and H02 are negative; H03 reveals an important methodological finding about null construction
sensitivity that challenges the robustness of earlier manifold_distance results.

---

## Command Trace

```bash
conda run -n subproject40-topology python3 iterations/iter_0027/run_iter0027_screen.py
# → iterations/iter_0027/h01_pc1_identity.json
# → iterations/iter_0027/h02_trrust_polarity.json
# → iterations/iter_0027/h03_threshold_sweep.json
```

All runs used:
- `conda run -n subproject40-topology python3`
- embeddings: `.../subproject_38.../outputs/cycle1_main/layer_gene_embeddings.npy` — [12, 4803, 512]
- TRRUST: `.../single_cell_mechinterp/external/networks/trrust_human.tsv`
- Dorothea: `.../single_cell_mechinterp/external/networks/dorothea_human.tsv`
- 209 named genes (all vocab genes matching HGNC pattern)

---

## H01: PC1 Biological Identity at L11 (intrinsic_dimensionality, refinement)

**Method:** At L11 (the most-collapsed layer, PR=1.68, PC1 EV=67.6%), compute PCA and extract
PC1 loadings for 209 named genes. Test whether PC1 separates known biological categories:
TF vs non-TF (TRRUST+Dorothea), HLA family, AP1 family, hub degree (TRRUST target count).
Mann-Whitney U AUROC and Spearman correlation vs hub degree. Also ran per-layer to check progression.

**Results:**

| Feature | N (positive) | AUROC | p-value | Interpretation |
|---------|-------------|-------|---------|----------------|
| TF vs non-TF | 69 TF / 140 non-TF | 0.492 | 0.86 | No separation |
| HLA vs non-HLA | 6 HLA genes | 0.433 | 0.58 | No separation |
| AP1 vs non-AP1 | 4 AP1 genes | 0.381 | 0.42 | No separation |
| Hub degree (n targets) | n=69 TFs | rho=−0.051 | 0.68 | No correlation |

PC1 EV ratio at L11: **67.6%** (i.e., the first component captures 2/3 of all variance).

**Top PC1 genes (high end):** PBK, COL1A1, TLR10, ITGB2, NEAT1, ALDH1A1, LCK, BIRC3, POU1F1, JUN
**Bottom PC1 genes (low end):** PAX5, PTGS2, NR4A3, TNF, HLA-DPB1, TBXAS1, KLF6, HLA-A, NCOA3, FOS

**Notable observation:** JUN (AP1 activator, top PC1) vs FOS (AP1 activator, bottom PC1) — two AP1
heterodimer partners land on opposite ends of the dominant axis. Similarly: LCK (T-cell kinase, top)
vs HLA-A / HLA-DPB1 (antigen presentation, bottom). This suggests the collapsed axis may encode a
T-cell effector vs antigen-presentation axis, but this is descriptive rather than statistically established.

TF AUROC per layer on PC1: `[0.459, 0.454, 0.460, 0.533, 0.526, 0.507, 0.502, 0.483, 0.480, 0.481, 0.468, 0.492]` — no significant trend.

**Decision: NEGATIVE** — PC1 at L11 does not preferentially separate TF-hood, HLA family, AP1 family,
or hub degree. The dominant axis of the nearly-1D collapse is not explained by these simple binary
features.

---

## H02: TRRUST Activation/Repression Polarity (module_structure, new_method)

**Method:** Loaded all TRRUST named-gene pairs (116 activation, 64 repression). Tested whether
activation vs repression targets are separated on PC1 at L11. Also tested per-layer and via L2 distance.
Within-TF: for 10 TFs with >=2 activation and >=2 repression named targets, computed delta = mean_act_PC1 − mean_rep_PC1.

**Results:**

Overall polarity on PC1 at L11:
- AUROC = **0.514**, p = 0.794 → **not significant**
- mean_act_PC1 = −0.316, mean_rep_PC1 = −0.631 (small difference, not significant)

Per-layer polarity AUROC:
`[0.523, 0.527, 0.512, 0.498, 0.492, 0.499, 0.509, 0.511, 0.509, 0.501, 0.505, 0.514]`
→ All ~0.50, no layer shows significant directional encoding.

L2 distance test (L11): mean_act = 7.64 vs mean_rep = 9.04, p = 0.178 → not significant.

Within-TF deltas (selected):
| TF | n_act | n_rep | delta (act−rep PC1) |
|----|-------|-------|---------------------|
| ETS1 | 10 | 2 | +10.85 |
| CIITA | 6 | 3 | −3.68 |
| FOS | 2 | 4 | −5.37 |
| JUN | 17 | 4 | −3.79 |
| HIF1A | 9 | 3 | −2.27 |
| ATM | 3 | 2 | +0.03 |
| EGR1 | 7 | 6 | +0.60 |

ETS1 shows a large positive delta (its activation targets are much higher on PC1 than its repression targets), but with only 2 repression targets this is unreliable. No consistent directionality across TFs.

**Decision: NEGATIVE** — TRRUST activation/repression directionality is not encoded in the PC1 space at any layer. Individual TFs show variable within-TF deltas but no consistent direction.

---

## H03: Dorothea Confidence Threshold Sensitivity Sweep (null_sensitivity, new_method)

**Method:** Used Dorothea confidence tiers as a score proxy (A=1.0, B=0.75, C=0.5, D=0.25). Swept
5 thresholds (D+, C+, B+, B+ upper, A-only) and for each computed AUROC of regulatory pairs being
closer than 500 random null pairs drawn from all 209 named genes. Compared with PR collapse curve from iter_0026.

**Results:**

| Threshold | Label | n_pairs | Best AUROC | Best Layer | PR at best layer |
|-----------|-------|---------|------------|------------|-----------------|
| 0.25 | D+ (all) | 1183 | 0.446 | L3 | 10.1 |
| 0.50 | C+ | 372 | 0.451 | L8 | 4.6 |
| 0.625 | B+/midpoint | 270 | 0.444 | L8 | 4.6 |
| 0.75 | B+ | 270 | 0.444 | L8 | 4.6 |
| 1.00 | A-only | 210 | 0.449 | L8 | 4.6 |

**Critical finding:** All thresholds yield AUROC < 0.5, meaning Dorothea regulatory pairs are
FARTHER from each other than random null pairs when using an all-gene null.

This directly contradicts iter_0026 H02 which reported AUROC = 0.658 for A/B pairs vs null.
The discrepancy is explained by the null construction: iter_0026 used a null drawn from the same
TF-target gene pool as the positive pairs. The current test uses null pairs from ALL 209 named genes.

If hub TFs (JUN, EGR1, HIF1A — each with 10-20 named targets) are central in embedding space,
then the TF-target pairs in Dorothea will generally be among the more central genes. A same-pool
null controls for gene centrality; an all-gene null does not. This implies the iter_0026 AUROC=0.658
was capturing "hub genes are central" rather than "regulatory pairs are specifically proximate."

Spearman(threshold, best_layer): rho=0.71, p=0.18 (marginally higher-conf pairs peak later, non-significant).

**Decision: NEGATIVE with IMPORTANT METHODOLOGICAL NOTE** — Dorothea threshold sensitivity reveals
that the earlier positive regulatory proximity signal was null-sensitive. With proper population-level
null, no threshold shows above-chance proximity. This challenges the specificity of manifold_distance
findings across iter_0022–0026.

---

## Overall Assessment

| H | Family | Decision | Key Metric | Novelty |
|---|--------|----------|------------|---------|
| H01 | intrinsic_dimensionality | **negative** | All features AUROC ~0.43–0.49 | refinement |
| H02 | module_structure | **negative** | Polarity AUROC=0.514, p=0.79 | new_method |
| H03 | null_sensitivity | **negative** (informative) | Best AUROC=0.451 < 0.5 | new_method |

All three negative. H03 is the most informative: it identifies a methodological fragility in prior
manifold_distance results. The PC1 identity question (H01) remains open — the dominant axis is not
captured by simple binary features, and the JUN-top / FOS-bottom / LCK-top / HLA-bottom pattern
deserves follow-up with a proper continuous biological annotation.

---

## Artifacts Generated

- `iterations/iter_0027/h01_pc1_identity.json` (machine-readable per-layer TF AUROC + gene loadings)
- `iterations/iter_0027/h02_trrust_polarity.json` (per-layer polarity AUROC + within-TF deltas)
- `iterations/iter_0027/h03_threshold_sweep.json` (per-threshold AUROC + PR alignment)
- `iterations/iter_0027/run_iter0027_screen.py` (full reproducible script)
