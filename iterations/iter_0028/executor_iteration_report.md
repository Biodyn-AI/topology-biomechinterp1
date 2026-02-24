# Executor Iteration Report — iter_0028

## Overview
Three hypotheses tested targeting: (H01) hub gene centrality in embedding space, (H02) PC1 as T-cell vs APC biological axis, (H03) kNN graph spectral gap across layers. All tests run to completion.

## Command trace

```bash
cd "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0028"
conda run -n subproject40-topology python run_iter0028_screen.py
```

Artifacts generated:
- `h01_hub_centrality.json`
- `h02_pc1_celltype_axis.json`
- `h03_spectral_gap.json`

---

## H01: Hub Gene Embedding Centrality

**Hypothesis:** High STRING-degree hub genes cluster nearer to the centroid of scGPT embedding space at deeper layers.

**Method:** For 201 named genes with ≥1 STRING edge (score≥0.4), compute mean L2 distance to the layer centroid. Test Spearman ρ(STRING_degree, −dist_to_centroid) per layer. Null: 500 permutations of degree labels.

**Results:**

| Layer | rho | p-value | null_p |
|-------|-----|---------|--------|
| L0  | +0.009 | 0.896 | 0.898 |
| L7  | +0.134 | 0.058 | 0.044 |
| L8  | +0.177 | 0.012 | 0.010 |
| L9  | +0.199 | 0.005 | 0.008 |
| L10 | **+0.225** | **0.00135** | **0.000** |
| L11 | +0.178 | 0.012 | 0.014 |

- Monotonic increase from L0 (rho=+0.009) through L10 (rho=+0.225).
- At L10: null_p=0.000 (0 out of 500 permutations exceed |rho|=0.225).
- Top most-central genes at L11: ETS1 (degree=53), IRF8 (61), TGFBR2 (36), REL (43), ITGB2 (57).

**Interpretation:** Hub genes (high network connectivity) progressively migrate toward the embedding centroid in deeper scGPT layers. This converts the iter_0027 null-construction confound into a characterizable geometric phenomenon. The signal is real, growing with depth, and survives permutation testing.

**Decision:** PROMISING — positive, survives null.

---

## H02: PC1 as T-cell vs APC Axis

**Hypothesis:** PC1 at L11 encodes a T-cell effector (JUN/LCK/PRF1/CD8A) vs APC (HLA-A/B/C/DRA/DRB1/FOS) axis.

**Method:** Curated T-cell markers in vocab: PRF1, LCK, IFNG, CD8A, JUN, JUNB (n=6). APC markers: HLA-A, HLA-B, HLA-C, HLA-DRA, HLA-DRB1, FOS, FOSB (n=7). Fit PCA per layer, compute AUROC and Spearman ρ(PC1, ±1 label) for axis genes.

**Results:**

| Layer | rho | p | AUROC | PC1_var |
|-------|-----|---|-------|---------|
| L0  | +0.412 | 0.161 | 0.738 | 19.1% |
| L1  | +0.454 | 0.120 | **0.762** | 24.6% |
| L7  | +0.000 | 1.000 | 0.500 | 41.6% |
| L11 | −0.371 | 0.212 | 0.714 | **76.7%** |

- L11 detail: T-cell markers (PRF1=+2.34, LCK=+2.21, JUN=+2.20) and APC markers (HLA-A=+2.26, HLA-B=+2.27, FOS=+2.33, FOSB=+2.39) all cluster at very similar PC1 values (+2.1 to +2.4). NO separation.
- The brainstormer's JUN-top vs FOS-bottom pattern does NOT hold at L11; both families are crowded at the same PC1 position.
- L1 AUROC=0.762 but n=13 axis genes and p=0.12 (not significant). Early layers show weak axis structure that collapses by L7.

**Decision:** NEGATIVE — PC1 does not encode a T-cell/APC biological axis at any layer. The L11 "axis" is a geometric artifact where both gene families collapse to the same PC1 extremum due to dimensional collapse (PC1_var=76.7%).

---

## H03: kNN Graph Spectral Gap (new family: graph_topology spectral)

**Hypothesis:** Spectral gap (lambda_2 of normalized graph Laplacian) of the kNN gene graph changes monotonically across scGPT layers, reflecting changing community structure.

**Method:** k=10 nearest neighbor graph on 209 named gene embeddings per layer. Normalized Laplacian eigenvalue spectrum. Spectral gap = smallest non-zero eigenvalue (Fiedler value). 100-permutation null at L11.

**Results:**

| Layer | Spectral Gap | Entropy | n_components |
|-------|-------------|---------|--------------|
| L0  | 0.1614 | 2.8238 | 2 |
| L3  | 0.1333 | 2.8126 | 2 |
| L6  | 0.0940 | 2.8008 | 2 |
| L9  | 0.0600 | 2.7936 | 2 |
| L11 | **0.0361** | **2.7911** | 2 |

- Spectral gap trend: rho=−0.993, p=1.3×10⁻¹⁰ (near-perfect monotonic decrease)
- Entropy trend: rho=−0.986, p=4.1×10⁻⁹
- A decreasing spectral gap = the kNN graph becomes more modular / community-structured with depth
- Null: permuted gene positions give identical spectral gap (null_p=0.190) — this is expected since spectral gap depends only on geometry, not gene labels

**Interpretation:** The kNN graph becomes progressively more modular with depth. This aligns with the PR (participation ratio) collapse seen in iter_0026 H03: dimensional collapse at deeper layers means genes cluster into tighter geometric communities with weaker inter-community connections. The spectral gap trend is a genuine geometric descriptor of increasing modularity.

The null test (permuted gene labels → same gap) confirms this is a property of the embedding geometry, not gene identity. A proper null (random Gaussian embeddings) would be needed to assess whether this exceeds baseline.

**Decision:** PROMISING — monotonic trend is striking (rho=−0.993), mechanistically coherent with PR collapse, but null testing is incomplete. Need random-embedding null to confirm.

---

## Summary table

| ID | Name | Family | Result | Decision |
|----|------|--------|--------|----------|
| H01 | Hub Gene Centrality | manifold_distance | rho=+0.225 at L10, null_p=0.000 | promising |
| H02 | PC1 T-cell vs APC Axis | intrinsic_dimensionality | AUROC=0.714 at L11 but no axis separation | negative |
| H03 | kNN Spectral Gap | graph_topology | gap rho=−0.993, p=1.3×10⁻¹⁰ | promising |
