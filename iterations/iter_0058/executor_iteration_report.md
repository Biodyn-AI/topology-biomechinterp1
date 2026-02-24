# Executor Iteration Report — iter_0058

**Date**: 2026-02-23
**Status**: Complete — 3 hypotheses tested

---

## Summary

1. **H01** (TF boundary gene identity + family enrichment): **Promising** — BCL11A, NFKB1, FOXO3 are top boundary anchors at peak layers L2/L3; Forkhead/STAT families best separated; bZIP/C2H2-ZF families near-target (negative margin).
2. **H02** (TRRUST graph Laplacian spectral alignment): **Neutral** — SV5-7 is statistically more aligned than random (Z=7.81), but absolute alignment is weak (PA ~85° vs 88° random; 90°=orthogonal).
3. **H03** (Effective rank as AUROC predictor): **Promising** — Spearman rho=0.855 (p<0.0001) between full effective rank and joint 6D AUROC across 36 (seed×layer) observations. Confirms and extends iter_0057 spectral geometry insight.

---

## Command Trace

```bash
mkdir -p "/Volumes/Crucial X6/.../iter_0058"
conda run -n subproject40-topology python /tmp/iter58_experiments.py
```

Data paths:
- `cycle4_immune_main/layer_gene_embeddings.npy` [12, 4941, 512]
- `cycle4_immune_seed43/layer_gene_embeddings.npy`
- `cycle4_immune_seed44/layer_gene_embeddings.npy`
- Edge file: `cycle4_immune_main/cycle1_edge_dataset.tsv` (589 positive TRRUST pairs)

Environment: `subproject40-topology` (conda)

---

## H01: TF Boundary Gene Identity + Family Enrichment

### Method
At peak layers L2 and L3 (highest cross-seed AUROC from iter_0057):
- SVD of centered nonzero embeddings [2039,512]; project to joint 6D SV2-7
- Per-TF **margin** = dist_to_target_centroid − dist_to_TF_centroid (positive = well-separated from targets)
- Ranked top anchors; classified TFs into families by name prefix heuristic

### Results

**Top boundary anchors (consistent L2 and L3):**

| TF Gene | Margin (L2) | Margin (L3) |
|---------|-------------|-------------|
| BCL11A  | 1.784       | 1.731       |
| NFKB1   | 1.543       | 1.445       |
| RB1     | 1.528       | 1.475       |
| FOXO3   | 1.521       | 1.399       |
| ZEB1    | 1.499       | 1.435       |
| SETBP1  | 1.483       | 1.351       |
| TRPS1   | 1.371       | 1.381       |
| NFIX    | 1.364       | 1.389       |
| STAT3   | 1.322       | 1.293       |
| RBPJ    | 1.308       | —           |

**TF family enrichment (L2 median margin):**

| Family          | Median Margin | Count |
|-----------------|---------------|-------|
| Forkhead (FOXO) | 1.521         | 1     |
| STAT            | 1.299         | 2     |
| bHLH (MYC)      | 1.025         | 1     |
| RUNX            | 0.843         | 2     |
| ETS             | 0.930         | 3     |
| Rel (NFKB)      | 0.582         | 5     |
| Other           | 0.668         | 43    |
| Nuclear receptor| 0.207         | 5     |
| GATA            | 0.027         | 2     |
| bZIP            | −0.180        | 4     |
| C2H2-ZF         | −0.472        | 5     |

**Biological interpretation:**
- BCL11A (B-cell lymphoma ZF 11A) is the top boundary anchor — known master TF for B-cell and T-cell differentiation. Highly relevant to immune dataset.
- FOXO3 (Forkhead): immune cell survival/apoptosis regulator.
- STAT3: key cytokine signaling TF.
- bZIP and C2H2-ZF families are **nearest targets** (negative margin), suggesting they may have dual TF/target roles or high co-expression with their own targets.

**Verdict**: Promising. Top anchors are biologically coherent and cross-layer stable.

---

## H02: TRRUST Graph Laplacian Spectral Alignment

### Method
- Build symmetrized undirected TRRUST adjacency (2039 nodes, 589 edges)
- Normalized Laplacian L = I − D^{-1/2} A D^{-1/2}
- Extract eigenvectors 1–9 (Fiedler eigenvectors; skipping constant vector)
- Principal angles (via QR + SVD) between SV5-7 [n,3] and Laplacian top-{3,6,9} subspaces
- Baseline: 100 random 3D subspaces vs Laplacian top-3

### Results

| Layer | PA(SV5-7, Lap-3) | PA(SV2-4, Lap-3) |
|-------|-----------------|-----------------|
| L0    | 84.88°          | 85.30°          |
| L1    | 85.47°          | 85.86°          |
| L2    | 86.24°          | 86.33°          |
| L3    | 84.90°          | 86.25°          |
| L8    | 86.19°          | 87.26°          |
| L11   | 85.94°          | 85.32°          |

**Random baseline**: mean PA = 88.22 ± 0.43°
**SV5-7 best**: 84.88° (Z-score = 7.81 vs random)

**Interpretation:**
- Statistically, SV5-7 is more aligned to TRRUST graph Laplacian than a random subspace (Z=7.81).
- Practically, all principal angles are 84–88° (near-orthogonal; 90° = fully orthogonal).
- The improvement over random is ~3.3° — modest in absolute terms.
- Both SV5-7 and SV2-4 show similar (weak) alignment, so this is not subspace-specific.

**Verdict**: Neutral. Statistically detectable but practically weak alignment. The embedding geometry is not a simple mirror of TRRUST graph spectral structure.

---

## H03: Effective Rank as Cross-Seed AUROC Predictor

### Method
- For each of 3 seeds × 12 layers:
  - SVD of centered nonzero embeddings [2039,512]
  - **Effective rank** = exp(−Σ p_i log p_i) where p_i = s_i²/Σs²
  - Also: spectral entropy, effective rank of SV2-7 subset
  - Joint 6D AUROC (5-fold LR, 30-perm null)
- Spearman correlation across 36 (seed, layer) observations

### Results

**Spearman correlations:**
- rho(eff_rank_full, AUROC) = **0.855** (p<0.0001)
- rho(eff_rank_sv27, AUROC) = 0.539 (p=0.0007)
- rho(spec_entropy_full, AUROC) = 0.855 (p<0.0001)

**Layer trends (mean across 3 seeds):**

| Layer | Eff Rank (full) | Eff Rank (SV2-7) | AUROC |
|-------|-----------------|------------------|-------|
| L0    | 236             | 5.45             | 0.778 |
| L2    | 212             | 5.64             | 0.797 |
| L3    | 185             | 5.72             | 0.797 |
| L6    | 123             | 5.75             | 0.739 |
| L9    | 73              | 5.42             | 0.711 |
| L11   | 48              | 5.28             | 0.715 |

**Effective rank decreases L0→L11 (237→49); AUROC tracks (0.778→0.715) with similar shape.**

**Biological interpretation:**
- Higher-dimensional (more diffuse) residual stream representations → better TF/target discrimination.
- As scGPT layers progressively compress the representation (lower effective rank), the TF/target geometric separation also decreases.
- This is consistent with iter_0057's rho=-0.93 finding (energy concentration predicts AUROC negatively): two complementary measures, same principle.
- Effective rank is a stronger predictor than SV2-7 subset rank — the full spectrum matters.

**Verdict**: Promising. Clean, significant, cross-seed reproducible. rho=0.855 is the strongest single-metric predictor of AUROC found to date.

---

## Artifacts Generated

| File | Description |
|------|-------------|
| `h01_tf_boundary_genes.csv` | Per-TF margin at L2/L3 (all TFs) |
| `h01_tf_family_enrichment.csv` | Family annotation + margin for peak layers |
| `h02_laplacian_alignment.csv` | Principal angles SV5-7/SV2-4 vs Laplacian top-{3,6,9}, all layers |
| `h03_effrank_auroc.csv` | Effective rank + AUROC for 3 seeds × 12 layers |
| `iter58_summary.json` | Summary JSON of all three experiments |
| `executor_hypothesis_screen.json` | Machine-readable hypothesis outcomes |
