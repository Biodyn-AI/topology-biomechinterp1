# Brainstormer Structured Feedback — iter_0037

## Gate Status
`passed_min_research_gate: true` — all checks passed.

## Signal Assessment

### What Worked
- **B-cell z=7.55** at L2 is the strongest single result across 37 iterations. The 5-gene panel (MS4A1, CD19, CD79A, BLK, PRDM1) produces robust, reproducible clustering.
- **Sub-lineage geometric split** (germinal center TFs near centroid; plasma TFs far) is a novel biological finding. BATF/SPIB/BACH2 at 86–96th pctile vs IRF4 at 20th pctile directly maps the B→plasma differentiation axis in embedding space.
- **T-cell permutation null** (emp_p=0.96) definitively closes the T-cell direction and validates B-cell specificity.
- **Centroid distance matrix** gives interpretable geometry: B-cell is an isolated cluster; DC and plasma are both far, but for different reasons.

### Key Biological Interpretation
The scGPT embedding geometry captures the B-cell germinal-center identity program as a coherent geometric manifold. Late differentiation TFs (IRF4, IRF8) exit this manifold, consistent with their context-switching role across cell types. This is a mechanistically interpretable result.

### What Was Not Tested (Gaps)
1. **Cross-model validation** — no Geneformer or scFoundation comparison yet. This is the single largest gap for the paper.
2. **Sub-lineage centroid quantification** — the GC-TF cluster vs plasma-TF cluster centroid distance has not been explicitly measured.
3. **Layer-resolved differentiation axis** — does the GC↔plasma split grow from L0→L11 (suggesting progressive specialization encoding)?
4. **Broader vocabulary scan** — PAX5, EBF1, BCL6, IKZF1 may be in the full 4803-gene vocabulary but were not searched.
5. **Topological structure** — persistent homology of the B-cell cluster has not been tested; all results are Euclidean.
6. **Manifold intrinsic structure** — local dimensionality, curvature, geodesics within the B-cell cluster not examined.

## Stale / Exhausted Directions
- **intrinsic_dimensionality SVD AUROC**: confirmed negative ≥3 iterations — retire.
- **TRRUST/GO proximity on random gene pairs**: confirmed negative — retire.
- **STRING PPI Euclidean proximity** (iter_0031 H02): confirmed negative — retire.
- **graph_topology Laplacian spectral gap**: one test, no biological anchor, no follow-up — retire.
- **module_structure GO BP**: negative — retire.
- **T-cell clustering**: definitively negative (4 iterations now) — retire as primary target; keep only as null control.
- **DC clustering**: negative, likely due to heterogeneity — no rescue potential with current gene set — deprioritize.

## Strategic Pivot Needed
The manifold_distance family has produced strong positive results but is becoming repetitive (same B-cell panel, same L2 metric). The next phase must:
1. **Cross-validate** the finding (different model, same protocol).
2. **Deepen mechanistically** (differentiation axis, layer progression, topological structure).
3. **Broaden biologically** (other cell types with larger marker sets, immune sub-programs).
4. **Add topology** — the persistent homology angle has never been seriously tested on the confirmed B-cell cluster.
