# E5 — Multiple-testing audit across the 183 paper hypotheses

## Global counts

- Total hypotheses (iter 1..63 inclusive): **183** (paper claim: 183 — match: yes)
- Hypotheses with extractable primary $p$: **112** (61.2%)

## Direction breakdown (paper's own classification)
```
direction
positive        117
negative         37
mixed            20
inconclusive      9
```

## Multiple-testing thresholds
- Bonferroni (paper denom $N=183$): $\alpha = 0.05/183 = 2.732e-04$
- Bonferroni (extractable denom $N=112$): $\alpha = 0.05/112 = 4.464e-04$
- BH-FDR at $q=0.05$ across the extractable set.

## Global pass counts
- BH-FDR $q\le0.05$: **70 / 112**
- Bonferroni (paper denom): **43 / 112**
- Bonferroni (strict denom): **45 / 112**

## Per-family

| family | n_hypotheses | n_extractable_p | n_bh_pass_within | n_bonferroni_pass_within | n_bh_pass_global | n_bonferroni_pass_paperdenom |
|---|---|---|---|---|---|---|
| module_structure | 33 | 33 | 21 | 13 | 21 | 12 |
| manifold_distance | 29 | 29 | 19 | 15 | 19 | 10 |
| intrinsic_dimensionality | 24 | 24 | 16 | 14 | 16 | 10 |
| graph_topology | 9 | 9 | 3 | 3 | 4 | 3 |
| null_sensitivity | 7 | 7 | 5 | 5 | 5 | 3 |
| topology_stability | 6 | 6 | 2 | 2 | 2 | 2 |
| persistent_homology | 4 | 4 | 3 | 3 | 3 | 3 |

## Headline-finding audit

| headline | claim | iter | hypothesis_id | primary_p | bh_q | bonferroni_pass | found |
|---|---|---|---|---|---|---|---|
| SV1_extracellular_enrichment | SV1 separates extracellular vs cytosolic; extracellular OR=6.37, p=2.6e-4 | iter_0006 | H02 | 0.0003 | 0.000764 | False | True |
| SV1_secretory_pathway_layered | Mitochondria/ER/extracellular sequence across layers (paper §2.2) | iter_0008 | H02 | 0 | 0 | True | True |
| SV2_PPI_copole_string07 | STRING≥0.7 PPI co-pole on SV2 at 12/12 layers | iter_0010 | H02 | 0 | 0 | True | True |
| SV2_TRRUST_copole | TRRUST TF–target co-pole on SV2; activation 12/12, repression 1/12 | iter_0011 | H01 | 0 | 0 | True | True |
| TF_vs_target_AUROC_joint | Joint SV2-SV7 TF/target AUROC=0.744 mean (max 0.789 at L3) | iter_0056 | H01 | 0 | 0 | True | True |
| celltype_marker_AUROC_851 | Cell-type marker AUROC=0.851 (12/12 layers p<1e-6) | iter_0023 | H01 | 0.004 | 0.00723 | False | True |
| edge_AUROC_decay | Edge-level AUROC peak 0.602 (L0); L0–L8 perm_p≤0.045 | iter_0062 | H01 | 0.045 | 0.0681 | False | True |
| string_confidence_gradient | STRING confidence gradient on SV2: rho=1.000, p=1.4e-24 | iter_0015 | H02 | 1.4e-24 | 7.47e-24 | True | True |
| bcell_attractor_convergence | BATF/BACH2 rank decay across L0→L11 (ρ=-0.972, p<1e-4) | iter_0042 | H03 | 0.0001 | 0.00026 | True | True |
| celltype_marker_AUROC_853_initial | Cell-type marker AUROC=0.853 (12/12 layers p<1e-6) | iter_0022 | H02 | 1e-06 | 3.61e-06 | True | True |
| TF_target_AUROC_cross_seed | Joint SV2-7 cross-seed: 0.744/0.753/0.757; p_perm=0.000 | iter_0057 | H01 | 0 | 0 | True | True |

## Notes / caveats
- p-values were extracted by regex from the executor's free-form `result_value`. Where the executor reported empirical/permutation p-values, those were preferred over GO Fisher p-values.
- 'Extractable' means at least one numeric p-value was parseable from the result_value string. Hypotheses with `direction='negative'` or `'mixed'` may have no headline p-value, so missingness is informative.
- Bonferroni with the paper's denominator (N=183) is the strictest defensible global correction, and is the threshold against which abstract claims are audited.
- Within-family BH/Bonferroni is reported because hypotheses within the same family are conceptually related and external readers may prefer family-wise correction over global.