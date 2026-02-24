# Next Steps — iter_0057

## Key Findings
1. **H01 confirmed**: Joint 6D (SV2-7) TF/target AUROC=0.751 is reproducible across 3 seeds. This is now a cross-validated, paper-ready claim.
2. **H02 negative**: Basis permutation ruled out. SV subspace misalignment (29°) is real but doesn't explain seed43 anomalies.
3. **H03 structural**: SV energy spreading (lower joint fraction → higher AUROC) is a strong correlation (rho=-0.93). This describes the geometry but needs interpretation.

## Priority Actions for iter_0058

### High priority
1. **Cross-seed gene-level consistency** (module_structure, new_method):
   - Do the *same* TF genes rank in the top decile of SV5-7 discriminative score across all 3 seeds?
   - Compute per-gene LR coefficients at L0 and L3 for each seed; Spearman correlation across seeds.
   - Positive result = genes are stable markers, not just seed-specific artifacts.

2. **Graph Laplacian spectral alignment with SV5-7** (new family: cross_model_alignment or module_structure):
   - Compute TRRUST regulatory graph Laplacian, extract leading eigenvectors.
   - Test whether graph spectral directions align with SV5-7 directions (principal angles).
   - This is potentially a paper-defining result: biological network topology encoded in embedding geometry.

3. **OOD generalization across cell types** (split_robustness):
   - Train joint 6D classifier on main seed (immune cells), test on a different cell type dataset.
   - Tests whether TF/target geometry is cell-type specific or universal.

### Medium priority
4. **H03 cross-seed validation**: Does rho(sv_joint_frac, AUROC) = -0.93 hold for seed43 and seed44?
5. **Persistent homology on 6D TF vs target subspace** (persistent_homology, rescue with 6D input instead of 3D):
   - Prior attempts used 3D; 6D provides richer topology signal. Previously retired but 6D is materially different.

## Retired directions
- Basis permutation (H02 negative): definitively ruled out
- Subspace rotation → AUROC (iter_0056 H03): negative

## Hypothesis retirement status
- `module_structure` (SV2-4 only): retired after iter_0051 H01 negative
- `persistent_homology` (3D): two negatives; eligible for 6D rescue
- Basis permutation (specific to seed43): retired now
