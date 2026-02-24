# Executor Next Steps — iter_0027

## Key Findings to Build On

1. **H01 negative but descriptive**: PC1 at L11 (67.6% EV) is not captured by TF-hood, HLA, AP1,
   or hub degree. Observed pattern: LCK/JUN at top, FOS/HLA-A/PAX5 at bottom — suggests cell-state
   polarity (T-cell effector vs B-cell/antigen presentation) but not confirmed statistically.

2. **H02 negative**: No population-level directional encoding of Act vs Rep TRRUST pairs. Individual
   TFs show variable deltas (ETS1 +10.85 but n=2). The Act/Rep question may require larger data.

3. **H03 important negative**: Dorothea threshold sweep shows AUROC < 0.5 with population-level null,
   contradicting earlier AUROC=0.658 (iter_0026). The discrepancy = null-construction sensitivity.
   Prior "manifold_distance positive" findings likely partly reflect hub-gene centrality confound.

## Retirement Candidates
- **Dorothea as primary regulatory proximity signal**: Retire until a hub-corrected null is designed.
  (Two consecutive negatives with changed method = retire per policy.)
- **Simple TF/HLA/AP1 binary features on PC1**: Confirmed negative after directed test.

## Recommended Hypotheses for iter_0028

### Priority 1 (high-reward, fills a gap exposed by H03 methodological critique):
**Hub-corrected regulatory proximity test**
- Build matched null: for each regulatory pair (TF→target), sample null pair (TF→random_target)
  using the SAME TF. This controls for TF centrality entirely.
- Test: within each TF, are its known targets closer than its randomly-sampled non-targets?
- Family: manifold_distance, new_method (rescue of retired direction with materially changed null)

### Priority 2 (new family, leverages existing data):
**Cell-type gene axis validation at L11**
- Follow up on the descriptive pattern from H01: LCK/ITGB2 (T-cell markers) at positive PC1,
  HLA-A/HLA-DPB1/PAX5 (B-cell/antigen-presentation) at negative PC1.
- Build cell-type label for 50+ of the 209 genes from literature.
- Test: does PC1 at L11 separate T-cell vs B-cell associated genes?
- Family: module_structure, new_method (cell-type axis at L11)

### Priority 3 (topology screening):
**Persistent homology betti curves at intermediate layers (L3, L6, L9)**
- Prior persistent homology tests (iter_0020, 0021) focused on L8. Now that we know PR collapses
  monotonically, test if Betti curves also change monotonically or show a phase transition at L10.
- Family: persistent_homology, new_method (multi-layer betti comparison)

## Methodological Notes
- For all future regulatory proximity tests: use within-TF null (same TF, different targets)
  instead of all-gene random null.
- Hub-corrected test is the cleanest way to rescue the manifold_distance hypothesis family.
- The dimensional collapse finding (PR 21→1.68) remains robust and is the strongest structural finding.
