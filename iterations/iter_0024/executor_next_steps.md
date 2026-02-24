# Next Steps — iter_0024

## Summary of iter_0024 outcomes

- **H01 POSITIVE**: Dorothea A/B confidence → AUROC=0.671 at L8, 12/12 layers sig, perm_p=0. Independent regulatory database replicates embedding proximity signal.
- **H02 NEGATIVE**: TF activation degree does not predict geometric centrality (0/12 sig).
- **H03 POSITIVE**: GO CC > MF > BP for embedding proximity signal. CC: Spearman=0.106 at L8, all 12 layers sig. BP shows layer-deepening pattern.

## Prioritized next hypotheses (iter_0025)

### H-A (HIGH PRIORITY): Multi-predictor joint model — STRING + Dorothea + GO CC
- Build a multivariate linear model: L2_distance ~ STRING_score + Dorothea_conf + GO_CC_Jaccard
- Test if all three predictors contribute independently (partial R², VIF)
- This would consolidate the three main positive signals into a unified quantitative claim
- Family: manifold_distance / module_structure hybrid
- Method: OLS/ridge regression with all three predictors + partial Spearman

### H-B (HIGH PRIORITY): GO CC subcellular compartment analysis
- Identify which specific CC terms drive the signal (nucleus, cytoplasm, membrane, mitochondria?)
- Genes sharing nucleus terms vs cytoplasm terms — which cluster more tightly?
- Could reveal compartment-specific geometric organization in scGPT manifold
- Family: module_structure (new_method)
- Method: Extract top CC terms, group genes by compartment, per-compartment AUROC

### H-C (MEDIUM): JUN activation cluster test
- JUN has 17 activation targets in TRRUST (top activator)
- Test if JUN's activation targets form a compact geometric cluster vs. random gene sets
- Family: module_structure (new_method — single-TF cluster analysis)
- Method: Mean pairwise L2 of JUN targets vs permutation control (1000 random same-size sets)

### H-D (NEW FAMILY): Chromosomal proximity negative control
- Test: do genes on the same chromosome/chromosome arm cluster in embedding space?
- Expected: NEGATIVE (chromosomal proximity ≠ functional proximity)
- This would confirm the signal is biologically functional, not genomic artifact
- Family: null_sensitivity (new_method)
- Method: Map named genes to chromosomes (mygene), compute AUROC same-chr vs different-chr

## Retired directions (from iter_0024)
- TF activation hub centrality (H02): first negative — monitor

## Most pressing open questions
1. Do STRING, Dorothea, GO CC signals survive when controlling for each other?
2. Which subcellular compartments drive GO CC proximity signal?
3. Is the TRRUST activation/repression asymmetry replicable in Dorothea if direction info were available?
