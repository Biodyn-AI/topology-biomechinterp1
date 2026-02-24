# Executor Next Steps — iter_0062

## Summary of iter_0062 outcomes

| Hypothesis | Direction | Decision |
|------------|-----------|----------|
| H01: edge AUROC cycle4 (735 pairs) | Positive | Promising |
| H02: AUROC mechanism (sparsity null) | Positive | Promising |
| H03: Cross-cycle Procrustes | Negative (design flaw) | Negative |

## Priority actions for iter_0063

### P1 (High): LR-based cross-cycle Procrustes transfer [H03 redesign]
**Rationale**: H03 failed because cosine AUROC is rotation-invariant. Redesign:
1. At cycle1 peak layer (L2 or L3): SVD → SV5-7, label genes as TF=1 or target-only=0
2. Train LR on cycle1 TF/target labels
3. Procrustes-align cycle4 SV5-7 → cycle1 SV5-7 coordinate frame
4. Apply cycle1 LR to cycle4-aligned embeddings
5. Evaluate: does the cycle1 classifier generalize to cycle4?
**Expected genes**: 80 shared, sufficient for classifier evaluation (but not training)
**Note**: This tests cross-condition geometric invariance of TF/target geometry.

### P2 (High): Edge-level AUROC with L-norm metric and SV6-10 exploration [H01 extension]
**Rationale**: H01 uses cosine similarity in SV5-7. Alternative metrics and subspaces untested:
- Try negative L2 distance as scoring function
- Try SV6-10 (expanding the biological subspace)
- Try SV1 alone (dominant axis)
**Expected**: Different subspace may show stronger or layer-shifted AUROC peak.

### P3 (Medium): Effective rank vs AUROC correlation [H02 extension]
**Rationale**: Sparsity (nz_count) is constant, but effective rank varies by layer. Test:
- Compute eff_rank at each layer (from SV fraction-explained)
- Spearman(AUROC, eff_rank) to check if another confound explains the decline
**Data**: Already computed in prior iterations (iter_0058/H03), just cross-correlate.

### P4 (New family): Persistent homology of TRRUST edges in SV5-7 at L0 vs L9
**Rationale**: AUROC signal strongest at L0, absent at L9. Does topological structure (Betti-1) change correspondingly?
- Compute Vietoris-Rips Betti-0/1 on (TF+target) gene embedding neighborhoods at L0 vs L9
- Compare TF-only vs target-only topological profiles
- Use giotto-tda or ripser

### P5 (New): kNN graph assortativity by TF/target label at each layer
**Rationale**: If TF-target cosine similarity is highest at L0, the kNN graph should show higher TF-target assortativity at L0. Test with networkx assortativity coefficient.

## Retired directions (do not re-enter without rescue rationale)
- FFL geometric ordering (permutation null ~0.50 by symmetry)
- L8 CKA boundary (consecutive CKA uniform)
- Dual-role gene ambiguity (all target-like)
- Topology stability via principal angles (negative, iter_0056/H03, iter_0057/H02)

## Key open questions
1. Why does SV2-4 AUROC fall below 0.50 (sub-null)? Is there active repulsion of TF-target pairs in that subspace?
2. cycle4_immune edge AUROC peak is at L0, not L2-L3 as in TF/target classification. Why?
3. Is the L10-L11 cross-seed CKA drop related to the AUROC collapse at L9?
