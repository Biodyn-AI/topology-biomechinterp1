# Next Steps — iter_0054

## Top priorities for iter_0055

### 1. Cross-seed replication of H01 + H03 (HIGH PRIORITY)
- Run H01 (LR AUROC) and H03 (layer trajectory) on cycle4_immune_seed43 and cycle4_immune_seed44
- Confirm: AUROC > 0.65 and peak trajectory at L9 in both seeds
- This would make H01+H03 the strongest validated findings in the project

### 2. SV2-4 control trajectory for H03 (CHEAP)
- Compute the same directionality trajectory (combined magnitude + Cohen's d) using SV2-4 subspace at all 12 layers
- Expectation from iter_0053 H01: SV2-4 should show opposite pattern — weak at L0, strong at L8
- Creates a complete "two-regime" picture: SV5-7 amplifies early→late; SV2-4 stronger late

### 3. Expanded feature set for H01 (CHEAP)
- Test logistic regression with SV1-7 combined (7D) vs SV5-7 alone
- Test SV1-7 at L8 for TF classification
- Tests whether combining early and late subspaces improves classification

### 4. Partial Spearman for H02 (MEDIUM)
- Control for TF/target binary label when testing depth-SV correlation
- Isolates "cascade position" signal from "TF vs target" binary signal
- If BFS depth independently predicts SV coordinates beyond TF/target label, that is a new finding

## Hypothesis families to explore further
- **cross_model_alignment**: Test if Geneformer shows analogous directionality pattern
- **null_sensitivity**: Feature-shuffle null for H01 (permute SV coordinates per gene, not labels)
- **topology_stability**: Bootstrap H03 magnitude at each layer (is L9 robustly the peak?)

## Retired directions (do not revisit without rescue rationale)
- TwoNN intrinsic dimensionality
- kNN graph topology (cell-type purity)
- SV1 alone as standalone test
- Betti loops (PH on small gene sets)
- Continuous degree→SV1 correlation

## Key open questions
1. Does the directionality trajectory peak at L9 cross-seed?
2. Does SV2-4 show a complementary (later-layer amplification) trajectory?
3. Can we directly classify TF vs target-only with >0.70 AUROC using more features?
4. Is the BFS depth effect independent of TF/target binary classification?
