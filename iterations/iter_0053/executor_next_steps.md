# Next Steps — iter_0053 → iter_0054

## Confirmed findings to build on

1. **SV5-7 L0-L2 regulatory proximity** (iter_0052 H03, fully validated this iteration):
   - rbc_resid=0.147 at L0, 95%CI=[0.084, 0.199], 100/100 bootstrap positive
   - Spectral specificity confirmed: SV5-7 = early layers only; SV2-4 = L8 only
   - This is now project's best-validated geometric finding

2. **TF→Target directional asymmetry** (H02 this iteration):
   - All 3 SV5-7 axes show consistent directional displacement at L0 (per-SV p<0.04)
   - Dramatically stronger at L8 (p_combined<0.001)
   - NEW: embedding geometry encodes regulatory direction, not just proximity

3. **SV1 = circuit identity proxy at L0-L1** (H03 this iteration, confirming iter_0051):
   - rbc=0.143 circuit vs non-circuit at L0 (p<0.001)

## Top 3 hypotheses for iter_0054

### Priority 1: TF role vs target role classification from SV5-7 (new_family)
- For each gene appearing in TRRUST, assign label: TF (has out-degree), target-only (only in-degree), or both
- Train logistic regression on SV5-7 coordinates at L0 to predict TF vs target-only role
- Evaluate AUROC
- Rationale: H02 shows systematic directional asymmetry; test if it's classifiable

### Priority 2: SV5-7 vs SV2-4 orthogonality in regulatory annotation (new_method)
- Test whether SV5-7 and SV2-4 encode DIFFERENT subsets of TRRUST edges
- Split positive edges by layer depth of maximum signal → compare annotation overlap
- Does SV5-7 capture early regulatory priming while SV2-4 captures mature context?
- Rationale: spectral decay shows clean separation; verify by edge-level analysis

### Priority 3: Cross-seed replication of SV5-7 L0 signal (split_robustness)
- Available: cycle4_immune_seed43, cycle4_immune_seed44 embeddings + edge datasets
- Replicate SV5-7 L0 rbc_resid in both alternative seeds
- Rationale: bootstrap validates single-run stability; seed replication validates data independence

## Retired directions (do not revisit)
- H1 Betti loops for regulatory circuit detection (retired iter_0052)
- TwoNN intrinsic dimensionality (retired iter_0047)
- kNN purity for cell-type markers (retired iter_0045)
- SV2-4 module structure without co-expression correction (retired iter_0051)

## SV spectral landscape summary (from H01)
| SV range | Best layer | rbc_resid | p |
|---------|-----------|-----------|---|
| SV2-4 | L8 | 0.083 | 0.0016 |
| SV5-7 | L0 | 0.148 | <0.001 |
| SV8-10 | L8 | 0.065 | 0.010 |
| SV11-14 | L8 | 0.078 | 0.003 |
| SV15-20 | L8 | 0.061 | 0.015 |
| SV21-30 | L0 | 0.071 | 0.006 |

Early-layer (L0): SV5-7 dominant, SV21-30 marginal signal
Deep-layer (L8): SV2-4 dominant, SV8-20 secondary signal
