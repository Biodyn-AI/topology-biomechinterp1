# Executor Next Steps — iter_0056 → iter_0057

## What was learned this iteration

1. **Joint SV2-7 AUROC is layer-stable at 0.72+ through L0-L8** (9/12 layers; all p_perm=0.000). The joint classifier consistently outperforms both individual subspaces, confirming complementary encoding. Ready for Claim 50.

2. **L9 crossover replicates in 2/3 seeds** (main + seed44 both show L9 crossover; seed43 does not). Claim 49 is partially validated but requires seed43 investigation before promotion.

3. **Subspace rotation angle does not predict AUROC trajectory** (Spearman rho=-0.27, p=0.42). The L2→L3 rotation is large (69.7°) but does not reduce AUROC. This direction is retired.

---

## Priority actions for iter_0057

### H01 (new, high-priority): Cross-seed validation of joint SV2-7 classifier
- Run the same joint 6D LR classifier on seed43 and seed44 embeddings
- Check if mean AUROC ≥ 0.72 at 9/12 layers holds in all three seeds
- Expected output: `h01_joint_auroc_crossseed.csv`
- If confirmed → Claim 50 is fully cross-seed validated

### H02 (refinement): Seed43 crossover anatomy
- For seed43, plot sv24_mag and sv57_mag trajectories across 12 layers
- Determine: does sv57_mag approach sv24_mag at L9 (near-crossover) or does the pattern differ qualitatively?
- Check if seed43 nonzero gene set has fewer valid TRRUST pairs that affect directionality
- Expected output: `h02_seed43_anatomy.csv`

### H03 (new, novel): Gene-level anatomy of SV5-7 at L8 minimum
- At L8 (SV57 AUROC minimum): which specific TF genes have the lowest SV5-7 projections?
- Rank TF genes by their L8 SV5-7 magnitude; identify top/bottom 10
- Compare L0 vs L8 vs L11 gene rankings — is there a systematic reordering?
- Novel family: intrinsic_dimensionality × module_structure
- Expected output: `h03_gene_ranking_L8.csv`

### H04 (new, novel): Wavelet/multi-scale analysis of AUROC trajectory
- Model the SV5-7 AUROC trajectory [12 values] as a 1D signal
- Decompose with DWT (pywt) and identify dominant frequency components
- Test if the mid-depth dip (L4-L8) corresponds to a specific spatial frequency
- Compare with SV2-4 trajectory as control
- Novel family: topology_stability (temporal regularity)
- Expected output: `h04_auroc_wavelet.csv`

---

## Retired directions

- **H03 (subspace rotation angle)**: negative — rotation angle does not predict AUROC. Maximum principal angle is uninformative as a scalar summary.
- **Persistent homology** (iter_0052 H02 negative, iter_0051 H03 inconclusive): retired.

---

## Paper update notes

- Add Claim 50: "The joint SV2-7 (6D) subspace provides a layer-stable TF/target classifier with mean AUROC 0.744 (null~0.505, all p_perm=0.00), exceeding either individual subspace at 11/12 layers."
- Update Claim 49 status: "2/3-seed validated; seed43 shows no crossover — claim conditioned on main and seed44."
- Note H03 negative: subspace rotation angle does not account for AUROC trajectory.
