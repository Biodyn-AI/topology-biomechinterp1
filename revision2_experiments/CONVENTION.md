# Declared analysis convention (round 2)

Every round-2 number is computed under exactly this convention. Where round 1 used a different
choice, the difference is recorded and a sensitivity analysis is reported.

| choice | round-2 canonical | round-1 (as deposited) | why |
|---|---|---|---|
| embedding matrix | `cycle1_main` (12, 4803, 512) | mixed: `cycle1_main` in revision code, `cycle4_immune_main` (4941) in the loop | 4,803 is what Methods 4.1 states; the revision code uses it |
| gene set | curated 209 from the cycle1 edge list, minus zero-norm at L0 -> **187** | 209 (E3) or "195" (manuscript, unreconstructible) | the manuscript's own stated OOV rule, applied consistently |
| SVD | on the full 4,803 x 512 layer matrix, mean-centred, then restricted to the gene set | full (E3) or 209-subset (loop `iter_0010`) | the axes should be properties of the model's representation, not of a curated subset; subset-SVD reported as sensitivity |
| poles | K = int(0.27 n) | K = 52 fixed (loop) or 0.27n (E3) | a fraction is comparable across gene sets of different size |
| TRRUST pairs | **deduplicated** distinct (TF, target) | not deduplicated (326 = row count) | row counts weight each pair by its number of supporting PMIDs, an undisclosed literature-bias confound |
| STRING pairs | distinct pairs, combined score >= 0.7 | same | -- |
| null | gene-label shuffle, as Methods 4.4 states | pair resampling (loop) despite the Methods text | tests whether the *axis* is special, which is the claim being made |
| permutations | **N = 100,000** (co-pole), **N = 10,000** (AUROC) | 50-1,000 | Bonferroni at alpha = 2.73e-4 requires N > 3,663 for the p-floor to clear it |
| empirical p | (b+1)/(N+1) | b/N, reported as 0.0 when b = 0 | a permutation p of exactly 0 does not exist |
| TF/target score | leave-one-out mean cosine to the TF set | self-similarity (= 1) included in each TF's own mean | self-inclusion adds 1/\|TF\| to every TF and nothing to any target: a systematic upward bias of +0.016 to +0.029 AUROC |
| GO enrichment | Fisher exact (analytic, not resolution-limited) + permutation cross-check | Fisher exact | -- |

**Cross-seed replication** is the one analysis that cannot use `cycle1_main`: the replication seeds
exist only in the `cycle4_immune_{main,seed43,seed44}` family (4,941-gene vocabulary). It is therefore
reported separately, on its own vocabulary, and labelled as such.
