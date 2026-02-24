# Executor Next Steps — iter_0019

## Current Status

Three hypotheses tested. All three produced quantitative results:
- H01: Cross-model alignment confirmed (GF MW_p=7.8e-127, scGPT–GF Spearman_abs=0.446)
- H02: Multi-axis composite enrichment 2.18x (monotonic gradient)
- H03: Repression=Activation=~2x attention enrichment; STRING has NO attention enrichment

## Top Priorities for iter_0020

### 1. Contextualized Geneformer Gene Embeddings (CRITICAL)
The iter_0019 H01 test used static Geneformer token embeddings (pre-contextualization). The cycle12 bootstrap data used contextualized embeddings (AUROC delta=0.039). The true cross-model test requires contextualized GF gene embeddings at the cell level, aggregated across cells. This requires running Geneformer inference on the lung dataset. If compute is available, this is the highest-priority next step — it would complete the cross-model alignment story.

### 2. Multi-Axis Composite P@k with Null Controls (HIGH)
H02 shows monotonic enrichment up to 2.18x. Next: add a shuffle null for the 3-axis count (shuffle gene labels within each axis independently), to confirm that the monotonic gradient is not an artifact of pole definition. Also: test quintile stratification within the 3-axis stratum (ranking by magnitude sum of SV2+SV3+SV4 projections).

### 3. Attention-SVD Joint Score (NEW FAMILY - attention geometry)
H03 establishes: attention=regulatory, SVD=PPI. Test whether combining both scores gives additive predictive power for TF-target prediction. For each named gene pair: compute score = (SVD co-pole count) + (attention > threshold). Does this outperform either alone for TRRUST pairs?

### 4. STRING Quintile × Multi-Axis (REFINEMENT)
From iter_0016/0017: STRING quintile gradient z-scores reach ~5 at high-confidence edges for SV2 alone. Re-run with composite 3-axis score: expect even steeper quintile gradient.

### 5. Geneformer Contextualized Fallback — Use Existing cycle12 Edge Data
If Geneformer inference is too slow, use cycle12 geneformer_feature_metrics.csv which already has per-edge Geneformer scores. Map named gene pairs to cycle12 edge dataset and compute correlation of Geneformer centered_cosine score with scGPT SV2 co-pole membership.

## Directions to Avoid (Retired)
- SV3/SV4 signed direction — negative in iter_0016
- GO sweeps — negative in iter_0011
- TRRUST signed direction in SVD geometry — negative in iter_0013/0017
- SV2 1D distance framing — weak (1.2x) in iter_0018
- Random Gaussian null — confirmed negative in iter_0018 (no need to rerun)
