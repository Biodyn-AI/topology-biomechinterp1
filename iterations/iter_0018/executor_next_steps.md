# Next Steps — iter_0018 → iter_0019

## Top Priority (iter_0019)

### 1. Attention geometry: TRRUST repression + per-layer test (high priority)
The H03 result (TRRUST activation MW_p=9.9e-09 via attention) is novel and needs follow-up:
- Test TRRUST repression pairs separately (are they also co-attended, or less so than activation?)
- If per-layer attention data exists, test whether attention encoding is layer-specific
- Test whether the TRRUST activation attention signal survives a gene-label shuffle null
This strengthens the mechanistic dissociation claim (attention=regulatory, SVD=PPI).

### 2. Precision@k using co-pole membership score (novel framing, graph_topology)
H02 showed 1.2× enrichment using SV2 distance. A better predictor:
- For each pair (i,j), compute co-pole rate across 12 layers × 3 axes (SV2+SV3+SV4)
- Score = fraction of (layer, axis) combinations where both genes are in same top-K or bottom-K pole
- Re-run precision@k with this composite score → expect much higher enrichment (~3-5×)

### 3. Geneformer cross-model validation (critical for paper, cross_model_alignment family)
Still the biggest gap (17 iterations, 0 cross-model results).
- Use `cycle12_geneformer_lung_bootstrap/geneformer_edge_dataset.tsv` and processed.h5ad
  Geneformer mean embeddings (need to check if available or compute from token embeddings)
- Fallback: use the pre-computed centered_cosine per-edge similarities from Geneformer
  and test if they show a STRING confidence gradient similar to scGPT SV2

### 4. SV2 co-pole membership as a predictor + out-of-sample hold-out test
Reserve 20% of STRING pairs as test set, train on the rest, measure co-pole enrichment in held-out set. This is a clean out-of-sample benchmark.

## Secondary

- Bootstrap stability of TRRUST attention signal (re-run 500× bootstrap of TRRUST pairs within attention)
- Test attention with TRRUST repression (expect weaker than activation, interesting if positive)
- Synthesize mechanistic model: scGPT uses attention for regulatory detection, residual stream for PPI organization

## Retired directions (maintain)
- TRRUST signed direction on SV3/SV4 (iter_0017 H03, 0/12 significant)
- SV2/SV3 GO enrichment sweeps (done, biologically characterized)
- Axis independence (done, paper-ready)

## Paper status
- Sections needed: (1) Random null confirms model-specificity, (2) Attention co-occurrence encodes regulatory pairs (new finding), (3) SVD vs Attention mechanical dissociation
- Most important missing: Geneformer replication for cross-model validity claim
