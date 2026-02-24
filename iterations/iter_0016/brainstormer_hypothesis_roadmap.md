# Brainstormer Hypothesis Roadmap — iter_0016 → iter_0017

---

## Retire / Deprioritize

| Direction | Reason | Action |
|---|---|---|
| Single-gene module permutation tests | Multiple negatives; gene-pair level is where signal lives | `retire_now` |
| TRRUST repression anti-pole (SV2) | Tested iter_0013; definitively negative | `retire_now` |
| Bootstrap CI for act/rep differential | Underpowered at n=64 repression pairs; no new data available | `retire_now` |
| GO BP enrichment in SV2 poles | Null result iter_0011 H03 | `retire_now` |
| SV5 axis characterization (primary) | Only 6/12 layers significant, weakest axis; deprioritize unless SV4 proves interesting | `deprioritize` |
| Quintile-binned confidence gradient for SV4 alone | Already planned; run as cheap completion task, not major hypothesis | `low_priority_completion` |

---

## New Hypothesis Portfolio

### H-A: SV4 STRING Confidence Quintile Gradient
**Hypothesis**: SV4 (like SV2, SV3) encodes a monotonic STRING confidence gradient across quintile-binned pairs.
**Test**: Same protocol as H01 (iter_0016), swap SV4_IDX=3. 5 quintiles, Spearman rho.
**Expected signal**: rho > 0.7, Q5 mean_z > 2.5 (weaker than SV3 per SV4's 7/12 significant layers baseline).
**Null**: Shuffled quintile labels.
**Value**: medium | **Cost**: low
**Family**: graph_topology

### H-B: Axis-Dominant Pair GO Composition (Mechanistic Orthogonality)
**Hypothesis**: The three axis-dominant pair subsets (SV2=1547, SV3=850, SV4=695 pairs) differ systematically in GO term composition, explaining why axes are orthogonal.
**Test**: For each pair subset, extract all member genes, run Fisher exact GO enrichment (591 terms), compare top hits across subsets. Compute Jaccard of top-20 GO hits between subsets.
**Expected signal**: SV2 pairs enriched for secretory/EV; SV3 pairs enriched for kinase/immune; SV4 pairs enriched for endosome/apoptosis. Jaccard < 0.3 between any two subsets.
**Null**: Random gene sets of matching size.
**Value**: high | **Cost**: low
**Family**: module_structure

### H-C: Geneformer SV2 PPI Co-pole (Cross-Model Conservation)
**Hypothesis**: Geneformer residual embeddings encode STRING PPI proximity on their leading SVD axis with z > 3.
**Test**: Load Geneformer residual embedding tensor for matched gene set. Compute SV2 (index 1) co-pole z-score for STRING pairs (score ≥ 0.4, N=3092) across available layers. Compare to scGPT SV2 z=5+.
**Expected signal**: Geneformer SV2 mean z > 2.0, at least 6/N layers significant.
**Null**: Gene-label shuffles (N=300).
**Value**: high | **Cost**: medium (depends on artifact availability)
**Family**: cross_model

### H-D: SV2+SV3 Joint 2D Co-pole (Intersection Specificity)
**Hypothesis**: Pairs co-localized on *both* SV2 and SV3 simultaneously (joint 2D pole) are more enriched for high-confidence PPIs than pairs co-localized on either axis alone.
**Test**: For each STRING pair, compute joint 2D co-pole indicator (co-pole on SV2 AND SV3). Compare mean STRING score for joint vs. single-axis vs. neither. Test z-score at each layer.
**Expected signal**: Joint 2D pairs have higher mean STRING score and higher OR for score ≥ 0.8 than single-axis pairs.
**Null**: Single-axis co-pole rate as baseline.
**Value**: high | **Cost**: low
**Family**: manifold_geometry

### H-E: Activation vs. Repression Geometry in SV3/SV4
**Hypothesis**: The signed regulatory geometry (activation > repression co-pole) found for SV2 generalizes to SV3 and/or SV4.
**Test**: Rerun TRRUST activation vs. repression co-pole test (SV3_IDX=2, SV4_IDX=3) across 12 layers. Compare layer-significant counts vs. SV2 (12/12 activation).
**Expected signal**: SV3 shows activation co-pole at ≥ 6/12 layers; SV4 possibly fewer.
**Null**: Gene-label shuffles (N=500).
**Value**: medium | **Cost**: low
**Family**: module_structure

### H-F: STRING Triangle (3-Clique) Geometry
**Hypothesis**: Sets of three mutually interacting STRING proteins concentrate in the same SV2 pole more than pairwise interactions predict (network motif beyond-pair effect).
**Test**: Enumerate all STRING triangles (all three edges score ≥ 0.5) in the named gene set. For each triangle, compute fraction of members in same pole. Compare to expected from pairwise copole rate.
**Expected signal**: Observed fraction of "all-same-pole" triangles > binomial expectation; OR > 2 vs. random gene triples.
**Null**: Random gene triples of same size.
**Value**: medium | **Cost**: medium (enumeration step)
**Family**: graph_topology

### H-G: Intrinsic Dimension of PPI vs. Non-PPI Gene Subsets
**Hypothesis**: Genes involved in ≥ 3 STRING interactions have lower intrinsic dimensionality in residual space than non-interacting genes, reflecting their tighter manifold organization.
**Test**: Split genes into high-degree (≥ 3 PPI partners in named set) and low-degree (0–1 partners). Compute TwoNN intrinsic dimension per layer per group. Test difference vs. shuffle null.
**Expected signal**: High-degree genes: ID 1–2 units lower; effect largest at layers 0–6 (where PPI copole is strongest).
**Null**: Random subsets matching size.
**Value**: medium | **Cost**: low
**Family**: manifold_geometry

### H-H: Per-Pair STRING Score Direct Regression (vs. Copole Z)
**Hypothesis**: STRING score predicts per-pair copole z-score in a continuous linear regression (not just quintile binning).
**Test**: For each STRING pair, compute mean copole z across 12 layers. Regress on STRING score (0.4–1.0) with OLS. Report R², slope, and permutation p.
**Expected signal**: R² > 0.05, slope > 0, permutation p < 0.01.
**Null**: Permuted STRING scores assigned to same pairs.
**Value**: medium | **Cost**: low
**Family**: graph_topology

### H-I: Axis-Dominant Assignment Stability Across Layers
**Hypothesis**: Each STRING pair's dominant axis (SV2/SV3/SV4) is stable across all 12 layers vs. being layer-specific.
**Test**: For each pair, compute per-layer dominant axis (which of SV2/SV3/SV4 has highest absolute copole z). Compute assignment entropy per pair across layers (stable = entropy ≈ 0, unstable = entropy ≈ log(3)). Compare entropy distribution to shuffled-axis null.
**Expected signal**: Most pairs have entropy < 0.5 bits (dominant axis consistent across layers).
**Null**: Randomly permuted per-layer axis assignments.
**Value**: medium | **Cost**: low
**Family**: module_structure

### H-J: Persistent Homology in SV2+SV3+SV4 Projection
**Hypothesis**: The 3D subspace defined by SV2, SV3, SV4 projections has lower H1 persistence (fewer loops) than a shuffled 3D control, reflecting tight manifold organization.
**Test**: Project all 209 named genes onto (SV2, SV3, SV4) per layer. Compute H1 barcodes (ripser). Compare max H1 persistence and bar count to feature-shuffle null (N=20 replicates).
**Expected signal**: Real H1 persistence significantly lower than shuffled (denser, less loop-like); z < −2 at ≥ 8/12 layers.
**Null**: Feature-shuffle nulls (same 3D projection but shuffled gene labels).
**Value**: medium | **Cost**: medium
**Family**: topology

### H-K: SV2 Copole Stability Under Subsampling (Robustness Audit)
**Hypothesis**: SV2 PPI copole z-score is robust to gene subsampling: computing on random 100-gene subsets still yields z > 2.
**Test**: Bootstrap: draw 100 genes from named set (N_boot=100), recompute SV2 copole z for STRING pairs within subset. Report distribution of z-scores across bootstraps.
**Expected signal**: Median z > 2.0, 5th percentile z > 1.0 (robust effect).
**Null**: Same bootstrap but with STRING score shuffled.
**Value**: medium | **Cost**: low
**Family**: module_structure (robustness)

### H-L: APP Layer-Dependent Pole Transition (Mechanistic Probe)
**Hypothesis**: APP switches SV4 poles between layers 7 (bottom/endosome) and 11 (top/synapse), reflecting a model-internal representation of APP's dual endosomal-synaptic function.
**Test**: Track APP's SV4 raw projection sign and rank across all 12 layers. Confirm pole at L7 (bottom) vs. L11 (top). Compare to control genes with stable assignments. Cross-reference gene list: does PSEN1, PSEN2, BACE1 (endosomal APP processing partners) track with APP at L7?
**Expected signal**: APP flips sign between L7–L10; APP's endosomal processing partners (PSEN1/PSEN2) co-localize at L7 bottom pole.
**Null**: Random genes with similar SV4 variance.
**Value**: medium | **Cost**: low
**Family**: module_structure (mechanistic)

### H-M: SV2 Copole Enrichment at STRING Functional Network Edges Only
**Hypothesis**: STRING "experimentally determined" channel (exp score ≥ 0.4) shows stronger copole enrichment than "co-expression" or "text mining" channels, confirming physical interaction geometry not annotation-proximity.
**Test**: Separate STRING edges by evidence channel (exp, coexp, textmine). Run SV2 copole test per channel subset. Compare mean z across 12 layers.
**Expected signal**: Experimental channel: mean z > 3; text-mining channel: mean z < 2.
**Null**: Gene-label shuffles (N=300).
**Value**: high | **Cost**: medium (STRING channel data download needed)
**Family**: graph_topology

---

## Top 3 for Immediate Execution

### #1 — High-Probability Discovery: H-B (Axis-Dominant Pair GO Composition)
**Rationale**: Zero new infrastructure; uses data already computed in H02 (axis-dominant pair labels). Will mechanistically explain why axes are orthogonal — the single most important missing piece for the paper's multi-axis claim. Expected positive result based on SV4 H03 biology.
**Concise plan**: Extract gene members for each of the 3 axis-dominant pair subsets. Run Fisher exact GO enrichment (591 terms). Report top 5 GO hits per subset, Jaccard between subsets.

### #2 — High-Risk/High-Reward: H-C (Geneformer Cross-Model Conservation)
**Rationale**: If Geneformer artifacts are available, a single positive result here transforms the paper from "scGPT-specific" to "biological LLM universal." Maximum impact per experiment cost.
**Concise plan**: `find` Geneformer .npy residual embeddings. If found: align gene index to named gene set, compute SVD, run SV copole test for STRING pairs (K=52, N=300 shuffles), report z-scores vs. scGPT baseline.

### #3 — Cheap Broad Screen: H-A + H-E combined (SV4 gradient + SV3/SV4 signed regulation)
**Rationale**: Both are near-zero marginal cost extensions of already-working scripts. H-A extends H01 script to SV4_IDX=3 (5 lines). H-E extends TRRUST copole to SV3/SV4 (loop over axis indices). Together they complete the axis-by-axis characterization and close the signed-regulation story.
**Concise plan**: Modify iter_0016 H01 script: add SV4_IDX=3 loop. Modify iter_0011/iter_0012 TRRUST script: add axis_idx in [2, 3] loop.
