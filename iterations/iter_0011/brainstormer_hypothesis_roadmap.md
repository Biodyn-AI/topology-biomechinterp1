# Brainstormer Hypothesis Roadmap — iter_0011 → iter_0012+

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| GO Biological Process → SV2 enrichment | **RETIRE NOW** | 0/591 terms significant at 3 layers. Fully exhausted. |
| Persistent homology rewiring-null | **RETIRED** (iter_0010) | Persistent negative, no rescue path. |
| GO BP cosine clustering | **RETIRED** (iter_0010) | Same data, same null outcome. |
| kNN graph transitivity | **RETIRED** (iter_0002) | Negative, no follow-up signal. |
| Mito L3→L4 pole-flip | **DEPRIORITIZE** | Falsified mechanistically; rank-based artifact explained. No further flip work needed. |

---

## New Hypothesis Portfolio

### H-A: Activation vs Repression co-pole differential — bootstrap test
- **Hypothesis**: The act_copole_rate − rep_copole_rate difference is statistically significant across layers (not just visually apparent).
- **Test**: Per-layer bootstrap (N=10,000 resamples of edge pairs with replacement). Compute Δ = obs_act − obs_rep. Report 95% CI, one-sided p against Δ=0. Also compute standardized effect size vs respective nulls.
- **Expected signal**: Δ > 0 at most layers with non-overlapping CIs; meta-analytic pooling gives FDR-controlled global p.
- **Null**: Δ = 0 (no difference between edge types).
- **Value**: high | **Cost**: low

### H-B: STRING PPI co-pole in SV2 (cross-graph convergent validation)
- **Hypothesis**: Physically interacting gene pairs (STRING high-confidence PPI, ≥700 score) are also co-localized in SV2 poles above null.
- **Test**: Map STRING human PPI to 209-gene set. Compute co-pole rate for PPI pairs vs N=1000 gene-label shuffles. Test at all 12 layers. Compare effect size to TRRUST.
- **Expected signal**: emp_p < 0.05 at multiple layers; effect size < TRRUST (PPI is less specific to regulatory geometry) but present.
- **Null**: Gene-label shuffle preserving degree distribution.
- **Value**: high | **Cost**: low

### H-C: SV1 axis identity — what does the first principal component encode?
- **Hypothesis**: SV1 encodes a biologically interpretable axis distinct from SV2's localization axis — likely cell-type, proliferation, or lineage.
- **Test**: Compute SV1 projections for 209 genes. Run same CC/MF enrichment screen as iter_0006+ but for SV1 poles. Check if TRRUST co-pole also holds for SV1. Compare CC enrichment ORs for SV1 vs SV2.
- **Expected signal**: SV1 encodes different CC categories than SV2 (e.g., nucleus vs cytoplasm), or shows no CC enrichment and instead tracks lineage markers.
- **Null**: Gene-label shuffles.
- **Value**: high | **Cost**: low

### H-D: Multi-SV joint geometry (SV1×SV2 joint clustering of regulatory pairs)
- **Hypothesis**: TF-target pairs cluster more tightly in joint SV1×SV2 space than in either axis alone — regulatory structure lives in a 2D subspace.
- **Test**: Project 209 genes into (SV1, SV2) 2D space. Compute Euclidean distance between TF-target pairs vs shuffled pairs. Also test: does TRRUST edge density increase in SV1×SV2 quadrant-sharing vs pole-sharing? Use KDE-based density test.
- **Expected signal**: TF-target Euclidean distance distribution is left-shifted vs null in 2D; quadrant co-occupancy enriched vs 1D pole co-occupancy.
- **Null**: Label shuffle, N=1000.
- **Value**: high | **Cost**: medium

### H-E: SV3/SV4/SV5 — functional enrichment screen of higher SVD components
- **Hypothesis**: SV3, SV4, or SV5 of the embedding matrix encode additional functional categories not captured by SV2 (e.g., stress response, cell cycle, immune activation).
- **Test**: For each of SV3–SV5, compute poles (top-K=40). Run Fisher exact for GO CC and GO MF (not BP, which is exhausted for SV2). Empirical null N=500.
- **Expected signal**: At least one higher SV encodes a distinct, interpretable CC or MF enrichment (e.g., nuclear envelope, ribosome, plasma membrane).
- **Null**: No enrichment above null in SV3–SV5.
- **Value**: medium | **Cost**: low

### H-F: TF hub-degree vs co-pole tightness
- **Hypothesis**: TFs with more targets in the 209-gene set show tighter co-pole clustering with their targets (higher TF-specific co-pole rate).
- **Test**: For each TF in TRRUST with ≥3 targets in named set, compute TF-specific co-pole rate (fraction of its targets in same pole as TF). Correlate with TF out-degree. Test at layers 4 and 11.
- **Expected signal**: Spearman ρ > 0 between TF degree and TF co-pole rate; high-degree TFs (NFKB1, TP53, SP1) show co-pole rates >0.4.
- **Null**: Degree-matched random TF-gene assignments.
- **Value**: medium | **Cost**: low

### H-G: Repression-specific geometry at layer 4 — what is special?
- **Hypothesis**: Layer 4 is architecturally special for repression geometry: either SV2 variance peaks at L4, or repressor TFs specifically cluster with targets at L4 before dispersing.
- **Test**: Compare SV2 variance (explained variance fraction) across 12 layers for repression pairs only. Check if repressor TFs at L4 are more pole-concentrated (smaller spread in SV2 projection). Identify which specific repressor pairs drive the L4 signal.
- **Expected signal**: SV2 variance fraction peaks near L4; repressor pairs that are co-pole at L4 shift to opposite poles by L7+.
- **Null**: L4 has same SV2 variance as surrounding layers.
- **Value**: medium | **Cost**: low

### H-H: HIF3A multi-flip characterization
- **Hypothesis**: HIF3A's multi-flip behavior (L0→1, L8→9, L10→11) reflects its unique position in the regulatory network as a hypoxia-response integrator with bidirectional targets.
- **Test**: Extract HIF3A embedding trajectory (L2 norms, SV1/SV2 projections at all layers). Map HIF3A's TRRUST regulators and targets in 209-gene set. Compare trajectory to stable mito genes (CCR7, PMAIP1 as top-pole stable; BNIP3, PDHA1 as stable bottom-pole). Test: does HIF3A flip correlate with changes in its TF regulator pole assignments?
- **Expected signal**: HIF3A regulators also flip at layers 8 and 10; HIF3A L2 norm spikes at flip layers.
- **Null**: Stable mito gene controls show no correlated regulator flips.
- **Value**: medium | **Cost**: low

### H-I: GO Molecular Function enrichment in SV2 poles
- **Hypothesis**: GO MF terms (not BP, which was negative) are enriched in SV2 poles — e.g., DNA-binding TFs cluster in one pole, structural proteins in another.
- **Test**: Filter GO MF terms to 3–50 named genes. Run same pipeline as H03 but with MF ontology. Empirical null N=200.
- **Expected signal**: DNA binding (GO:0003677) and transcription factor activity enriched in one SV2 pole; cytoskeletal/structural MF terms enriched in the other.
- **Null**: 0 MF terms significant (as with BP).
- **Value**: medium | **Cost**: low

### H-J: STRING functional coupling (co-expression + text-mining) co-pole
- **Hypothesis**: STRING gene pairs with high combined text-mining + co-expression scores show SV2 co-pole enrichment distinct from PPI.
- **Test**: Filter STRING to text-mining + coexpression channel scores ≥400 (separate from PPI channel). Map to 209 genes. Compare to PPI co-pole effect size.
- **Expected signal**: Text-mining pairs show similar or stronger co-pole signal than PPI, suggesting semantic co-mention correlates with embedding geometry.
- **Null**: Label shuffles.
- **Value**: medium | **Cost**: low

### H-K: Cross-layer co-pole consistency per TF-target pair
- **Hypothesis**: Some TF-target pairs are co-pole at all 12 layers ("stable co-pole pairs"), while others are layer-specific. Stable co-pole pairs are more enriched for known direct binding (TRRUST high-confidence).
- **Test**: For each of 333 TRRUST pairs, compute fraction of layers in which the pair is co-pole. Identify pairs with 12/12, 11/12, 10/12 co-pole consistency. Test if stable co-pole pairs are enriched for activation mode vs repression.
- **Expected signal**: ~30–50 pairs are co-pole at 10+/12 layers; these are predominantly activation edges.
- **Null**: Shuffled-label pairs, same K.
- **Value**: medium | **Cost**: low

### H-L: SV2 pole stability vs expression variance correlation
- **Hypothesis**: Genes that are stably in the same SV2 pole across all 12 layers have lower expression variance across cell types in the training data, suggesting the embedding geometry reflects expression stability.
- **Test**: For each of 209 named genes, compute flip frequency across 12 layers. Correlate with any available expression variance or gene dispersion metric from scGPT pretraining context. If no direct metric, use proxy: mito vs non-mito, housekeeping marker genes vs regulatory genes.
- **Expected signal**: Low-flip genes are more often housekeeping / constitutively expressed; high-flip genes are regulatory and context-sensitive.
- **Null**: No correlation between flip frequency and gene stability class.
- **Value**: medium | **Cost**: medium

### H-M: SV2 asymmetry — top vs bottom pole size and gene composition
- **Hypothesis**: The SV2 axis is compositionally asymmetric: the top and bottom poles encode distinct functional modules, not just opposite ends of the same axis.
- **Test**: At layers 4 and 11, compute SV2 projection distribution for all 209 genes. Compare gene composition of top-K=52 vs bottom-K=52 by: GO CC distribution, TF vs non-TF ratio, TRRUST edge density within each pole. Test asymmetry with Fisher exact for each category.
- **Expected signal**: One pole enriched for TFs + nucleus; other enriched for cytoplasmic structural genes. TF-TF pair density higher in one pole.
- **Null**: Top/bottom poles are functionally symmetric.
- **Value**: medium | **Cost**: low

### H-N: SV2 geometry across scGPT model sizes — is this a scale-independent property?
- **Hypothesis**: The SV2 compartment and regulatory geometry is present in scGPT regardless of model size, implying it emerges from pretraining data structure not architecture scale.
- **Test**: If subproject_38 has embeddings from a smaller/larger scGPT variant, run SV2 co-pole for TRRUST pairs and CC enrichment. If only one model available, compare first 6 vs last 6 layers as a proxy for "shallow" vs "deep" processing.
- **Expected signal**: Same qualitative result in shallow-layer subset as full 12-layer; robustness across depth implies early-layer encoding.
- **Null**: SV2 co-pole only emerges in deep layers of full model.
- **Value**: high | **Cost**: medium

### H-O: Activation network sub-graph clustering coefficient in SV2 space
- **Hypothesis**: The activation-only TRRUST sub-graph forms a denser geometric cluster in SV2 space than a random sub-graph of equal size — not just co-pole, but spatially concentrated.
- **Test**: For activation pairs, compute mean pairwise SV2 projection distance (not just pole assignment). Compare to N=1000 random same-size sub-graphs from 209 genes. Test whether the spatial spread (variance) of activation-pair genes in SV2 is smaller than null.
- **Expected signal**: Activation-pair gene set has smaller SV2 spread than random same-size gene sets at most layers; spatial concentration increases toward layer 4.
- **Null**: Same spread as random gene sets.
- **Value**: high | **Cost**: low

---

## Top 3 for Immediate Execution

### Candidate 1 — HIGH PROBABILITY DISCOVERY
**H-A: Activation vs Repression bootstrap differential**

This converts the visually striking but statistically informal act/rep difference into a publication-grade claim with CIs. The data is already in hand (h01_trrust_copole_12layers.json has obs_act and obs_rep per layer). Cost: trivially low. Probability of confirming: >90% given the raw numbers already show consistent Δ. This completes the signed regulatory geometry finding.

### Candidate 2 — HIGH RISK / HIGH REWARD
**H-D: Multi-SV joint geometry (SV1×SV2 2D regulatory clustering)**

If TF-target pairs cluster in 2D embedding subspace beyond what either axis alone captures, this reframes the story: the regulatory network geometry is multi-dimensional, not just SV2. This would be a significant structural finding. Risk: the 2D signal may not exceed 1D; the gain over SV2 alone may be marginal.

### Candidate 3 — CHEAP BROAD SCREEN
**H-B: STRING PPI co-pole**

Reuses existing infrastructure (just swap edge list). If STRING PPI pairs also show SV2 co-pole enrichment, we have a cross-graph-type convergent result: both regulatory and physical interactions are geometrically organized. Tests in ~5 minutes of compute at 12 layers with N=1000 shuffles. Low cost, high interpretability if positive.

---

## Execution Order for iter_0012

1. H-A (bootstrap act/rep test) — completes H01 line, no new data needed
2. H-B (STRING PPI co-pole) — new edge data, existing pipeline
3. H-C (SV1 identity screen) — expand the SVD story before going multi-SV

Run H-D (multi-SV) in iter_0013 once SV1 identity is characterized.
