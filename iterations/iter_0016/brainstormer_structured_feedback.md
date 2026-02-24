# Brainstormer Structured Feedback — iter_0016

## Gate Status
`passed_min_research_gate: true`. All three hypotheses yielded positive, quantitative results. No recovery action needed.

## What This Iteration Established

### H01 — SV3 STRING confidence gradient (rho=0.90, p=0.037)
Confirms monotonic confidence encoding is not SV2-exclusive. The Q2/Q3 tie (1.447 vs 1.397) is a minor wrinkle; the Q4→Q5 jump (1.906→3.813) is the real signal. SV3 encodes high-confidence PPIs with the same directional logic as SV2, just with less perfect monotonicity.

**Assessment**: Solid confirmation. Extends the confidence-gradient story to a second axis. Marginal statistical significance (p=0.037) due to n=5 quintile points — this is a structural limitation of the test design, not evidence of weakness.

### H02 — SV2/SV3/SV4 axis independence (max r=0.247)
Near-zero inter-axis copole correlations confirm the multi-axis framing is not redundancy. SV2 dominance (1547/3092 pair-assignments) is expected given its higher mean z-scores; SV3+SV4 together cover half the pairs independently. The Jaccard overlaps (0.13–0.18) indicate partial but not full gene membership overlap — axes are geometrically distinct but biologically related.

**Assessment**: Clean result. Max r=0.247 (SV3–SV4) is low enough to call these independent axes. The dominant-axis pair count framing is publication-worthy.

### H03 — SV4 GO biology profile
Endosome/apoptosis (bottom pole) vs. antiviral/T-cell/synapse (top pole) is a coherent biology. The contrast with SV2 (secretory/EV) and SV3 (kinase/immune/MHC-II) is real. APP appearing in both endosome (L7 bottom) and synapse (L11 top) is interesting — APP is an endosomal processing substrate with synaptic function; the layer-dependent pole shift merits tracking.

**Assessment**: Biologically interpretable SV4 profile. The endosome→synapse transition for APP across layers could be a hook for the paper's narrative.

## Key Narrative Position After iter_0016

The project now has a **4-axis multi-dimensional PPI geometry story**:
- SV2: secretory/EV + PPI (rho=1.00 confidence gradient, z~10 mean, 12/12 layers)
- SV3: kinase/immune/MHC-II + PPI (rho=0.90 confidence gradient, z~7.2 mean, 12/12 layers)
- SV4: endosome/apoptosis/antiviral + PPI (7/12 layers, weaker)
- SV5: partial PPI signal (6/12 layers)

The axes are near-orthogonal (max r=0.247). This is the central organizing structure for the paper.

## What's Missing for a Complete Story

1. **SV4 confidence gradient** — straightforward extension of H01 to SV4_IDX=3. Expected to be weaker than SV3 but important for completeness.
2. **SV3 GO biology at all layers** — SV3 was characterized in iter_0014 but the full comparison across SV2/SV3/SV4 using the same protocol has not been done.
3. **Axis-dominant pair biological labels** — the 3 clusters (SV2=1547, SV3=850, SV4=695) have different biological contents; showing this mechanistically explains the orthogonality.
4. **Cross-model validation** — Geneformer SV2 copole is the single highest-upside experiment remaining.
5. **Continuous confidence score vs. copole z** — the quintile approach throws away within-quintile variance; a direct regression (per-pair STRING score → copole z) would be cleaner.

## Stale Directions
- Single-gene permutation tests in module_structure family: retired (previously confirmed)
- TRRUST TF→target directional (repression anti-pole): retired (previously confirmed)
- GO BP enrichment in SV2 poles: retired (iter_0011 H03 null result)
- Bootstrap CI for act/rep differential: underpowered, no rescue value

## Directions to Open

1. **SV4/SV5 full confidence gradient** — completes the axis-by-axis gradient characterization.
2. **Axis-dominant GO composition** — would mechanistically explain axis orthogonality.
3. **Geneformer cross-model** — conservation test; highest-upside.
4. **Direct STRING score regression** — cleaner than quintile binning.
5. **Activation vs. repression geometry in SV3/SV4** — signed regulation may not be SV2-exclusive.
6. **Layer-resolved copole rate trajectory** — per-axis layer profiles could reveal when each axis "activates."
7. **PH (persistent homology) on multi-axis subspace** — test topological structure in 2D/3D SV2+SV3+SV4 projection.
8. **SV2–SV3 2D copole** — pairs co-localized on both axes simultaneously; could be more specific than per-axis.
9. **Network motif geometry** — triangles/cliques in STRING: do triangle members concentrate in the same SV pole more than pairs?
10. **Gene module assignment stability across layers** — does axis-dominant assignment flip across layers?
