# Executor Next Steps — iter_0043

## Retire / Negative

1. **PRDM1/IRF4 anti-convergence** (H01): Negative. PRDM1 converges with rho=-1.000. GC repressor doesn't diverge geometrically. Anti-convergence hypothesis retired for this circuit.
2. **BCL6 divergence** (H02): Negative. BCL6 converges (rho=-0.993). Rank anomaly at L2 is a secondary curiosity.

## New Positive

3. **Effector CD8 ID compression vs pan-T expansion** (H03): Strong dissociation at L7. Effector CD8 compresses (ΔID=-45.8), general T-cell expands (ΔID=+26.2). New hypothesis family extension needed.

## Recommended next hypotheses (iter_0044)

### Priority 1 (HIGH): T-cell functional compression extended panel
- Test naive T markers (CCR7, LEF1, TCF7, SELL) — prediction: expand like general_t
- Test NK markers (NCAM1, KLRD1, FCGR3A, GNLY) — prediction: compress like effector CD8
- Test gamma-delta T (TRGC1, TRGC2, TRDV2) if in vocab
- Generates a full "T-cell ID landscape" figure with clear functional interpretation

### Priority 2 (HIGH): GC circuit unity implications — cross-lineage comparison
- The entire GC circuit converges toward B-cell centroid at rho≈-1.
- Test if myeloid circuit genes (SPI1/PU.1, CEBPA, CEBPB, IRF8) do the same within myeloid context
- Test if T-cell circuit genes (TBX21/T-bet, GATA3, RORC) converge within T-cell context
- If true: this is a **universal attractor property** — regulatory circuits always consolidate near their lineage centroid

### Priority 3 (MEDIUM): BCL6 rank anomaly at L2 cross-correlation
- BCL6 rank jumps from 654 (L0) → 1311 (L2) while metabolic cluster is densest at L0-L4
- Test correlation: BCL6 rank trajectory vs BCL6 cluster size (from h02_metabolic_isolation.json)
- Hypothesis: high metabolic cluster density at early layers causes BCL6 to be ranked as "far" from B-cell centroid in relative terms even as absolute distance decreases

### Priority 4 (NEW FAMILY): Null sensitivity — B-cell marker set label shuffle
- Shuffle which 4 genes are called "B-cell centroid" (random draw from vocab)
- Recompute rho for GC circuit members
- If rho still ≈ -1.0 for random centroid, the convergence is trivial (all genes converge everywhere)
- This is a **critical control** for the circuit-unity finding

## Retired directions (cumulative)
- Cross-model alignment (vocab mismatch, iter_0040)
- GC-plasma subspace angles (iter_0040)
- Chromosomal proximity (iter_0040)
- PC2/PC3 axes (iter_0040)
- PRDM1/IRF4 anti-convergence (iter_0043: negative)
- BCL6 divergence from B-cell centroid (iter_0043: negative)
