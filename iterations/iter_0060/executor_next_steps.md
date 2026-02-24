# Executor Next Steps — iter_0060

## Completed This Iteration
- H01: Signed displacement projection AUROC = 0.563 (max) — directional geometry does NOT improve over scalar distance. Brainstormer's prediction of >0.62 falsified.
- H02: TF family-stratified margin trajectory reveals bHLH (HIF1A) sign-reversal at L8-L11, bZIP consistently target-like; SV5-7 AUROC increases across layers for TRRUST gene subset.
- H03: FFL triangle geometry — t_mean > 0 significantly at L0-L6 (p<0.0001), signal collapses at L8-L11. Only 22 FFLs, insufficient power for strong conclusions.

## Retire Directions
- **Signed displacement projection → AUROC > 0.62**: falsified. Max = 0.563, same as scalar distance. Retire this specific sub-branch.

## Priority Hypotheses for iter_0061

### High Priority
1. **FFL geometry with expanded motif set** (rescue rationale: prior N=22 too small; method change: include ALL TF→TF edges even if TF-TF pair is not in pos_edges source-to-source):
   - The FFL t_mean signal at L0-L6 is promising but underpowered
   - Expand by: using any gene where TRRUST TF-gene-TF chain can be inferred from neg_edges too
   - Or: use STRING/co-expression to find more TF→TF edges, then TRRUST for TF→target
   - Target: N≥50 FFLs for adequate power

2. **HIF1A single-gene deep-layer geometry** (new_family: biological_anchor):
   - HIF1A shows dramatic flip from target-like (L0-L7) to TF-like (L8-L11) in SV5-7
   - Test: does HIF1A's neighborhood in the SV5-7 embedding change composition (who are its nearest neighbors)?
   - Hypothesis: at L8-L11, HIF1A's nearest neighbors shift from bZIP-like target genes to HIF1A-regulated targets (VEGFA, EPO, etc.)

3. **Within-layer seed AUROC vs eff_rank** (carried over from iter_0059 brainstormer suggestion):
   - Partial_r showed eff_rank is a layer confound, but what about WITHIN single layers?
   - For each of 12 layers, compute seed-to-seed AUROC differences (3 seeds) and eff_rank differences
   - N=3 per layer is small, but test across all 12 layers jointly (N=36 layer×seed pairs)

4. **SV2-4 triangle geometry vs SV5-7** (refinement of H03):
   - Current FFL test uses SV5-7; try SV2-4 and SV2-7 combined
   - If SV2-4 gives stronger t_mean signal, it would identify the geometrically meaningful subspace for regulatory circuit encoding

### Medium Priority
5. **TF family enrichment in FFL motifs** (analytical, quick):
   - Are the 22 FFLs enriched for specific TF families?
   - Is STAT4→TBX21 or RUNX1→JUN known regulatory circuit? Annotate against KEGG/literature.

6. **Cross-layer t_mean trajectory for FFL** (new_method):
   - The FFL t_mean drops sharply at L8. Does this correlate with bHLH flip?
   - Are the same FFLs driving the t_mean collapse (decompose by triplet identity)?

## Critical Open Questions
1. Why does the FFL geometric ordering (t_mean>0) collapse at L8 exactly when bHLH AUROC flips?
2. Is the increasing SV5-7 AUROC across layers (L0: 0.570 → L11: 0.739) a real trend or an artifact of only using TRRUST-subset genes?
3. What drives the persistent boundary geometry (mean AUROC ~0.73-0.76 at L2-L3 in prior work)? We know it's not hub TFs (iter_0059) and not eff_rank (iter_0059), and not directional geometry per se (iter_0060).
