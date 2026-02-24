# Executor Next Steps — iter_0014

## Top priorities for iter_0015

### P1: SV3 biological axis — drug/perturbation anchor (new_family, high-value)
H01 established SV3 encodes kinase/signaling (early) → DNA binding (late) vs stress/immune axes. Next: Map SV3 pole genes to perturbation datasets (CRISPRi, drug sensitivity). Test if SV3 projection correlates with known kinase inhibitor sensitivity or inflammatory perturbation response. Could be the bridge to clinical relevance.

### P2: SV3 community detection — do kinase and immune communities cluster separately? (new_method)
SV3 shows polarization by biology. Test: compute kNN graph on SV3 projections alone (or SV2+SV3 joint), run Louvain/modularity. Do identified communities match TRRUST/STRING subnetworks? Compare community membership to GO annotations from H01.

### P3: Joint SV2+SV3 geometry — 2D PPI embedding analysis (new_method)
SV2 and SV3 are both significant PPI axes. Plot 209 genes in SV2-SV3 plane for a key layer (e.g., L3 where SV3 peaks). Annotate by GO category. Test if PPI edges are shorter in this 2D space than in random 2D projections. This produces a direct visualization of the PPI geometry claim.

### P4: Non-hub cross-validation at score>=0.5 threshold (refinement of H03)
H03 used score>=0.4. Test score>=0.5 as intermediate to confirm robustness of non-hub finding. Low effort, tightens the control argument.

## Retired directions (do not re-run without rescue rationale)
- Repression anti-pole (iter_0013 H02): 0/12 sig, z=-1.41 mean. Definitive negative.
- Module structure with TRRUST (iter_0012 H01, iter_0013 H02): 2+ negatives/inconclusive.

## Key evidence to cite in paper (from this iteration)
1. SV3 biological axis characterization (H01): 12/12 layers significant, depth-progressive GO transitions.
2. Hub-degree control resolution (H03): non-hub mean z=3.11 (12/12 sig) > hub z=2.60 (9/12 sig). PPI geometry is not a hub artifact.
