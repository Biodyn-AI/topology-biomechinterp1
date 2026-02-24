# Brainstormer Structured Feedback: iter_0044

## Gate Status
`passed_min_research_gate: true` — three experiments ran, all artifacts present. Research gate passed.

## Outcome Summary
Three hypotheses tested; all negative or inconclusive. This is the third consecutive all-negative iteration (iter_0042 H01/H02 BCL6, iter_0043 H01 BCL6, iter_0044 H01/H03). The BCL6 divergence family and cross-lineage ID compression law are now definitively retired.

## What iter_0044 Actually Told Us

**H01 (cross-lineage ID compression):** The 5-point TwoNN is geometrically measuring a 5-simplex, not local ID within the full gene manifold. The methodological mismatch with iter_0043 is now documented. This is a useful negative because it clarifies what iter_0043's ΔID=-45.8 was actually measuring: TwoNN over ~4941 genes with local neighborhood extraction. Future ID work must use the full manifold.

**H02 (TCR circuit attractor):** The rank improvement (L0=92 → L11=41, Δ=-50.7) is directionally interesting but failed null significance (p=0.275). The CD247 confound (in both circuit and centroid) inflates precision@20 to 0.67 at L0 before any convergence has occurred. With CD247 removed from the centroid and only CD28+LCK as circuit genes, this could be a real signal. Not retired — deserves one clean retest.

**H03 (BCL6 metabolic isolation):** BCL6 B-cell ranks and BCL6 metabolic ranks co-move (rho=0.5). No directed metabolic attractor. This is the final retirement of the BCL6 story (also retired in iter_0042 and iter_0043).

## Methodological Insight Worth Recording
The iter_0044 TwoNN discrepancy (ΔID ~0.5–3 vs iter_0043's 45.8) reveals a critical distinction:
- **5-point TwoNN in isolation**: measures pairwise distance ratios in a 5-point set; meaningless for ID estimation
- **Full-manifold TwoNN with local subset**: measures how many dimensions are needed to explain local neighborhoods in the ~4941-gene space

All future ID work must use the full manifold. The B-cell ID reduction finding (iter_0042, rho=-0.951) is robust — it was computed correctly.

## Positive Result Inventory (Unchanged)
1. B-cell PC1 negative pole, stable all 12 layers (iter_0033)
2. B-cell kNN community, OR=10–16 from L2 onward (iter_0034)
3. GC-TF attractor onset Layer 3, BATF/BACH2/PAX5 (iter_0040)
4. GC attractor is lineage-unique — T-cell/myeloid TFs are pre-wired from L0, GC-TFs show delayed convergence (iter_0041)
5. B-cell intrinsic dimensionality reduction: 8.16→5.09, rho=-0.951, breakpoint at L3 (iter_0042)
6. GC circuit unity: PRDM1 (repressor) converges to B-cell centroid at same rate as activators BATF/BACH2 (iter_0043)
7. H1 persistence elevated vs null in scGPT lung (iter_0003)
8. kNN clustering coefficient elevated vs null, all layers (iter_0002)

## Direction Assessment

| Direction | Status | Reason |
|---|---|---|
| BCL6 divergence / metabolic isolation | **RETIRE** | Three negatives across iter_0042–0044 |
| Cross-lineage ID compression (5-point TwoNN) | **RETIRE** | Methodologically invalid; wrong computation |
| TCR circuit attractor | **Rescue once** | Has directional signal; CD247 confound must be fixed |
| GC circuit unity extension (other lineages) | **New, promising** | Positive result in B-cell; test T-cell/myeloid |
| Plasma cell divergence from GC attractor | **New, high probability** | PRDM1 is now in GC neighborhood; when does it diverge? |
| Full-manifold TwoNN ID trajectory | **Priority** | Foundational; validates iter_0042/0043 and enables new ID tests |
| kNN graph community structure (4941 genes) | **Untested, positive prior** | iter_0036 positive on 195-gene set; scale up |
| Persistent homology on immune data | **Untested in this dataset** | Only done on lung (iter_0003); import to cycle4_immune |
| Ollivier-Ricci curvature transition at L3 | **New, high risk/reward** | Layer 3 is established structural transition; curvature could formalize it |
| Cross-model B-cell axis (scGPT → Geneformer) | **Untested, high value** | Is the B-cell PC1 axis model-invariant? |
