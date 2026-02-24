# Brainstormer Structured Feedback — iter_0024

## Gate Status
`passed_min_research_gate: true` — all checks passed, paper updated, artifacts present.

## Result Assessment

### H01: Dorothea Confidence Stratification — STRONG POSITIVE
- AUROC=0.671 at L8, 12/12 layers significant, 0/500 permutations exceed real signal.
- Key interpretive upgrade: Dorothea A/B signal is independent of STRING (different database, regulatory directionality not PPI). This is a 4th independent biological anchor for scGPT embedding geometry.
- Next logical step: partial out STRING score to confirm independence, then claim regulatory confidence as a geometric dimension.
- This result is publication-ready pending the joint-predictor test.

### H02: TF Activation Hub Centrality — NEGATIVE (retire)
- Spearman(act_degree, proximity) = -0.047, p=0.74. 0/12 layers sig.
- The pairwise activation proximity signal (iter_0023 H03) does not generalize to hub centrality.
- Interpretation: TF regulatory function is encoded locally (pair-level), not globally (degree-centrality). No rescue attempt warranted — this is a genuine null.

### H03: GO Ontology Comparison (CC > MF > BP) — STRONG POSITIVE
- CC Spearman=0.106 at L8 (peak 0.124 at L5), 12/12 layers all three ontologies sig.
- Key finding: subcellular localization (GO CC) is the single best GO predictor of scGPT proximity, stronger than molecular function or biological process.
- BP shows layer-deepening (0.050 → 0.083) while CC is early-dominant; this dissociation hints at different encoding timescales in the transformer stack.
- This is directly actionable: GO CC should be the primary ontology predictor in all multi-predictor models going forward.

## Portfolio Assessment

### Confirmed Multi-Anchor Picture (as of iter_0024)
The scGPT embedding geometry is now anchored by:
1. STRING PPI proximity (AUROC=0.614, iter_0022) — physical interaction
2. TRRUST-exclusive TF-target proximity (AUROC=0.573, activation: 0.640, iter_0023) — regulatory co-program
3. GO BP Jaccard proximity (Spearman=-0.077, iter_0023) — functional process
4. GO CC Jaccard proximity (Spearman=0.106, iter_0024) — subcellular localization ← new
5. Dorothea A/B confidence pairs (AUROC=0.671, iter_0024) — regulatory confidence ← new
6. Cell-type marker clustering (AUROC=0.851, iter_0022/0023) — identity structure

### Key Open Questions
1. Are all anchors independent after controlling for each other? (multi-predictor regression pending)
2. Which specific GO CC terms (nucleus vs cytoplasm vs membrane vs mitochondria) dominate?
3. Does the BP layer-deepening pattern reflect a mechanistically different encoding than CC early dominance?
4. Can we extend activation proximity asymmetry to Dorothea (which lacks activation/repression directionality in the current data)?
5. What is the topological structure of the multi-anchor distance space (manifold dimensionality per anchor)?

## Directions to Retire/Deprioritize
- **TF activation hub centrality** (H02 this iter): retire now — clear null, no rescue signal
- **Pure SV pole tests without biological anchoring**: largely exhausted through iter_0013-0017; avoid revisiting unless new axis found
- **GO BP enrichment in SV poles**: already retired (iter_0011 H03) — do not revisit
- **Repression anti-pole hypothesis**: retired iter_0013 H02 — confirmed null

## Highest-Value Unexplored Territory
1. Multi-predictor joint model (consolidates all 5+ anchors)
2. GO CC compartment decomposition (nucleus/cytoplasm/membrane/mito)
3. Single-TF target cluster analysis (JUN, TP53 as archetypes)
4. Chromosomal proximity negative control (validate functional specificity)
5. Geneformer GO CC alignment (does cross-model consistency hold for CC?)
6. Topological persistence of GO CC clusters vs STRING clusters (comparative PH)
7. Cell-type marker expansion to epithelial/myeloid with GO CC integration
