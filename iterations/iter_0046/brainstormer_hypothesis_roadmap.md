# Brainstormer Hypothesis Roadmap — iter_0046

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| Lineage centroid orthogonality (cosine-based) | SV1 at 93.4% variance flattens all cosine similarities to ~0.99; z < 0.3 across all layers | `retire_now` |
| Any cosine-centroid test without SV1 subtraction | Same contamination; will always appear null | `retire_now` |
| GO BP enrichment in raw SV2 poles (lung, iter_0011 approach) | Explored; move to testing whether residual-SV1 subspace carries independent GO structure | `deprioritize` |
| Cross-model feature-effect alignment via low-dim summaries | Inconclusive since iter_0003; only revisit with matched residual tensors from both models | `deprioritize` |

---

## New Hypothesis Portfolio

### FAMILY: Dominant Direction (SV1) Biology

**H-A: SV1 biological identity**
- Hypothesis: The SV1 direction (93.4% of variance at L11) corresponds to a known biological axis — expression level, gene length, or GC content — rather than a functional regulatory axis.
- Test: Project all nonzero gene embeddings onto SV1 at L11; correlate SV1 score with: (1) mean expression level from scRNA reference, (2) gene length, (3) number of GO annotations, (4) TRRUST out-degree. Spearman ρ and Fisher-exact enrichment for top/bottom quartiles.
- Expected signal if true: |ρ| > 0.4 with expression level or annotation density; if false (SV1 is functionally structured), significant GO/pathway enrichment for top vs bottom SV1 quartile.
- Null: permutation of gene-SV1 score labels; random gene sets matched for size.
- Value: HIGH | Cost: LOW

**H-B: Residual-SV1 subspace ID and biology**
- Hypothesis: After projecting out SV1, the remaining ~18 local nonlinear dimensions are biologically structured (carry GO enrichment and regulatory co-membership signal).
- Test: At each layer, subtract the SV1 component from all gene embeddings (project out rank-1 direction). Re-run TwoNN ID estimate on residuals. Re-run TRRUST activation co-pole test on SV2-of-residuals. Re-run GO compartment enrichment.
- Expected signal: Residual ID is stable across layers (no compression after removing SV1); biological enrichment persists in residual SV2 but at lower effect size.
- Null: Feature-shuffled residuals at same rank-1 projection.
- Value: HIGH | Cost: MEDIUM

**H-C: SV1 trajectory across layers**
- Hypothesis: SV1 direction changes its biological character across layers — early layers encode expression noise, late layers encode a specific functional axis.
- Test: At each of 12 layers, compute SV1 projection for all nonzero genes; run GO enrichment on top/bottom quartile. Track which GO terms enter/exit significance across layers.
- Expected signal: GO terms enriched in early SV1 (e.g., highly-expressed housekeeping genes) differ from late SV1 terms.
- Null: Random gene sets; annotation-density controlled.
- Value: HIGH | Cost: MEDIUM

---

### FAMILY: Zero-Norm Sparsification

**H-D: Zero-norm gene identity and biology**
- Hypothesis: The 2902 zero-norm genes at L11 are biologically distinguishable from the 2039 nonzero genes — they are systematically less connected, less annotated, or from specific functional categories.
- Test: Compare zero vs nonzero gene sets at L11: GO enrichment (Fisher exact), STRING degree distribution, TRRUST in/out-degree, median expression in immune reference. Test across layers 0/3/6/9/11.
- Expected signal: Zero genes are enriched for low-annotation, low-connectivity categories; OR they correspond to tissue-specific genes not relevant to the immune training context.
- Null: Random split of all 4941 genes into same-sized sets.
- Value: HIGH | Cost: LOW

**H-E: Zero-norm emergence tracking**
- Hypothesis: Genes go to zero in a layer-coordinated manner within biological pathways — members of the same GO term or STRING cluster go to zero at the same layer transition.
- Test: For each gene, find the first layer where norm drops below threshold (0.01). Compute pairwise concordance of zero-onset layers for STRING-connected pairs vs random pairs. Test with permutation null.
- Expected signal: PPI partners have more concordant zero-onset layers than random pairs (z > 2).
- Null: Same-degree random gene pairs; gene-label permutation.
- Value: MEDIUM | Cost: LOW

**H-F: Zero-fraction drives ID compression**
- Hypothesis: The monotone ID compression is partly artifactual — as more genes go to zero, TwoNN samples fewer true points and the effective ID drops mechanically.
- Test: Re-run TwoNN ID estimation at each layer restricted to nonzero-norm genes only. Compare ID curve (nonzero-only) to full curve. If ID still compresses monotonically among nonzero genes, compression is structural not artifactual.
- Expected signal: If hypothesis true, nonzero-restricted ID is flat; if false, nonzero-restricted ID still compresses from ~33 to ~19-ish.
- Null: No explicit null needed — this is a diagnostic decomposition.
- Value: HIGH | Cost: LOW (critical validity check)

---

### FAMILY: Cross-Tissue and Cross-Model Generalization

**H-G: ID compression in cycle1 (tissue generalization)**
- Hypothesis: TwoNN ID compression from ~33 to ~19 is universal across tissue contexts, not specific to cycle4_immune.
- Test: Run TwoNN ID at 12 layers on cycle1 embeddings (different tissue). Compare shape and magnitude of compression curve to cycle4 results. Bootstrap CI at L0 and L11.
- Expected signal: Same monotone compression, L11 ID in range 15–22, non-overlapping CI with L0.
- Null: Feature-shuffled embeddings at each layer.
- Value: MEDIUM | Cost: LOW

**H-H: SVD spectral collapse cross-tissue**
- Hypothesis: The 14× effective rank collapse is specific to the immune context, where one dominant lineage axis dominates; other tissue contexts may show less extreme collapse.
- Test: Compute SVD effective rank at 12 layers for cycle1. Compare L11 eff_rank and top-1 variance fraction to cycle4 values (1.64, 93.4%).
- Expected signal: Different tissue shows eff_rank at L11 > 2 and top-1 var < 85%, suggesting immune context has unusually strong dominant direction.
- Null: Feature-shuffled cycle1 embeddings at L11.
- Value: MEDIUM | Cost: LOW

**H-I: Signed regulatory geometry cross-tissue (immune context)**
- Hypothesis: The activation-specific TRRUST co-pole signal (iter_0011, found in lung/scGPT) replicates in the immune embedding (cycle4_immune_main) on SV2 of the immune model.
- Test: Extract SV2 from cycle4_immune_main at all 12 layers. Compute co-pole rate for TRRUST activation vs repression pairs (using immune-expressed genes). Permutation null.
- Expected signal: Activation pairs co-pole at z > 2 in ≥ 8/12 layers; repression pairs do not.
- Null: Gene-label shuffle null (500 replicates).
- Value: HIGH | Cost: MEDIUM

---

### FAMILY: Topology and Manifold Structure

**H-J: Persistent homology on nonzero-gene subspace**
- Hypothesis: Prior H1 persistence results were contaminated by zero-norm gene points collapsing topology. Restricting to nonzero genes at each layer will show cleaner/stronger H1 persistence signal.
- Test: At each of 12 layers, run Ripser H1 persistence on nonzero-norm genes only (n ≤ 2039). Compare mean H1 lifetime to feature-shuffled null. Compare to full-set results from iter_0003.
- Expected signal: Cleaner separation between real and null; possibly different layer-dependence of H1 signal.
- Null: Feature-shuffled embeddings restricted to same nonzero-gene set.
- Value: MEDIUM | Cost: MEDIUM

**H-K: ID compression phase transition layer**
- Hypothesis: ID compression is not smooth — there is a specific layer (L6–L8 based on SVD data) where the compression rate accelerates, representing a functional processing transition.
- Test: Fit piecewise-linear model to TwoNN ID-vs-layer curve across 12 layers (all 3 seeds). Detect breakpoint via CUSUM or RSS minimization. Compare to SVD eff_rank breakpoint (visually, L8–L10 shows rapid collapse). Test whether breakpoint layers coincide.
- Expected signal: Both TwoNN and SVD show same breakpoint layer (±1 layer).
- Null: Smooth monotone fit (linear or exponential); bootstrap CI on breakpoint.
- Value: MEDIUM | Cost: LOW

**H-L: Ollivier-Ricci curvature across layers**
- Hypothesis: The kNN graph curvature (Ollivier-Ricci) changes systematically across layers, with early layers showing more positive curvature (clustered structure) and late layers showing negative curvature (tree-like or radial structure consistent with dominant axis).
- Test: Compute Ollivier-Ricci curvature for kNN graph (k=10) of gene embeddings at each layer. Mean edge curvature ± null (feature-shuffled). Use GraphRicciCurvature package.
- Expected signal: Curvature decreases (more negative) as layers increase, reflecting transition from locally clustered to globally radial structure.
- Null: Feature-shuffled embeddings at each layer.
- Value: MEDIUM | Cost: MEDIUM

**H-M: STRING PPI co-pole in immune SV2**
- Hypothesis: STRING PPI partners (high-confidence ≥ 700) co-localize in SV2 poles in cycle4_immune_main at layer 11, extending the TRRUST regulatory co-pole result to physical interactions.
- Test: Extract top-500 SV2 and bottom-500 SV2 genes at L11. Test STRING edge co-pole rate vs permutation null. Stratify by STRING channel (coexpression vs experimental).
- Expected signal: Experimental PPI partners show higher co-pole rate than coexpression partners; both exceed null.
- Null: Degree-matched random gene pairs.
- Value: HIGH | Cost: LOW

**H-N: SV1-corrected centroid test**
- Hypothesis: After subtracting SV1 from all embeddings, lineage centroids (B/T/Myeloid) become distinguishable from random gene set centroids via cosine similarity.
- Test: Project out SV1 from all embeddings at each layer; recompute centroid cosine similarities between lineage centroids; recompute null z-scores.
- Expected signal: z-scores rise above 2.0 at ≥ 1 layer after SV1 removal; lineage centroids become more separated.
- Null: Same random gene set baselines, also SV1-corrected.
- Value: MEDIUM | Cost: LOW (direct recovery of H03)

---

## Top 3 for Immediate Execution

### #1: High-probability discovery candidate
**H-F + H-D combined: Zero-norm sparsification analysis**

Rationale: Two cheap tests that directly address whether the core ID compression finding is valid (H-F) and reveal a completely unexplored structural feature (H-D). 2902/4941 genes being zero at L11 is a dramatic structural property that could reframe all prior geometry results. If H-F shows compression survives restriction to nonzero genes, the finding is robustly validated. If H-D shows zero genes are biologically distinct, it's a novel discovery.

Implementation: Single script: (1) track zero-norm fraction per layer; (2) characterize zero vs nonzero gene sets at L11 by GO/STRING; (3) re-run TwoNN on nonzero-only genes.

Expected runtime: ~20 min. Artifacts: zero_norm_analysis.json.

---

### #2: High-risk / high-reward candidate
**H-A: SV1 biological identity**

Rationale: SV1 explains 93.4% of embedding variance at L11. If it encodes expression level or annotation density, it's a confound that partially explains prior results. If it encodes a genuine functional axis (e.g., immune activation state), it's the most interpretable finding of the project. This test is decisive either way.

Implementation: Correlate SV1 gene scores with expression level, GO annotation count, TRRUST degree. Run GO enrichment on top/bottom SV1 quartiles. Need: reference expression data for immune context (cycle4 training data gene means, or public PBMC expression).

Expected runtime: ~30 min (if expression reference available), ~1h if need to load scRNA data. Artifacts: sv1_identity.json.

---

### #3: Cheap broad-screen candidate
**H-K + H-G: ID compression curve shape + cycle1 generalization**

Rationale: Two cheap experiments: (1) fit breakpoint model to existing TwoNN data (requires no new GPU runs, just analysis of h01 artifact); (2) re-run TwoNN on cycle1 embeddings (same script as H01, different data path). Together they characterize the compression curve shape and test tissue generalization in a single executor run.

Implementation: Load h01_id_crossseed_bootstrap.json; fit piecewise-linear to ID-vs-layer; then run TwoNN on cycle1 layer embeddings. Artifacts: id_breakpoint.json, cycle1_id.json.

Expected runtime: ~15 min for breakpoint, ~20 min for cycle1 TwoNN.
