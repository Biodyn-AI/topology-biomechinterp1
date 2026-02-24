# Brainstormer Hypothesis Roadmap — iter_0012 → iter_0013+

---

## Retire / Deprioritize

| Direction | Reason | Status |
|---|---|---|
| Bootstrap CI on act-rep differential | Underpowered at n=64; directional qualitative result already captured by per-set empirical p-values | `retire_now` |
| GO BP enrichment in SV2 poles | 0/591 terms pass emp_p<0.05; null distribution covers observed | `retire_now` |
| Feature-column shuffle null | Degenerate for L2 norms (confirmed iter_0007) | `retire_now` (done) |
| Drift TF enrichment as primary result | Fails gene-label-shuffle specificity | `retire_now` (done) |
| Higher-order PH (H2 Betti) on full embeddings | Ripser on 4803×512 is computationally prohibitive; topological signal in H1 already established | `deprioritize` |

---

## New Hypothesis Portfolio

### H-NEW-01 — STRING hub degree confounder control
**Hypothesis**: The STRING PPI co-pole enrichment (z=3.3–6.5) survives after controlling for hub degree, confirming it reflects geometric organization rather than hub-clustering artifacts.
**Test**: For each of 12 layers, compute degree-matched null: for each gene, select random replacement gene with degree within ±2 (within the 209-gene induced subgraph); run co-pole enrichment on degree-matched null (N=500). Compare observed co-pole rate to degree-matched null vs. uniform null. Also: exclude top-10 highest-degree genes and rerun.
**Expected signal if true**: Co-pole enrichment z-score remains significant (z>2) after degree control; removing hubs does not eliminate the effect.
**Null/control**: Degree-matched shuffle; hub-removed subset.
**Value**: high | **Cost**: low (uses cached STRING + existing embeddings)

---

### H-NEW-02 — STRING channel decomposition (experimental vs. co-expression channel)
**Hypothesis**: The SV2 PPI co-pole enrichment is driven by experimentally validated physical interactions (STRING experimental channel), not co-expression links.
**Test**: Parse STRING v12.0 per-channel scores from cached `string_ppi_named_genes.json`. Filter separately: (a) experimental ≥ 0.4, (b) coexpression ≥ 0.4, (c) database ≥ 0.4. Run co-pole test on each subset across 12 layers. Compare z-score profiles.
**Expected signal if true**: Experimental channel shows high z-scores; co-expression channel shows reduced/absent enrichment.
**Null/control**: Uniform gene-label shuffle (same as H03).
**Value**: high | **Cost**: low (cached data, new filtering step only)

---

### H-NEW-03 — SV axis specificity of PPI co-pole (SV1 vs SV2 vs SV3)
**Hypothesis**: STRING PPI co-pole enrichment is SV2-specific; SV1 and SV3 show significantly lower enrichment.
**Test**: Run STRING co-pole test (score≥0.7, 1022 pairs, N=500 shuffles) using poles defined on SV1 and SV3 projections (top/bottom K=52 each), across 12 layers. Compare z-score profiles for SV1 vs SV2 vs SV3.
**Expected signal if true**: SV2 z-scores (3.3–6.5) significantly exceed SV1 and SV3 at most layers; SV2 is the primary PPI axis.
**Null/control**: Gene-label shuffle; compare same-K poles on SV1/SV3.
**Value**: high | **Cost**: low (reuses cached STRING + existing SVD projections)

---

### H-NEW-04 — Repression anti-pole (cross-pole) geometry
**Hypothesis**: Repression TF-target pairs are enriched in *opposite* SV2 poles (one gene top-K, one gene bottom-K), whereas activation pairs co-localize in the same pole.
**Test**: Compute cross-pole rate (gene_i in top-52, gene_j in bottom-52 OR vice versa) for activation and repression pairs separately at each of 12 layers vs. gene-label-shuffle null. Also run signed SV2 displacement test: mean(SV2_i × SV2_j) should be positive for activation (same sign), negative for repression (opposite sign).
**Expected signal if true**: Repression cross-pole rate significantly elevated vs. null; signed product test confirms act > 0, rep < 0 at majority of layers.
**Null/control**: Gene-label shuffle null for cross-pole rate; sign test against zero for product statistic.
**Value**: high | **Cost**: low

---

### H-NEW-05 — Continuous signed displacement test for act vs. rep
**Hypothesis**: Activation pairs have same-sign SV2 projections (both positive or both negative) more than expected; repression pairs have opposite-sign projections more than expected.
**Test**: For each pair (i,j), compute sign(SV2_i) × sign(SV2_j). Mean over activation pairs (expect +1 for co-polar); mean over repression pairs (expect −1 for anti-polar). Test against null distribution (N=500 shuffles). Also test continuous product SV2_i × SV2_j mean.
**Expected signal if true**: Activation mean product > 0 (p<0.05); repression mean product < 0 (p<0.05); significant act vs. rep difference at majority of layers.
**Null/control**: Gene-label shuffle.
**Value**: high | **Cost**: low

---

### H-NEW-06 — CORUM protein complex membership in SV2 poles
**Hypothesis**: Genes belonging to the same annotated protein complex (CORUM database) are co-localized in SV2 poles beyond PPI co-pole enrichment.
**Test**: Download CORUM human complexes. For each complex with 3+ named genes, compute co-pole rate (within same SV2 pole) vs. null. Also compare to STRING PPI co-pole at same K; CORUM complex should show higher z-scores if physical complex co-localization is a stronger predictor than generic PPI.
**Expected signal if true**: CORUM co-pole z-scores > STRING PPI z-scores; protein complexes are more tightly localized in SV2 than general PPIs.
**Null/control**: Gene-label shuffle; CORUM random subsets of matched size.
**Value**: high | **Cost**: medium (CORUM download required)

---

### H-NEW-07 — Cross-model SV2 alignment: scGPT vs. Geneformer
**Hypothesis**: The SV2 gene ordering (projection of matched genes) is preserved across scGPT and Geneformer, demonstrating that PPI/regulatory geometry is model-independent.
**Test**: Extract matched gene SV2 projections from scGPT layer-11 and Geneformer final layer for the 209 named genes (or the intersection). Compute Spearman r(scGPT SV2, Geneformer SV2). Run permutation null (N=500 gene-label shuffles). Also test whether Geneformer SV2 alone shows PPI co-pole enrichment.
**Expected signal if true**: Spearman r > 0.5 (empirical p < 0.01); Geneformer SV2 also shows PPI co-pole enrichment. Would be the strongest result across all iterations.
**Null/control**: Gene-label shuffle null on Spearman r.
**Value**: high | **Cost**: high (requires Geneformer embeddings to be loaded/computed)

---

### H-NEW-08 — SV2 as a PPI link predictor (held-out edge prediction)
**Hypothesis**: SV2 inter-gene distance predicts held-out STRING PPI edges (AUROC > 0.6 vs. random edge classifier).
**Test**: Hold out 20% of STRING high-confidence pairs (score≥0.7) as test set. Train threshold on SV2 pairwise distance using 80% training pairs. Compute AUROC on test set vs. matched non-edges. Repeat for SV1 and SV3. Also test multi-axis distance (SV1+SV2+SV3 Euclidean in 3D projection space).
**Expected signal if true**: SV2 AUROC > 0.6 for PPI prediction; SV2 outperforms SV1 and SV3; 3D distance outperforms any single axis.
**Null/control**: Random edge classifier (AUROC=0.5); degree-matched random edge classifier.
**Value**: high | **Cost**: medium

---

### H-NEW-09 — Layer-resolved STRING z-score peak decomposition
**Hypothesis**: The STRING PPI co-pole z-score peak at L1 (z=6.52) and L11 (z=5.45), with trough at L9-10, reflects distinct biological processes at early vs. late layers.
**Test**: At L1 (peak) and L9 (trough), run STRING co-pole enrichment stratified by STRING channel (experimental, co-expression, database) and also by interaction GO category (if available). Compare which interaction types drive the L1 peak vs. which drop at L9.
**Expected signal if true**: L1 peak is dominated by experimental/complex interactions (physical core preserved from embedding); L9 trough correlates with ER transient seen in iter_0009 (layer-specific geometric reorganization temporarily disrupts PPI co-localization).
**Null/control**: Same gene-label shuffle; compare channel-stratified z-scores at L1 vs. L9.
**Value**: medium | **Cost**: low

---

### H-NEW-10 — SV2 co-expression module alignment
**Hypothesis**: SV2 poles align with scRNA-seq co-expression gene modules (Leiden clustering on Pearson correlation) more than expected by chance.
**Test**: Compute Pearson expression correlation matrix for 209 named genes across cells in the dataset. Apply Leiden clustering (resolution ~0.5) to get co-expression modules. For each module, compute what fraction of its genes are in the same SV2 pole. Test against gene-label-shuffle null.
**Expected signal if true**: Leiden co-expression modules have higher same-pole rate than random; specific modules map to top vs. bottom SV2 pole cleanly.
**Null/control**: Gene-label shuffle; module size-matched random gene sets.
**Value**: medium | **Cost**: medium

---

### H-NEW-11 — Multi-axis (SV1+SV2+SV3) 3D proximity for PPI pairs
**Hypothesis**: PPI partners are co-localized in the 3D subspace spanned by (SV1, SV2, SV3) projections more than in any single axis, suggesting the full top-3 spectral subspace is the biologically organized space.
**Test**: For each layer, compute mean Euclidean distance of PPI pairs in (SV1, SV2, SV3) 3D projection. Compare to null (gene-label shuffle). Also compare z-scores: single-axis SV2 vs. 3D distance. Run for TRRUST activation pairs too.
**Expected signal if true**: 3D Euclidean distance shows higher z-scores than SV2 alone; SV1+SV2+SV3 together predict PPI proximity better than any single axis.
**Null/control**: Gene-label shuffle; compare z-scores across axis configurations.
**Value**: medium | **Cost**: low

---

### H-NEW-12 — Persistent homology on SV2-projected regulatory subgraph
**Hypothesis**: The SV2 1D projection of the 209-gene regulatory subgraph has higher H0 persistence (fewer disconnected components at larger distance) than random same-size subgraphs, indicating connected regulatory clusters in 1D SV2 space.
**Test**: Project 209 named genes onto SV2. Compute Vietoris-Rips H0 persistence on the 1D point cloud (SV2 values). Compare to N=200 gene-label-shuffle nulls. Also test H1 for loops in SV2 × SV3 2D space.
**Expected signal if true**: H0 persistence shows fewer long gaps (more compact clusters) than null; clustering structure in 1D SV2 reflects biological modules.
**Null/control**: Gene-label shuffle.
**Value**: medium | **Cost**: low (1D persistence is fast)

---

### H-NEW-13 — TF family-resolved SV2 co-pole enrichment
**Hypothesis**: Certain TF families (IRF, KLF, RUNX, bHLH) show stronger within-family SV2 co-pole enrichment than others, revealing TF family specificity of geometric co-localization.
**Test**: Annotate 209 named genes with TF family (TFclass database or manual HGNC). For each family with ≥4 members in the set, compute within-family co-pole rate vs. null. Report family-resolved z-scores.
**Expected signal if true**: IRF/RUNX families (known to interact physically) show high z-scores; KLF family (diverse targets) shows lower z-scores.
**Null/control**: Gene-label shuffle; family-size-matched random gene sets.
**Value**: medium | **Cost**: low

---

### H-NEW-14 — SV2 pole composition: annotate with known co-activator/co-repressor complexes
**Hypothesis**: The top SV2 pole is enriched for known co-activator complex genes; the bottom pole is enriched for co-repressor complex genes.
**Test**: Collect known co-activator complexes (e.g., SAGA, Mediator) and co-repressor complexes (e.g., NuRD, CoREST, Sin3A) gene members. Fisher exact test for enrichment in top vs. bottom SV2 pole at layers 7-11.
**Expected signal if true**: Co-activators enriched in one pole, co-repressors in other, separating two poles functionally by their regulatory directionality.
**Null/control**: Gene-label shuffle; complex-member count-matched random sets.
**Value**: medium | **Cost**: low

---

## Top 3 for Immediate Execution

### #1 — High-probability discovery: H-NEW-01 + H-NEW-02 combined (STRING confound resolution)
**Why**: The STRING PPI co-pole finding (H03, z=3.3–6.5) is the most compelling result but needs confound control before it can be stated cleanly. Hub degree control and channel decomposition both use already-cached data (`string_ppi_named_genes.json`). If both controls confirm the result, H03 advances to primary claim status. If hub degree explains the enrichment, that changes the interpretation substantially. This must be done before building further on H03.
**Deliverable**: Degree-matched null z-scores vs. uniform null z-scores at 12 layers; channel-stratified co-pole z-scores (experimental vs. co-expression vs. database).

### #2 — High-risk/high-reward: H-NEW-04 + H-NEW-05 combined (repression anti-pole signed geometry)
**Why**: Iter_0011 found activation co-poles (12/12 layers) but repression non-co-poles (1/12). The natural follow-up is whether repression pairs are *anti*-co-polar (opposite poles). If true, this establishes that SV2 encodes the *sign* of regulation as spatial direction, not just co-regulation magnitude. This would be the strongest mechanistic claim in the project. If negative (repression is neither co-polar nor anti-polar), it means repression pairs are simply dispersed in SV2. Both outcomes are informative. The test is cheap (re-uses all existing data).
**Deliverable**: Cross-pole rate and signed product statistic for repression pairs across 12 layers with null controls.

### #3 — Cheap broad-screen: H-NEW-03 (SV axis specificity of PPI co-pole)
**Why**: If STRING PPI co-pole is SV2-specific (SV1 and SV3 show z<2), this strongly supports SV2 as the dedicated "interaction geometry axis" and gives a clean mechanistic claim. If SV1 also shows PPI co-pole (which it might, since SV1 explains 93% of variance), the claim is weaker. Takes ~20 lines of code on existing data. Results directly interpret the scope and uniqueness of the SV2 PPI geometry finding.
**Deliverable**: SV1 vs. SV2 vs. SV3 z-score profiles for STRING PPI co-pole across 12 layers.
