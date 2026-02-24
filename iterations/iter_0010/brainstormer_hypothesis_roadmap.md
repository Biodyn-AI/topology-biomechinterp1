# Brainstormer Hypothesis Roadmap — iter_0010 → iter_0011

---

## Retire / Deprioritize

| Direction | Action | Reason |
|-----------|--------|--------|
| Persistent homology rewiring-null | **RETIRE NOW** | 4 consecutive negatives, no new angle available |
| GO BP cosine clustering | **RETIRE NOW** | Structurally underpowered, SVD pole approach supersedes it |
| kNN transitivity | **RETIRE NOW** | Negative result at maximum statistical separation |
| PCA-only TRRUST co-target clustering | **DEPRIORITIZE** | SVD pole approach is a strict improvement; revisit only if SVD approach fails |
| SV4+ without biological anchor hypothesis | **DEFER** | Only add if SV1–SV3 biology is fully characterized |

---

## New Hypothesis Portfolio

### H-A: TRRUST layer profile + edge-type stratification
**Hypothesis**: TF-target co-pole signal peaks at a specific intermediate layer (not L11) and activation edges show stronger co-pole enrichment than repression edges, implying signed regulatory geometry.
**Test**: Compute SV1/SV2 co-pole rates for TRRUST pairs across all 12 layers. Stratify 326 pairs by edge type (Activation/Repression). Fisher test within each stratum. Null: N=1000 random pairs.
**Expected signal**: Co-pole rate curve peaks at L8–L11; activation pairs >repression pairs in co-pole rate; emp_p<0.01 for activation.
**Null/control**: Edge-type-stratified random pairs; compare activation vs repression effect sizes.
**Value**: HIGH | **Cost**: LOW (extension of H02, same infrastructure)

### H-B: STRING PPI co-pole test
**Hypothesis**: Physically interacting gene pairs (STRING high-confidence PPI) co-localize in SVD poles more than random pairs, providing convergent evidence from physical (not regulatory) networks.
**Test**: Download/use STRING human PPI at confidence≥700. Filter to 209-gene set. Compute co-pole rates for SV1 and SV2 at L11. Null: N=1000 random pairs.
**Expected signal**: PPI co-pole rate > null, emp_p<0.05 for at least SV2.
**Null/control**: Random same-size pairs; compare STRING vs TRRUST effect sizes.
**Value**: HIGH | **Cost**: LOW-MEDIUM (need STRING file; same test code)

### H-C: SV2 cytoskeletal gene identity + cross-layer pole stability
**Hypothesis**: The 7 cytoskeletal genes driving SV2 OR=20.35 occupy the same pole across all 12 layers with near-perfect stability (Spearman r ≈ 1.0 for their SV2 rank positions).
**Test**: Extract SV2 ranks of all cytoskeletal-annotated genes at each layer. Compute Spearman r between layers. Identify gene names. Visualize rank drift.
**Expected signal**: Cytoskeletal genes clustered in top or bottom K consistently; Spearman r >0.9 between adjacent layers.
**Null/control**: Same Spearman r for randomly selected 7 genes.
**Value**: MEDIUM-HIGH | **Cost**: LOW (analytical only, same data)

### H-D: SV2 L8 decomposition — which compartments drive the variance peak
**Hypothesis**: The SV2 variance peak at L8 (9.9%) is driven by the mitochondrial transient, not by the cytoskeletal signal (which is stable). At L8, mito genes dominate SV2 pole; before/after L8, cytoskeletal genes dominate.
**Test**: At L7, L8, L9: compute SV2 top-K and bottom-K. For each, compute compartment enrichments. Track which compartment shows largest OR shift at L8. Compare mito vs cytoskeleton OR at these three layers.
**Expected signal**: Mito OR peaks at L8; cytoskeletal OR is stable across L7–L9; overlap between mito-SV2 genes and cytoskeletal-SV2 genes is low.
**Null/control**: Same shuffled null (N=200) as H01.
**Value**: HIGH | **Cost**: LOW (subset of H01 data, targeted re-analysis)

### H-E: Pole-flip gene tracking (mitochondrial polarity inversion)
**Hypothesis**: A specific subset of mitochondrial genes switch SV2 poles between L3 and L4 (the observed mito polarity flip), while non-mitochondrial genes show stable SV2 pole membership.
**Test**: For each layer pair (L0→L1, L1→L2, ..., L10→L11), compute SV2 top-K and bottom-K. Track which genes switch pole. Identify peak pole-switching layer. Test whether mito-annotated genes are enriched in switchers.
**Expected signal**: Peak switcher activity at L3→L4 transition; mito genes overrepresented in switchers (OR>>1 by Fisher test).
**Null/control**: Random gene-to-compartment relabeling.
**Value**: HIGH | **Cost**: MEDIUM (requires per-layer gene tracking, ~30 min compute)

### H-F: Reactome pathway co-pole test
**Hypothesis**: Genes sharing a Reactome metabolic or signaling pathway co-localize in SVD poles more than random gene pairs.
**Test**: Download Reactome gene sets (GMT file). For each pathway with ≥5 genes in the 209-gene set, compute SV2 co-pole rate. Null: N=1000 random same-size gene sets.
**Expected signal**: Mitochondrial/ETC pathways show strong co-pole enrichment; general/housekeeping pathways show weaker signal.
**Null/control**: Random gene sets; FDR correction over pathways.
**Value**: HIGH | **Cost**: MEDIUM (need Reactome GMT, new pathway loop)

### H-G: Bootstrap SV2 pole stability (robustness check)
**Hypothesis**: SV2 compartment enrichments are stable under gene-level bootstrap resampling — removing any single gene does not abolish the cytoskeletal signal.
**Test**: 200 bootstrap resamplings of 190/209 genes (leave-10-out). At each, compute SV2 top-K, bottom-K, Fisher OR for cytoskeleton and mito. Distribution of ORs.
**Expected signal**: Cytoskeletal OR ≥10 in >95% of bootstraps; mito OR ≥3 in >80% of bootstraps at early layers.
**Null/control**: Bootstrap distribution width vs gene count (robustness quantification).
**Value**: MEDIUM | **Cost**: LOW-MEDIUM (~15 min, N=200 bootstraps × 12 layers × 2 compartments)

### H-H: Cross-model SV2 alignment: scGPT vs Geneformer
**Hypothesis**: The SV2 axis of scGPT (cytoskeletal axis) aligns with a corresponding axis in Geneformer gene embeddings, demonstrating model-independent biological encoding in sub-leading SVD directions.
**Test**: If Geneformer embeddings available: compute matched-gene SVD for both models. Test Procrustes/CKA alignment between scGPT SV2 and Geneformer SV1–SV4. If not: use iter_0003 Geneformer feature-effect vectors aligned to scGPT SV2 weights.
**Expected signal**: CKA(scGPT SV2, Geneformer SV2) > CKA(scGPT SV2, Geneformer noise axis).
**Null/control**: CKA with permuted Geneformer axes.
**Value**: HIGH | **Cost**: MEDIUM-HIGH (depends on Geneformer data availability)

### H-I: Higher-order SVD structure — SV4/SV5 scan
**Hypothesis**: SV4 and SV5 encode additional separable biological programs (e.g., secretory pathway, lysosomal, or cell-cycle related genes), extending the multi-axis compartment decomposition.
**Test**: Same scan as H01 for SV4 and SV5. 12 layers × 9 compartments × N=200 shuffles.
**Expected signal**: ≥1 compartment with emp_p<0.05 for SV4 or SV5 at any layer.
**Null/control**: Same shuffle null as H01.
**Value**: MEDIUM | **Cost**: LOW (identical infrastructure to H01)

### H-J: GO Biological Process (not just compartment) enrichment in SVD poles
**Hypothesis**: SVD poles encode not just subcellular localization but also functional gene programs (e.g., "oxidative phosphorylation", "actin cytoskeleton organization") as GO BP terms.
**Test**: Repeat H01-style scan but with GO BP gene sets (filter to sets with ≥5 genes in 209-gene pool). Test top-50 most populated BP terms for SV1/SV2 at L11. Fisher + N=200 shuffles.
**Expected signal**: Oxidative phosphorylation and cytoskeletal organization BP terms significant in SV2.
**Null/control**: Same shuffle null; FDR correction.
**Value**: HIGH | **Cost**: MEDIUM (need GO BP gene sets; ~20 min)

### H-K: SV2 projection as a continuous "cytoskeletal score" — correlation with known cytoskeletal regulators
**Hypothesis**: The SV2 projection value of each gene (continuous, not binary pole) correlates with expert-curated cytoskeletal gene importance (e.g., F-actin network genes have most extreme SV2 scores).
**Test**: At L11, compute SV2 projection for all 209 genes. Rank. Cross-reference with CytoSkDb or GO cytoskeleton sub-terms. Spearman r between SV2 rank and "cytoskeletal centrality" (number of cytoskeletal GO annotations per gene).
**Expected signal**: Top SV2 genes have high cytoskeletal annotation counts; Spearman r > 0.3, p<0.01.
**Null/control**: Same test for SV1 projection; Spearman r with non-cytoskeletal annotation count.
**Value**: MEDIUM | **Cost**: LOW (no new compute, just re-analysis of existing SVD output)

### H-L: Regulatory hub test — do high-degree TFs sit in specific SV2 poles?
**Hypothesis**: TRRUST TFs with high out-degree (many targets in the 209-gene set) are concentrated in one SV2 pole, while their target genes cluster in the same or opposite pole.
**Test**: For TFs with ≥10 targets in the 209-gene set, test whether TF is in top-K or bottom-K SV2 pole at L11. For each TF, compute the fraction of its targets also in the same pole. Null: random same-size gene subsets.
**Expected signal**: High-degree TFs preferentially in one pole; their targets in the same pole more than chance.
**Null/control**: Random gene subsets matched for set size; also test SV1.
**Value**: MEDIUM | **Cost**: LOW (analytical, same existing data)

### H-M: Layer-resolved co-pole rate curve for TRRUST + Spearman with spectral ratio
**Hypothesis**: The TRRUST co-pole enrichment profile across 12 layers tracks the SV2 variance fraction profile (both peak near L8), suggesting that regulatory geometry is encoded most strongly when SV2 is largest.
**Test**: Compute TRRUST SV2 co-pole rate at all 12 layers. Compute Spearman r between co-pole rate vector and SV2 var% vector from H03. Also plot co-pole rate vs layer.
**Expected signal**: Spearman r > 0.6 between co-pole rate and SV2 var% across layers.
**Null/control**: Same co-pole rate from random pairs × 12 layers.
**Value**: HIGH | **Cost**: LOW (extend H02 to all layers, ~5 min additional)

### H-N: Compartment-resolved co-pole test — are TF-target pairs enriched within specific compartments?
**Hypothesis**: The TRRUST co-pole signal is compartment-specific — TF-target pairs where both genes are cytoskeletal (or both mitochondrial) show stronger co-pole enrichment than pairs crossing compartment boundaries.
**Test**: At L11, stratify 326 TRRUST pairs by whether both genes share the same GO compartment annotation. Compute co-pole rate within and across compartment classes. Fisher test.
**Expected signal**: Within-compartment TF-target pairs show higher co-pole rate; across-compartment pairs show weaker or null effect.
**Null/control**: Random pairs stratified by same compartment-sharing criteria.
**Value**: HIGH | **Cost**: LOW (analytical stratification of existing H02 data)

---

## Top 3 for Immediate Execution

### Candidate 1 — HIGH PROBABILITY DISCOVERY
**H-A: TRRUST layer profile + edge-type stratification**
- Direct extension of the strongest result (H02, emp_p=0.000).
- Can be run in <30 min with existing TRRUST data and H02 infrastructure.
- If activation > repression in co-pole rate: first signed regulatory geometry finding in LLM gene embeddings.
- Pairs naturally with H-M (layer profile vs spectral ratio), which can run in the same script.

### Candidate 2 — HIGH RISK / HIGH REWARD
**H-E: Pole-flip gene tracking (mitochondrial polarity inversion)**
- The SV2 mito polarity flip (top L1–L3 → bottom L4+) is mechanistically unexplained.
- Identifying exactly which genes cross poles, and whether there is a coordinated flip event, could reveal how scGPT builds its mitochondrial representation over depth.
- Risk: may find only 2–3 genes flipping, too small for statistics. Reward: a clean per-gene trajectory story.
- Could also test whether flip genes are enriched in import machinery (TOMM/TIMM complexes) vs structural components.

### Candidate 3 — CHEAP BROAD SCREEN
**H-J: GO Biological Process enrichment in SVD poles**
- Compartment-level biology is confirmed (SV1–SV3). BP-level biology tests whether *functional programs* are also encoded — a qualitatively different claim.
- Low cost (use existing H01 infrastructure, substitute GO BP sets for GO CC sets).
- A positive result (e.g., "actin cytoskeleton organization" enriched in SV2) strengthens the cytoskeletal axis interpretation from compartment annotation to functional process annotation.
- Can run alongside H-A in the same iteration.
