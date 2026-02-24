# Brainstormer Hypothesis Roadmap — iter_0009 → iter_0010

---

## Retire / Deprioritize

| Direction | Status | Reason |
|-----------|--------|--------|
| Geodesic/rewiring-null PH survival | **retire_now** | 4 consecutive iterations with 0/24 significant results across all null formulations (degree-preserving, metric-matched, edge-length-constrained, bridge-conditioned). No rescue candidate remains. |
| GO BP cosine-distance clustering (209-gene vocab) | **retire_now** | Structurally underpowered; z~-1.1 consistently. Signal lives in SVD projections not cosine clustering. |
| TRRUST co-target Euclidean clustering | **deprioritize** | Two attempts inconclusive. Rescue possible if re-run in SVD projection space; not worth raw-Euclidean retries. |

---

## New Hypothesis Portfolio

### A. SV2/SV3 Axis Characterization (Direct Extensions of H01-H03)

**A1 — SV2 Full Layer Scan**
- **Hypothesis:** SV2's EV enrichment at the bottom pole and IL-4/immune enrichment at the top pole are both present across all 12 layers, potentially strengthening in depth.
- **Test:** Replicate H02 design (12-layer × 8-compartment Fisher scan, N=200 gene-label shuffles) using SV2 instead of SV1 projections. Add GO:0032753 (IL-4 production regulation) to the 8 compartments.
- **Expected signal:** EV persistent (as SV2 bottom), IL-4/immune persistent (SV2 top); both may strengthen at L11.
- **Null/control:** Gene-label shuffle N=200 per cell.
- **Value:** high | **Cost:** low (same code, change SV index)

**A2 — SV3 Characterization Scan**
- **Hypothesis:** The third SVD axis encodes a biologically distinct program not captured by SV1 (secretory) or SV2 (EV/immune), potentially metabolic or cell-cycle related.
- **Test:** At layer 11, project 209 named genes onto SV3. Fisher enrichment for 15 GO terms spanning major programs (mitochondrion, cell cycle, translation, cytoskeleton, DNA repair, lipid metabolism, etc.) with N=200 shuffles.
- **Expected signal:** At least one GO term enriched in SV3 top or bottom pole with OR>3.
- **Null/control:** Gene-label shuffle.
- **Value:** high | **Cost:** low

**A3 — SV2 Cross-Layer Rank Stability**
- **Hypothesis:** SV2's gene rank is less stable across layers than SV1's (0.929), reflecting a more context-dependent immune/EV axis.
- **Test:** Replicate H03 design (Spearman r between adjacent layer SV projections, N=200 row-permutation null) for SV2 and SV3. Report mean adjacent r and variance explained per SV.
- **Expected signal:** SV2 mean r < SV1 mean r; SV3 mean r lowest, indicating increasing axis instability for higher-order components.
- **Null/control:** Row-permutation null N=200.
- **Value:** high | **Cost:** low (same code, extend loop over SVs)

---

### B. L10→L11 Transition Mechanism

**B1 — Attention Weight Axis Alignment at L10→L11**
- **Hypothesis:** The L10→L11 SV1 rank discontinuity (r=0.663) is driven by a specific set of attention heads in layer 10 that restructure gene relationships. These heads should show abnormally high gene-pair distance rewiring relative to other layers.
- **Test:** Extract attention weights from layer 10 (mean over heads and cells). Compute gene-pair weight matrix. Correlate with SV1 projection change (|Δproj_L10→L11| per gene). Null: permute gene-attention assignments.
- **Expected signal:** Genes with large SV1 projection change between L10-L11 have systematically higher mean attention from layer-10 heads.
- **Null/control:** Gene index permutation.
- **Value:** high | **Cost:** medium (requires attention extraction from model)

**B2 — Mito Transience at L2-4 and L10: Single or Dual Mechanism?**
- **Hypothesis:** The mito enrichment signal on SV1 at L2-4 and its resurgence at L10 (H02 result) reflect two distinct mechanisms: early embedding initialization vs. late residual restructuring at the L10-L11 boundary.
- **Test:** Extract the specific mito-enriched gene set from SV1 top pole at L2 vs L10. Compute Jaccard overlap between the two sets. Compare to bootstrap null (random same-size gene subsets). If overlap > expected → same mechanism; if near-zero → distinct.
- **Expected signal:** Jaccard overlap significantly lower than expected, indicating layer-specific mito gene subsets.
- **Null/control:** Bootstrap null over random same-size gene sets.
- **Value:** medium | **Cost:** low

---

### C. Biological Anchoring: Confounder Tests and Extended Annotation

**C1 — Annotation Density Confounder (SV1 position ~ GO annotation count)**
- **Hypothesis:** Genes with more GO annotations are systematically placed at the SV1 top pole, confounding the ER lumen/secreted enrichment.
- **Test:** Compute Spearman correlation between per-gene GO annotation count and SV1 projection value at layer 11. Compare against permuted null (N=1000).
- **Expected signal:** If |r| < 0.2 with p > 0.05, confounder is absent. If significant positive r, confounder is present and all SVD enrichments need covariate-adjusted reanalysis.
- **Null/control:** Permutation of SV1 projection values.
- **Value:** high (critical validity check) | **Cost:** low

**C2 — TRRUST TF→Target SVD Projection Alignment**
- **Hypothesis:** TRRUST TF-target pairs are more co-projected onto SV1 (both at same pole) than random gene pairs of matched size.
- **Test:** For each TRRUST TF with ≥5 targets, compare mean |SV1_TF - SV1_target_mean| to null (random same-size gene sets). Also test SV2. Null N=200.
- **Expected signal:** TF-target pairs share SV1 pole more than random (z < -1.96) for secretory-pathway TFs (XBP1, ATF6, CREB3).
- **Null/control:** Gene-label shuffle.
- **Value:** high | **Cost:** low

**C3 — STRING High-Confidence PPI Pairs: SV1/SV2 Co-Projection**
- **Hypothesis:** Protein-protein interaction partners (STRING score > 700) tend to be co-projected onto the same SV1/SV2 pole.
- **Test:** Extract STRING edges among 209 named genes at confidence > 700. For each edge, compute signed SV1 product (same sign = co-pole). Compare mean signed product to null (random edge permutation).
- **Expected signal:** Mean signed SV1 product > 0 for interacting pairs vs null.
- **Null/control:** Random edge permutation among 209 genes.
- **Value:** high | **Cost:** low

---

### D. Cross-Model and Cross-Domain

**D1 — Geneformer SV1 Axis Recovery**
- **Hypothesis:** Geneformer's residual stream SVD, on the same 209 named genes, reveals a similar primary axis (ER lumen / secretory vs. immune) with Spearman r > 0.5 against scGPT SV1 projection rankings.
- **Test:** Load Geneformer embeddings (same gene set), compute SVD, project 209 genes onto Geneformer SV1. Compute Spearman r vs scGPT L11 SV1 projections. Null: permute gene labels.
- **Expected signal:** r > 0.5, emp_p < 0.05.
- **Null/control:** Gene-label permutation N=1000.
- **Value:** high | **Cost:** medium (requires Geneformer embedding access)

**D2 — Cross-Cell-Context SVD Stability (lung vs immune vs external-lung)**
- **Hypothesis:** The SV1 secretory axis and SV2 EV/immune axis are context-invariant (lung, immune, external-lung gene embeddings all give similar SV1/SV2 projections for the 209 shared genes).
- **Test:** Load all three context embeddings at layer 11. Compute SVD per context, project 209 named genes. Compute pairwise Spearman r of SV1 projections across contexts (3×3 matrix). Null: gene-label permutation.
- **Expected signal:** All off-diagonal r > 0.7, emp_p < 0.05.
- **Null/control:** Gene-label permutation N=500.
- **Value:** high | **Cost:** medium

---

### E. Spectral/Geometric Novelty

**E1 — Singular Value Gap Profile as Layer Function**
- **Hypothesis:** The SV1/SV2 singular value ratio (measuring spectral dominance of the primary axis) increases monotonically with depth, mirroring the variance-explained trajectory found in H03.
- **Test:** At each of 12 layers, compute SVD of 209-gene mean-centered embedding. Record SV1/SV2 ratio and effective rank. Correlate with layer depth. Null: feature-shuffle N=50.
- **Expected signal:** Monotonically increasing SV1/SV2 ratio; L10-L11 transition coincides with ratio jump.
- **Null/control:** Feature-shuffle destroys monotonicity.
- **Value:** medium | **Cost:** low

**E2 — Intrinsic Dimension of SVD-Filtered Subspace**
- **Hypothesis:** The manifold intrinsic dimension within the top-K SVD subspace (K=3 or 5) is lower than in the full space, indicating that the biological signal is concentrated in a near-rank-1 manifold within the SVD subspace.
- **Test:** Project 209 genes onto top-3 SVD components at each layer. Apply TwoNN estimator. Compare to TwoNN in full PCA-20 space (existing iter_0004 data).
- **Expected signal:** TwoNN ID in SVD-3 subspace < TwoNN in full PCA-20.
- **Null/control:** Feature-shuffle null N=50.
- **Value:** medium | **Cost:** low

**E3 — GO Term Centroid Displacement Across Layers (directional drift)**
- **Hypothesis:** GO term gene centroids undergo directional (not isotropic) drift across layers — secretory/ER centroids drift toward the SV1 top pole, while immune centroids drift toward SV2 top.
- **Test:** For 8 GO terms from H02, compute the centroid (mean projection onto SV1) at each of 12 layers. Compute monotonicity (Kendall τ) of centroid drift vs layer depth. Null: feature-shuffle N=200.
- **Expected signal:** ER lumen centroid has significant positive Kendall τ (monotonic approach to SV1 top). Immune centroid has significant drift on SV2.
- **Null/control:** Feature-shuffle.
- **Value:** medium | **Cost:** low

---

## Top 3 for Immediate Execution

### Slot 1: High-Probability Discovery Candidate
**A1 + A3 (combined as one experiment): SV2/SV3 Layer Scan + Stability**
- Run H02-style 12-layer × 9-compartment scan on SV2 and SV3 projections.
- Simultaneously run H03-style Spearman stability for SV2 and SV3.
- Both extend proven working code with trivial parameter changes.
- Expected outcome: SV2 EV persistent early, IL-4 strengthens late; SV3 reveals a new axis; SV2 less stable than SV1.
- **Directly answers the most pressing open question from H01-H03.**

### Slot 2: High-Risk / High-Reward
**C2 + C3 combined: Regulatory/PPI Structure in SVD Projection Space**
- Test TRRUST TF→target co-pole enrichment (C2) and STRING high-confidence PPI co-pole enrichment (C3) in one script.
- If positive: establishes that SVD axes encode regulatory and physical interaction networks — a strong mechanistic claim.
- If negative: rules out simple co-regulation explanation for the SV1 axis structure.
- New null design (gene-label shuffle on edges, not random groups) is a methodological advance over the failed iter_0004 TRRUST test.

### Slot 3: Cheap Broad Screen
**C1 + E1 combined: Annotation Density Confounder + SV Ratio Profile**
- C1 (annotation density confounder check) is a <10-line addition to existing analysis scripts.
- E1 (SV1/SV2 ratio per layer) is another <20-line addition.
- Together they (a) validate that prior SVD results are not annotation-density artifacts, and (b) characterize the spectral dominance profile across depth.
- Cheap, fast, high-validity-impact.
