# Hypothesis Roadmap: Post-iter_0051

---

## Retire / Deprioritize

### RETIRE NOW: SV2-4 Regulatory Proximity Claim
- **Reason**: Definitively falsified by iter_0051 H01. Co-expression explains 11/12 layers of TRRUST proximity in SV2-4.
- **What survives**: SV2-4 encodes co-expression geometry. TF enrichment in SV2-4 is an artifact of TFs being co-expression hubs.
- **Paper update required**: Replace "SV2-4 encodes regulatory proximity" → "SV2-4 encodes co-expression structure enriched for TF hub genes; regulatory-specific topology is confounded by co-expression at 11/12 layers."

### RETIRE NOW: 0-dim PH Discrete Clustering in SV2-4
- **Reason**: H03 showed fragmentation integral not significant (p=0.116). Circuit genes are compact but not topologically isolated.
- **Rescue**: Only via H1 homology (loops), which asks a different topological question. Retire the discrete-island framing.

### DEPRIORITIZE: Uncontrolled SV2-4 TRRUST Proximity Tests
- **Reason**: Any test without co-expression regression is uninformative given H01 result.
- **Status**: `retire_now` for all variants without co-expression control.

### DEPRIORITIZE: Full-512d PH
- Already flagged in iter_0050. Signal lives in spectral subspaces.

---

## New Hypothesis Portfolio

### H-N: Housekeeping Gene Enrichment in SV1-High Population
**Hypothesis**: The genes with highest SV1 loading (the non-regulatory "background" population identified in H02) are enriched for constitutively expressed housekeeping genes relative to a baseline gene list.
**Test**: Use Eisenberg & Levanon 2013 housekeeping gene list (~3804 human genes). Map to ENSG IDs or symbols. Fisher exact test: housekeeping overlap in top-20% SV1 genes vs. bottom-20% at L0 and L11. Also test ribosomal gene subset as control.
**Expected signal**: OR > 2.0 for housekeeping enrichment in SV1-high (p < 0.01). If ribosomal genes are also enriched, SV1 axis = expression ubiquity axis.
**Null**: No enrichment (OR ≈ 1.0). If confirmed, would suggest the SV1 axis reflects something else (e.g., expression level rather than cell-type specificity).
**Value**: high | **Cost**: low (all embeddings exist; need gene name mapping + housekeeping list)
**Addresses**: definitive biological identity of SV1 axis; resolves potential circular-reasoning concern in H02

---

### H-O: SV5-7 Regulatory Signal After SV2-4 Residualization
**Hypothesis**: The TF enrichment signal in SV5-7 (rbc=0.152, iter_0050 H03) is partially orthogonal to co-expression and survives residualization relative to SV2-4 co-expression structure.
**Test**: At each of 12 layers: (1) compute gene projections onto SV5-7 [2039×3]; (2) regress out the SV2-4 component via OLS (SV5-7_proj ~ SV2-4_proj); (3) run Mann-Whitney on residual projection magnitudes for TF vs. non-TF genes. Report rbc per layer.
**Expected signal**: If SV5-7 encodes an orthogonal regulatory dimension, TF enrichment survives with rbc > 0.05 in ≥3 layers. If SV5-7 TF enrichment is co-expression confounded like SV2-4, rbc collapses to ≈0.
**Null**: rbc < 0.02 in all layers after residualization.
**Value**: high | **Cost**: low (extends existing SVD code)
**Addresses**: the one remaining open question on spectral regulatory signal

---

### H-P: H1 Persistent Homology (Loops) on Circuit Genes at L8
**Hypothesis**: Circuit genes (TF+target, n=295) form topologically non-trivial loop structures (Betti-1 ≥ 1 survivors at significant filtration radius) in SV2-4 space at L8 that exceed what random gene sets produce.
**Test**: At L8, compute pairwise SV2-4 distances for circuit genes (n=295) and 10× bootstrap-sampled non-circuit control sets (n=295 each). Compute H1 Betti curves (ripser or gudhi). Compare: (a) max Betti-1 count, (b) total H1 persistence (area under curve), (c) longest surviving loop birth/death interval. Permutation p-value.
**Expected signal**: If regulatory gene pairs form cyclic co-expression motifs (A→B→C→A), these should manifest as 1-cycles. rbc for loop persistence > 0.1 in circuit vs. control.
**Null**: Betti-1 statistics indistinguishable from size-matched random sets.
**Value**: high | **Cost**: medium (ripser on 295×295 matrix: fast; bootstrap: <10 min)
**Addresses**: topological structure beyond compactness; novel finding if positive

---

### H-Q: TRRUST Module Internal Cohesion (Within-Regulon Compactness)
**Hypothesis**: Genes sharing a TF regulator (same regulon) are more compact in SV2-4 space at L8 than random gene sets of equal size, even after controlling for all-TRRUST average proximity.
**Test**: For each TF with ≥5 named targets in the 2039-gene set: compute mean pairwise SV2-4 distance within that regulon. Compare to 1000 random gene sets of same size drawn from non-circuit genes. Report empirical p-value per TF, and meta-p across TFs.
**Expected signal**: ≥60% of regulons show within-module compactness p < 0.05. If co-expression explains TF proximity, the strongest regulons should show the strongest within-module clustering.
**Null**: Within-regulon compactness not different from random. Would mean TF-target proximity is global artifact, not per-TF structure.
**Value**: high | **Cost**: low (all data present; simple pairwise computation)

---

### H-R: SV1 Signed Delta Gene Identity (L11 vs. L0 Movers)
**Hypothesis**: Genes whose SV1 loading shifts most from L0 to L11 (high |SV1_L11 − SV1_L0|) are enriched for cell-type-specific immune markers (L0-dominant movers) vs. housekeeping genes (L11-dominant movers).
**Test**: For each gene, compute delta = SV1_loading_L11 − SV1_loading_L0. Rank by signed delta. Top-50 positive delta = "L11 gainers"; top-50 negative delta = "L0 gainers." Map to gene symbols. Overlap with: (a) PanglaoDB immune cell markers, (b) Eisenberg housekeeping list, (c) TRRUST TF list.
**Expected signal**: L0 gainers = immune cell-type markers (variable, context-dependent); L11 gainers = housekeeping (constitutive). This would mechanistically explain the SV1 rotation as the model shifting from cell-type-specific to constitutive expression structure.
**Null**: No differential enrichment; gainers at both poles are random.
**Value**: high | **Cost**: low (uses h01_sv1_vectors.npy + gene name mapping)

---

### H-S: Cross-Layer Geodesic Path Length Predicts TF Status
**Hypothesis**: TF genes traverse longer total paths through embedding space across 12 layers than non-TF genes (sum of consecutive layer L2 distances), reflecting greater representational transformation.
**Test**: For each of 2039 genes: path_length = Σ_{l=0}^{10} ||emb[l+1] − emb[l]||₂. Mann-Whitney: TF genes (n=73 named TFs) vs. non-TF genes. Spearman correlation: path_length vs. SV2-4 projection magnitude at L8.
**Expected signal**: TF path length > non-TF (one-sided p < 0.05, rbc > 0.1). If TFs are co-expression hubs, they may also be the most "transformed" across layers.
**Null**: No path length difference by TF status.
**Value**: medium | **Cost**: low (simple L2 computation across layers)

---

### H-T: GO Biological Process Enrichment of SV1-Low Gene Set
**Hypothesis**: The SV1-low gene set (bottom 20% at L0 and L11 — enriched for TFs/targets) shows specific GO-BP enrichment for transcriptional regulation, immune response, and signal transduction, distinguishing them from the SV1-high anonymous population.
**Test**: At L0 and L11: extract bottom-20% SV1 genes with known names (within the ~295 named genes). Run GO-BP hypergeometric enrichment (gseapy) against all named genes as background. Report top-10 enriched terms per layer.
**Expected signal**: Strong enrichment for transcription factor activity, cytokine signaling, immune regulation. Should match known biology of TRRUST circuit genes.
**Null**: No GO enrichment beyond baseline.
**Value**: medium | **Cost**: low (gseapy or scipy hypergeometric; small gene set)

---

### H-U: SV1 Axis Reproducibility Under Gene Subsampling (Stability Test)
**Hypothesis**: The SV1 biological identity (TF depletion in top-20%) is stable under 80% gene subsampling, validating that the result is not driven by specific outlier genes.
**Test**: 100× bootstrap: sample 80% of 2039 genes, recompute SVD, align SV1 sign (procrustes), compute Fisher OR for TF depletion in top-20% SV1 at L0. Report distribution of ORs across bootstrap samples.
**Expected signal**: Median OR < 0.2 with narrow CI — confirms robustness. If OR is unstable (wide CI), the H02 result may be driven by a small gene subset.
**Null**: OR varies widely across bootstraps (CI spans 1.0).
**Value**: medium | **Cost**: low (fast resampling; validates H02)

---

### H-V: L8 Specialization — Attention-Pattern Proxy via Residual Magnitude
**Hypothesis**: Layer 8 produces a qualitatively different transformation of TF gene embeddings (relative to adjacent layers) that explains why L8 is the only layer retaining residual regulatory signal after co-expression regression.
**Test**: For each gene at each layer, compute residual magnitude = ||emb[l] − emb[l-1]|| (layer-to-layer change). At L7→L8 transition, compare residual magnitudes for TF vs. non-TF genes (Mann-Whitney). Also compare the L7→L8 jump size to L6→L7 and L8→L9 for TF genes specifically.
**Expected signal**: TF genes show larger residual magnitudes at L7→L8 than non-TF genes or than adjacent layer transitions, consistent with L8 being a functionally critical attention layer for regulatory gene processing.
**Null**: No TF-specific spike at L7→L8.
**Value**: medium | **Cost**: low (embedding differences already computable from layer_gene_embeddings.npy)

---

### H-W: Spectral Scan SV8-SV15 for Regulatory Signal (Extended Decay Curve)
**Hypothesis**: Regulatory proximity signal (TRRUST rbc) decays monotonically through SV groups (SV2-4 > SV5-7 > SV8-10 > SV11-15), following a geometric (log-linear) decay that quantifies how TF regulatory structure is distributed across spectral dimensions.
**Test**: At L8, compute pairwise distances in SV8-10, SV11-13, SV14-16 projections. Mann-Whitney TRRUST vs. negative pairs. Combine with existing rbc values for SV2-4 (rbc≈0.156) and SV5-7 (rbc≈0.152). Fit log-linear regression: rbc ~ log(SV_group_rank).
**Expected signal**: Monotone decay with R² > 0.7 and significant negative slope. Would imply regulatory information is spread across all spectral dimensions in a power-law pattern.
**Null**: Non-monotone (e.g., flat after SV5-7) or increase at higher SVs.
**Value**: medium | **Cost**: low (extends prior SVD code)

---

### H-X: STRING Co-expression Validation of SV2-4 Co-expression Claim
**Hypothesis**: STRING high-confidence co-expression edges (score ≥700) show SV2-4 proximity (rbc) ≥ TRRUST proximity, quantitatively confirming that co-expression network structure (not regulatory topology) drives the SV2-4 geometry.
**Test**: Download STRING v12 human network. Filter: score ≥700, co-expression channel. Map to the 2039-gene set. Compute SV2-4 pairwise distances for STRING edges vs. non-STRING pairs at each layer. Compare rbc(STRING) vs. rbc(TRRUST) per layer.
**Expected signal**: rbc(STRING) ≥ rbc(TRRUST) in ≥8/12 layers, providing the strongest external validation that SV2-4 is a co-expression encoder.
**Null**: rbc(TRRUST) > rbc(STRING) — would require revising the co-expression narrative.
**Value**: high | **Cost**: medium (STRING download, parsing, mapping)

---

### H-Y: 0-dim PH at Layers L5-L6 (Replication Test)
**Hypothesis**: If circuit genes show topological clustering (0-dim PH fragmentation integral), it appears at L5 or L6 (where rbc was elevated: L6=0.127, L5=0.113) rather than only at L8.
**Test**: Replicate H03 method at L5 and L6. Extract SV2-4 of circuit vs. matched non-circuit genes. Compute fragmentation integral and bootstrap null (n=500). Report ratio and p-value.
**Expected signal**: p < 0.05 at ≥1 of {L5, L6}. If the H03 mixed result at L8 is due to layer choice, testing L5/L6 may reveal topological cluster structure.
**Null**: p > 0.05 at both layers.
**Value**: low | **Cost**: low (directly reuses H03 code with layer parameter change)

---

### H-Z: Regulatory Edge Prediction AUC: SV1 × SV2-4 Cross-Term vs. SV2-4 Alone
**Hypothesis**: A logistic regression edge predictor using the product (SV1_loading_difference × SV2-4_distance) achieves higher TRRUST edge prediction AUC than SV2-4 distance alone, indicating the regulatory membership axis (SV1) improves regulatory edge prediction when combined with co-expression geometry.
**Test**: For each TRRUST pair (pos n=589, neg n=1523): compute (a) SV2-4 distance, (b) |SV1_loading_i − SV1_loading_j| at L8. Train logistic regression with these features. Cross-validate AUC. Compare to SV2-4-only baseline AUC.
**Expected signal**: AUC gain > 0.02 (e.g., 0.58 → 0.60+) when SV1 loading difference is included. Mechanistic: pairs where both genes are low-SV1 (regulatory genes) AND close in SV2-4 are most likely true TRRUST edges.
**Null**: AUC unchanged or lower with SV1 cross-term.
**Value**: medium | **Cost**: low (sklearn logistic regression; data all present)

---

## Top 3 for Immediate Execution

### #1 — High-Probability Discovery: H-N (Housekeeping Gene Enrichment in SV1-High)
**Rationale**: H02 established a strong positive (OR=0.108 for TF depletion in SV1-high). The natural biological interpretation is that SV1-high genes are housekeeping genes. This test directly confirms or refutes that interpretation using an independent, widely-cited gene list (Eisenberg housekeeping set). All embedding data exists. Gene name mapping is already established (295 named genes). The result is binary: either SV1 = housekeeping axis (strong narrative for paper) or it's not (requires alternative interpretation). Zero ambiguity, low cost. **Expected runtime**: <5 min. If also checking ribosomal genes and running GO-BP simultaneously (H-T), the full SV1 biological identity story can be wrapped up in a single execution.

### #2 — High-Risk/High-Reward: H-P (H1 Betti Curves on Circuit Genes at L8)
**Rationale**: H03 showed circuit genes are compact (0-dim confirmed) but not topologically isolated (0-dim cluster topology rejected). The logical next step is H1 homology — do regulatory feedback loops leave geometric traces as 1-cycles in the SV2-4 representation? If yes, this is a genuinely novel topological finding about how LLMs encode regulatory circuit structure. If no, we cleanly retire all PH-based regulatory circuit topology claims. ripser can handle 295×295 distance matrix in seconds. Bootstrap permutation of 10 non-circuit control sets adds <5 min. **Expected runtime**: <15 min.

### #3 — Cheap Broad Screen: H-O (SV5-7 Regulatory Signal After SV2-4 Residualization)
**Rationale**: SV5-7 showed rbc=0.152 for TF enrichment in iter_0050 H03, but co-expression confound was not tested. H01 showed SV2-4 is co-expression confounded. Whether SV5-7 is orthogonally regulatory is the critical open question for whether any spectral axis in these embeddings encodes regulatory topology beyond co-expression. The computation is minimal (extend existing SVD code with one residualization step). Result gates the entire "spectral regulatory axes" research thread: if SV5-7 survives, there is a genuine regulatory signal somewhere in the spectrum; if not, the SV1 biological identity story becomes the primary positive finding. **Expected runtime**: <5 min.
