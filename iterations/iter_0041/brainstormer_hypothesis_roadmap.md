# iter_0041 Brainstormer Hypothesis Roadmap

---

## Retire / Deprioritize

| Direction | Reason | Decision |
|-----------|--------|----------|
| GC-plasma subspace angle | Near-orthogonal throughout, no signal path remaining | `retire_now` |
| Functional annotation → proximity (TRRUST/GO Jaccard) | AUROC ≈ 0.50 across all layers, confirmed null (iter_0035) | `retire_now` |
| PC2/PC3 cell-type axes at L11 | All null, underpowered (iter_0034 H02) | `retire_now` |
| Cross-model Geneformer layer-wise alignment | Missing artifact (layer-wise Geneformer tensors); not rescuable without new data extraction | `retire_now` (rescue only if Geneformer layer embeddings become available) |
| Chromosomal proximity | Confirmed null control (iter_0025 H03) | `retire_now` |
| BCL6 as GC attractor member | Refuted — BCL6 lives in metabolic neighborhood throughout all layers | `retire_now` (replace with BCL6-metabolic-isolation claim) |

---

## New Hypothesis Portfolio

### H-A: TwoNN ID Onset Alignment (Inflection at L3?)
**Hypothesis**: The TwoNN intrinsic dimensionality decrease for B-cell markers shows a change-point or inflection near L3, coinciding with the GC-TF attractor onset.
**Test**: Compute dID/dL per layer; fit piecewise-linear model or change-point detector (e.g., PELT); test whether the best-fit break is at L3 vs. other layers.
**Expected signal**: Change-point at L3±1; before L3 flat/slow ID decline, after L3 accelerated decline.
**Null**: Break point is random within L0-L11 (permute layer order and refit).
**Value**: high | **Cost**: low

### H-B: PAX5 L0 Pre-wiring Permutation Null
**Hypothesis**: PAX5's B-cell surface receptor neighborhood at L0 is statistically significant vs. random — B-cell neighbor content is higher than expected by chance.
**Test**: Shuffle gene-to-embedding assignment 500 times; recompute PAX5's k=20 neighbor B-cell content at L0; compare to observed value.
**Expected signal**: Observed B-cell content (1-2/20) > 95th percentile of shuffle distribution.
**Null**: Permuted B-cell content ≥ observed.
**Value**: high | **Cost**: low

### H-C: BCL6 GO Enrichment of Metabolic Neighborhood
**Hypothesis**: BCL6's top-20 neighbors are enriched in metabolic stress GO terms (glycolysis, glutamine metabolism, NAD synthesis) relative to background.
**Test**: Map BCL6 top-20 neighbors per layer to GO BP/MF terms using mygene. Fisher exact test for metabolic process GO terms vs. 195-gene background.
**Expected signal**: Significant enrichment (p<0.05) in "metabolic process", "response to nutrient levels", or specific glycolysis/NAD terms.
**Null**: GO term frequency matches background rate.
**Value**: high | **Cost**: low

### H-D: T-cell Attractor Replication with Cycle1 Vocab
**Hypothesis**: T-cell TFs (FOXP3, GATA3, TBX21) are available in cycle1 vocab; with full T-cell TF coverage, a stronger T-cell attractor signal may emerge, or the pre-wiring characterization can be finalized.
**Test**: Load cycle1 embeddings; check availability of FOXP3, GATA3, TBX21; if present, run same rank/centroid convergence test as H01. Compare rho and onset layer.
**Expected signal**: Either (a) T-cell attractor rho improves with more TFs and confirms pre-wiring characterization, or (b) reveals T-cell sub-structure not visible with RUNX3 only.
**Null**: rho remains near -0.543 (single-TF result).
**Value**: high | **Cost**: medium

### H-E: B-cell ID Bootstrap Confidence Intervals
**Hypothesis**: The TwoNN ID decrease for B-cell markers is robust to subsampling within the marker set.
**Test**: Bootstrap resample the 4 B-cell markers with replacement 200 times; compute TwoNN ID at each layer for each bootstrap; derive 95% CI per layer.
**Expected signal**: CIs narrow at deeper layers (consistent with compression), wider at L0 (noisier estimates).
**Null**: CIs overlap with T-cell/myeloid ID range at all layers.
**Value**: medium | **Cost**: low

### H-F: Metabolic Isolation as a General Principle (Other Context-Dependent TFs)
**Hypothesis**: BCL6 is not unique — other context-dependent TFs (STAT3, HIF1A, MYC) that serve multiple lineages also occupy metabolic/stress neighborhoods rather than lineage-specific neighborhoods.
**Test**: For STAT3, HIF1A, MYC: find k=20 neighbors across layers; count B-cell/T-cell/myeloid content; compare to master TFs like PAX5 (B-cell) and SPI1 (myeloid).
**Expected signal**: Context-dependent TFs have <1 lineage-specific neighbor (like BCL6); lineage-specific master TFs have ≥2.
**Null**: Context-dependent TFs are indistinguishable from master TFs in neighborhood content.
**Value**: high | **Cost**: low

### H-G: GC-TF Attractor Onset Precision (Layer 3 Characterization)
**Hypothesis**: The L3 onset of GC-TF convergence is a sharp transition (significant between L2 and L3) rather than a gradual drift.
**Test**: For each of BATF, BACH2, PAX5: compute rank near B-cell centroid at all layers; test L2 vs L3 rank difference using paired t-test (n=3 TFs); compute effect size.
**Expected signal**: Significant drop in mean rank between L2 and L3 (p<0.05 per TF, consistent direction).
**Null**: No significant L2-L3 transition; gradient is smooth.
**Value**: medium | **Cost**: low

### H-H: Persistent Homology of B-cell Submanifold (H0 Betti count vs depth)
**Hypothesis**: The B-cell submanifold (5 core markers) simplifies topologically with depth — H0 component count decreases from multiple connected components to one as the cluster converges.
**Test**: Run Ripser on B-cell marker embedding vectors (k=5 points) at each layer; compute H0 lifetime and H0 count. Compare to T-cell (k=10) and myeloid (k=6) at same layers.
**Expected signal**: H0 count reduces from 2-3 components (L0) to 1 (L6+) for B-cell; no such reduction for T-cell.
**Null**: All cell types show 1 component throughout (already clustered), or H0 count doesn't change.
**Value**: medium | **Cost**: low

### H-I: Attention Weight Concentration on GC-TF vs BCL6 at L3
**Hypothesis**: At the L3 onset layer, GC-TFs (BATF, BACH2, PAX5) receive elevated attention from B-cell marker genes, while BCL6 does not — mechanistically explaining the differential convergence.
**Test**: Load attention scores tensor (if layer-resolved attention exists); extract attention weights from B-cell markers (CD19, MS4A1, CD79A) toward GC-TFs and BCL6 at each layer. Compare normalized attention at L2 vs L3.
**Expected signal**: GC-TF attention from B-cell markers increases sharply at L3; BCL6 attention remains low.
**Null**: Attention is uniform; no TF-specific enrichment at L3.
**Value**: high | **Cost**: medium (depends on attention artifact availability)

### H-J: Centroid Distance Trajectory for PAX5 vs BCL6 (Full Layer Profile)
**Hypothesis**: PAX5 maintains consistently low rank near B-cell centroid throughout all 12 layers (flat or slightly improving), while BCL6 shows no convergence and possibly diverges.
**Test**: Plot rank near B-cell centroid at all 12 layers for PAX5 and BCL6 separately. Fit linear trend; test slope significance.
**Expected signal**: PAX5 slope ≈ 0 or negative (consistently close); BCL6 slope ≈ 0 or positive (consistently far or diverging).
**Null**: Both show similar trajectories.
**Value**: medium | **Cost**: low

### H-K: Myeloid Geometric Compression Absence (Formal Test)
**Hypothesis**: The myeloid TwoNN ID shows a slight *increase* (rho=+0.699 from H03), which is a distinct inverse-compression pattern. Test if this is significant and characterize its biological meaning.
**Test**: Formalize the myeloid ID increase: compute bootstrap CI on the slope, test against null slope=0. If significant, characterize L0 vs L11 ID change as an "expansion" and check if myeloid markers become more dispersed (anti-clustering).
**Expected signal**: myeloid ID slope is positive and significant; myeloid precision@k DECREASES with depth.
**Null**: myeloid ID change not distinguishable from flat.
**Value**: medium | **Cost**: low

### H-L: Multi-lineage Centroid Separation Trajectory (Normalized)
**Hypothesis**: Normalized centroid separation (iter_0035 method) applied across B-cell/T-cell/myeloid shows that ONLY B-cell normalized separation increases monotonically with depth.
**Test**: Replicate iter_0035 normalized centroid separation for T-cell and myeloid panels; compare rho values.
**Expected signal**: B-cell rho ≈ +0.97 (established); T-cell and myeloid rho ≈ 0 or negative.
**Null**: All lineages show similar normalized separation trends.
**Value**: medium | **Cost**: low

### H-M: Cycle4 vs Cycle1 Embedding Consistency for PAX5/BCL6
**Hypothesis**: The PAX5 pre-wiring and BCL6 metabolic isolation findings replicate in cycle1 embeddings (different training cycle).
**Test**: Load cycle1 embeddings; find PAX5 and BCL6 neighborhoods at L0 and L11; count B-cell content.
**Expected signal**: PAX5 B-cell neighbor count ≥ 1 at L0; BCL6 B-cell neighbor count = 0.
**Null**: Results differ between cycles (training-cycle artifact).
**Value**: medium | **Cost**: low

### H-N: Functional Neighborhood Stability Across Layers for GC-TFs
**Hypothesis**: GC-TFs (BATF, BACH2) show increasing B-cell neighborhood content across layers (progressive incorporation into B-cell geometric context) mirroring their centroid convergence.
**Test**: For BATF and BACH2: compute k=20 neighbor B-cell content at each of L0, L3, L6, L9, L11; track B-cell neighbor count.
**Expected signal**: B-cell content increases from 0-1 at L0 to 2-4 at L11.
**Null**: No systematic change in B-cell neighborhood content.
**Value**: high | **Cost**: low

---

## Top 3 for Immediate Execution

### #1 — High-probability discovery candidate
**H-A: TwoNN ID Onset Alignment (Inflection at L3)**

The B-cell TwoNN ID decrease (rho=-0.951) combined with the GC-TF attractor onset at L3 creates a testable prediction: the ID change should accelerate at L3. A quantitative change-point test directly linking two independently discovered signals is a publishable landmark finding. Low cost, uses existing H03 data with one additional analysis pass.

### #2 — High-risk/high-reward candidate
**H-F: Metabolic Isolation as General Principle (STAT3, HIF1A, MYC neighborhoods)**

If context-dependent multi-lineage TFs broadly occupy metabolic neighborhoods rather than lineage-specific ones, this becomes a general mechanistic principle: scGPT organizes TFs by their regulatory specificity, with promiscuous TFs gravitating toward metabolic hubs. This would substantially elevate the BCL6 finding from an anecdote to a principle. Risk: the pattern may not generalize; metabolic neighborhood for BCL6 may be idiosyncratic.

### #3 — Cheap broad-screen candidate
**H-N + H-C combined: GC-TF Neighborhood Trajectory + BCL6 GO Annotation**

Run k=20 neighbor tracking for BATF, BACH2 across all layers (H-N) AND BCL6 GO enrichment (H-C) in a single script. Both are low-cost extensions of the H02 neighborhood analysis already implemented. Together they complete the neighborhood characterization story: GC-TFs converge into B-cell neighborhood (as expected from attractor claim), BCL6 is biologically annotated as metabolic (GO confirmation). Deliverable: two clean supplementary-quality figures + GO table.
