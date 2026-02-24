# Brainstormer Hypothesis Roadmap — iter_0002

## Retire / Deprioritize

| Hypothesis | Status | Reason | Effort to Retest |
|------------|--------|--------|------------------|
| H05: Distance-permutation null | RETIRE_NOW | 0/12 layer-tests significant across domains; over-adversarial | High (already done) |
| H06–H09: Rewiring-null (immune) | RETIRE_NOW | Metric matching + edge-length binning failed; no positive signal after multi-iteration refinement | High (diagnostic exhausted) |
| H11: Bridge-conditioned attribution | RETIRE_NOW | Split-confounded strata make attribution unidentifiable; pooled delta contradicts hypothesis | High (unfixable without protocol redesign) |
| H10: Metric calibration shift | RETIRE_NOW | Geodesic–Euclidean shift = +0.180; too small to improve conclusions | Medium (done; no rescue) |
| H12: Quantile-constrained rewiring | RETIRE_NOW | Edge-length binning preserved; same negative result as unconstrained rewiring | High (already done) |

**Closure narrative:** These branches exhausted calibration space within a single-null paradigm. Feature-shuffle null remains the positive anchor; stronger nulls uniformly fail, suggesting signal is *real but fragile* (depends on feature relationships, not just distance/connectivity). Next step: test whether fragile signal maps to biological ground truth (TRRUST, GO, cell ontology).

---

## New Hypothesis Portfolio

### **TIER 1: BIOLOGICAL ANCHORING (High priority)**

#### **N01: TRRUST Co-Regulatory Network Alignment**
- **Hypothesis:** Gene pairs with high TRRUST co-regulatory edges show shorter geodesic distances in scGPT residual embeddings than non-regulatory pairs (degree-preserving permutation null).
- **Test design:**
  - Collect TRRUST edge list (TF-target pairs) for relevant cell types
  - For each seed-layer combination: compute pairwise geodesic distances in scGPT embeddings
  - Compute mean geodesic distance for (i) TRRUST-connected pairs, (ii) random-node-shuffled null
  - Measure effect size (mean difference / std error) and Fisher-combined p-value across seeds-layers
  - Control: randomize TF labels but preserve edge structure
- **Expected signal:** Observed geodesic distance << null (Δ < 0, z < -3)
- **Null:** Random rewiring of TF-target edges; preserve in-degree/out-degree
- **Value:** HIGH (directly links topology to gene regulatory biology)
- **Cost:** MEDIUM (requires TRRUST curation, but geodesic computation already available)

#### **N02: Gene Ontology Functional Clustering**
- **Hypothesis:** Genes sharing GO terms (e.g., "immune response", "cell cycle") cluster tightly in kNN neighborhoods (kNN CC elevation, local homogeneity).
- **Test design:**
  - Retrieve GO annotations for lung/immune cell-type genes
  - For each layer: compute mean kNN CC within GO-term subsets vs. between-term pairs
  - Measure within-term CC elevation relative to between-term CC
  - Use Fisher-combined p-value (Wilcoxon signed-rank per layer) across seeds
  - Control: shuffle GO labels but preserve gene list
- **Expected signal:** Within-term CC >> between-term CC (Δ > +0.05, z > +2)
- **Null:** Permuted GO assignment
- **Value:** HIGH (validates CC signal is functionally coherent, not statistical artifact)
- **Cost:** LOW (requires GO API access; CC metrics already computed)

#### **N03: STRING Network Path-Length Coupling**
- **Hypothesis:** Genes connected by short STRING PPI paths cluster in residual embeddings (kNN cliques); mean clique size increases with STRING-path proximity.
- **Test design:**
  - Build weighted STRING network for protein interactions
  - For each layer: compute shortest STRING path distance for all gene pairs
  - Stratify kNN neighborhoods by path distance bins (distance 1, 2, 3, ≥4)
  - Measure within-neighborhood edge density as function of string path distance
  - Test correlation (Spearman ρ) between path distance and edge density
  - Permutation null: randomize STRING edges while preserving degree
- **Expected signal:** Shorter paths → higher edge density (ρ > +0.3, Fisher p<0.05)
- **Null:** Permuted STRING topology
- **Value:** HIGH (connects graph topology to protein-level biology)
- **Cost:** MEDIUM (requires STRING download + path computation; ρ test is fast)

#### **N04: Cell Ontology Type Enrichment in Neighborhoods**
- **Hypothesis:** In immune/lung embeddings, kNN neighborhoods are enriched for specific cell types (macrophage, T cell, epithelial) beyond random chance.
- **Test design:**
  - Map genes to cell type of origin (from metadata)
  - For each layer: compute purity = fraction of same-type genes in kNN neighborhood
  - Measure mean purity vs. null (random neighbor assignment)
  - Use Fisher-combined p-value across seeds-layers
  - Control: shuffle cell-type labels but preserve kNN structure
- **Expected signal:** Purity >> random (Δ > +0.1, z > +3)
- **Null:** Shuffled cell-type assignment
- **Value:** HIGH (directly validates that embeddings preserve cell-identity structure)
- **Cost:** LOW (only requires metadata; kNN neighborhoods already computed)

---

### **TIER 2: ALTERNATIVE TOPOLOGY FAMILIES (Medium-high priority)**

#### **N05: Betti-0 (Connected Components) Persistence**
- **Hypothesis:** Betti-0 lifetime (connected-component lifetime in Ripser filtration) is elevated in real embeddings vs. feature-shuffle null.
- **Test design:**
  - Use same scGPT embeddings (3 seeds, 12 layers)
  - Run Ripser with Rips filtration on PCA(20)-reduced data
  - Extract H0 death-birth differences (lifetime of components) → aggregate as sum/mean
  - Compare against feature-shuffle null (20 replicates)
  - Compute Fisher-combined p-value across seeds-layers
- **Expected signal:** If topology is robust, H0 lifetime should be >0 and elevated (Δ > +2, z > +2)
- **Null:** Feature-shuffle null (same as H1 protocol)
- **Value:** MEDIUM (tests whether topology family matters; H1 might be specific instance)
- **Cost:** LOW (Ripser already integrated; minimal code change)

#### **N06: Persistent Cluster Morphology (H2 / Betti-2)**
- **Hypothesis:** Betti-2 lifetime (persistent 2-dimensional voids in filtration) is lower in real embeddings vs. null, indicating sparser, less-cavity-dense topology.
- **Test design:**
  - Run Ripser H2 computation on same setup as H1
  - Extract H2 lifetime sum/mean
  - Compare against feature-shuffle null
  - Compute Fisher-combined p across seeds-layers
- **Expected signal:** H2 lifetime lower in real (negative delta), consistent with sparsity interpretation
- **Null:** Feature-shuffle
- **Value:** MEDIUM (complementary to H1; tests if signal is loop-specific or general topology property)
- **Cost:** LOW (Ripser H2 already available)

#### **N07: Filtration-Order Robustness (Vietoris-Rips vs. Alpha Complex)**
- **Hypothesis:** H1 persistence signal is invariant to choice of filtration (Rips vs. Alpha complex); signal is topological property, not filtration artifact.
- **Test design:**
  - Recompute H1 persistence using Alpha complex (Delaunay-based) on same scGPT data
  - Compare alpha-H1 vs. feature-shuffle null (same 20 replicates)
  - Compare alpha-H1 effect size to Rips-H1 effect size (should be similar)
  - Compute Spearman ρ(Rips_H1, Alpha_H1) across layer-seeds
- **Expected signal:** ρ > +0.7 and both Rips/Alpha significantly positive
- **Null:** Feature-shuffle
- **Value:** MEDIUM (validates that signal is topological, not filtration-dependent)
- **Cost:** MEDIUM (requires Alpha complex library; adds computational time)

#### **N08: Persistent Homology Under Geodesic vs. Euclidean Distance**
- **Hypothesis:** H1 persistence computed using geodesic distance on data manifold is stronger than Euclidean H1, indicating signal is manifold-intrinsic.
- **Test design:**
  - Compute manifold geodesics (diffusion distance, kNN geodesics) for each layer
  - Run Ripser on geodesic distance matrix (instead of Euclidean)
  - Compare geodesic-H1 vs. feature-shuffle null
  - Compute (geodesic-H1 z - euclidean-H1 z) delta across layers
- **Expected signal:** Geodesic-H1 is more strongly positive (Δz > +2) and/or lower p-values
- **Null:** Feature-shuffle (same null applies to both metrics)
- **Value:** MEDIUM (tests whether signal is distance-metric dependent; informs geometric interpretation)
- **Cost:** MEDIUM (geodesic computation is fast for kNN; Ripser rerun required)

---

### **TIER 3: MANIFOLD GEOMETRY & CROSS-MODEL (Medium priority)**

#### **N09: Intrinsic Dimension Layer Consistency**
- **Hypothesis:** Layers with higher H1 persistence have higher participation-ratio intrinsic dimension (robust across seeds and domains).
- **Test design:**
  - Compute participation-ratio dimension (Skutil) for each layer-seed in lung/immune/external-lung
  - Recompute H1 effect sizes (delta vs. feature-shuffle null) for all domains
  - Compute Spearman ρ(H1_delta, Participation_Ratio) per seed; Fisher-combine across seeds
  - Test if ρ is positive and significant within each domain
  - Extend to all 12 layers (not just top/weak), compute across-domain correlation
- **Expected signal:** ρ > +0.3, Fisher p<0.05 in all domains (differs from iter_0004, which was mixed)
- **Null:** Permute H1 deltas across layers; recompute ρ
- **Value:** MEDIUM (mechanistic insight: does topology couple to manifold dimensionality?)
- **Cost:** LOW (metrics already computed or easily computable)

#### **N10: Geneformer Clustering Coefficient Cross-Validation**
- **Hypothesis:** Geneformer embeddings show same kNN CC elevation vs. feature-shuffle null as scGPT; signal is model-agnostic.
- **Test design:**
  - Load Geneformer residual embeddings (if available in project; otherwise skip or fetch)
  - Run same kNN CC screen as H01 (k=10, PCA-20, 300 genes, 15 null replicates) on 2–3 Geneformer seeds/layers
  - Compare mean z-score to scGPT baseline (mean z=+27.64)
  - Compute Spearman ρ(scGPT_CC, Geneformer_CC) across layers
- **Expected signal:** Geneformer z > +5 (significant CC elevation) and ρ > +0.6 (consistent layer-ordering)
- **Null:** Feature-shuffle (same as iter_0002)
- **Value:** HIGH (validates signal is not scGPT-specific; critical for generalization claim)
- **Cost:** MEDIUM (depends on Geneformer artifact availability; if available, test is quick)

#### **N11: CKA Cross-Model Alignment (with Permutation Null)**
- **Hypothesis:** CKA (Centered Kernel Alignment) between scGPT and Geneformer full-layer residual embeddings is >0.6 and significantly higher than random-seed null.
- **Test design:**
  - Procure matched-gene Geneformer residual tensors (shape: layers × genes × dims)
  - For each layer: compute CKA(scGPT_layer, Geneformer_layer) on full gene set
  - Use permutation null: randomize gene order in Geneformer, recompute CKA (100 replicates)
  - Compute Fisher-combined p-value and effect-size z-score
  - Report mean CKA, min/max, count of layers with CKA > 0.6
- **Expected signal:** Mean CKA > 0.6, all p<0.05, effect z > +2
- **Null:** Permuted gene order
- **Value:** HIGH (tests if model-invariant structure exists; critical for biological validity)
- **Cost:** MEDIUM–HIGH (requires Geneformer residuals; CKA computation is fast)

#### **N12: Procrustes Alignment with Geodesic Matching**
- **Hypothesis:** After optimal Procrustes rotation, scGPT-Geneformer geodesic structures (pairwise distances) are highly correlated (ρ > +0.7).
- **Test design:**
  - Apply Procrustes to align scGPT/Geneformer in matched PCA(20) space
  - Compute pairwise geodesic distances within each model (kNN geodesic, 500 gene pairs sample)
  - Compute Spearman ρ between scGPT-geodesics and Geneformer-geodesics (after alignment)
  - Test if ρ is significantly higher than random-gene-permutation null
  - Report per-layer ρ values
- **Expected signal:** Mean ρ > +0.6, min ρ > +0.4, all p<0.05
- **Null:** Permuted gene pairing
- **Value:** MEDIUM (tests if manifold structure is aligned across models)
- **Cost:** MEDIUM (Procrustes + geodesic computation)

---

### **TIER 4: ALGORITHMIC SIGNATURES & ADVERSARIAL (Lower priority, exploratory)**

#### **N13: Layer-Wise Gradient Flow Signature**
- **Hypothesis:** Layers with high H1 persistence show distinct gradient-flow properties during training (loss-weighted contribution per layer, gradient variance, condition number).
- **Test design:**
  - Re-run model training with gradient tracking enabled
  - Compute layer-wise gradient variance, gradient-to-weight ratio, and condition number of Jacobian for each layer
  - Correlate layer-wise gradient metrics with H1 persistence effect sizes from iter_0003–0004
  - Permutation null: randomize layer ordering, recompute correlations
- **Expected signal:** Positive ρ between H1-delta and gradient metrics (ρ > +0.3, p<0.05)
- **Null:** Permuted layer ordering
- **Value:** LOW–MEDIUM (mechanistic insight into why topology emerges; requires training access)
- **Cost:** HIGH (requires model access + gradient instrumentation)

#### **N14: Adversarial-Sampling Null (Reverse-Rank Embedding)**
- **Hypothesis:** H1 persistence is robust even under adversarial null: reverse-ranking genes by feature correlation (to create anti-correlated, distance-preserving null).
- **Test design:**
  - For each layer: rank genes by correlation with top PC
  - Create "reverse-rank" null by permuting genes to anti-rank while preserving pairwise distances
  - Compute H1 persistence on reverse-rank null (20 replicates)
  - Compare observed H1 to reverse-rank null (expect observed > null)
  - Compute Fisher-combined p-value
- **Expected signal:** Observed H1 > reverse-rank null (positive delta), z > +3, Fisher p<0.05
- **Null:** Reverse-rank embedding
- **Value:** MEDIUM (if feature-shuffle succeeds but reverse-rank fails, suggests signal depends on feature correlations, not just distance preservation)
- **Cost:** MEDIUM (requires reverse-ranking algorithm; Ripser rerun)

#### **N15: Cell Differentiation Pathway Activation Coupling**
- **Hypothesis:** Genes in the same cell differentiation pathway (e.g., hematopoiesis) show shorter geodesic distances in residual embeddings; embedding geometry reflects developmental trajectories.
- **Test design:**
  - Curate cell differentiation pathways (KEGG, Reactome) for immune/lung
  - For each pathway: compute mean within-pathway geodesic distance vs. between-pathway distance
  - Measure effect size and Fisher-combined p across seeds-layers
  - Correlate pathway distance reduction with degree of pathway activation in cell state
- **Expected signal:** Within-pathway distance << between-pathway distance (z > +3); significant correlation with activation
- **Null:** Permuted pathway assignment
- **Value:** MEDIUM (links topology to functional cell state; biological validation)
- **Cost:** MEDIUM (requires pathway curation + activation scoring)

#### **N16: Layer-Depth Mechanistic Interpretation (Attention-Structure Correspondence)**
- **Hypothesis:** kNN CC and H1 persistence elevation are layer-invariant because embedding structure is independent of transformer depth; topology is learned by shallow layers and preserved downstream.
- **Test design:**
  - Examine attention weights across layers (if available or re-extract)
  - Compute attentional entropy, head-specialization metrics per layer
  - Correlate attention metrics with H1 persistence effect sizes
  - Test whether early-layer attention structure predicts late-layer topology
- **Expected signal:** Weak/non-significant correlation between depth and H1-delta; strong correlation between early-layer attention and CC elevation
- **Null:** Random layer permutation
- **Value:** LOW (exploratory mechanistic insight; does not directly test topology robustness)
- **Cost:** HIGH (requires attention analysis; large exploratory space)

---

## Summary Table

| ID | Name | Family | Value | Cost | Key Risk |
|---|---|---|---|---|---|
| N01 | TRRUST Co-Regulatory Network | Bio Anchoring | HIGH | MED | Incomplete TRRUST coverage for cell type |
| N02 | GO Functional Clustering | Bio Anchoring | HIGH | LOW | GO annotations may not reflect context |
| N03 | STRING Path Coupling | Bio Anchoring | HIGH | MED | Path distance may be weak proxy for regulation |
| N04 | Cell Ontology Enrichment | Bio Anchoring | HIGH | LOW | Cell-type assignment noise |
| N05 | Betti-0 Persistence | Alt. Topology | MED | LOW | May show same signal as H1 (redundant) |
| N06 | Betti-2 (H2) Persistence | Alt. Topology | MED | LOW | May be empty/sparse (negative result likely) |
| N07 | Filtration Robustness | Alt. Topology | MED | MED | Alpha complex may be expensive; expect consistency |
| N08 | Geodesic vs. Euclidean H1 | Alt. Topology | MED | MED | Geodesic may marginally improve; not transformative |
| N09 | Intrinsic Dimension Coupling | Manifold Geom | MED | LOW | Expected to match iter_0004 (mixed result) |
| N10 | Geneformer CC Cross-Val | Cross-Model | HIGH | MED | Geneformer artifacts may be unavailable |
| N11 | CKA with Permutation Null | Cross-Model | HIGH | MED–HIGH | Requires Geneformer residuals; high computational cost |
| N12 | Procrustes Geodesic Alignment | Cross-Model | MED | MED | May be post-hoc; depends on N11 success |
| N13 | Gradient-Flow Signature | Mechanistic | LOW–MED | HIGH | Requires model re-training access |
| N14 | Reverse-Rank Adversarial Null | Null Paradigm | MED | MED | May show surprising fragility; useful diagnostic |
| N15 | Differentiation Pathway Activation | Bio Anchoring | MED | MED | Requires pathway curation + activation scoring |
| N16 | Attention-Structure Correspondence | Mechanistic | LOW | HIGH | Large exploratory space; low specificity |

---

## Top 3 for Immediate Execution

### **CANDIDATE 1: HIGH-PROBABILITY DISCOVERY (N01 + N02 + N04 bundle)**
**"Biological Anchoring Trio"**
- **Rationale:** All three require only data-curation + geodesic/neighborhood computation (already available). If any succeed, directly validates that topological signal maps to biological structure; this is the critical missing link in the narrative.
- **Expected ROI:** HIGH. Success (>2/3 tests positive) immediately promotes entire topology family to "biologically grounded." Failure (0/3) indicates signal is statistical artifact.
- **Sequence:** N04 first (fast, low-cost), then N02 (moderate cost), then N01 (highest cost). Bail if N04 fails.
- **Estimated effort:** 3–4 hours of experimentation; 2 hours curation.
- **Output:** 3 new hypothesis screen JSON files + biological interpretation narrative.

### **CANDIDATE 2: HIGH-RISK / HIGH-REWARD (N10 + N11)**
**"Cross-Model Consistency Validation"**
- **Rationale:** If Geneformer residuals are available, CC and CKA tests directly answer "is signal model-agnostic?" If both succeed, topology family becomes major finding (general feature of transformer residual geometry). If both fail, entire project needs re-direction.
- **Expected ROI:** Transformative if successful; if unsuccessful, identifies that signal is model-specific quirk.
- **Sequence:** N10 first (quick, low-cost, diagnostic). If positive, proceed to N11. If negative, stop.
- **Estimated effort:** 1–2 hours if Geneformer available; 6+ hours if re-curation needed.
- **Output:** 2 hypothesis screen JSON files + model-consistency narrative.

### **CANDIDATE 3: CHEAP BROAD-SCREEN (N05 + N14)**
**"Topology Family Robustness Screen"**
- **Rationale:** Both are low-cost re-runs of existing code (Ripser H0, reverse-rank null). Together they test (i) whether signal is H1-specific or general topology property, and (ii) whether stronger adversarial null breaks signal. Low cost, high information density.
- **Expected ROI:** HIGH. Positive H0 confirms topology is robust across homology families. Positive reverse-rank null indicates signal is fundamental, not feature-correlation artifact. Either finding is publishable.
- **Sequence:** N05 first (1 hour), then N14 (2–3 hours). Both are parallel-capable if compute permits.
- **Estimated effort:** 3–4 hours total.
- **Output:** 2 hypothesis screen JSON files + robustness narrative.

---

## Execution Priority
1. **Immediate (next iteration):** Candidate 1 (N01+N02+N04) — Critical missing validation.
2. **Contingent on artifact availability:** Candidate 2 (N10+N11) — Decisive for generalization.
3. **Parallel or fallback:** Candidate 3 (N05+N14) — Cheap insurance against topology-family specificity.

If Candidate 1 shows >2/3 success, elevate N03 (STRING) to next iteration. If Candidate 2 succeeds, proceed immediately to N12 (Procrustes alignment) and Candidate 3.

