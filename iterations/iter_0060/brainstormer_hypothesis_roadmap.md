# Brainstormer Hypothesis Roadmap — iter_0060
Date: 2026-02-23

---

## Retire / Deprioritize

| Direction | Reason | Status |
|-----------|--------|--------|
| Signed displacement projection (SV5-7) for edge AUROC | Two iterations; ceiling at 0.563/0.565; method exhausted | `retire_now` |
| Scalar pairwise distance (SV5-7) for edge AUROC | Same ceiling 0.565; original finding from iter_0059 | `retire_now` |
| Hub-TF degree confound hypothesis | Tested iter_0059, not explanatory | `retired` |
| eff_rank layer confound for edge prediction | Tested iter_0059, partial_r showed no residual | `retired` |
| SV5-7 binary betweenness permutation as FFL metric | NS at p=0.30; replace with continuous t_mean | `deprioritize` (use t_mean instead) |

---

## New Hypothesis Portfolio

### Cluster A: FFL / Circuit Geometry (Expand H03)

**A1 — FFL geometry with expanded motif set and multiple SV subspaces**
- Hypothesis: The FFL t_mean>0 signal (intermediate TF geometrically between master TF and target) is real but underpowered at N=22. Expanding via STRING TF-TF regulatory links to N≥50 FFLs and testing across SV1-4, SV2-7, SV5-7 will reveal which subspace encodes regulatory hierarchy most strongly.
- Test: Build directed graph from TRRUST (288 edges) + STRING (score≥700, TF-TF pairs only). Enumerate all A→B→C FFLs where A,B are TFs. For each subspace {SV1-4, SV2-7, SV5-7}: compute t_mean per layer, one-sample t-test vs 0. Primary metric: t_mean with 95% CI.
- Expected signal: t_mean>0 at L0-L6 in the best subspace, permutation p<0.05 with N≥50.
- Null: t_mean~0 at all layers; all subspaces equivalent.
- Value: **high** | Cost: **medium**

**A2 — Regulatory cascade chain geometry (beyond FFLs to length-3+ paths)**
- Hypothesis: If FFLs show geometric ordering (B between A and C), longer regulatory chains (A→B→C→D) should show monotone geometric progression along the A→D axis.
- Test: Enumerate directed paths of length 3 and 4 in TRRUST+STRING graph. For each path [v1,v2,v3,v4], compute projection parameters t2=(v2-v1)·(v4-v1)/|v4-v1|² and t3=(v3-v1)·(v4-v1)/|v4-v1|². Test ordering: t2<t3 (monotone chain). One-sample t-test on (t3-t2). Permutation test at peak layer.
- Expected signal: t3>t2 significantly, indicating monotone geometric embedding of regulatory cascades.
- Null: t3-t2 ~ 0, no ordering.
- Value: **high** | Cost: **medium**

**A3 — FFL geometry per triplet identity across layers (decomposition)**
- Hypothesis: The t_mean collapse at L8 is driven by specific FFL triplets, not uniform across all 22. Identifying which triplets lose geometric ordering at L8 will reveal which regulatory circuits are disrupted by the deep-layer representation shift.
- Test: For each of 22 FFLs, compute t_value per layer. Cluster triplets by t_value trajectory. Annotate clusters by TF family. Compare L0-L7 vs L8-L11 t_value per triplet. Fisher's exact: are bZIP or bHLH triplets enriched in the L8-collapse cluster?
- Expected signal: STAT4/RUNX1 triplets maintain ordering; bZIP-containing triplets collapse.
- Null: No family-stratified difference in collapse pattern.
- Value: **medium** | Cost: **low**

---

### Cluster B: L8 Transition Characterization

**B1 — Global intrinsic dimension change at L8**
- Hypothesis: L8 marks a structural transition in the full embedding manifold's intrinsic dimension, explaining both the bHLH flip and FFL collapse.
- Test: For each of 12 layers, compute two-NN intrinsic dimension estimator (Facco et al. 2017) on full nonzero gene embeddings [~1972, 512]. Plot d_intrinsic vs layer. Test: is there a significant jump or drop at L8?
- Expected signal: d_intrinsic drops or changes regime at L8 (more compressed/collapsed manifold in deep layers).
- Null: d_intrinsic changes monotonically or not at L8.
- Value: **high** | Cost: **low**

**B2 — TF/target intraclass variance trajectory across layers**
- Hypothesis: TF class tightens (lower intraclass variance in SV5-7) as layers increase, explaining the monotone AUROC increase; the rate changes at L8.
- Test: For each layer, compute mean within-class pairwise distance (TF vs target) in SV5-7. Compute ratio (between_class / within_class). Plot across layers. Test for inflection at L8.
- Expected signal: Ratio increases monotonically, with a rate change at L8.
- Null: Ratio flat or decreasing.
- Value: **medium** | Cost: **low**

**B3 — Layer-to-layer representation similarity (CKA) trajectory**
- Hypothesis: The L8 transition manifests as a discrete jump in representation similarity between consecutive layers (low CKA between L7-L8 compared to adjacent pairs).
- Test: Compute linear CKA between all pairs of consecutive layers {L0-L1, L1-L2, ..., L10-L11} using full gene embedding matrix. Plot similarity vs layer transition. Test: is L7-L8 CKA an outlier?
- Expected signal: L7→L8 CKA drop significantly larger than mean adjacent-layer drop.
- Null: CKA decreases monotonically across layers (as expected for deeper networks).
- Value: **high** | Cost: **low**

---

### Cluster C: HIF1A Biological Anchor

**C1 — HIF1A nearest-neighbor composition before/after L8**
- Hypothesis: HIF1A's 20 nearest neighbors in SV5-7 shift from TRRUST target genes at L0-L7 to HIF1A-regulated targets (VEGFA, EPO, CA9, LDHA, etc.) at L8-L11, reflecting a representation transition from "regulated gene" to "master regulator."
- Test: For each of 12 layers, find HIF1A's top-20 nearest neighbors in SV5-7. Annotate: (a) are they TFs or targets in TRRUST? (b) are they known HIF1A targets (HIF1A target gene list from MSigDB/GSEA)? Compute fraction in each category per layer.
- Expected signal: Fraction of known HIF1A targets in top-20 NN increases dramatically at L8-L11.
- Null: Neighbor composition stable across layers.
- Value: **high** | Cost: **low**

**C2 — HIF1A vs other bZIP TFs: trajectory divergence as biological signal**
- Hypothesis: HIF1A's trajectory diverges from bZIP TFs (which remain target-like) because HIF1A is both a target gene (regulated by oxygen/VHL) and a master TF. The model learns its dual role. bHLH family has insufficient representation (n=1) — test whether other genes with dual TF/target roles also show sign-reversal.
- Test: Identify all genes in the dataset that appear as BOTH TF source (row in TRRUST) AND target (column in TRRUST). Compute margin trajectory for these "dual-role" genes across layers. Compare to single-role TFs and single-role targets.
- Expected signal: Dual-role genes show intermediate or sign-reversing margin trajectories.
- Null: Dual-role gene margins follow TF class throughout.
- Value: **high** | Cost: **low**

---

### Cluster D: Topology / PH

**D1 — Persistent homology on TRRUST subgraph embedding (Vietoris-Rips, per layer)**
- Hypothesis: The TF regulatory subgraph embedded in SV5-7 has topological features (connected components, loops) that track layer — H0 betti decreases (merging clusters), H1 betti peaks at L2-L3 (matches AUROC peak), collapses in deep layers.
- Test: For each layer, extract SV5-7 coordinates of the 58 TRRUST TF genes. Compute Vietoris-Rips PH up to H1. Record birth/death of features. Compare H0 and H1 betti numbers across layers.
- Expected signal: H1 betti number peaks at L2-L3 where TF/target AUROC peaks; decreases in deep layers.
- Null: PH stable or monotone across layers.
- Value: **medium** | Cost: **medium**

**D2 — Local curvature (second-order geometry) at TF genes vs target genes**
- Hypothesis: TF genes occupy regions of higher manifold curvature (more curved local geometry) than target genes, providing a geometric explanation for the ~0.73 AUROC boundary.
- Test: For each gene, estimate local curvature using k-NN residual variance: fit PCA to k=20 nearest neighbors, curvature = 1 - explained_variance_ratio(1st PC). Compare mean curvature of TF genes vs target genes per layer. One-sided t-test.
- Expected signal: TF genes have significantly different curvature from target genes, explaining the TF/target boundary.
- Null: No curvature difference; boundary is linear.
- Value: **medium** | Cost: **low**

---

### Cluster E: Cross-Model / Biological Alignment

**E1 — STRING co-expression module membership vs SV cluster membership**
- Hypothesis: STRING functional clusters (co-expression modules, score≥700) correspond to geometric clusters in SV1-4, revealing that the model represents biological co-function as geometric proximity.
- Test: Retrieve STRING co-expression edges for the 4803 genes. Run Leiden clustering on STRING graph. For each SV subspace (SV1-4, SV5-7), run k-means (k=10). Compute adjusted mutual information (AMI) between STRING clusters and geometric clusters per layer. Compare to random shuffle.
- Expected signal: AMI is significantly above chance at mid layers (L3-L7), tracking biological co-function.
- Null: AMI near 0; geometric clusters not biologically interpretable.
- Value: **medium** | Cost: **medium**

**E2 — GO term enrichment gradient along SV1 axis**
- Hypothesis: The first singular vector of the embedding (SV1) captures a biological gradient (e.g., metabolic vs. immune, housekeeping vs. stress-response), which can be identified via GO enrichment at the top/bottom extremes of SV1.
- Test: For each layer, extract SV1 coordinates. Rank genes by SV1 value. Run GSEA/hypergeometric test on top-200 and bottom-200 genes against GO Biological Process. Report top enriched terms per layer.
- Expected signal: SV1 captures a stable biological axis (metabolic vs immune) that differs between early and late layers.
- Null: No coherent GO enrichment; SV1 is noise.
- Value: **medium** | Cost: **low**

---

### Cluster F: Mechanistic Motifs

**F1 — Autoregulatory loop geometry (TF→TF self or TF→itself via intermediate)**
- Hypothesis: Genes that appear as both source and target in TRRUST (autoregulated or near-autoregulated TFs) occupy geometrically central positions in the SV embedding.
- Test: Find TRRUST TFs where the TF appears in its own target list (direct autoregulation) or where TF→X→TF loops exist. Compute their centroid distance from TF class centroid and target class centroid. Compare to single-role TFs.
- Expected signal: Autoregulated TFs are significantly closer to the TF/target boundary (intermediate positions).
- Null: Autoregulated TFs cluster with all TFs.
- Value: **medium** | Cost: **low**

---

## Top 3 for Immediate Execution

### 1. High-Probability Discovery Candidate
**B3 — Layer-to-layer CKA trajectory**
- Hypothesis: L7→L8 CKA is an outlier drop, confirming L8 as a structural boundary.
- Rationale: This is fast (compute CKA for 11 consecutive layer pairs), directly tests the L8 boundary hypothesis that emerges from both H02 and H03, and would provide a model-level mechanistic anchor for the multiple L8 phenomena. If confirmed, this is a clean, novel result.
- Cost: **low** (CKA is a simple matrix operation on already-loaded embeddings)

### 2. High-Risk / High-Reward Candidate
**C2 — Dual-role gene margin trajectory**
- Hypothesis: Genes appearing as both TF source and target in TRRUST show sign-reversing margins across layers.
- Rationale: If true, this links the HIF1A finding (n=1) to a generalizable principle: the model separates regulatory roles dynamically across layers. This would be a novel mechanistic claim with direct biological interpretability. Risk: few dual-role genes may exist in the dataset (may have same n=1 problem).
- Cost: **low** (identify dual-role genes, extract margins from H02 artifact)

### 3. Cheap Broad-Screen Candidate
**A1 — FFL geometry with expanded motif set and SV subspace sweep**
- Hypothesis: N≥50 FFLs in the best SV subspace will give significant t_mean with permutation p<0.05.
- Rationale: Directly rescues H03 with the two identified weaknesses fixed (N and subspace). The signal is promising (p<0.0001 at L2 for t_mean). The test is already coded; just needs more FFLs (STRING) and a subspace loop.
- Cost: **medium** (needs STRING edge download and subspace loop, but core logic exists)
