# Next Steps — iter_0025

## Portfolio Status

### Strong positives to build on
1. **Multi-predictor validation (H01):** STRING + Dorothea + GO CC each independently significant (partial R², VIF<2). This is publishable evidence of biologically specific manifold organization.
2. **Layer-resolved timeline (H02):** GO_CC peaks L5, GO_BP peaks L7, suggesting functional biology encodes progressively deeper. Peak span = 7 layers across 5 anchors.
3. **Genomic null (H03):** AUROC=0.515 confirms all prior functional signals are NOT genomic artifacts.

### Retired directions
- TF activation hub centrality (H02 iter024): clean null
- Null sensitivity (repeated family): consistently negative

---

## Recommended Hypotheses for iter_0026

### H-A (HIGH PRIORITY): Cross-model validation of joint model
**Family:** cross_model_alignment | **Novelty:** new_method
**Rationale:** H01 showed STRING+Dorothea+GO_CC each independently predict scGPT geometry. Does Geneformer show the same partial R² structure? If both models show the same predictor ranking, this validates model-independent biological encoding.
**Method:** Load Geneformer word embeddings (20275×1152 from HuggingFace cache). Identify named gene subset. Compute pairwise L2 distances. Run same partial-Spearman R² model (STRING + Dorothea + GO_CC). Compare ranking vs scGPT.
**Cost:** Medium. Geneformer embeddings previously loaded in iter_0019.

### H-B (HIGH PRIORITY): Layer-specific biological content test
**Family:** module_structure | **Novelty:** new_method
**Rationale:** H02 shows GO_CC peaks early (L5), GO_BP peaks late (L7). Test if the "excess" signal in early layers for GO_CC is specifically explained by cellular compartment structure not shared with GO_BP.
**Method:** At each layer, partial out GO_CC from GO_BP Spearman and vice versa. Measure if early-layer GO_CC advantage (L0-L4 vs L8-L11) exceeds what's expected from their inter-correlation. Test statistical significance.
**Cost:** Low. Uses existing artifacts + small new computation.

### H-C (MEDIUM PRIORITY): Functional module membership test
**Family:** module_structure | **Novelty:** new_method
**Rationale:** STRING defines network communities (clusters of highly-interconnected genes). Test if genes in the same STRING community are geometrically closer than cross-community pairs.
**Method:** Build kNN graph on STRING pairs (score≥0.5). Use Louvain community detection. For each pair, assign same-community vs cross-community label. AUROC at each layer. Null: label shuffle.
**Cost:** Medium.

### H-D (MEDIUM PRIORITY): Partial R² trajectory — STRING vs regulatory anchors across layers
**Family:** manifold_distance | **Novelty:** new_method
**Rationale:** H01 shows STRING partial R² increases from L0 (0.00457) to L8 (0.00962) then decreases to L11 (0.00854). Dorothea partial R² also peaks mid-to-late. Test if this trajectory is significantly non-flat (quadratic peak) vs uniform.
**Method:** Fit quadratic regression to 12-point partial R² trajectory for each predictor. Test if quadratic coefficient is significant. Identify optimal encoding layer per predictor.
**Cost:** Low. Uses only h01_multi_predictor_joint_model.json.

### H-E (LOW COST SCREEN): Gene expression correlation control
**Family:** null_sensitivity | **Novelty:** new_method
**Rationale:** Co-expressed genes might be artifactually closer in scGPT embeddings due to shared cell-type context, not functional relationship. Test if STRING AUROC survives after controlling for expression correlation.
**Method:** From lung processed.h5ad, compute pairwise Pearson correlation of gene expression profiles (all cells). Add this as a 5th predictor in the H01 joint model. Check if STRING partial R² changes substantially.
**Cost:** Medium (requires loading processed.h5ad).

---

## Paper Update Priority
The H01 multi-predictor results are the keystone paper finding. The next iteration should extend this with:
1. Cross-model validation (H-A above)
2. Formal statistical test of the STRING trajectory (H-D above)
