# Executor Next Steps: iter_0035

## Key Findings This Iteration

1. **H01 (mixed):** B-cell community purity is NOT monotonically increasing across layers (rho=0.08, p=0.80). However, B-cell geometry is strong at L2–L3 (OR=12–16) before weakening mid-layers and returning at L11. The story is not "late-layer crystallization" but "early-embedding geometry + L11 binary collapse."

2. **H02 (negative):** PC2 and PC3 show no significant T-cell or myeloid enrichment (all p>0.14) on the 195-gene set. The positive PC1 pole is not organized by lineage type in lower PCs.

3. **H03 (positive, strongest control):** B-cell PC1 signal survives bootstrap null (empirical p=0.000, z=-3.06) and L2-norm regression (residual AUROC=0.214, p=0.002). The B-cell geometry is structural, not an expression-level artifact.

---

## Retired Directions (do not revisit without rescue rationale)
- PC2/PC3 T-cell/myeloid axes (195-gene set) — 2 negatives, too few myeloid markers
- STRING PPI L2 distance — retired (iter_0031)
- Dorothea sign-split test — retired (iter_0033)

---

## Highest-Priority Hypotheses for iter_0035

### A. Early-Layer B-cell Geometry Validation (high-priority, new_method)
**Family:** module_structure / topology_stability
**Rationale:** H01 shows OR=16 at L2, OR=12 at L3 — stronger than L11! Test whether the L2 community partition already captures B-cell identity. Run Fisher test on L2 partition specifically; also test ARI between L2 and L11 partitions to see if the same genes are being grouped together from the start.

### B. Confound Replication at L2 (high-priority, refinement of H03)
**Family:** null_sensitivity
**Rationale:** H03 proved the B-cell signal structural at L11. Replicate the bootstrap confound test at L2 (OR=16), the layer with the strongest signal. If it also passes there, the finding is layer-agnostic.

### C. Cross-Layer B-cell Gene Membership Stability (novel, split_robustness)
**Family:** topology_stability
**Rationale:** The B-cell community fluctuates across layers (identity shuffling). Track each gene's layer-wise community membership: which genes are consistently in the B-cell-enriched community across all layers? This produces a "core B-cell geometry gene set."

### D. Myeloid Signal in Full Vocab (novel, new_method)
**Family:** intrinsic_dimensionality
**Rationale:** H02 negative for 195-gene set (only 6 myeloid markers in vocab). The full 4803-gene vocab likely has 50+ myeloid markers. Test myeloid enrichment on PC1-PC3 using the full vocab projection.

### E. TF Regulatory Specificity Test (novel, module_structure)
**Family:** module_structure
**Rationale:** TFs are not enriched in either community (iter_0033). But TFs that specifically regulate B-cell genes (PRDM1, IRF4, PAX5, EBF1 targets) should be in Community 1. Test whether B-cell-regulatory TFs (with TRRUST links to B-cell markers) are in Community 1.

---

## Recommended Portfolio for iter_0035
1. **A (early-layer community validation)** — new_method, extends strongest finding
2. **C (cross-layer gene membership stability)** — novel, topology_stability
3. **D (full-vocab myeloid PC axis)** — new_method, intrinsic_dimensionality
