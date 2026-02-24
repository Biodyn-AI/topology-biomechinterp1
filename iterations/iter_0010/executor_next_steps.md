# Next Steps — iter_0010 → iter_0011

## Priority Queue (ordered by expected information gain)

### P1 — TRRUST co-pole × 12 layers + edge-type stratification (new_method, family=module_structure)
**Why**: H02 found emp_p=0.000 for SV2 co-pole at layer 11 with 326 TRRUST pairs. Need to:
- Track co-pole rate across all 12 layers to find peak layer.
- Stratify by edge type (Activation vs Repression) — does one direction dominate?
- Test STRING PPI pairs for the same co-pole enrichment (convergent evidence from physical interactions).
**Target**: If activation edges show stronger co-pole than repression, it implies signed regulatory geometry.

### P2 — SV2 cytoskeleton gene characterization + Spearman stability (refinement, family=intrinsic_dimensionality)
**Why**: OR=20.35 (emp_p=0.000) for cytoskeleton in SV2 is the strongest compartment signal seen. Only 7 cytoskeletal genes in our set.
- Which 7 genes? Are they consistently in the same SV2 pole across all 12 layers?
- Compute Spearman r for SV2 rank stability (like H03 of iter_0009 did for SV1).
- Cross-check with TRRUST: are cytoskeletal genes TF targets of the same regulators?

### P3 — Cross-model SV alignment: scGPT SV2 vs Geneformer SV2 (new_family, family=cross_model_alignment)
**Why**: iter_0004 showed CKA alignment for SV1. SV2 alignment across models would establish that compartment structure is model-independent.
- If Geneformer data is available: compute Geneformer SVD, align SV2 axes, measure CKA or Procrustes distance.
- If Geneformer unavailable: fallback = bootstrap SV2 compartment stability (resample genes).

### P4 — SV2 spectral peak at L8: mito transient mechanistic test (refinement, family=module_structure)
**Why**: SV2 variance peaks at L8 (9.9%), coincident with mito transient (iter_0009). But L8 may also be where cytoskeletal signal peaks.
- Decompose SV2 variance at L8 vs L11: which compartments drive the L8 peak?
- Rank cytoskeletal genes at L7, L8, L9 to see if SV2 represents mito OR cytoskeleton at each layer.

## Retired Directions (do not revisit without rescue rationale)
- Persistent homology rewiring-null (4 consecutive negatives)
- GO BP cosine clustering (structurally underpowered)

## Data Availability Notes
- TRRUST: available at `/Volumes/Crucial X6/.../single_cell_mechinterp/external/networks/trrust_human.tsv`
- scGPT embeddings: available at `layer_gene_embeddings.npy [12, 4803, 512]`
- Geneformer: check for availability in `subproject_38` outputs
- STRING PPI: check for local file or download
