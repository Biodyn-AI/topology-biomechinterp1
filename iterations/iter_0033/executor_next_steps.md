# Executor Next Steps: iter_0034

## Summary of iter_0033 key findings
- Community 1 = B-cell (all 9 B-cell markers; OR=10.60, p=0.008)
- Community 0 = non-B-cell (T-cell/myeloid/other)
- PC1 negative pole = B-cell, stable across all 12 layers (p<0.01 at every layer)
- L11 community partition stable at k=15–30 (ARI=0.77–0.90)
- Community partition is qualitatively different at L11 (2 communities) vs L0–L10 (3–4 communities)

## Priority 1 (HIGH): Refine B-cell vs T-cell/myeloid community identity
**Rationale:** Community 1 = B-cell confirmed. What is Community 0? Is it T-cell + myeloid mixed, or is there substructure?
- Expand Community 0 analysis: test T-cell (n=8), myeloid (n=6), NK cell markers in Community 0
- Test if k=15 partition (3 communities) shows T-cell/myeloid separation within Community 0

## Priority 2 (HIGH): Cross-model validation — Geneformer PC1 B-cell polarity
**Rationale:** We have a clean finding (B-cell at negative PC1 pole, stable across layers). Now test if Geneformer embeddings show the same axis.
- Load Geneformer embeddings if available
- Compute PC1 on 195 in-vocab named genes, test B-cell marker enrichment at one pole
- Compute CKA or cosine alignment between scGPT and Geneformer PC1 directions

## Priority 3 (MEDIUM): Community dynamics — when does B-cell community emerge?
**Rationale:** L11 uniquely has 2 communities; L0–L10 have 3–4. The B-cell community presumably consolidates progressively.
- Track B-cell community purity across layers: at each layer, find the community with highest B-cell enrichment and compute its Fisher p-value and OR
- Test: does B-cell community become purer monotonically?

## Priority 4 (LOWER): Regulatory geometry × community membership
**Rationale:** Dorothea high-conf AUROC decays L0→L11 (iter_0032). Now test: are TF-target pairs more likely to be in the same community at each layer?
- For 205 Dorothea high-conf pairs in 195-gene vocab: compute fraction of pairs that are community co-members per layer
- Compare to random pairs (Fisher exact or AUROC over community co-membership)

## Retirements
- STRING PPI → L2 distance: definitively retired (multiple negative results)
- Dorothea activation vs repression: blocked (no sign/mor data in available files)
- intrinsic_dimensionality pure PR: finding established; no new hypothesis

## Artifacts to produce in iter_0034
- `h01_community0_subcommunity.json` — k=15 subcommunity within Community 0
- `h02_bcell_community_emergence.json` — B-cell community purity across layers
- `h03_regulatory_community_co_membership.json` — TF-target community co-membership per layer
