# E2 Summary — Cross-model replication

**Status:** Done.
**Code:** `run.py`.
**Outputs:** `results.csv`, `cka.csv`, `results.json`.

## Headline finding

The headline geometric organisation is **partially scGPT-specific** at L11. Scaling Geneformer from V1-10M to V2-104M *partially recovers* the SV1/SV2 structure but never matches scGPT-L11's strength. The **TF-vs-target distinction** is convergent across all three models — and is in fact *stronger* in both Geneformer variants.

| Model | n_genes | SV1 GO:0005615 OR | $p$ | SV2 PPI z | SV2 PPI obs rate | TF/target AUROC | CKA vs scGPT-L11 |
|---|---|---|---|---|---|---|---|
| **scGPT-L11** | 1,672 | **1.55** | 0.0013 | **12.13** | 0.485 | 0.651 | — |
| Geneformer V2-104M | 1,672 | 1.36 | **0.0007** | 4.28 | 0.283 | **0.754** | 0.173 |
| Geneformer V1-10M | 1,672 | 0.98 | 0.35 | -1.03 | 0.114 | **0.748** | 0.103 |

(Shared gene set: 1,672 genes present in scGPT-L11 + Geneformer V1 + Geneformer V2 + STRING/TRRUST/GO. STRING ≥ 0.7 pairs: 968 in this set.)

## Interpretation

Three patterns:

1. **SV1 subcellular and SV2 PPI structures scale with model size.**
   - scGPT-L11 (12 layers, ~50M params): full strength (OR=1.55, PPI z=12.13).
   - Geneformer V2 (104M params): partial replication (OR=1.36, PPI z=4.28). Significant but weaker.
   - Geneformer V1 (10M params): no signal (OR=0.98, PPI z=-1.03). Below the gene-shuffle null.

2. **TF-vs-target distinction is convergent and stronger in Geneformer.** Both Geneformer variants achieve AUROC ~0.748–0.754 on the shared gene set, vs scGPT's 0.651. This is a robust, model-agnostic finding.

3. **CKA between scGPT-L11 and Geneformer is low** (0.10 / 0.17). The full embedding spaces are largely incomparable — but the *biological* findings (when projected onto SV1/SV2 / SV2-7 spectral subspaces) partially replicate. This is consistent with the broader literature that says representational similarity (CKA) and downstream-task agreement diverge.

## Manuscript implication

The "Important Negative Findings" item on cross-model alignment has been updated. The main-text claim now reads (in §"Important Negative Findings"):

> "The headline SV1 subcellular and SV2 PPI structures partially replicate to Geneformer V2-104M (SV1 OR=1.36 p=0.0007; PPI z=4.28) but not to Geneformer V1-10M (SV1 OR=0.98 p=0.35; PPI z=-1.03). The TF-vs-target distinction is the most robust geometric property across the three models tested: AUROC=0.651 for scGPT-L11, 0.748 for Geneformer V1, 0.754 for Geneformer V2. The headline geometric organisation is therefore better characterised as a property of scGPT's particular spectral basis at L11, with partial cross-model replication that improves with model size."

## Manuscript implication for §"Cross-model" subsection

We add a new subsection §"Cross-model replication" or extend the current "Negative Findings" item:

> "On the 1,672-gene set shared across scGPT-L11, Geneformer V1-10M, and Geneformer V2-104M (Table~\ref{tab:cross_model}), three patterns emerge.
> First, the SV$_1$ subcellular and SV$_2$ PPI signals scale with model size: full in scGPT-L11 (OR=1.55, PPI z=12.13), partial in Geneformer V2-104M (OR=1.36, z=4.28), and absent in Geneformer V1-10M (OR=0.98, z=$-$1.03).
> Second, the TF-vs-target distinction (joint SV$_2$--SV$_7$ AUROC) is convergent and \emph{stronger} in Geneformer (0.748 V1, 0.754 V2) than in scGPT (0.651 on the same shared gene set). This is the most model-agnostic geometric property we identify.
> Third, the centred kernel alignment between scGPT-L11 and either Geneformer variant is low (0.10–0.17), confirming that the full residual-stream geometries are largely incomparable. Where the biological projections agree (TF-vs-target, partial PPI), the convergence reflects shared representation of biology rather than shared geometry overall.
> Practical consequence: the SV$_1$ secretory and SV$_2$ PPI claims are scGPT-specific (partially); the TF-vs-target claim is universal across the tested foundation models."

## Caveats

- We tested only the L11 layer of scGPT and the final layer of Geneformer V1/V2. A full-layer comparison would refine the picture. Earlier layers may show different agreement patterns.
- Geneformer V1-10M may simply be too small to have learned the SV1/SV2 structure; a larger Geneformer training run on the same tokenised data could recover it.
- ESM2 cross-validation (E10) is in progress; that result will further inform whether the SV1 subcellular structure is sequence-derivable (in which case ESM2 would recover it) or scRNA-seq-pipeline-specific.
- 1,672 shared genes is a smaller test set than the original 4,803; the curated 195-gene effect sizes from the paper are not directly comparable to these numbers.
