# Executor Iteration Report: iter_0052

## Summary

Three hypotheses tested per brainstormer guidance:
- **H01** (Housekeeping enrichment in SV1-high): **INCONCLUSIVE** — only 7 housekeeping genes in the 2039-gene nonzero set, insufficient statistical power. Trend: 4 HK genes in SV1-low at L0 (DNAJB1, HSPA1A, HSPA1B, UBB), 0 in SV1-high (p=0.062 one-sided), suggesting HK genes may be SV1-low, opposite direction. At L5, 4 in SV1-high. Dataset coverage too sparse for firm conclusions.
- **H02** (H1 Betti curves on circuit genes at L8): **NEGATIVE** — Circuit genes have fewer H1 1-cycles (108 vs 119) and shorter mean lifetime (0.00143 vs 0.00178) than matched non-circuit genes. MW p=0.953 (circuit NOT greater). No evidence for regulatory feedback loops as geometric 1-cycles.
- **H03** (SV5-7 regulatory signal after co-expression residualization): **POSITIVE** — SV5-7 axes retain significant TRRUST proximity after regressing out embedding cosine similarity at layers L0-L2 (rbc=0.148, 0.119, 0.083 at L0/L1/L2; all p<0.01). Signal disappears at L3+. This contrasts with SV2-4, which only retained signal at L8 (rbc=0.083). Different spectral subspaces encode regulatory proximity at different processing stages.

---

## Command Trace

```bash
# Write experiment script v2
cat > /tmp/iter52_experiments_v2.py << 'PYEOF'
# ... full script ...
PYEOF

# Run all three experiments
conda run -n subproject40-topology python /tmp/iter52_experiments_v2.py
```

### Data sources
- Embeddings: `/Volumes/Crucial X6/.../cycle4_immune_main/layer_gene_embeddings.npy` — shape [12, 4941, 512]
- Gene names from: `TSP14.h5ad` (AnnData var_names, 4941 genes)
- Edge dataset: `cycle1_edge_dataset.tsv` — 589 TRRUST positive / 1523 negative pairs (nonzero-gene subset)
- SV1 loadings: `iter_0049/h01_sv1_vectors.npy` — [12, 2039]
- Housekeeping gene set: manually curated from Eisenberg 2003 (238 genes); 7 overlap with nonzero embedding set

---

## H01: Housekeeping Gene Enrichment in SV1-high

**Hypothesis**: SV1-high genes (large first singular vector loadings) represent "background"/housekeeping genes, explaining their depletion of TFs/targets.

**Method**:
1. Load gene names from TSP14.h5ad (4941 genes, 2039 nonzero).
2. Compile 238-gene housekeeping reference set from Eisenberg 2003 categories (ribosomal proteins, glycolytic enzymes, chaperones, proteasome, ubiquitin, RNA processing, etc).
3. Split nonzero genes into top/bottom 20% by SV1 loading at L0, L5, L11.
4. Fisher exact test: HK enrichment in SV1-high vs SV1-low.

**Results**:

| Layer | n_high | HK_in_high | n_low | HK_in_low | OR | p_greater |
|-------|--------|-----------|-------|-----------|-----|-----------|
| L0 | 408 | 0 | 408 | 4 | 0.000 | 1.000 |
| L5 | 408 | 4 | 408 | 0 | inf | 0.062 |
| L11 | 408 | 1 | 408 | 2 | 0.499 | 0.876 |

- Only 7 HK genes in nonzero set (PSMC4, UBB, HSPA1A, ALDOB, DNAJB1, HMBS, HSPA1B).
- At L0: 4 HK genes in SV1-**low** (DNAJB1, HSPA1A, HSPA1B, UBB), 0 in SV1-high — direction suggests HK genes may be low-SV1, not high-SV1.
- At L5: pattern reverses (4 in SV1-high).
- **Conclusion**: INCONCLUSIVE. Dataset has 98% unnamed genes; named housekeeping gene coverage is insufficient. The direction at L0 (HK in low-SV1) is consistent with HK genes having regulatory annotations, which would place them in SV1-low per the H02 iter_0051 finding.

**Decision**: inconclusive

---

## H02: H1 Betti Curves (Loop Topology) on Circuit Genes at L8

**Hypothesis**: Regulatory feedback circuits (TF→target→TF loops) leave traces as geometric 1-cycles (H1 homology) in the SV2-4 subspace at L8.

**Method**:
1. At L8, SVD of centered nonzero embeddings [2039, 512]; project to SV2-4 (3D).
2. Circuit genes: 295 nonzero TRRUST-annotated TFs/targets.
3. Matched non-circuit control: 295 randomly sampled non-circuit nonzero genes.
4. `ripser` (Vietoris-Rips) on each set, maxdim=1.
5. Compare H1 bar counts and mean lifetime (MW test).

**Results**:
- Circuit H1 bars: **108**, mean lifetime: 0.001429
- Non-circuit H1 bars: **119**, mean lifetime: 0.001778
- Total H1 lifetime: circuit=0.1543, non-circuit=0.2115
- MW test (circuit > non-circuit): p=0.953
- rbc = 0.129 (AGAINST circuit having more loops)

**Conclusion**: NEGATIVE. Circuit genes form *fewer* and *shorter* 1-cycles than non-circuit genes. H1 topology does not differentiate regulatory circuits from background in this subspace. Non-circuit genes may be more scattered/loopy due to their heterogeneous identity.

**Decision**: negative

---

## H03: SV5-7 Regulatory Signal After Co-Expression Residualization

**Hypothesis**: SV5-7 subspace (beyond SV2-4) encodes TRRUST regulatory proximity independent of co-expression.

**Method**:
1. At each of 12 layers: SVD of centered nonzero embeddings [2039, 512].
2. Project genes to SV5-7 (indices 4,5,6 of U matrix).
3. Compute pairwise Euclidean distances in SV5-7 for TRRUST pos (n=589) and neg (n=1523) pairs.
4. Also compute full-embedding cosine similarity (co-expression proxy) per pair.
5. OLS regression: regress cosine similarity out of SV5-7 distances.
6. Mann-Whitney test on residualized distances (TRRUST < neg).

**Results**:

| Layer | rbc_raw | p_raw | rbc_residual | p_residual |
|-------|---------|-------|-------------|-----------|
| L0 | **0.205** | <0.001 | **0.148** | <0.001 |
| L1 | **0.186** | <0.001 | **0.119** | <0.001 |
| L2 | **0.163** | <0.001 | **0.083** | 0.0015 |
| L3 | 0.107 | <0.001 | 0.022 | 0.218 |
| L4 | 0.105 | <0.001 | 0.006 | 0.420 |
| L5 | 0.048 | 0.042 | -0.079 | 0.998 |
| L6 | 0.008 | 0.393 | -0.095 | 1.000 |
| L7-L11 | negative | >0.5 | negative | >0.99 |

**Key finding**: SV5-7 retains co-expression-independent regulatory proximity signal at L0-L2 (rbc=0.148/0.119/0.083). This contrasts with SV2-4, which lost signal everywhere except L8. SV5-7 appears to encode regulatory relationships in early layers before being washed out by deeper processing.

**Decision**: promising

---

## Artifacts Generated
- `h01_housekeeping_sv1.json` — HK enrichment by layer (H01)
- `h02_h1_betti_circuit_L8.json` — Ripser H0/H1 diagrams for circuit vs non-circuit (H02)
- `h03_sv57_regulatory_signal.csv` — SV5-7 raw and residualized rbc/p by layer (H03)
