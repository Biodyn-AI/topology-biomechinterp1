# Executor Iteration Report — iter_0010

## Summary
Three hypotheses tested. All three positive. Key new finding: TRRUST TF-target pairs co-localize in SVD poles (SV2 emp_p=0.000), establishing the first regulatory-network → geometry link. SV2 also strongly encodes cytoskeleton (OR=20.35) and mitochondrion (OR=6.36) compartments across multiple layers. Zero annotation-density confounders found, strengthening prior positive claims.

---

## Command Trace

```bash
conda run -n subproject40-topology python \
  /Volumes/Crucial\ X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/iterations/iter_0010/run_iter0010_screen.py
```

Runtime: ~8 min (H01: 2 SVs × 12 layers × 9 compartments × 200 shuffles each).

Data inputs:
- `layer_gene_embeddings.npy` [12, 4803, 512] (scGPT residual stream)
- `cycle1_edge_dataset.tsv` (209 named genes)
- `gene2go_all.pkl` (GO annotations)
- `trrust_human.tsv` (9396 TF-target edges)

---

## H01: SV2/SV3 12-layer × 9-compartment scan

**Hypothesis**: Higher SVD components (SV2, SV3) of the 209-gene layer-wise embedding also encode GO compartment structure, as SV1 did (iter_0009).

**Method**: For each of 12 layers and 2 SVs (SV2, SV3), compute SVD of 209-gene mean-centered embedding. Take top-K=52 and bottom-K=52 genes. For each of 9 GO compartments, Fisher exact OR vs empirical null (N=200 gene-label shuffles; test both poles, take best).

**Results**:
- **61 significant cells (emp_p<0.05)** across SV2+SV3 × 12 layers × 9 compartments (total 216 cells).
- **SV2 cytoskeleton**: OR=20.35, emp_p=0.000, stable across L1–L11 (7 cytoskeletal genes all in one pole).
- **SV2 mitochondrion**: OR=6.36, emp_p=0.000, present at L1–L3 (top pole) then shifts to bottom pole L4+.
- **SV2 extracellular_vesicle**: OR=4.2–5.4, emp_p=0.000 at L4–L11.
- **SV3**: 15 significant cells; strongest hit nucleus OR=4.3 at L3, plasma_membrane OR=2.9 at L4.
- **Zero annotation density confounders** (Spearman r between n_ann and OR per layer: 0/18 tests p<0.05).

**Decision**: PROMISING (positive, null-controlled, multi-layer stable).

---

## H02: TRRUST TF-target co-pole test in SVD space

**Hypothesis**: TF-target gene pairs from TRRUST co-localize in SVD poles more than random gene pairs — i.e., regulatory interaction graphs are encoded geometrically.

**Method**: At layer 11, compute SVD of 209-gene embedding. For SV1 and SV2, define top-K=52 and bottom-K=52 poles. 326 TRRUST pairs have both genes in our 209-gene set. Measure co-pole rate (both genes in same pole). Null: N=1000 random same-size gene pairs.

**Results**:
- **SV1**: obs co-pole rate=0.163, null mean=0.123 ± 0.027, **emp_p=0.016** (1-sided).
- **SV2**: obs co-pole rate=0.206, null mean=0.122 ± 0.027, **emp_p=0.000** (1-sided).
- SV2 effect is larger and more significant than SV1.
- 326 pairs from TRRUST overlap our 209-gene dataset.

**Decision**: PROMISING (positive, null-controlled, SV2 result is highly significant at emp_p=0.000).

---

## H03: Spectral ratio profile + annotation density confounder check

**Hypothesis**: SV1/SV2/SV3 variance fractions evolve systematically across layers; enrichment results are not driven by annotation density.

**Method**: At each of 12 layers, compute SVD of 209-gene embedding and track var(SV1)/var(total), var(SV2)/var(total), var(SV3)/var(total), and SV1/SV2 magnitude ratio. Also: Spearman r between n_ann and best OR across layers per compartment × SV combination.

**Results — Spectral profile**:

| Layer | SV1 var% | SV2 var% | SV3 var% | SV1/SV2 ratio |
|-------|----------|----------|----------|----------------|
| 0     | 19.1%    | 5.7%     | 3.3%     | 1.83           |
| 3     | 29.6%    | 6.6%     | 3.3%     | 2.12           |
| 7     | 41.6%    | 8.5%     | 3.9%     | 2.21           |
| 10    | 69.0%    | 7.1%     | 3.1%     | 3.13           |
| 11    | 76.7%    | 6.4%     | 2.2%     | 3.46           |

- SV1 variance grows monotonically 19%→77%, confirming continued compression of biological variation onto SV1.
- SV2 peaks at L8 (9.9%) then drops; this coincides with the mito transient observed in iter_0009.
- SV3 stays flat ~3-4% throughout.
- SV1/SV2 ratio jumps notably at L9–L11 (2.51→3.13→3.46), matching the SV1 stability break at L10→L11 found in iter_0009.

**Confounder check**: 0 of 18 tested (SV2/SV3 × 9 compartments) showed Spearman r between n_ann and OR with p<0.05. Enrichment results are NOT confounded by annotation density.

**Decision**: PROMISING / NEUTRAL (validation pass; spectral profile adds new structural insight).

---

## Artifacts Generated

| File | Type | Contents |
|------|------|----------|
| `h01_sv2sv3_layer_compartment_scan.json` | JSON | Full 2×12×9 OR/emp_p table |
| `h02_trrust_copole_result.json` | JSON | TF-target co-pole test results |
| `h02_trrust_copole_null_sv1.npy` | NPY | N=1000 null co-pole rates for SV1 |
| `h02_trrust_copole_null_sv2.npy` | NPY | N=1000 null co-pole rates for SV2 |
| `h03_spectral_profile.csv` | CSV | SV1/SV2/SV3 var fractions × 12 layers |
| `h03_annotation_density_confounder.json` | JSON | Spearman r(n_ann, OR) per SV×compartment |
| `iter0010_results.json` | JSON | Master results summary |

---

## Key Quantitative Findings

1. **SV2 cytoskeleton**: OR=20.35, emp_p=0.000 — strongest compartment signal seen in any SV/layer.
2. **SV2 TF-target co-pole**: emp_p=0.000 — regulatory network topology is encoded in SV2 geometry.
3. **SV1 TF-target co-pole**: emp_p=0.016 — also significant but weaker.
4. **61 of 216 SV2/SV3 compartment tests significant** — both higher SVDs are biologically structured.
5. **Zero confounder hits** — annotation-density does not explain enrichment results.
6. **SV2 peaks at L8** (9.9% var) — temporally coincides with mito/ER processing transient.
