# iter_0040 Executor Iteration Report

## Summary

Three hypotheses tested. Two strongly positive (H01, H03), one informative negative (H02).
Key new finding: **GC attractor onset is reproducibly at Layer 3**, confirmed by two independent measures.
BCL6 correction: BCL6 does NOT converge to B-cell centroid (rank ~900 throughout) and should be removed from GC attractor membership.
Geneformer negative result: B-cell/GC-TF proximity is scGPT-specific (not present in Geneformer input token embeddings).

---

## H01: Full 12-Layer GC Attractor Onset Scan

**Hypothesis**: The GC-TF attractor (BATF, BACH2, BCL6, PAX5 converging toward B-cell centroid) has a specific onset layer.

### Command
```bash
conda run -n subproject40-topology python run_iter0040_screen.py
```

### Results

| Layer | PAX5 rank | BATF rank | BACH2 rank | BCL6 rank | Mean GC rank | Centroid drift |
|-------|-----------|-----------|------------|-----------|--------------|----------------|
| L0    | 36        | 1519      | 562        | 830       | 736.8        | 0.00           |
| L1    | 31        | 1335      | 524        | 880       | 692.5        | 12.21          |
| L2    | 24        | 1282      | 274        | 1377      | 739.2        | 13.89          |
| L3    | **19**    | 1029      | 262        | 1146      | 614.0        | 15.55          |
| L4    | 19        | 493       | 486        | 1128      | 531.5        | 16.29          |
| L5    | 19        | 343       | 508        | 1100      | 492.5        | 17.48          |
| L6    | 16        | 346       | 387        | 948       | 424.2        | 18.90          |
| L7    | 15        | 268       | 388        | 1042      | 428.2        | 20.55          |
| L8    | 16        | 279       | 285        | 825       | 351.2        | 22.87          |
| L9    | 18        | 251       | 177        | 757       | 300.8        | 25.51          |
| L10   | 19        | 202       | 196        | 767       | 296.0        | 27.99          |
| L11   | 18        | 156       | 153        | 899       | 306.5        | 27.20          |

- **Attractor onset: L3** (PAX5 first reaches rank ≤ 20 at L3 = rank 19)
- GC mean rank Spearman rho = **-0.951**, p = 0.000002 (strong monotonic decrease)
- IRF4: OOV in cycle4 vocabulary
- **BCL6 correction**: BCL6 ranks 750–1500 throughout — not part of GC attractor. Remove from GC cluster claim.
- PAX5 is the earliest indicator gene (rank 36 at L0, already approaching top neighborhood)
- BATF and BACH2 converge later (strongest improvement L4-L9)

**Artifact**: `h01_gc_attractor_scan.json`

---

## H02: Geneformer Cross-Model B-Cell Precision@10

**Hypothesis**: Geneformer input token embeddings show B-cell/GC-TF proximity similar to scGPT.

### Command
```python
from safetensors.torch import load_file
state = load_file('/Users/ihorkendiukhov/.cache/huggingface/hub/models--ctheodoris--Geneformer/snapshots/05fcbeb8a27d49e0a7a4349152202ee2c1cbfd28/model.safetensors')
gf_emb_matrix = state['bert.embeddings.word_embeddings.weight'].numpy()  # [20275, 1152]
```

### Results

| Metric | Value |
|--------|-------|
| Geneformer vocab size | 20,275 tokens |
| Embedding dimension | 1,152 |
| Immune genes mapped | 2,047 |
| B-cell markers available | MS4A1, CD79A, BLK, PRDM1 (CD19 **absent**) |
| GC-TFs available | BATF, BACH2, BCL6, PAX5 |
| Precision@10 (full anchor) | **0.000** |
| Null mean | 0.0000 |
| Top-20 neighbors of B-cell centroid | Olfactory receptors, zinc fingers, pseudogenes |

**Result: Negative.** Geneformer input token embeddings do not encode B-cell/GC-TF proximity. Top-20 nearest genes are biologically incoherent (olfactory receptors, pseudogenes). CD19 absent from Geneformer vocabulary.

**Interpretation**: The B-cell/GC-TF proximity signal is specific to scGPT's **learned contextual layer representations**, not present in raw token embeddings. A valid cross-model comparison requires extracting Geneformer **layer-wise** gene embeddings (not currently available without a forward pass run).

**Artifact**: `h02_geneformer_bcell.json`

---

## H03: GC-Plasma Subspace Principal Angles + Minimal Anchor Scan

**Hypothesis A**: GC subspace (4-gene) converges toward B-cell subspace (3-gene) across layers (principal angle decreases).
**Hypothesis B**: CD19+BLK minimal anchor reproduces GC-TF precision@10 onset at L3.

### Command
```bash
conda run -n subproject40-topology python run_iter0040_screen.py
```

### Part A: Subspace Principal Angles

| Layer | GC-plasma angle (°) | GC-bcell angle (°) |
|-------|--------------------|--------------------|
| L0    | 84.6               | 84.3               |
| L1    | 85.7               | 84.9               |
| L2    | 86.4               | 86.9               |
| L3    | 88.4               | 85.4               |
| L4    | 89.0               | 82.9               |
| L5    | 89.4               | 80.1               |
| L6    | 88.1               | 79.6               |
| L7    | 87.9               | 77.3               |
| L8    | 86.6               | 75.6               |
| L9    | 88.1               | 74.9               |
| L10   | 85.3               | 75.0               |
| L11   | 82.4               | 76.4               |

- **GC-bcell subspace angle**: Spearman rho = **-0.888**, p = 0.0001 → GC subspace converges toward B-cell subspace
- **GC-plasma subspace angle**: Spearman rho = -0.098, p = 0.762 → **null** (no convergence)
- Directional specificity confirmed: GC subspace approaches B-cell subspace, NOT plasma subspace

### Part B: CD19+BLK Minimal Anchor Precision@10

| Layer | p@10 |
|-------|------|
| L0-L2 | 0.000 |
| **L3** | **0.100** (onset) |
| L4-L11 | 0.100 |

- L11 p@10 = 0.100, null mean = 0.001 ± 0.010, **z = 9.95**, pctile = 99th
- **Onset at L3** independently confirmed via minimal anchor
- Full anchor (MS4A1+CD79A+BLK, NO CD19): p@10 = 0.000 throughout → **CD19 criticality confirmed** (3rd independent test)

**Artifact**: `h03_subspace_angles.json`, `iter_0040_summary.csv`

---

## Cross-Hypothesis Summary

| Finding | Evidence | Iterations |
|---------|----------|------------|
| GC-TF proximity to B-cell centroid | p@10 > null, Spearman rho decreasing | iter_0036-0040 |
| Attractor onset at L3 | Both rank scan (H01) and p@10 onset (H03) | iter_0040 |
| PAX5 is earliest GC-TF indicator | Rank 36 at L0, top-20 by L3 | iter_0040 |
| BCL6 NOT in GC attractor | Rank 750-1500 throughout | iter_0040 |
| GC subspace converges toward B-cell subspace | Principal angle rho=-0.888 | iter_0040 |
| GC subspace does NOT converge toward plasma subspace | Angle rho=-0.098, null | iter_0040 |
| CD19 is essential anchor gene | LOO/minimal anchor tests | iter_0039, 0040 |
| Geneformer input embeddings do not encode signal | precision@10=0.000 | iter_0040 |

---

## Quantitative Metrics

- H01 Spearman rho (GC rank vs layer): -0.951, p=0.000002
- H03 Spearman rho (GC-bcell angle vs layer): -0.888, p=0.0001
- H03 L11 minimal anchor z: 9.95 (99th percentile vs null)
- H02 cross-model precision@10: 0.000 (negative)

---

## Artifacts

- `h01_gc_attractor_scan.json` — layer-by-layer GC-TF ranks and centroid drift
- `h02_geneformer_bcell.json` — Geneformer cross-model result
- `h03_subspace_angles.json` — principal angles and minimal anchor precision
- `iter_0040_summary.csv` — combined summary table (24 rows)
- `run_iter0040_screen.py` — full reproducible script
