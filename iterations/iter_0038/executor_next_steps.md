# Executor Next Steps — iter_0038

## Completed this iteration
- H01: GC/plasma proximity divergence across layers — PROMISING (plasma z drifts +0.55 to positive by L11)
- H02: Directional drift — INCONCLUSIVE (GC-plasma orthogonality at L11 is novel but drift doesn't align with plasma axis)
- H03: BLOCKED by vocabulary limitations

## Priority next actions (iter_0039)

### 1. Cross-model validation — HIGH PRIORITY
The B-cell specificity signal (z=7.55) and GC-TF geometry (BATF/SPIB/BACH2 at 86-96th pctile) need to be replicated in Geneformer. This is the main remaining gap for the paper.
- Load Geneformer embeddings (if available)
- Run same precision@10 analysis for B-cell markers
- Compare GC-TF proximity pattern

### 2. Extended vocabulary scan — MEDIUM
Use `cycle2_maxgenes1024` embeddings to check if PAX5/EBF1/BCL6 are included.
- If yes: master B-cell TF proximity test (H03 rescue)
- Also check for larger NK/myeloid marker sets

### 3. Differentiation axis geometry — MEDIUM
H02 found GC-plasma angle becomes 94° at L11. Two follow-ups:
- What direction DOES the B-cell centroid drift toward? Find nearest gene to the drift endpoint.
- Test if the GC-TF subspace and plasma-TF subspace are truly orthogonal (SVD on 2-gene covariance).

### 4. Retirement status
- intrinsic_dimensionality: retired (repeated negatives)
- module_structure: retired (repeated inconclusives)
- STRING PPI Euclidean: retired
- T-cell as primary target: retired
- manifold_distance (B-cell): KEEP — still most productive direction

## Vocabulary limitation note
The cycle1_main vocabulary (195 in-vocab genes) lacks PAX5, EBF1, BCL6, and most NK/myeloid markers.
Switching to cycle2_maxgenes1024 should be done in iter_0039 for H03-type screens.
