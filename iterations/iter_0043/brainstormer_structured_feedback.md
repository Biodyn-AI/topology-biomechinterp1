# Brainstormer Structured Feedback — iter_0043

## Gate status
`passed_min_research_gate: true` — 3 hypotheses screened, artifacts present, paper updated.

---

## Findings assessment

### H01: GC repressor anti-convergence (PRDM1) — NEGATIVE
The prediction that PRDM1 diverges from B-cell centroid is cleanly falsified. PRDM1 rho=-1.000, tightest pair at L11 is BACH2-PRDM1 (3.940). This is the third consecutive negative in the anti-convergence family (H02_iter0042 also negative). **Retire the repressor-divergence hypothesis family.**

The surprise positive finding — "GC circuit unity" — is genuinely interesting: activators and repressors co-locate in the same embedding region. This implies scGPT encodes the GC regulatory circuit as a coherent spatial cluster rather than encoding regulatory logic (activation vs repression) as geometric opposition. This is a meaningful mechanistic statement worth following up.

### H02: BCL6 divergence from B-cell centroid — NEGATIVE
Confirmed for the second time (H02_iter0042 lineage). BCL6 rho=-0.993 (converges). The rank anomaly at L2 (654→1311) is an interesting secondary signal co-occurring with peak metabolic cluster density. This is worth one targeted test but does not rescue the divergence hypothesis.

### H03: T-cell subtype ID stratification — POSITIVE (clean)
Effector CD8 ΔID=-45.8 vs general T-cell ΔID=+26.2. Dissociation = 72 units. This is the cleanest new positive in several iterations. The interpretation is biologically compelling: functional specialization → low dimensionality; pan-lineage identity → high dimensionality. This direction has strong generalization potential.

---

## Direction health

| Direction | Status | Reason |
|-----------|--------|--------|
| GC repressor anti-convergence | **RETIRE** | 3+ consecutive negatives, clearly falsified |
| BCL6 B-cell divergence | **RETIRE** | 2 explicit negatives, concept exhausted |
| BCL6 metabolic cluster (specificity) | **Saturating** — one more test max | Already confirmed ~90x enrichment; diminishing returns |
| B-cell GC attractor / ID compression | **Consolidating** — extend, don't repeat | Core finding well-established; extend to other lineages |
| Effector CD8 compression | **NEW LEAD — expand aggressively** | One clean positive; strong generalization prediction |
| T-cell activation circuit attractor | **NEW — explore** | Untested analog to B-cell GC attractor |
| PRDM1/plasma paradox | **NEW — one test** | PRDM1 at L11 near GC cluster; plasma markers (JCHAIN/SDC1) may diverge |

---

## Key open mechanistic questions

1. Is "functional specialization → ID compression" a general law, or specific to effector CD8?
2. Does the T-cell lineage have an attractor analog (activation circuit: CD28/LAT/LCK/ZAP70)?
3. What is the geometry inside the GC circuit unity cluster — spherical or elongated? Do activators and repressors at least form distinct sub-clusters within their shared neighborhood?
4. Does the B-cell GC attractor pattern generalize to the plasma cell program (JCHAIN/SDC1)?
5. Does BCL6 rank anomaly at L2 quantitatively track with BCL6 metabolic cluster density?
