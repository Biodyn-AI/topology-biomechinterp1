# Brainstormer Structured Feedback: iter_0051

## Gate Status
`passed_min_research_gate: true` — 3 hypotheses tested, all with artifacts.

---

## Result Assessment

### H01: Residual TF Signal After Co-expression Regression — NEGATIVE
**Assessment**: This is a definitive result that closes a long-running open question. SV2-4 encodes co-expression structure; regulatory proximity is explained by the co-expression confound at 11/12 layers. The single surviving layer (L8, rbc=0.083) is small but real. Key implication: the paper must revise any claim that SV2-4 encodes TF-target regulatory proximity specifically. Closing this is progress — it sharpens the framing.

**What survives**: L8 is anomalous. Something specific about L8 produces regulatory-adjacent structure that partially survives co-expression regression. This is worth investigating mechanistically.

### H02: SV1-High Gene Identity — POSITIVE (strong)
**Assessment**: OR=0.108 at L0 is a large effect. The SV1 axis cleanly separates non-regulatory "background" genes from TF/target-enriched genes. This is the strongest positive result across the last several iterations and anchors a coherent biological narrative: the dominant spectral axis of gene embeddings is a regulatory membership axis. The layer-by-layer inversion of SV1 therefore reflects a reorganization of how the model positions regulatory vs. non-regulatory genes across depth.

**Next logical step**: Confirm housekeeping gene identity of SV1-high genes (Eisenberg list or ribosomal gene check). If housekeeping = high-SV1, the paper gets a clean biological interpretation: SV1 = housekeeping vs. regulatory axis.

**What's missing**: The executor correctly identified that 295/2039 named genes is a limitation. The unnamed gene mass (1744 genes) is exactly the SV1-high population — anonymous genes without TRRUST annotations. This is self-consistent but could be circularity if "unnamed = non-regulatory by definition." Need to test with an independent non-TRRUST gene set (housekeeping list) to rule out definitional circularity.

### H03: 0-dim PH on Circuit Genes at L8 — MIXED
**Assessment**: The proximity signal is real (pairwise KS p<1e-6, persistence KS p=0.007) but 0-dim PH does not detect additional topological structure beyond compactness. Circuit genes are in a tighter cloud, not in discrete islands. This is consistent with SV2-4 encoding continuous co-expression similarity rather than discrete regulatory modules. The mixed result is scientifically clean: it closes the "discrete topological clustering" question without falsifying the proximity finding.

**What was not tested**: H1 homology (loops). If regulatory circuits form cyclic structures (e.g., feedback loops are geometrically encoded as circles in representation space), H1 would detect them where H0 cannot. This is the highest-value next test.

---

## Direction Assessment

| Direction | Status | Rescue Potential |
|-----------|--------|-----------------|
| SV2-4 regulatory proximity (uncontrolled) | RETIRE | None — co-expression fully explains it |
| SV2-4 as co-expression encoder | KEEP | Core validated claim |
| SV1 biological identity | ACTIVE | Strong — housekeeping test pending |
| 0-dim PH on circuit genes | CLOSE for 0-dim | Rescue via H1 homology |
| SV5-7 regulatory signal (independent of co-expression) | UNTESTED | High potential — H-C from iter_0050 never executed |
| L8 layer specialization | NEW | Mechanistic question — why does L8 show residual regulatory structure? |

---

## Momentum Assessment
The iter_0051 work produced two definitive closures (H01 negative, H03 mixed) and one strong positive (H02). The biological picture is converging: SV1 = regulatory membership axis; SV2-4 = co-expression geometry. The next critical open question is whether any spectral direction beyond SV2-4 carries regulatory signal independent of co-expression — specifically SV5-7 (rbc=0.152 from iter_0050, untested for co-expression independence). This is the highest-value untested hypothesis from prior planning.

---

## Research Gate Concerns: None
All artifacts present, paper updated, code trace recorded. Proceed normally.
