# E9 / A6 — the automated screen: provenance, claims ledger, multiplicity — independently verified

**Size.** 409 hypothesis records in 137 JSON files (iterations 1–138; iter_0001, iter_0003, iter_0139 have no
JSON; 3 records belong to an aborted duplicate run). Iterations 1–63 hold the **183** records the preprint
(Kendiukhov 2026, arXiv:2602.22247) reports. The JSON uses 8 families within the 183 (15 overall); the
preprint's "13 families" is its own regrouping (`S1_family_mapping_json_x_preprint.csv`).

**Matrices actually used (S2_hypotheses.csv).** Of the 183: **115 on the 4,803-gene Tabula Sapiens KIDNEY
matrix** (`cycle1_main`, described in the preprint as immune), 67 on genuine immune (`cycle4_immune_main`),
1 on a Gaussian null, 2 partly on Geneformer. Kidney-derived preprint claims: SV1 secretory/compartment axis,
STRING/PPI co-pole results, cell-type and gene-family clustering, GC–plasma angle, B-cell precision@10.
Immune-derived: effective rank 14.4×, SV1 93.4%, TwoNN, TF/target AUROC 0.744, edge AUROC 0.602, BATF/BACH2.

**Three procedural defects, each identifiable from the records.**
- **D1 tissue.** As above.
- **D2 never-observed genes in the SVD.** SV1 of the zero-row matrix correlates −0.9998 with the observed
  indicator. OR 6.37 (iter_0006) reproduces exactly there (21 vs 5 genes; all 22 never-observed named genes in the
  "extracellular" pole); on the observed-gene SVD: **OR 1.00, p = 0.60**. 14.4-fold and 93.4% come from iter_0046 on
  IMMUNE using a 2,000-row subsample with ~59% zero rows (reproduced: 14.3-fold; 53.5% → 93.1%). **The screen itself
  flagged this at iter_0047** ("dominated by zero vectors"; true collapse 4.9×); the preprint kept 14.4× and 93.4%.
- **D3 label misalignment** in 36 of the 183 (iterations 27, 29–39): scripts dropped blank lines of `gene_list.txt`
  before enumerating, so all 209 labels sat on rows 0–208 (0 of 209 on their own row). The preprint's "195 in-vocabulary /
  14 OOV" is this artefact; truth is 187 observed / 22 never observed. On true rows the GC–plasma angle trend reverses
  (77° → 94° as executed; 62.7° → 55.0° on true rows); B-cell neighbour clustering survives (z 6.3–8.3).

**Other mis-statements against the records.** Marker-AUROC depth slope has the wrong sign (record −0.036; preprint
"+0.036, increases with depth"); SV2–4 co-expression residual significant at 1/12 layers, not 0/12; "co-expression" in
iterations 50–53 is full-embedding cosine; edge AUROC 0.6017 reproduces with 1,523 negatives (not 2,351); the immune
label set behind 0.744 is 73 TFs / 222 targets (preprint: 51 / 144); STRING gradient p = 1.4e-24 is scipy's
t-approximation at ρ = 1, n = 5 (exact p = 0.0167).

**Bands were not prespecified (band_provenance.md, verbatim quotes).** SV2–4 assembled from single-axis STRING co-pole
tests on the kidney 209-gene subset (22 zero rows), iterations 13–19; under the brainstormer's stated criterion (z > 3 at
≥ 6/12 layers) **SV4 passed by z = 3.00003 and SV5 failed (2/12)**, and the executor ran at STRING ≥ 0.4 instead of the
specified ≥ 0.7. SV5–7 = "the next three" (iter_0049). SV2–7 formed after complementary layer profiles were seen
(iter_0055); its "≥ 0.72 at 10+ layers" threshold relaxed after 9/12 (iter_0056). The screen's own Bonferroni window scan
(iter_0063) found no surviving window.

**Crosswalk (geometry only; legacy_band_mapping.json).** kidney_A legacy full-matrix SVk ↔ observed SV(k−1) in 94.7% of
layer × direction cases (mean |cos| at exactly k−1: 0.956); legacy SV1 has no counterpart. Bands: SV2–4 → SV1–3
(cos² 0.990), SV5–7 → SV4–6 (0.938), SV2–7 → SV1–6 (0.966). On immune the screen used the observed basis (|cos| 1.000).

**Multiplicity.** 183 records: 112 extractable p, 18 stored as 0, Bonferroni 43 → 25 after flooring, BH 70 — reproduced
exactly by script. Headline findings surviving Bonferroni: 8 → 5 (floors) → 4 (exact gradient p) → **3** (dropping the
across-layer Stouffer p for TF/target). All 409: 130 extractable, 28 Bonferroni (0.05/130), 79 BH; excluding the aborted
run 127 / 27 / 76.

**Ledger.** 61 claims (`E9_claims_ledger.csv`): 32 recomputed (25 reproduce exactly or approximately, 2 partly, 5 do not);
27 provenance-only; 2 not recomputable (their iteration folders were overwritten by an earlier loop).
