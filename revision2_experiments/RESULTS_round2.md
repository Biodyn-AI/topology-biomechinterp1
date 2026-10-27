# Round-2 canonical rebuild — results

All numbers computed under `CONVENTION.md`. Sources: `rebuild_headline.py` ->
`headline_rebuilt.csv` / `rebuild_headline.log`; `sensitivity.py` -> `sensitivity_results.txt`.
Permutations: N = 100,000 (co-pole), N = 10,000 (AUROC). Empirical p = (b+1)/(N+1).

## 1. Every headline claim, before and after

| claim | manuscript | round-2 canonical | verdict |
|---|---|---|---|
| Compression, "effective rank" | 23.6 -> 1.6, **14.4x** | effective rank 383 -> 133, **2.9x** | **mislabelled statistic** (see 2) |
| Participation ratio | 58 -> 9.5, 6.1x | 59.5 -> 10.3, **5.8x** | **reproduces** |
| SV1 variance fraction (full vocab) | 53.7% -> 93.4% | 57.2% -> 93.2% | **reproduces** |
| Compression is learned (feature shuffle) | rank rises 17.6x | 133 -> 361 under shuffle | **reproduces** (direction and force) |
| Depth trend | Spearman rho = -1.000 | rho = -1.000, p ~ 0 | **reproduces** |
| SV1 extracellular @L11 | OR = 6.37, p = 2.6e-4 | **OR = 0.46, p = 0.13** (direction reversed) | **fails** |
| SV1 extracellular, best layer | -- | L8: OR = 2.41, p = 0.031 | weak, does not survive Bonferroni |
| SV2 TRRUST co-pole | 0.206 vs 0.122, "significant at all 12 layers", survives Bonferroni | 0.220 vs 0.141, best p = 1.2e-2; **5/12 at p<0.05, 0/12 at Bonferroni** | **fails Bonferroni** |
| SV2 STRING co-pole | 0.226 vs 0.124, 12/12 layers | 0.246 vs 0.141 best, 8/12 at p<0.05, **3/12 Bonferroni** | partially survives |
| SV3 STRING co-pole | 0.198, 12/12 | 0.297 best, 11/12 at p<0.05, **9/12 Bonferroni** | **survives, strongest axis** |
| SV4 STRING co-pole | 0.158, 7/12 | 3/12 at p<0.05, 2/12 Bonferroni | weak |
| TF/target AUROC SV2-SV7 | mean **0.744**, max 0.789 @L3 | mean **0.612**, max 0.640 @L11 (null 0.496) | **effect halves** |
| TF/target significance | "p<0.01 all 12 layers, surviving Bonferroni" | 12/12 at p<0.05; **0/12 individually at Bonferroni**; Stouffer across layers **Z = 8.69, p = 1.8e-18** | **survives, but only in combination** |

## 2. Why "14.4-fold" is wrong

Methods 4.7 defines effective rank as the exponential of the Shannon entropy of the normalised
singular-value spectrum. Computed that way it falls **2.9x** (383 -> 133) on the full matrix and 1.3x
on the curated set — never 14.4x. The manuscript's "23.6 -> 1.6" is closely matched by the
**participation ratio on the curated 209-gene set** (21.2 -> 1.68, 12.6x), and its separately reported
participation ratio "58 -> 9.5" is matched by the **participation ratio on the 187-gene set**
(59.5 -> 10.3, 5.8x). The paper reports the same quantity twice on two different gene sets, calling
one of them by the other's name.

## 3. Why the TF/target AUROC drops from 0.744 to 0.612

Two separable causes, both convention choices:

| convention | AUROC (mean over 12 layers) |
|---|---|
| subset-SVD + self-similarity included (round 1) | ~0.76 |
| subset-SVD, leave-one-out | 0.738 |
| **full-SVD, leave-one-out (canonical)** | **0.612** |

* **Self-inclusion** (+0.02): each TF's own similarity to itself (= 1) enters its own score and no
  target's, biasing TFs upward by 1/\|TF\|.
* **Subset-SVD** (+0.13): fitting the singular vectors to the same 187 genes that are then classified
  in that basis. The axes are fitted to the evaluation set. Computing the SVD on the full 4,803-gene
  matrix — so the axes are a property of the representation rather than of the test set — and then
  restricting removes this, and is the only non-circular version of the test.

The claim survives in direction and significance (null 0.496; Stouffer p = 1.8e-18) but at roughly
half the reported effect size.

## 4. What the permutation upgrade changes

At N = 100,000 the p-floor is 1.0e-5, comfortably below alpha = 2.73e-4, so Bonferroni statements are
now meaningful. Three co-pole tests reach the floor genuinely (b = 0 of 100,000). Under the round-1
counts (N = 50-1,000, floor 2.0e-2 to 1.0e-3) **no permutation test could have cleared that threshold
at all**, which is the arithmetic error documented in `FINDING_01_reproducibility.md`.

The honest Bonferroni tally over the headline set is now:

* **survives**: SV3 STRING co-pole (9/12 layers), SV2 STRING co-pole (3/12), SV4 STRING (2/12),
  spectral depth trend (analytic), TF/target AUROC via Stouffer combination
* **fails**: SV1 extracellular enrichment (p = 0.13 at L11, reversed direction),
  SV2 TRRUST co-pole (best p = 1.2e-2)

## 5. Net effect on the paper's thesis

The thesis — *what survives compression is what is fixed below the level of transcription;
the relational fact does not survive any interpretable low-rank summary* — **is not damaged, and in
one respect is strengthened**:

* the **protein-interaction** axis is the robust one (SV3 especially), surviving a 100,000-permutation
  gene-label-shuffle null at 9 of 12 layers;
* the **TF/target** distinction survives but at half the claimed effect;
* the **subcellular-localisation** claim, which the abstract leads on, **does not survive** in the
  canonical basis and must be retired or heavily reframed;
* the **regulatory-edge** claim was already conceded as weak and remains so.

The paper that follows from this is narrower and more defensible than the current draft: compartment
drops out, physical interaction and regulatory class survive, and the relational fact still does not.
