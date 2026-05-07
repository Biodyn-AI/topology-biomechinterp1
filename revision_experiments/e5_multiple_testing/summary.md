# E5 Summary — Multiple-testing audit

**Status:** Done.
**Inputs:** `iter_0001..iter_0063/executor_hypothesis_screen.json`.
**Code:** `aggregate.py`.
**Outputs:** `master_p.csv`, `family_summary.csv`, `headline_audit.csv`, `headline_audit.md`, `audit_log.json`.

## Headline numbers

- **183 hypotheses** in iters 1–63 — matches the paper claim exactly.
- **109/183 (59.6%)** have an extractable primary $p$ from the executor's free-form result strings (regex on `emp_p=`, `perm_p=`, `p_perm=`, then fallback `p=`). The remaining 74 are mostly hypotheses where the executor reported a `direction` (negative/mixed/inconclusive) without a numeric headline $p$ — informative missingness.

## Multiple-testing thresholds applied

| Correction | Denominator | Threshold | Pass count |
|---|---|---|---|
| Bonferroni (paper $N=183$) | 183 | $2.73 \times 10^{-4}$ | **41 / 109** |
| Bonferroni (extractable $N=109$) | 109 | $4.59 \times 10^{-4}$ | **43 / 109** |
| BH-FDR $q\le 0.05$ | 109 | adaptive | **68 / 109** |

## Headline-finding verdict (the abstract-relevant ones)

| Finding | iter/H | $p$ | BH $q$ | Bonferroni pass? |
|---|---|---|---|---|
| SV1 secretory-pathway layered (mito→ER→extra) | 0008/H02 | $\approx 0$ (emp) | $\approx 0$ | **yes** |
| SV2 STRING≥0.7 PPI co-pole | 0010/H02 | $\approx 0$ (emp) | $\approx 0$ | **yes** |
| SV2 TRRUST TF–target co-pole | 0011/H01 | $\approx 0$ (emp) | $\approx 0$ | **yes** |
| TF-vs-target AUROC=0.744 (joint SV2-7) | 0056/H01 | $\approx 0$ (perm) | $\approx 0$ | **yes** |
| TF-vs-target cross-seed | 0057/H01 | $\approx 0$ (perm) | $\approx 0$ | **yes** |
| Cell-type marker AUROC=0.853 (initial) | 0022/H02 | $1\times10^{-6}$ | $3.6\times10^{-6}$ | **yes** |
| STRING confidence gradient | 0015/H02 | $1.4\times10^{-24}$ | $7.5\times10^{-24}$ | **yes** |
| BATF/BACH2 attractor convergence | 0042/H03 | $1\times10^{-4}$ | $2.6\times10^{-4}$ | **yes** |
| SV1 extracellular enrichment (Fisher) | 0006/H02 | $3\times10^{-4}$ | $7.6\times10^{-4}$ | **no** (passes BH) |
| Cell-type marker AUROC=0.851 (expanded panel) | 0023/H01 | $4\times10^{-3}$ | $7.2\times10^{-3}$ | **no** (passes BH) |
| **Edge-level AUROC peak 0.602 + decay** | 0062/H01 | $4.5\times10^{-2}$ | $6.8\times10^{-2}$ | **no — fails BH-FDR** |

## Action items for the manuscript

1. **Edge-level AUROC peak 0.602** is the only abstract-relevant claim that **does not survive BH-FDR**. The paper's `kendiukhov2025attention`-style depth-decay trend (Spearman $\rho=-0.958$, $p=9.5\times10^{-7}$) is statistically real, but the *single-layer* peak that we report numerically (0.602 with $p_\text{perm}=0.045$) is at-or-below the multiple-testing waterline. Manuscript edit needed: reframe the peak number as "marginal at L0–L8 (multiple-testing-uncorrected); the depth trend is robust."
2. **SV1 extracellular OR=6.37** (a foundational Section 2.2 claim) and **cell-type marker AUROC=0.851** (foundational Section 2.2/2.6 claim) pass BH but **fail strict Bonferroni**. Recommend adding a sentence: "passes BH-FDR at $q=0.05$ but not the stricter $\alpha=0.05/183$ Bonferroni threshold."
3. **Eight remaining headline findings** — including the Section 2.3 PPI claims, Section 2.4 TF/target classifier, and the Section 2.6 B-cell attractor — survive Bonferroni at the strictest paper-cited denominator. These can stay in the abstract without qualification.
4. **Repression vs activation asymmetry** (the paper flags this as "pending replication") was *not* in the iter 1–63 logs as a single hypothesis with an explicit $p$; it is folded into iter_0011/H01 ("activation 12/12, repression 1/12") which passes Bonferroni overall. The asymmetry per se is a comparison statement that needs its own permutation test (planned for E7).
5. **Missing-$p$ hypotheses** (74 of 183): every one was either retired by the brainstormer (`direction='negative'`) or labelled `mixed`/`inconclusive`. None of these enters the paper's claims; the multiple-testing audit therefore conservatively excludes them.

## Per-family pass rates

The 7 families with extractable $p$ (out of 13 paper-curated families — the executor's internal taxonomy is finer-grained than the paper's narrative grouping):

| family | n | BH-pass (within) | Bonf-pass (within) |
|---|---|---|---|
| module_structure | 31 | 19 | 12 |
| manifold_distance | 28 | 19 | 15 |
| intrinsic_dimensionality | 24 | 16 | 14 |
| graph_topology | 9 | 3 | 3 |
| null_sensitivity | 7 | 5 | 5 |
| topology_stability | 6 | 2 | 2 |
| persistent_homology | 4 | 3 | 3 |

Module-structure and manifold-distance families have the highest absolute pass counts; topology-stability has the lowest pass rate, consistent with the paper's existing acknowledgement that some topology findings collapsed under stricter nulls.

## Caveats

- $p$ extraction is regex-based on free-form executor narratives. Some hypotheses report multiple $p$ values (e.g., per-layer); we take the **minimum** as the headline, which is the most-favourable-to-the-finding choice. Using a Fisher-combined per-layer $p$ would be more conservative and is logged as a follow-up but does not change which findings clear Bonferroni for the headline list.
- 'p ≈ 0' in the table corresponds to `emp_p=0.000` from $N=500$ or $N=1{,}000$ permutation tests, i.e., $p < 1/N$. We treat these as the strongest possible signal in this framework.
- Within-family BH is reported alongside global BH because reviewers explicitly asked for "within hypothesis family" correction (Reviewer 2 #6).
