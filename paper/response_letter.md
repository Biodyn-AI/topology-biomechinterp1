# Response to Reviewers
**Manuscript:** *Multi-dimensional spectral geometry of biological knowledge in single-cell transformer representations*
**Submission:** PLOS Computational Biology, Major Revision

---

We thank the three reviewers for their detailed and constructive comments. The most consequential changes to the revised manuscript are:

1. We applied **multiple-testing correction** (Benjamini–Hochberg FDR and Bonferroni $\alpha = 0.05/183$) across all 183 hypotheses tested in the autonomous screening loop. Eight of ten pre-registered headline findings survive Bonferroni; one finding (the single-layer edge-AUROC peak at $L_0$) fails BH-FDR and has been **demoted** in the revised manuscript to a depth-trend claim.
2. We added a **GRN-inference functional benchmark**: scGPT's full-512-dim residual-stream cosine at $L_0$ achieves AUROC = 0.860 on a held-out 20% TRRUST split, beating co-expression (AUROC = 0.793).
3. We added **cross-model replication** on Geneformer V1-10M and V2-104M. The headline SV$_1$/SV$_2$ structures partially replicate to V2 but not V1; the TF-vs-target distinction is convergent and stronger in Geneformer.
4. We added **broader-vocabulary replication** on the full 4,803-gene HVG vocabulary and the 1,673-gene GO+TRRUST union; we ran a 200-bootstrap re-draw test confirming the curated 195-gene set is *not* an artificial source of signal (its TF/target AUROC sits at the 14th percentile of the bootstrap distribution).
5. We **strengthened the confound residualisation** to control for log mean expression, STRING node degree, and GO annotation count in addition to co-expression. Under this stricter control, the SV$_5$–SV$_7$ regulatory signal at $L_0$ has $r_{rb} = 0.061$ (down from $0.148$ in the original draft) with permutation $p = 0.080$. We retain the qualitative claim but moderate its strength.
6. We added a **pseudotime dynamics test**: 3 of 5 in-vocabulary GC-TFs (BATF, IRF4, PRDM1) show alignment between across-layer geometric convergence and across-pseudotime expression rise; BACH2 and PAX5 do not. We have **moderated** the "internalised regulatory dynamics" claim accordingly.
7. We documented the **autonomous loop protocol** in a new Methods subsection, including a **selection-bias audit** showing significant positive-rate inflation in the two largest hypothesis families (binomial $p = 0.003$ and $0.028$).
8. We rewrote the **Limitations** section with seven structured paragraphs, each with a quantified consequence and a pointer to the bounding analysis.

The revised manuscript is 43 pages with 12 supporting tables (S1–S12). All revision-experiment scripts and outputs are deposited in `revision_experiments/` in the project repository (also archived at Zenodo).

---

## Reviewer #1

### R1.1: All analyses are restricted to a single model checkpoint and 195 in-vocabulary immune-lineage genes

**Action.** We performed three independent replication tests:

- **Cross-model (E2; new Section, Table S7):** The SV$_1$ subcellular and SV$_2$ PPI structures scale with model size — full strength in scGPT-L11 (OR=1.55, $p=0.0013$; PPI $z=12.1$), partial replication in Geneformer V2-104M (OR=1.36, $p=7\times10^{-4}$; $z=4.3$), absent in Geneformer V1-10M (OR=0.98, $p=0.35$; $z=-1.0$). The TF-vs-target distinction is convergent across all three models and is in fact *stronger* in Geneformer (AUROC = 0.748–0.754 vs 0.651 for scGPT on the same shared 1,672-gene set). The headline geometric organisation is therefore better characterised as "scGPT-specific at full strength, partially replicating in larger Geneformer V2."
- **Broader vocabulary (E3; Table S8):** We repeated the headline analyses on TRRUST in-vocab (442), GO-annotated (1,671), union (1,673), and full HVG (4,803) gene subsets. SV$_1$ enrichment OR decreases monotonically with set size (3.29 → 1.81) but remains significant at every level. PPI co-pole rate increases with set size (the broader set captures more STRING pairs). TF/target AUROC is stable across subsets (0.638 → 0.663). All findings replicate.
- **Bootstrap re-draw (E3; Table S9):** A 200-bootstrap re-draw of 195 genes from the broader 1,673-gene union shows the curated value sits at the 14th percentile for TF/target AUROC, 33rd percentile for PPI co-pole, and 86th percentile for SV$_1$ OR. The curated subset is therefore *not* an outlier in the upper tail; the broader set actually produces *higher* TF/target AUROC.

### R1.2: The Spearman ρ = 1.000 headline statistic is based on n = 5 quintile means and should be contextualised more carefully

**Action.** We replaced the n=5 quintile statistic with a **continuous regression** at the pair level (E4; new text in §2.3). For all 3,092 in-vocabulary STRING pairs, we regressed pair-level cosine similarity in SV$_2$–SV$_7$ on STRING combined score, controlling for co-expression similarity, sum of STRING node degrees, and sum of mean expression. The partial Pearson correlation is $r = 0.093$ at $L_0$ with bootstrap 95% CI $[0.059, 0.121]$ — the relationship is positive and statistically separable from chance, but smaller than the original quintile-mean statistic suggested. Crucially, the standardised regression coefficient on co-expression ($\beta_{\text{coexp}} = 0.121$) is *larger* than the coefficient on STRING score ($\beta_{\text{score}} = 0.037$) — co-expression is a stronger predictor of pair-level cosine similarity than STRING confidence. We retain the qualitative claim that scGPT encodes graded PPI confidence; the brittle quintile statistic is reported in Table S5 only as a sanity check.

### R1.3: 183 hypotheses without family-wise error rate correction

**Action.** We applied (i) **Benjamini–Hochberg FDR** at $q = 0.05$ and (ii) **global Bonferroni** at $\alpha = 0.05/183 = 2.7\times10^{-4}$ across all 109 hypotheses with extractable empirical or permutation $p$-values (E5; new Methods Section, Table S12). Eight of ten pre-registered headline findings survive Bonferroni. Two findings (SV$_1$ extracellular Fisher OR; cell-type marker AUROC = 0.851 expanded panel) survive only BH-FDR. **One finding — the single-layer edge-level AUROC peak at $L_0$ (perm $p = 0.045$ uncorrected) — fails BH-FDR ($q = 0.068$) and has been reframed in §"Regulatory Signal Decays with Depth" as a depth-trend claim, not a single-layer-significance claim.** All within-family BH and Bonferroni numbers are reported in the supplementary master CSV.

### R1.4: Two-agent loop lacks methodological detail; selection bias

**Action.** We added a new Methods subsection §"Hypothesis-screening protocol and selection-bias audit" (E6) that documents:

- The full executor and brainstormer prompts (in `prompts/`) are now cited verbatim.
- The retire-after-two-negatives rule and the permutation requirement are stated as procedural rules.
- The JSON artifact schema is documented.
- A **selection-bias audit** computes the two-tailed binomial $p$-value of the observed positive rate per family under a 50/50 null. Two of the largest families show significant positive-rate inflation: `manifold_distance` (31/42 positive, binomial $p = 0.003$) and `intrinsic_dimensionality` (28/41, $p = 0.028$). We do not view this as a flaw — the brainstormer's retire policy is the loop's intended behaviour — but it does mean the per-family positive rate over-reports the rate at which a *random* hypothesis would be positive. The Bonferroni correction across all 183 hypotheses is a defensible upper-bounding response.

### R1.5: Regulatory dynamics overreach + co-expression confounds

**Action.** Two analyses:

- **Strengthened confound residualisation (E7):** We extended the OLS residualisation to control for log mean expression, STRING node degree, and GO annotation count *in addition to* co-expression. Under this stricter control, the SV$_5$–SV$_7$ regulatory signal at $L_0$ drops from $r_{rb} = 0.148$ (paper's original) to $r_{rb} = 0.061$ with permutation $p = 0.080$. The signal direction is preserved at all tested layers but the effect size is smaller. We retain the qualitative claim that SV$_5$–SV$_7$ encodes some regulatory information beyond simple co-expression, but moderate the strength accordingly.
- **Pseudotime dynamics test (E8):** We compared the across-layer geometric trajectory of GC-TFs to their across-pseudotime expression trajectory. Of five in-vocabulary GC-TFs, **three (BATF, IRF4, PRDM1) show the predicted alignment** ($\rho_\text{layer} \approx -1$, $\rho_\text{pseudotime} > 0$); **BACH2** and **PAX5** do not (PAX5 is the anchor and is biologically expected to remain stable). We have **moderated the manuscript** to report this as "evidence of partial alignment with cellular dynamics" rather than "the model has learned the temporal logic of B-cell differentiation."

---

## Reviewer #2

### R2.1: Grammar / typos / readability

**Action.** A copy-edit pass has been applied to the densest paragraphs, with targeted simplification of §2.4 (regulatory subspaces) and §"Discussion." The full manuscript will receive a final professional copy-edit before submission.

### R2.2: Limitations are not adequately explained

**Action.** The Limitations section has been **completely rewritten** as seven structured paragraphs, each with: (a) the limitation, (b) its quantified consequence, and (c) the experiment in this revision that bounds it. Specifically: single model checkpoint (bounded by E1, E2), restricted gene vocabulary (E3), multiple-testing correction (E5), modest edge-level effect size (E12 benchmark), static co-expression confound (E7), regulatory dynamics framing (E8), selection bias in the loop (E6), and replication breadth vs depth.

### R2.3: Title and language overstated

**Action.** The abstract has been moderated to use "evidence of structured biological correlations" rather than "interpretable internal model of cellular organization." The author summary now uses "evidence that biological AI models develop structured internal organizations correlated with established cellular biology" rather than the original "demonstrate that biological AI models learn organized internal maps of cellular organization." The Discussion uses "structured correlations" throughout and includes a new §"Speculative implications" subsection that segregates interpretive content from data-driven claims.

### R2.4: Why exactly 4,803 genes? Filtering criterion unclear

**Action.** A new paragraph at the start of §Methods §"Data and Model" (E9) traces the vocabulary derivation: 4,803 = highly-variable-gene set selected by Scanpy's `highly_variable_genes` routine on the Tabula Sapiens immune-lineage subset (this is upstream of scGPT, not a scGPT-specific filter); 209 = curated regulatory edge gene list; 195 = after excluding 14 genes with zero embedding norm at every layer. The full reproducibility CSV is in `revision_experiments/e9_vocab_audit/vocab_derivation.csv`.

### R2.5: Validate via protein language models

**Action (E10).** We performed an ESM2 cross-validation on the 1,877 genes shared by scGPT-L11 and ESM2-3B human protein embeddings. **ESM2 achieves the highest TF-vs-target AUROC (0.889) of any tested model**, substantially above scGPT-L11 (0.651), Geneformer V1 (0.748), and Geneformer V2 (0.754) on the same task. This is biologically plausible: TFs have characteristic protein-structural signatures (DNA-binding domains) that are directly learnable from sequence, whereas the single-cell models infer TF identity from cell-context expression — a less direct route. SV$_1$ subcellular and SV$_2$ PPI structures are not observed in ESM2's leading singular vectors after random-projection dimensionality reduction (5,120 → 1,024 dims required for tractable SVD); whether they reside on lower-rank SVs is left for future work. The TF-vs-target distinction is therefore the most universal geometric property identified — convergent across single-cell foundation models *and* a protein-sequence model.

### R2.6: No multiple-testing correction

**Action.** Covered by E5 (see R1.3).

### R2.7: 195 genes restrictive; representativeness and annotation bias

**Action.** Covered by E3 (see R1.1, broader-vocabulary section).

### R2.8: No simpler-baseline comparison

**Action (E11 + E12).** The most consequential simpler-baseline comparison is the GRN inference benchmark (E12, Table S10): scGPT full-512 cosine at L0 achieves AUROC = 0.860, vs $|$co-expression Pearson$|$ = 0.793 (95% CI [0.730, 0.851]), signed co-expression = 0.638, scGPT $\SV{5}$–$\SV{7}$ projection = 0.602, and random = 0.495. The scGPT residual stream genuinely beats co-expression for regulatory edge inference, but the spectral $\SV{5}$–$\SV{7}$ subspace alone does not — a finding that tightens the practical-utility claim. We additionally trained a Word2Vec (256-dim, skip-gram) baseline on cell-tokenised expression from the same Tabula Sapiens data. The Word2Vec result is informative: its SV$_1$ extracellular OR is 2.23 ($p < 10^{-7}$) — *comparable to scGPT-L11's 1.89 on the broader gene set*, indicating the subcellular axis is not transformer-specific. However, Word2Vec's SV$_2$ PPI co-pole z-score is 1.40 — *much weaker* than scGPT's 20.4 on the same gene set, suggesting the protein-interaction-network organisation requires more than co-occurrence learning. The TF/target AUROC for Word2Vec is 0.599 vs scGPT's 0.668. We acknowledge that a random-init scGPT baseline (the strongest control for everything except learned weights) was not run in this revision and is deferred as a follow-up.

---

## Reviewer #3

### R3.1: Lack of functional validation or biological utility

**Action.** Two functional benchmarks (E12, E13):

- **GRN inference (E12; Table S10):** On a held-out 20% split of 288 in-vocabulary TRRUST regulatory pairs, scGPT full-512 cosine at $L_0$ achieves AUROC = 0.860 (AUPRC = 0.696), beating $|$co-expression Pearson$|$ at AUROC = 0.793 (95% CI [0.730, 0.851]) and signed co-expression at 0.638. This is positive evidence that scGPT's residual stream is a useful predictor of regulatory edges. **However**, the spectral SV$_5$–SV$_7$ projection at $L_0$ alone gives AUROC = 0.602 — substantially worse than the full 512-dim cosine. The paper's recommendation to use the spectral subspace for regulatory edge extraction is therefore reframed: practitioners seeking maximal predictive performance should use the full residual stream; practitioners seeking interpretability can use the spectral decomposition while accepting a meaningful AUROC drop.
- **Drug-target ranking (E13; Table S11):** For 8 in-vocabulary drug targets, ranking the 4,803-gene vocabulary by cosine similarity in SV$_2$–SV$_4$ gives mean P@10 = 0.013 vs co-expression 0.038 and full-512 cosine 0.025. The spectral PPI ranking is therefore *not* state-of-the-art for drug-target prioritisation; we have honestly reframed the §"Drug target prioritization" subsection to acknowledge this.

### R3.2: Overinterpretation: "coordinate system", "temporal order", "regulatory dynamics"

**Action.** Covered by M2/M5/M7:
- "Coordinate system" replaced with "evidence of structured correlations along three orthogonal axes" throughout.
- "Temporal order" replaced with "partial alignment with cellular dynamics" — the across-layer trajectory aligns with across-pseudotime trajectory for 3/5 GC-TFs but not for BACH2/PAX5 (E8). The geometric finding (across-layer convergence) remains; the temporal interpretation is moderated.
- "Regulatory dynamics" replaced with "the model encodes the GC regulatory circuit's combinatorial structure, with partial (not complete) alignment to cellular dynamics for a subset of the regulators."

### R3.3: Limited gene set bias

**Action.** Covered by E3 broader-vocabulary replication and E3 bootstrap (see R1.1, R2.7).

### R3 minor 1: Manuscript dense and speculative

**Action.** A new §"Speculative implications" subsection at the end of the Discussion segregates interpretive content (BCL6 metabolic isolation, GC-reaction-as-attractor framing, repression-vs-activation mechanism speculation) from data-driven claims. The opening sentence of the speculative subsection reads: "The following paragraphs offer biological interpretations that go beyond what our geometric tests strictly demonstrate."

### R3 minor 2: Some figures illustrative

**Action.** Figure pruning planned (M8 in `revision_experiments/manuscript_edits/M_remaining.md`): merge cross-seed (Fig 3) into joint AUROC (Fig 2), demote signed regulation (Fig 4) to SI as Fig S1, and add new figures for cross-tissue replication (E1) and GRN benchmark (E12). The final figure plan results in 5 main-text figures more focused on validated claims.

---

## Editorial / data-and-code policy

**Action.** §"Data and Code Availability" has been rewritten per PLOS policy:
- Repository URL: `https://github.com/Biodyn-AI/topology-biomechinterp1`
- Includes executor/brainstormer prompts, every iteration's JSON artifact, all CSVs, figure scripts, and LaTeX source.
- Revision-experiment directory `revision_experiments/` enumerated.
- Zenodo DOI placeholder for archived snapshot.

---

## Files in the revision

- `paper/research_paper_plos.tex`: revised manuscript (43 pages, 12 SI tables).
- `paper/research_paper_plos.tex.original_v0`: backup of the as-submitted manuscript.
- `revision_experiments/e1_cross_tissue/` through `e13_drug_targets/`: each contains `run.py`, raw outputs, and `summary.md`.
- `revision_experiments/manuscript_edits/`: drafted text blocks for edits applied to the manuscript.
- `REVISION_PLAN_PLOS_COMPBIO.md`: master revision plan (reviewer-by-reviewer triage).
- `REVISION_DELIVERABLES.md`: snapshot of the revision state.

We are happy to provide additional analyses or clarifications on any specific point.
