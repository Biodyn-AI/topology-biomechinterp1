# E10 — Cross-context variance, kidney "non-transfer", detection, lineage TFs, transfer, stability (sensitivity)

Written by the orchestrator from the implementer's findings with the verifier's Corrections 1–3 and wording fixes
applied. CONVENTION3 hash verified. Verifier: no result-changing bug; headline numbers reproduced exactly from raw
inputs (V1–V8).

## Setup
16 contexts × 12 layers; V-centred embeddings, observed-gene SVD. Statistics: S1 STRING-physical edge cosine AUROC vs
one Maslov–Sneppen rewired non-pair set; S2 detection-stratified LOO-cosine AUROC, Lambert TF vs non-TF (= H2 primary
statistic); S2b TRRUST TF vs target-only (full dim; SV2–SV7 band); S3 within-TF stratified AUROC of signed cosine on
h3_pairs (= H3 primary statistic). Joint Poisson gene bootstrap (B = 1,000) shared across contexts; S3 two-way TF ×
target. A first bootstrap removed only one copy of a resampled TF (replicates +0.05); fixed, old outputs archived in
`partial/biased_v1`. S1 uses a node bootstrap (SE within 5% of delete-d subsampling) and one rewiring draw (draw SD
0.002–0.005, negligible).

## Q1 — tissue dependence: NOT CLAIMED for any statistic
Rule (stated in the script docstring; timing not verifiable from files): D = τ²_tissue − pooled within-tissue floor
(resample + donor + HVG axis) must have 95% lower bound > 0 (parametric random-effects bootstrap, R = 10,000).
Layer-mean D (`varcomp.csv`): S1 −0.00047 [−0.0023, 0.0037]; S2 −0.0014 [−0.0031, 0.0027]; S2b_full −0.00012
[−0.0015, 0.0008]; S2b_band −0.0025 [−0.0061, 0.0050]; S3 +0.0007 [−0.0027, 0.0038]. No per-layer lower bound > 0.
Components (layer-mean τ, AUROC units): largest is the HVG axis (kidney_A vs kidney_B, same donor): S1 0.075, S2 0.073,
S3 0.052; between-tissue S1 0.020, S2 0, S3 0.033; checkpoints on identical cells S1 0.023, S2 0.027, S2b_full 0.043;
immune cell resamples S2 0.655 / 0.609 / 0.672 (I² 0.67, p_Q 0.034). **Not equivalence:** Q-profile upper bounds on
tissue τ are S1 0.086, S2 0.034, S2b_full 0.029, S2b_band 0.050, S3 0.139; tissue is confounded with HVG axis /
extraction regime; k = 2–4 per group; datasets share donors. Exploratory: S2b_band at L2–L3 parametric p 2e-4 / 1e-4
(would pass Bonferroni over 60 descriptive tests) — an axis-confounded signal (immune 0.756 vs large-axis 0.55–0.62).

## Q2 — the kidney "non-transfer" (`kidneyA_vs_B.csv`)
Submitted computation reproduced exactly: kidney_B 0.493, 0.500, 0.529, 0.551 at L0–L3 (59 TFs / 108 targets);
kidney_A layer mean 0.612. (i) Only L0–L3 were reported; L4–L11 give 0.589–0.664. (ii) Early "chance" values are the
zero-row SVD artefact: observed-gene SVD on the same 167 genes gives 0.677, 0.665, 0.648, 0.665. (iii) Full-dim on
those genes 0.644 (kidney_A 0.681). Hanley–McNeil SE at 59/108 is 0.047. Matched re-evaluation (1,126 shared genes;
49 TRRUST TFs / 103 targets; 76 Lambert TFs; paired bootstrap B = 2,000), Δ = A − B: band own-V SVD +0.055 [−0.033,
0.134]; band SVD fitted on M −0.016 [−0.057, 0.024]; S2b_full +0.030 [−0.023, 0.088]; S2 own strata +0.011 [−0.035,
0.056]; shared strata −0.003 [−0.039, 0.035]. Mechanism: axes 4,803 vs 8,143 genes (2,766 shared); 68% of a kidney_B
cell's detected genes are not on kidney_A's axis; 191/192 kidney_B cells exceed the 512-gene cap vs 65/256; median
non-zero genes per cell 1,054 vs 436; median TF n_obs on M 24 vs 14 (rank ρ 0.97). Basis-instability test not
supported at layer mean (|Δ_band| − |Δ_full| = +0.025 [−0.053, 0.092]); mixed per layer.

## Q3 — detection confound flips sign (`detection.csv`)
log n_obs alone separates Lambert TFs from non-TFs at 0.65–0.72 on immune axes (TF median n_obs 40–66 vs 4–9) but at
0.450–0.461 (CIs below 0.5) in kidney_B, lung, lung_external and checkpoints; kidney_A 0.555. S2 restricted to n_obs ≥ 10:
lung 0.686 [0.638, 0.723], lung_external 0.674, kidney_B 0.637, immune 0.701 [0.636, 0.765] (unrestricted lung 0.654,
immune 0.655).

## Q4 — lineage TFs (`lineage_tfs.csv`)
Specificity-selected top-30 sets recover canonical regulators (immune BACH2, STAT4, ETS1, BCL11B, IKZF1, POU2F2, PAX5,
LEF1, IRF4, EBF1; kidney EMX2, HNF4A, PAX8, PAX2, HNF1B, HNF1A, SALL1; lung EPAS1, NKX2-1, ELF3, EHF, GATA6, FOXA2, FOXF1).
Detection-stratified separation from non-TFs (layer mean): kidney set in kidney_B 0.962 [0.94, 0.99]; lung set in lung /
lung_external 0.903 / 0.916; immune set in immune 0.709 (random-subset reference 0.60), in lung 0.911, kidney_B 0.873.
Set AUROC vs median n_obs ρ = −0.39 (descriptive only: rows not independent; some sets < 10 TFs).

## Q5 — transfer (`transfer.csv`)
L2 logistic TF classifier trained on immune genes (OOF 0.835, slightly optimistic) scores unseen genes at kidney_B 0.652
(in-context 0.784), lung 0.674 (0.786), lung_external 0.672 (0.798) — below in-context at 12/12 layers; kidney_A 0.871 vs
0.876 in-context (CI includes 0 at all layers). Shared genes 0.955–0.966 (embedding similarity, not generalisation).

## Q6 — stability (`stability.csv`)
k = 20 neighbour Jaccard (chance 0.004–0.028): cell resample 0.64–0.74; donor 0.59–0.76; cell type 0.53–0.77; HVG axis
0.55–0.61; tissue pairs 0.37–0.51; lung vs lung_external 0.44–0.48; checkpoints 0.03–0.12. RSA 0.72–0.95 except
checkpoints 0.33–0.54. TF > non-TF stability only on immune axes and is a detection effect after decile adjustment.
Power matching to immune's sizes (R = 100): S2/S2b gaps vanish (z ≤ 0.9); S1 and S3 stay below immune (gaps 0.040–0.045
and 0.056–0.062) but at 1.2–1.9 combined SE — **not resolved**.

## Interpretation
- Direct: no between-tissue excess over the within-tissue floor for any statistic; the kidney non-transfer is an artefact
  of layer selection, the zero-row SVD and small n; matched kidney extractions agree; detection's role reverses between
  small and large gene axes; TF read-outs transfer partially (fully to kidney_A).
- Inference: context dependence is driven mainly by the extraction (gene axis → each cell's inputs; 512-gene truncation;
  detection; observed set), not by tissue biology. (Post hoc regime grouping is exploratory only.)
- Hypothesis: lineage TFs group tightly wherever their cell populations contribute a coherent programme.
- Not shown: tissue invariance (τ_tissue up to 0.03–0.09, 0.14 for S3, is not excluded).
- Manuscript: withdraw "the TF/target distinction does not transfer across tissue".
