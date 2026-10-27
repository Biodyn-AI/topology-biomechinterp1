# E11 — Compartment axis re-test and ledger failures (sensitivity analysis; CONVENTION3 §5)

Written by the orchestrator from the implementer's and verifier's returned findings, with the verifier's fix and
corrections applied. CONVENTION3 hash verified in `results.json`.

## What was run
Observed-gene SVD (paper convention) in lung, lung_external, kidney_B, immune, kidney_A. Family: SV1–SV12 ×
12 layers × 2 primary terms (GO:0005615 extracellular space, GO:0005829 cytosol), sign-invariant AUROC, max-T over
N = 2,000 joint gene-label permutations; detection-stratified version (within log n_obs deciles); annotated-only
universe; gene-bootstrap CIs; secondary terms mitochondrion, ER lumen; Fisher test at the 25% poles. Ledger
reproduction of the screen's OR = 6.37 (iter_0006) and of the submitted retraction. Candidate "failed" hypotheses:
(i) no GO enrichment at SV2 poles, (ii) FFL intermediacy, (iii) persistent homology under rewiring,
(iv) effective rank vs AUROC.

## Direct evidence
- Best cells (extracellular AUROC [95% gene-bootstrap CI], max-T z; null-max 95th pct ≈ 4.6; p_adj = 1/2001 floor
  in every context) (`maxT.csv`, `axis_summary.csv`): lung L3 SV1 0.773 [0.754, 0.792], z 34.5; lung_external L3
  SV1 0.756 [0.736, 0.777], z 30.8; kidney_B L2 SV1 0.707 [0.681, 0.733], z 20.4; immune L0 SV2 0.760
  [0.732, 0.788], z 21.6. Lung L3 SV1 top-pole Fisher OR 5.52, p = 7.9e-81; opposite pole OR 0.13.
- Replication lung → lung_external: same direction (|r| = 0.985 on 4,535 matched genes), same orientation,
  interval excludes 0.5. Matches kidney_B (0.970) and immune L0 SV2 (0.920).
- Oriented AUROCs on the lung axis cell: cytosol 0.392, mitochondrion 0.400, ER lumen 0.667.
- Detection-stratified AUROC 0.776 / 0.767 / 0.700 / 0.761 (lung / lung_external / kidney_B / immune); |ρ| of the
  coordinate with log n_obs ≤ 0.25. Annotated-only universe unchanged.
- Depth: lung SV1 extracellular AUROC 0.730 (L0), 0.773 (L3), 0.611 (L11) (`sv1_by_layer.csv`). At L11 SV1:
  lung 0.611; lung_external 0.548 with p_adj 0.037 (0.049 under another seed; **borderline**); immune 0.585
  (p_adj 0.001–0.002); kidney_A 0.646 (floor); kidney_B null (p_adj 1.0).
- kidney_A (corrected by verifier): SV1 carries an **extracellular-vs-rest** signal that survives detection
  stratification (oriented stratified AUROC 0.80–0.81 at L0–L3, 0.73 at L11) and the annotated-only universe
  (L0 0.658), but the cytosol pole fails (within-decile cytosol AUROC ≈ 0.56; 0.50 annotated-only), the direction is
  not the lung axis (r = 0.12), and SV1 is detection-correlated (ρ to −0.80). It is not the
  extracellular-vs-cytosol axis of the other four contexts.
- Ledger: OR = 6.3677, p = 2.57e-4 (21 vs 5 genes) reproduces exactly under the screen's convention, on the
  zero-row SVD (SV1 r with observed indicator −0.9998); all 22 never-observed named genes sit in the
  "extracellular" quartile (63.6% extracellular vs 17.6% of observed named genes); dropping them gives OR 2.92,
  p = 0.072. The submitted retraction (L11 OR 0.46, p 0.13; L8 OR 2.41, p 0.031) reproduces exactly and was
  also computed on the zero-row SVD with 187 curated genes. Restricted to the 95–187 named genes under the
  paper convention, the test is underpowered and null (L11 OR 0.71–2.23, p ≥ 0.06).
- (i) GO at SV2 poles: the screen's null was degenerate (one bijection applied to pole sets and term set; emp p 1.0
  and null sd 0 by construction; 181 of 591 terms not BP). Under the paper convention 411/1,387 BP terms BH-sig in
  lung, 112/297 immune, 130/249 kidney_A (364/93/121 counting OR > 1 only). Not detection-stratified.
- (ii) FFL intermediacy: lung L8 SV2–4 frac(0 ≤ t ≤ 1) 0.681 vs role-preserving null 0.616 (max-T over 36 cells,
  p_adj 0.018); immune p_adj 0.066; lung_external does not replicate (0.089 role, 0.176 label). Mean t: label-null
  max-T p_adj 0.019 in lung_external, 0.077 lung, 0.13 immune; role-preserving null ≥ 0.20 everywhere. Counts 264
  (immune), 4,017 (lung), 4,056 (lung_external) FFLs.
- (iii) PH under rewiring: screen result only; the screen's seeded draws were 58.3% (31.7–82.2%) never-observed
  zero rows; design compares Euclidean with geodesic H1.
- (iv) ER vs AUROC: raw ρ 0.860 (lung), 0.944 (immune); ρ(layer, ER) = −1.00, so **the partial correlation on layer
  is not estimable within one context** (verifier fix; earlier "partial ρ ≈ 0" values were floating-point noise and
  must not be used). Immune-triplet within-layer ρ changes sign with definitions (0.114, 0.532, −0.424, −0.844).

## Verification
Independent SVD/AUROC/max-T (seed 777) reproduces every best cell; Fisher and r reproduce; ledger reproduced from
the iter_0006 source without E11 code; FFL counts from raw TRRUST; GO BP counts. gene2go_all.pkl is not
propagated (direct annotations only). One bug (iv partial Spearman) fixed; pre-fix copies in `verify/pre_fix/`.

## Interpretation
- Direct: an extracellular/ER-lumen vs cytosol/mitochondrion axis is the leading early-layer direction in two
  tissues and an independent lung atlas, not explained by detection, and fades with depth.
- Inference: the published OR = 6.37 and the submitted retraction were both statements about the zero-row
  indicator, not about the representation; the compartment axis is real but sits at early layers, not L11.
- Hypothesis: secretory-pathway expression programmes (co-regulated secreted/ER genes across cell types) shape early
  scGPT gene geometry; whether this is learned or inherited from input statistics is not tested here.

## Caveats
p_adj at the permutation floor bounds rather than measures the evidence; CIs conditional on the estimated
direction; GO direct-only; no cross-model/static-token baseline for the axis; candidate (i) not detection-stratified;
borderline cells (lung_external L11 SV1) not resolved at N = 2,000.

## Addendum — is the axis specific to the expression-trained residual stream? (`e11_static_baselines.py`)
Run by the orchestrator after the E11 workflow (sensitivity; no pass/fail). One common gene set per context
(V ∩ scGPT static vocabulary ∩ ESM2 file: lung 5,413, lung_external 5,306, kidney_B, immune); representations scGPT
residual L0/L3/L11, scGPT static input-token embedding, ESM2-15B per-protein embedding (5,120-d, the UCE file).
(a) own-SVD SV1 and best-of-SV1–12 sign-free AUROC; (b) out-of-fold L2 logistic probe on 256 PCs (5 outer folds,
inner-CV C), gene-bootstrap CI, detection-stratified AUROC, paired bootstrap Δ vs scGPT L3 (`static_baselines.csv`,
`static_baselines_paired.csv`).
- Linear probe, extracellular: scGPT 0.80–0.84 at L0, L3 and L11 alike (|Δ| ≤ 0.008 vs L3 everywhere, CIs include
  0); scGPT static token embedding the same (|Δ| ≤ 0.008); ESM2-15B higher by +0.062 to +0.079 (CIs 0.040–0.095,
  all exclude 0). Cytosol: scGPT 0.68–0.76, ESM2 +0.079 to +0.102. Detection-stratified AUROCs within 0.02 of raw.
- Unsupervised geometry: scGPT L0/L3 SV1 is the compartment axis (lung 0.746/0.779, lung_external 0.730/0.766,
  kidney_B 0.694/0.724); immune uses SV3–4 on this gene set. ESM2 SV1 0.54–0.56 everywhere; its best direction is
  SV8 (0.69–0.76). Static token embedding SV1 0.51–0.58, best SV2 (0.70–0.76).
- Interpretation. Direct: the compartment information is linearly present to the same degree at every scGPT layer and
  already in the model's input gene embedding; the protein-sequence model carries more of it. Inference: depth does
  not add or remove compartment information; what the early layers do is make it the leading axis of variation among
  genes (in ESM2 and in the static token table it sits on minor directions). Hypothesis: sequence features of the
  secretory pathway (signal peptides, transmembrane segments) are why ESM2 reads compartment better, and the
  expression context is what promotes compartment to the leading scGPT direction.
