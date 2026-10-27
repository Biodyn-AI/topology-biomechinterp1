# E7 — H4 rank vs full (confirmatory) — orchestrator copy of the generated SUMMARY.txt (verifier-audited: delete-d PPI CIs; k = 1 tie fix)

```
E7 - H4: does truncating to the top-k spectral directions improve recovery over the full representation?
Pre-registered primary (CONVENTION3 §4), sha256 3e3249d1d0a8ce1bc3cf64639fe233b016535456f82432a8bbfcacdd64a05c01 (matches FROZEN.txt: True).

SELECTION (immune_seed43, layer-mean, family T, ties -> smallest k; kstar.json written before any evaluation context was computed):
  PPI  k* =  64   layer-mean metric(k*) = 0.7047, metric(512) = 0.7021, selection gap +0.0026
  NODE k* = 512   layer-mean metric(k*) = 0.6088, metric(512) = 0.6088, selection gap +0.0000
  EDGE k* = 128   layer-mean metric(k*) = 0.6199, metric(512) = 0.6198, selection gap +0.0001

MAIN TABLE: layer-mean Delta = metric(k*) - metric(512) [95% CI], delta.csv (layer == 'mean')
  task    k*                        lung               lung_external                    kidney_B                      immune
  PPI     64  -0.0011 [-0.0024,+0.0003]  -0.0006 [-0.0019,+0.0007]  -0.0013 [-0.0026,-0.0000]  +0.0022 [-0.0018,+0.0064]
  (PPI, pre-registered sqrt(m/n) scaling)         [-0.0020,-0.0001]          [-0.0015,+0.0003]          [-0.0022,-0.0004]          [-0.0006,+0.0052]
  NODE   512  +0.0000 [+0.0000,+0.0000]  +0.0000 [+0.0000,+0.0000]  +0.0000 [+0.0000,+0.0000]  +0.0000 [+0.0000,+0.0000]
  EDGE   128  +0.0004 [-0.0022,+0.0026]  +0.0003 [-0.0021,+0.0026]  -0.0002 [-0.0027,+0.0023]  -0.0004 [-0.0070,+0.0061]
  PPI CIs in the table use the delete-d vertex-subsampling scaling sqrt(m/(n-m)) (CONVENTION3_ERRATA item 3, manuscript convention; delta.csv ci_lo_deleted/ci_hi_deleted = ci_lo_ms/ci_hi_ms); the pre-registered sqrt(m/n) intervals are kept in delta.csv ci_lo/ci_hi. NODE/EDGE: ci_lo_ms = ci_lo.
  full-dimension layer-mean metric(512): PPI: lung 0.655, lung_external 0.658, kidney_B 0.649, immune 0.700; NODE: lung 0.654, lung_external 0.644, kidney_B 0.634, immune 0.655; EDGE: lung 0.530, lung_external 0.530, kidney_B 0.541, immune 0.610
  per-layer CIs: delta.csv; CI methods: NODE gene bootstrap, EDGE two-way TF x target cluster bootstrap, PPI vertex subsampling (m = n/2) of the paired difference; B = 2,000.
  per-layer: number of layers (of 12) with CI above 0 / below 0 (PPI: delete-d; pre-registered in brackets):
    PPI: lung 0/6 [0/9]; lung_external 0/0 [0/2]; kidney_B 0/7 [0/10]; immune 0/0 [0/0]
    NODE: lung 0/0; lung_external 0/0; kidney_B 0/0; immune 0/0
    EDGE: lung 0/0; lung_external 0/0; kidney_B 0/0; immune 0/0

ANSWERS (frozen form):
  PPI: Rank reduction (k* = 64 of 512) does not improve PPI-partner recovery by layer-mean Delta [95% CI] = lung -0.001 [-0.002, +0.000]; lung_external -0.001 [-0.002, +0.001]; kidney_B -0.001 [-0.003, -0.000]; immune +0.002 [-0.002, +0.006] in 0 of 4 contexts (CI below 0 in 1 of 4). [PPI CIs: delete-d vertex-subsampling scaling per CONVENTION3_ERRATA item 3; with the pre-registered sqrt(m/n) scaling: lung [-0.002, -0.000]; lung_external [-0.002, +0.000]; kidney_B [-0.002, -0.000]; immune [-0.001, +0.005], CI above 0 in 0 of 4, below 0 in 2 of 4.]
  NODE: Rank reduction does not improve TF-identity recovery by Delta = 0 [0, 0] in 0 of 4 contexts: the selection context chose k* = 512 (the full representation), so Delta is identically zero; every truncation scored below full in immune_seed43.
  EDGE: Rank reduction (k* = 128 of 512) does not improve regulatory-edge recovery by layer-mean Delta [95% CI] = lung +0.000 [-0.002, +0.003]; lung_external +0.000 [-0.002, +0.003]; kidney_B -0.000 [-0.003, +0.002]; immune -0.000 [-0.007, +0.006] in 0 of 4 contexts.

DESCRIPTIVE joint 'any k beats full' (any_k_test.csv; N = 500 joint permutations across k and layers; statistic max over k<512 of the layer-mean gap):
  lung           PPI  observed max gap +0.0002 (k=256), null p = 0.848 (null q95 +0.0042);  any layer: +0.0007 (k=128, L0), p = 0.964
  lung           NODE observed max gap +0.0001 (k=128), null p = 0.772 (null q95 +0.0295);  any layer: +0.0030 (k=64, L4), p = 0.852
  lung           EDGE observed max gap +0.0004 (k=128), null p = 0.808 (null q95 +0.0069);  any layer: +0.0021 (k=128, L1), p = 0.934
  lung_external  PPI  observed max gap -0.0000 (k=256), null p = 0.946 (null q95 +0.0041);  any layer: +0.0003 (k=128, L8), p = 0.990
  lung_external  NODE observed max gap +0.0020 (k=64), null p = 0.445 (null q95 +0.0201);  any layer: +0.0051 (k=32, L4), p = 0.629
  lung_external  EDGE observed max gap +0.0003 (k=128), null p = 0.832 (null q95 +0.0075);  any layer: +0.0024 (k=128, L0), p = 0.880
  kidney_B       PPI  observed max gap -0.0000 (k=256), null p = 0.914 (null q95 +0.0042);  any layer: +0.0001 (k=256, L0), p = 0.996
  kidney_B       NODE observed max gap -0.0009 (k=256), null p = 0.990 (null q95 +0.0123);  any layer: +0.0037 (k=32, L2), p = 0.850
  kidney_B       EDGE observed max gap +0.0005 (k=256), null p = 0.786 (null q95 +0.0082);  any layer: +0.0014 (k=256, L2), p = 0.960
  immune         PPI  observed max gap +0.0022 (k=64), null p = 0.649 (null q95 +0.0164);  any layer: +0.0092 (k=32, L1), p = 0.707
  immune         NODE observed max gap -0.0012 (k=256), null p = 0.960 (null q95 +0.0312);  any layer: +0.0021 (k=128, L0), p = 0.972
  immune         EDGE observed max gap -0.0002 (k=64), null p = 0.960 (null q95 +0.0216);  any layer: +0.0039 (k=64, L1), p = 0.962
  Caveat: the label permutation erases the signal in every representation, so this p measures the gap against label noise only; it is not a test of H0 metric(k) <= metric(512).

POST-HOC curve shape (posthoc_best_k.csv; best top-k per context, layer-mean; descriptive):
  immune         EDGE best k=512 0.6096 vs full 0.6096 (+0.0000); k=10 0.5659, k=64 0.6094, k=128 0.6092
  immune         NODE best k=512 0.6553 vs full 0.6553 (+0.0000); k=10 0.5714, k=64 0.6433, k=128 0.6529
  immune         PPI  best k= 64 0.7020 vs full 0.6998 (+0.0022); k=10 0.6823, k=64 0.7020, k=128 0.7012
  immune_seed43  EDGE best k=128 0.6199 vs full 0.6198 (+0.0001); k=10 0.5832, k=64 0.6189, k=128 0.6199
  immune_seed43  NODE best k=512 0.6088 vs full 0.6088 (+0.0000); k=10 0.5539, k=64 0.5910, k=128 0.6028
  immune_seed43  PPI  best k= 64 0.7047 vs full 0.7021 (+0.0026); k=10 0.6844, k=64 0.7047, k=128 0.7045
  kidney_A       EDGE best k= 64 0.6350 vs full 0.6317 (+0.0034); k=10 0.6032, k=64 0.6350, k=128 0.6333
  kidney_A       NODE best k=256 0.7420 vs full 0.7417 (+0.0003); k=10 0.7003, k=64 0.7372, k=128 0.7417
  kidney_A       PPI  best k= 64 0.7703 vs full 0.7674 (+0.0029); k=10 0.7596, k=64 0.7703, k=128 0.7692
  kidney_B       EDGE best k=256 0.5417 vs full 0.5413 (+0.0005); k=10 0.5275, k=64 0.5403, k=128 0.5411
  kidney_B       NODE best k=512 0.6343 vs full 0.6343 (+0.0000); k=10 0.6006, k=64 0.6335, k=128 0.6333
  kidney_B       PPI  best k=512 0.6486 vs full 0.6486 (+0.0000); k=10 0.6343, k=64 0.6474, k=128 0.6484
  lung           EDGE best k=128 0.5304 vs full 0.5300 (+0.0004); k=10 0.5158, k=64 0.5289, k=128 0.5304
  lung           NODE best k=128 0.6542 vs full 0.6542 (+0.0001); k=10 0.6122, k=64 0.6537, k=128 0.6542
  lung           PPI  best k=256 0.6551 vs full 0.6549 (+0.0002); k=10 0.6408, k=64 0.6539, k=128 0.6549
  lung_external  EDGE best k=128 0.5303 vs full 0.5300 (+0.0003); k=10 0.5186, k=64 0.5298, k=128 0.5303
  lung_external  NODE best k= 64 0.6456 vs full 0.6436 (+0.0020); k=10 0.6254, k=64 0.6456, k=128 0.6448
  lung_external  PPI  best k=512 0.6579 vs full 0.6579 (+0.0000); k=10 0.6467, k=64 0.6573, k=128 0.6577

FAMILY R remove-top-j (remove_top_j.csv; layer-mean diff vs full):
  immune_seed43  PPI  R1:+0.0219  R2:+0.0123  R3:+0.0130  R5:+0.0037  R10:-0.0413
  immune_seed43  NODE R1:+0.0338  R2:+0.0923  R3:+0.1014  R5:+0.1170  R10:+0.1299
  immune_seed43  EDGE R1:-0.0004  R2:+0.0096  R3:+0.0141  R5:+0.0184  R10:+0.0045
  lung           PPI  R1:-0.0068  R2:-0.0099  R3:-0.0238  R5:-0.0270  R10:-0.0410
  lung           NODE R1:+0.0364  R2:+0.0602  R3:+0.0907  R5:+0.0991  R10:+0.0892
  lung           EDGE R1:+0.0093  R2:+0.0186  R3:+0.0327  R5:+0.0405  R10:+0.0170
  lung_external  PPI  R1:+0.0022  R2:-0.0225  R3:-0.0296  R5:-0.0364  R10:-0.0491
  lung_external  NODE R1:+0.0281  R2:+0.0982  R3:+0.0991  R5:+0.1072  R10:+0.0915
  lung_external  EDGE R1:+0.0108  R2:+0.0314  R3:+0.0375  R5:+0.0385  R10:+0.0143
  kidney_B       PPI  R1:+0.0012  R2:-0.0065  R3:-0.0195  R5:-0.0175  R10:-0.0329
  kidney_B       NODE R1:+0.0249  R2:+0.0698  R3:+0.0907  R5:+0.1082  R10:+0.1081
  kidney_B       EDGE R1:+0.0092  R2:+0.0186  R3:+0.0229  R5:+0.0234  R10:+0.0026
  immune         PPI  R1:+0.0267  R2:+0.0181  R3:+0.0147  R5:-0.0000  R10:-0.0453
  immune         NODE R1:+0.0426  R2:+0.0822  R3:+0.0862  R5:+0.0924  R10:+0.0994
  immune         EDGE R1:+0.0052  R2:+0.0200  R3:+0.0247  R5:+0.0258  R10:+0.0217
  kidney_A       PPI  R1:+0.0312  R2:+0.0116  R3:+0.0038  R5:-0.0238  R10:-0.0681
  kidney_A       NODE R1:+0.0056  R2:-0.0164  R3:-0.0022  R5:-0.0067  R10:-0.0008
  kidney_A       EDGE R1:+0.0033  R2:-0.0052  R3:-0.0089  R5:-0.0119  R10:-0.0251

FAMILY N screen bands vs same-width top-k (bands_vs_topk.csv; layer-mean value, diff vs partner):
  immune_seed43  PPI  nom SV2-4:0.636(+0.032 vs T3)  nom SV5-7:0.646(+0.042 vs T3)  nom SV2-7:0.685(+0.033 vs T6)  kA-rule SV1-3:0.604(+0.000 vs T3)  kA-rule SV4-6:0.630(+0.027 vs T3)  kA-rule SV1-6:0.652(+0.000 vs T6)
  immune_seed43  NODE nom SV2-4:0.520(+0.004 vs T3)  nom SV5-7:0.542(+0.026 vs T3)  nom SV2-7:0.561(+0.037 vs T6)  kA-rule SV1-3:0.516(+0.000 vs T3)  kA-rule SV4-6:0.493(-0.023 vs T3)  kA-rule SV1-6:0.523(+0.000 vs T6)
  immune_seed43  EDGE nom SV2-4:0.516(-0.016 vs T3)  nom SV5-7:0.538(+0.006 vs T3)  nom SV2-7:0.547(-0.012 vs T6)  kA-rule SV1-3:0.533(+0.000 vs T3)  kA-rule SV4-6:0.533(+0.000 vs T3)  kA-rule SV1-6:0.560(+0.000 vs T6)
  lung           PPI  nom SV2-4:0.596(-0.018 vs T3)  nom SV5-7:0.577(-0.036 vs T3)  nom SV2-7:0.619(-0.013 vs T6)  kA-rule SV1-3:0.614(+0.000 vs T3)  kA-rule SV4-6:0.578(-0.036 vs T3)  kA-rule SV1-6:0.632(+0.000 vs T6)
  lung           NODE nom SV2-4:0.549(+0.023 vs T3)  nom SV5-7:0.571(+0.045 vs T3)  nom SV2-7:0.595(+0.015 vs T6)  kA-rule SV1-3:0.526(+0.000 vs T3)  kA-rule SV4-6:0.577(+0.051 vs T3)  kA-rule SV1-6:0.580(+0.000 vs T6)
  lung           EDGE nom SV2-4:0.496(+0.000 vs T3)  nom SV5-7:0.541(+0.045 vs T3)  nom SV2-7:0.512(+0.007 vs T6)  kA-rule SV1-3:0.496(+0.000 vs T3)  kA-rule SV4-6:0.532(+0.036 vs T3)  kA-rule SV1-6:0.505(+0.000 vs T6)
  lung_external  PPI  nom SV2-4:0.612(-0.014 vs T3)  nom SV5-7:0.577(-0.049 vs T3)  nom SV2-7:0.634(-0.007 vs T6)  kA-rule SV1-3:0.626(+0.000 vs T3)  kA-rule SV4-6:0.575(-0.051 vs T3)  kA-rule SV1-6:0.642(+0.000 vs T6)
  lung_external  NODE nom SV2-4:0.560(+0.003 vs T3)  nom SV5-7:0.619(+0.063 vs T3)  nom SV2-7:0.619(+0.026 vs T6)  kA-rule SV1-3:0.556(+0.000 vs T3)  kA-rule SV4-6:0.564(+0.008 vs T3)  kA-rule SV1-6:0.592(+0.000 vs T6)
  lung_external  EDGE nom SV2-4:0.500(+0.008 vs T3)  nom SV5-7:0.543(+0.050 vs T3)  nom SV2-7:0.519(+0.011 vs T6)  kA-rule SV1-3:0.492(+0.000 vs T3)  kA-rule SV4-6:0.534(+0.042 vs T3)  kA-rule SV1-6:0.508(+0.000 vs T6)
  kidney_B       PPI  nom SV2-4:0.602(-0.009 vs T3)  nom SV5-7:0.571(-0.039 vs T3)  nom SV2-7:0.623(-0.005 vs T6)  kA-rule SV1-3:0.610(+0.000 vs T3)  kA-rule SV4-6:0.571(-0.039 vs T3)  kA-rule SV1-6:0.628(+0.000 vs T6)
  kidney_B       NODE nom SV2-4:0.550(+0.004 vs T3)  nom SV5-7:0.590(+0.044 vs T3)  nom SV2-7:0.597(+0.020 vs T6)  kA-rule SV1-3:0.547(+0.000 vs T3)  kA-rule SV4-6:0.546(-0.001 vs T3)  kA-rule SV1-6:0.578(+0.000 vs T6)
  kidney_B       EDGE nom SV2-4:0.511(+0.006 vs T3)  nom SV5-7:0.540(+0.036 vs T3)  nom SV2-7:0.531(+0.013 vs T6)  kA-rule SV1-3:0.504(+0.000 vs T3)  kA-rule SV4-6:0.533(+0.029 vs T3)  kA-rule SV1-6:0.519(+0.000 vs T6)
  immune         PPI  nom SV2-4:0.645(+0.036 vs T3)  nom SV5-7:0.647(+0.037 vs T3)  nom SV2-7:0.686(+0.031 vs T6)  kA-rule SV1-3:0.609(+0.000 vs T3)  kA-rule SV4-6:0.635(+0.026 vs T3)  kA-rule SV1-6:0.655(+0.000 vs T6)
  immune         NODE nom SV2-4:0.552(+0.039 vs T3)  nom SV5-7:0.477(-0.037 vs T3)  nom SV2-7:0.564(+0.019 vs T6)  kA-rule SV1-3:0.513(+0.000 vs T3)  kA-rule SV4-6:0.501(-0.012 vs T3)  kA-rule SV1-6:0.545(+0.000 vs T6)
  immune         EDGE nom SV2-4:0.511(-0.015 vs T3)  nom SV5-7:0.542(+0.016 vs T3)  nom SV2-7:0.540(-0.008 vs T6)  kA-rule SV1-3:0.527(+0.000 vs T3)  kA-rule SV4-6:0.535(+0.009 vs T3)  kA-rule SV1-6:0.548(+0.000 vs T6)
  kidney_A       PPI  nom SV2-4:0.734(+0.016 vs T3)  nom SV5-7:0.697(-0.021 vs T3)  nom SV2-7:0.770(+0.023 vs T6)  kA-rule SV1-3:0.718(+0.000 vs T3)  kA-rule SV4-6:0.699(-0.019 vs T3)  kA-rule SV1-6:0.747(+0.000 vs T6)
  kidney_A       NODE nom SV2-4:0.560(-0.060 vs T3)  nom SV5-7:0.625(+0.005 vs T3)  nom SV2-7:0.663(-0.014 vs T6)  kA-rule SV1-3:0.620(+0.000 vs T3)  kA-rule SV4-6:0.628(+0.008 vs T3)  kA-rule SV1-6:0.677(+0.000 vs T6)
  kidney_A       EDGE nom SV2-4:0.538(-0.030 vs T3)  nom SV5-7:0.532(-0.036 vs T3)  nom SV2-7:0.561(-0.015 vs T6)  kA-rule SV1-3:0.568(+0.000 vs T3)  kA-rule SV4-6:0.518(-0.050 vs T3)  kA-rule SV1-6:0.576(+0.000 vs T6)
  kidney_A rule: SV1-3 and SV1-6 ARE top-3 / top-6 (diff 0 by construction); immune: mapping = identity.

SENSITIVITY with paired CIs (sens_R_N_ci.csv; same resamples as the primary; layer-mean diff [95% CI]; PPI CIs delete-d scaled, columns ci_lo_ms/ci_hi_ms):
  lung           PPI  R1-T512:-0.007[-0.012,-0.002]  R2-T512:-0.010[-0.020,+0.000]  R3-T512:-0.024[-0.037,-0.010]  R5-T512:-0.027[-0.042,-0.013]  R10-T512:-0.041[-0.056,-0.026]  Nnom_SV2-4-T3:-0.018[-0.024,-0.012]  Nnom_SV5-7-T3:-0.036[-0.053,-0.019]  Nnom_SV2-7-T6:-0.013[-0.018,-0.008]  Nmap_SV5-7-T3:-0.036[-0.053,-0.018]
  lung           NODE R1-T512:+0.036[+0.027,+0.043]  R2-T512:+0.060[+0.047,+0.070]  R3-T512:+0.091[+0.074,+0.102]  R5-T512:+0.099[+0.080,+0.111]  R10-T512:+0.089[+0.063,+0.103]  Nnom_SV2-4-T3:+0.023[+0.008,+0.036]  Nnom_SV5-7-T3:+0.045[+0.012,+0.083]  Nnom_SV2-7-T6:+0.015[+0.006,+0.022]  Nmap_SV5-7-T3:+0.051[+0.017,+0.090]
  lung           EDGE R1-T512:+0.009[-0.007,+0.021]  R2-T512:+0.019[-0.007,+0.045]  R3-T512:+0.033[-0.009,+0.071]  R5-T512:+0.041[-0.006,+0.081]  R10-T512:+0.017[-0.034,+0.060]  Nnom_SV2-4-T3:+0.000[-0.020,+0.015]  Nnom_SV5-7-T3:+0.045[-0.009,+0.095]  Nnom_SV2-7-T6:+0.007[-0.011,+0.021]  Nmap_SV5-7-T3:+0.036[-0.015,+0.084]
  lung_external  PPI  R1-T512:+0.002[-0.004,+0.009]  R2-T512:-0.023[-0.036,-0.008]  R3-T512:-0.030[-0.044,-0.014]  R5-T512:-0.036[-0.051,-0.021]  R10-T512:-0.049[-0.065,-0.033]  Nnom_SV2-4-T3:-0.014[-0.022,-0.006]  Nnom_SV5-7-T3:-0.049[-0.067,-0.030]  Nnom_SV2-7-T6:-0.007[-0.014,-0.000]  Nmap_SV5-7-T3:-0.051[-0.069,-0.033]
  lung_external  NODE R1-T512:+0.028[+0.019,+0.037]  R2-T512:+0.098[+0.079,+0.108]  R3-T512:+0.099[+0.080,+0.112]  R5-T512:+0.107[+0.087,+0.122]  R10-T512:+0.092[+0.060,+0.105]  Nnom_SV2-4-T3:+0.003[-0.012,+0.020]  Nnom_SV5-7-T3:+0.063[+0.035,+0.103]  Nnom_SV2-7-T6:+0.026[+0.014,+0.043]  Nmap_SV5-7-T3:+0.008[-0.024,+0.046]
  lung_external  EDGE R1-T512:+0.011[-0.005,+0.026]  R2-T512:+0.031[-0.006,+0.067]  R3-T512:+0.037[-0.007,+0.075]  R5-T512:+0.039[-0.007,+0.078]  R10-T512:+0.014[-0.033,+0.056]  Nnom_SV2-4-T3:+0.008[-0.020,+0.040]  Nnom_SV5-7-T3:+0.050[-0.003,+0.104]  Nnom_SV2-7-T6:+0.011[-0.008,+0.029]  Nmap_SV5-7-T3:+0.042[-0.012,+0.091]
  kidney_B       PPI  R1-T512:+0.001[-0.005,+0.007]  R2-T512:-0.006[-0.018,+0.005]  R3-T512:-0.019[-0.033,-0.005]  R5-T512:-0.017[-0.033,-0.002]  R10-T512:-0.033[-0.048,-0.016]  Nnom_SV2-4-T3:-0.009[-0.017,-0.001]  Nnom_SV5-7-T3:-0.039[-0.055,-0.022]  Nnom_SV2-7-T6:-0.005[-0.011,+0.001]  Nmap_SV5-7-T3:-0.039[-0.056,-0.021]
  kidney_B       NODE R1-T512:+0.025[+0.015,+0.035]  R2-T512:+0.070[+0.050,+0.079]  R3-T512:+0.091[+0.065,+0.104]  R5-T512:+0.108[+0.084,+0.122]  R10-T512:+0.108[+0.077,+0.123]  Nnom_SV2-4-T3:+0.004[-0.010,+0.019]  Nnom_SV5-7-T3:+0.044[+0.012,+0.083]  Nnom_SV2-7-T6:+0.020[+0.005,+0.036]  Nmap_SV5-7-T3:-0.001[-0.036,+0.041]
  kidney_B       EDGE R1-T512:+0.009[-0.005,+0.024]  R2-T512:+0.019[-0.007,+0.043]  R3-T512:+0.023[-0.014,+0.060]  R5-T512:+0.023[-0.018,+0.067]  R10-T512:+0.003[-0.040,+0.047]  Nnom_SV2-4-T3:+0.006[-0.016,+0.027]  Nnom_SV5-7-T3:+0.036[-0.008,+0.078]  Nnom_SV2-7-T6:+0.013[-0.004,+0.029]  Nmap_SV5-7-T3:+0.029[-0.017,+0.082]
  immune         PPI  R1-T512:+0.027[+0.011,+0.042]  R2-T512:+0.018[-0.004,+0.041]  R3-T512:+0.015[-0.008,+0.038]  R5-T512:-0.000[-0.027,+0.027]  R10-T512:-0.045[-0.081,-0.011]  Nnom_SV2-4-T3:+0.036[+0.009,+0.063]  Nnom_SV5-7-T3:+0.037[-0.002,+0.075]  Nnom_SV2-7-T6:+0.031[+0.010,+0.051]  Nmap_SV5-7-T3:+0.026[-0.009,+0.063]
  immune         NODE R1-T512:+0.043[+0.006,+0.070]  R2-T512:+0.082[+0.024,+0.112]  R3-T512:+0.086[+0.021,+0.110]  R5-T512:+0.092[+0.035,+0.122]  R10-T512:+0.099[+0.037,+0.146]  Nnom_SV2-4-T3:+0.039[-0.015,+0.103]  Nnom_SV5-7-T3:-0.037[-0.094,+0.066]  Nnom_SV2-7-T6:+0.019[-0.026,+0.065]  Nmap_SV5-7-T3:-0.012[-0.073,+0.079]
  immune         EDGE R1-T512:+0.005[-0.019,+0.031]  R2-T512:+0.020[-0.025,+0.061]  R3-T512:+0.025[-0.017,+0.064]  R5-T512:+0.026[-0.021,+0.070]  R10-T512:+0.022[-0.058,+0.086]  Nnom_SV2-4-T3:-0.015[-0.064,+0.045]  Nnom_SV5-7-T3:+0.016[-0.062,+0.080]  Nnom_SV2-7-T6:-0.008[-0.049,+0.030]  Nmap_SV5-7-T3:+0.009[-0.062,+0.074]
  kidney_A       PPI  R1-T512:+0.031[+0.015,+0.046]  R2-T512:+0.012[-0.015,+0.037]  R3-T512:+0.004[-0.026,+0.032]  R5-T512:-0.024[-0.061,+0.012]  R10-T512:-0.068[-0.110,-0.028]  Nnom_SV2-4-T3:+0.016[-0.005,+0.038]  Nnom_SV5-7-T3:-0.021[-0.061,+0.016]  Nnom_SV2-7-T6:+0.023[+0.008,+0.038]  Nmap_SV5-7-T3:-0.019[-0.057,+0.019]
  kidney_A       NODE R1-T512:+0.006[-0.047,+0.027]  R2-T512:-0.016[-0.067,+0.028]  R3-T512:-0.002[-0.055,+0.044]  R5-T512:-0.007[-0.090,+0.044]  R10-T512:-0.001[-0.100,+0.059]  Nnom_SV2-4-T3:-0.060[-0.124,-0.024]  Nnom_SV5-7-T3:+0.005[-0.057,+0.078]  Nnom_SV2-7-T6:-0.014[-0.060,+0.013]  Nmap_SV5-7-T3:+0.008[-0.064,+0.082]
  kidney_A       EDGE R1-T512:+0.003[-0.046,+0.046]  R2-T512:-0.005[-0.095,+0.065]  R3-T512:-0.009[-0.109,+0.072]  R5-T512:-0.012[-0.109,+0.076]  R10-T512:-0.025[-0.116,+0.067]  Nnom_SV2-4-T3:-0.030[-0.099,+0.034]  Nnom_SV5-7-T3:-0.036[-0.186,+0.095]  Nnom_SV2-7-T6:-0.015[-0.068,+0.033]  Nmap_SV5-7-T3:-0.050[-0.180,+0.060]

STATIC MODELS on the lung Tier B gene set (static_curves.csv / static_summary.csv; descriptive):
  scgpt_static           PPI  d=  512 full 0.6609  best k=  64 0.6623 (+0.0014); k=10 0.6475
  scgpt_static           NODE d=  512 full 0.7512  best k= 512 0.7512 (+0.0000); k=10 0.6690
  scgpt_static           EDGE d=  512 full 0.5661  best k= 128 0.5706 (+0.0045); k=10 0.5514
  gf_v2_316m             PPI  d= 1152 full 0.6499  best k=   5 0.6541 (+0.0042); k=10 0.6539
  gf_v2_316m             NODE d= 1152 full 0.8209  best k=1152 0.8209 (+0.0000); k=10 0.6432
  gf_v2_316m             EDGE d= 1152 full 0.5610  best k=  16 0.5649 (+0.0039); k=10 0.5487
  gf_v1_10m              PPI  d=  256 full 0.5833  best k= 256 0.5833 (+0.0000); k=10 0.5328
  gf_v1_10m              NODE d=  256 full 0.7254  best k= 256 0.7254 (+0.0000); k=10 0.5439
  gf_v1_10m              EDGE d=  256 full 0.5177  best k= 256 0.5177 (+0.0000); k=10 0.4811
  gf_v2_104m             PPI  d=  768 full 0.6512  best k=   7 0.6679 (+0.0167); k=10 0.6628
  gf_v2_104m             NODE d=  768 full 0.8250  best k= 512 0.8262 (+0.0012); k=10 0.6522
  gf_v2_104m             EDGE d=  768 full 0.5624  best k=  64 0.5728 (+0.0104); k=10 0.5518
  aido_100m              PPI  d=  640 full 0.6050  best k= 256 0.6054 (+0.0004); k=10 0.5945
  aido_100m              NODE d=  640 full 0.6269  best k=  64 0.6489 (+0.0221); k=10 0.5808
  aido_100m              EDGE d=  640 full 0.5463  best k= 640 0.5463 (+0.0000); k=10 0.5290
  scfoundation           PPI  d=  512 full 0.5415  best k=  32 0.5415 (+0.0000); k=10 0.5410
  scfoundation           NODE d=  512 full 0.4535  best k=   5 0.4569 (+0.0033); k=10 0.4517
  scfoundation           EDGE d=  512 full 0.5024  best k=  32 0.5040 (+0.0015); k=10 0.5019
  tgpt                   PPI  d= 1024 full 0.5745  best k=1024 0.5745 (+0.0000); k=10 0.5656
  tgpt                   NODE d= 1024 full 0.6121  best k= 512 0.6132 (+0.0011); k=10 0.5515
  tgpt                   EDGE d= 1024 full 0.5189  best k= 512 0.5199 (+0.0010); k=10 0.5082
  esm2_uce               PPI  d= 5120 full 0.6385  best k= 512 0.6387 (+0.0002); k=10 0.5997
  esm2_uce               NODE d= 5120 full 0.8799  best k=  32 0.8817 (+0.0018); k=10 0.8560
  esm2_uce               EDGE d= 5120 full 0.4985  best k=   1 0.5080 (+0.0095); k=10 0.4909
  esm2_650m              PPI  d= 1280 full 0.6154  best k= 512 0.6155 (+0.0001); k=10 0.6017
  esm2_650m              NODE d= 1280 full 0.8461  best k=  64 0.8468 (+0.0006); k=10 0.8421
  esm2_650m              EDGE d= 1280 full 0.4966  best k=   1 0.5031 (+0.0066); k=10 0.4910
  seq_features           PPI  d=   32 full 0.6668  best k=  32 0.6668 (+0.0000); k=10 0.6624
  seq_features           NODE d=   32 full 0.8670  best k=  32 0.8670 (+0.0000); k=10 0.8635
  seq_features           EDGE d=   32 full 0.4987  best k=   5 0.4998 (+0.0011); k=10 0.4949
  hpa_coexp              PPI  d=  128 full 0.5951  best k=  16 0.5975 (+0.0024); k=10 0.5970
  hpa_coexp              NODE d=  128 full 0.5735  best k= 128 0.5735 (+0.0000); k=10 0.5578
  hpa_coexp              EDGE d=  128 full 0.5149  best k=  16 0.5156 (+0.0007); k=10 0.5146
  scgpt_lung_contextual  PPI  d=  512 full 0.6493  best k= 512 0.6493 (+0.0000); k=10 0.6352
  scgpt_lung_contextual  NODE d=  512 full 0.6572  best k= 128 0.6574 (+0.0002); k=10 0.6116
  scgpt_lung_contextual  EDGE d=  512 full 0.5296  best k= 128 0.5301 (+0.0006); k=10 0.5162

INTERPRETATION (question -> result -> interpretation -> implication; direct evidence vs inference vs hypothesis marked)
[revised by the independent E7 verifier, 2026-10-26: PPI intervals now delete-d scaled (CONVENTION3_ERRATA item 3);
k = 1 values corrected (exact +-1 cosines); small-k loss range corrected. See AUDIT below.]
Q (R2.6): does truncating scGPT gene embeddings to the top-k observed-gene SVD directions improve biological
recovery over the full 512-d representation?
Direct evidence (pre-registered H4). k* was fixed on immune_seed43 before any evaluation context was computed
(kstar.json: PPI 64, NODE 512, EDGE 128). In the four evaluation contexts the layer-mean Delta never has a CI above
zero: PPI -0.0011 [-0.0024,+0.0003] lung, -0.0006 [-0.0019,+0.0007] lung_external, -0.0013 [-0.0026,-0.0000]
kidney_B, +0.0022 [-0.0018,+0.0064] immune (delete-d scaling; pre-registered sqrt(m/n) intervals in delta.csv ci_lo/ci_hi:
[-0.0020,-0.0001], [-0.0015,+0.0003], [-0.0022,-0.0004], [-0.0006,+0.0052]); EDGE +0.0004, +0.0003, -0.0002, -0.0004 with
CIs spanning about +/-0.0025 (lung) to +/-0.0065 (immune); NODE Delta = 0 by construction because the selection context
chose the full representation. Every |Delta| is <= 0.0022 AUROC, an order of magnitude below the between-task and
between-context differences (full PPI 0.649-0.700, NODE 0.634-0.655, EDGE 0.530-0.610). The descriptive joint 'any k
beats full' test (any_k_test.csv) finds no k whose layer-mean gap exceeds label-permutation noise in any context or task
(layer-mean p >= 0.445; any-layer p >= 0.629).
Inference. Top-k curves rise and saturate by k ~ 64-128 (posthoc_best_k.csv: e.g. lung NODE k=10 0.612, k=64 0.654,
full 0.654; lung PPI k=10 0.641, full 0.655). Low-rank summaries lose clearly: at k = 10 the layer-mean loss vs full is
0.011-0.084 AUROC across the 12 evaluation context x task cells (largest immune NODE -0.084, immune EDGE -0.044, lung NODE
-0.042); at k = 3 (the width of the screen's bands) it is 0.032-0.142 (curves_<ctx>.csv, diff_vs_full). The information
the three tasks use is spread over many directions, not concentrated in a few leading ones: rank reduction is at best
neutral (k >= 64) and harmful at small k.
Sensitivity (not pre-registered; sens_R_N_ci.csv, PPI CIs delete-d). Removing the LEADING directions changes the picture
by task: dropping SV1..SV5 raises the H2 TF-identity statistic by +0.099 [+0.080,+0.111] lung, +0.107 [+0.087,+0.122]
lung_external, +0.108 [+0.084,+0.122] kidney_B, +0.092 [+0.035,+0.122] immune (not in kidney_A: -0.007
[-0.090,+0.044]), while it lowers PPI recovery in lung/lung_external/kidney_B (-0.027 [-0.042,-0.013], -0.036
[-0.051,-0.021], -0.017 [-0.033,-0.002]) and leaves EDGE within noise (+0.023 to +0.041, all CIs include 0).
Hypothesis (not tested here): the high-variance leading axes dominate cosine similarity with a gene-wide factor (e.g.
residual expression level / detection within decile, or a cell-composition axis) that carries partner co-location but
masks TF identity; TF identity sits in lower-variance directions. This is post hoc and was chosen among 5 j values x 3
tasks; it should be reported as an effect size, not as a confirmed finding, and the 'remove leading axes' operation is
not what R2.6 asked about.
Screen bands (family N, bands_vs_topk.csv / sens_R_N_ci.csv). Against same-width top-k, the screen's bands are not
uniformly better: in lung/lung_external/kidney_B they are worse for PPI (e.g. SV5-7 vs top-3 -0.036 [-0.053,-0.019]
lung) and better for NODE (SV5-7 +0.045 [+0.012,+0.083] lung), mirroring family R (bands exclude SV1). On immune, where
the bands were found, they beat same-width top-k for PPI by +0.031 to +0.037 (SV2-4 +0.036 [+0.009,+0.063], SV2-7 +0.031
[+0.010,+0.051]; SV5-7 +0.037 [-0.002,+0.075] includes 0 under delete-d scaling). Under the kidney_A submitted-convention
mapping, SV2-4 and SV2-7 become SV1-3 / SV1-6, i.e. exactly top-3 / top-6.
Static models (static_curves.csv, descriptive, lung Tier B 5,226 genes): the same shape holds; full vs best top-k
differs by <= 0.022 AUROC for every model and task (largest: AIDO.Cell NODE best k=64 +0.022, Geneformer V2-104M PPI
best k=7 +0.017), with no consistent k across models (the ESM2 EDGE 'best k = 1' rows are at AUROC ~0.50-0.51).
Implication for the manuscript: answer R2.6 with the three frozen sentences above; state that truncation to the
leading directions neither improves PPI, TF-identity nor edge recovery, that the full representation is the right
default, and that the earlier spectral-band results should not be read as evidence that a low-rank subspace is
privileged for biological recovery.
Limitations / failure modes: (1) the selection context (immune_seed43) shares its tissue and gene axis with immune,
so immune is not independent of selection (immune's PPI Delta is the only positive one); (2) Delta CIs are tiny because
k* = 64/128 representations are nearly the full space; the analysis is informative about 'no gain', not about small-k
losses (those are descriptive); (3) the PPI non-pair set is one rewiring; (4) NODE CIs keep detection strata fixed per
gene; all CIs condition on the fitted SVD (the basis is not re-estimated within resamples); (5) static-model SVDs are
fitted on the evaluated gene set, so 'full' there is the set-centred embedding; (6) the any-k p-values are against a
label-noise null, not a test of H0 metric(k) <= metric(512).
Files: kstar.json, delta.csv, any_k_test.csv, curves_<ctx>.csv (immune_seed43, lung, lung_external, kidney_B, immune,
kidney_A), posthoc_best_k.csv, remove_top_j.csv, bands_vs_topk.csv, sens_R_N_ci.csv, static_curves.csv,
static_summary.csv, results.json; code e7lib.py, e7_run.py, e7_sens.py, e7_static.py, e7_assemble.py; logs e7_*.log;
verifier: verify/e7_verify.py, verify/verify_*.json, verify/rerun*.log, verify/pre_audit_backup/ (pre-audit outputs).

AUDIT (independent adversarial verifier; details in results.json 'adversarial_verification' and verify/audit.json)
Checked and confirmed: primary matches CONVENTION3 §4 (statistic, universe, paired-bootstrap units, B = 2,000, 3 tasks x
4 contexts with CIs and no test, frozen-form sentence); k* chosen on immune_seed43 only, kstar.json (14:05:16) written
before any evaluation-context curve (first: kidney_B 14:07:56); observed-gene SVD (cache obs_idx == non-zero rows of the
canonical embedding at L0 and L11), V-centring, U_k S_k with SV1 = column 0, R_j = columns j..511, bands = columns a-1..z-1;
gene alignment (obs_genes == h5ad var[obs_idx]); joint permutations shared across all k and layers, within detection
strata (NODE) / within TF (EDGE); LOO excludes self under permuted labels; full-dimension NODE and EDGE values equal E3
h2_primary and E4a h3_primary bit-for-bit in all 5 contexts used; STRING edge counts equal E2 prep_meta.
Fixed (2 items): (A) PPI vertex-subsampling intervals: E7 reported only the pre-registered sqrt(m/n) scaling;
CONVENTION3_ERRATA item 3 requires delete-d scaling in the manuscript. Added ci_lo_deleted/ci_hi_deleted and
ci_lo_ms/ci_hi_ms to delta.csv and sens_R_N_ci.csv (pre-registered columns kept). Consequence: layer-mean PPI CI below 0
in 1 of 4 contexts (kidney_B) instead of 2 of 4; per-layer CIs below 0: lung 6/12 (was 9), lung_external 0/12 (was 2),
kidney_B 7/12 (was 10); the decision 'does not improve in 0 of 4' is unchanged. (B) k = 1 cosines: stats3.l2norm_rows
divides by norm + 1e-12, so 1-D cosines (exactly +-1, should tie) were ordered by |SV1|; e7lib.unit_rep now normalises
exactly. Only k = 1 values changed (all other curve values bit-identical; k* unchanged when re-derived from the new
immune_seed43 curves; delta npz untouched and equal to the new curves): k = 1 layer-mean PPI changed by -0.010 to -0.024,
EDGE by -0.012 to +0.003, NODE 0; any-k lung EDGE any-layer max changed from +0.0072 (k=1, L10, p = 0.457) to +0.0021
(k=128, L1, p = 0.934); static ESM2-UCE EDGE best-k (k=1) gap +0.0072 -> +0.0095, ESM2-650M +0.0094 -> +0.0066.
kstar.json is the untouched selection record; its k = 1 curve entries are the pre-fix values (PPI 0.5147, EDGE 0.5321;
corrected 0.5052 / 0.5202), which cannot affect the argmax.
Also corrected in this text: the 'small k loses 0.01-0.04' statement (actual k = 10 loss 0.011-0.084).
Independently reproduced from raw embeddings with own code (verify/e7_verify.py, sklearn AUROCs, own SVD/LOO/
resampling, different seeds): k* = 64/512/128 (all T curves k >= 2 equal to 1e-15); lung layer-mean Delta PPI -0.001051
and EDGE +0.000354; lung NODE R5 - full +0.0991; lung PPI Delta CI (own vertex subsets) pre-registered
[-0.00202,-0.00006], delete-d [-0.00242,+0.00034] vs E7 [-0.00199,-0.00008] / [-0.00237,+0.00032]; lung EDGE Delta CI
(own two-way bootstrap) [-0.0022,+0.0025] vs [-0.0022,+0.0026]; lung_external NODE any-k gap +0.00199 (k=64), own-seed
p = 0.409 vs 0.445 (MC SE ~0.022); static Geneformer V2-104M PPI k=7 0.66785 and full 0.65119 (E7 0.6679 / 0.6512).
```
