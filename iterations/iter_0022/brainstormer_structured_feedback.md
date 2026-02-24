# Brainstormer Structured Feedback — iter_0022

## Gate Status
`passed_min_research_gate: true`. All three hypotheses positive. Clean iteration.

---

## Signal Assessment

### H01 — STRING continuous AUROC + TRRUST-exclusive proximity
**Strong.** Continuous Spearman (rho=−0.093, N=3092, all 12 layers) replaces the underpowered quintile test. More important: TRRUST-exclusive pairs (N=141, NOT in STRING) retain proximity (effect=−0.030, AUROC=0.573, all 12 layers p<0.007). This is the key independence claim — regulatory program geometry is not a STRING artifact. This is now publish-worthy for the overlap-corrected TRRUST result.

**What's missing**: Direction asymmetry (activation vs repression) and TF hub structure are un-tested and would add mechanistic depth. The STRING decile test would pin down whether there's a nonlinear threshold or it's truly linear throughout.

### H02 — Cell-type marker gene cluster separation
**Transformative.** AUROC=0.853 is 40 AUROC points above STRING proximity. Effect deepens monotonically layer 0→11 (−0.164→−0.276). This is the headline finding. Three things remain open:
1. Is AUROC=0.853 specific to the chosen 14 marker genes, or does it hold across a broader cell-type vocabulary? Must test with 5-6 cell types (add macrophage, NK, myeloid, endothelial).
2. Contamination control: do non-marker named genes show any residual cell-type signal?
3. Cross-model: does Geneformer replicate this geometry? If yes, the finding is architecture-agnostic.

### H03 — Co-polarity bootstrap (all 12 layers)
**Confirmed and closeable.** Mean ratio=1.542, all 12 CIs exclude 1.0. No layer trend. This result is publication-ready as-is. The proxy entropy (rho=−0.434, p=0.16) is insufficient — true persistence entropy requires re-running ripser and storing H1 lifetime arrays. That's a separate, more expensive task.

---

## Methodological Note
The gene index bug (`.split()` vs `line.strip()`) was discovered and fixed within the iteration. Critical: verify all prior iterations used `line.strip()`. The executor confirmed consistency with iter_0020 stored results. No retroactive correction needed unless re-inspection reveals discrepancy.

---

## Paper Narrative Update Needed
The paper's primary claim should now be restructured around three tiers:
1. **Manifold existence** (kNN clustering, CKA layer-invariance, intrinsic dimensionality) — established.
2. **Biological anchoring** (STRING PPI, TRRUST regulatory) — established, STRING-independent.
3. **Cell-type identity encoding** (AUROC=0.853, layer-deepening) — NEW, strongest signal, should become the lead result.

The cell-type result reframes the paper from "geometry exists" to "geometry encodes cell-type identity progressively across transformer depth."
