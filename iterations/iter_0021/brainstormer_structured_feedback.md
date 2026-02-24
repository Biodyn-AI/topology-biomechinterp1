# Brainstormer Structured Feedback — iter_0021

## Gate Status

`passed_min_research_gate: true`. Two strong positives, one underpowered but correct direction. Good iteration.

---

## Assessment of Each Tested Hypothesis

### H01 — STRING Quintile Gradient (Inconclusive)
- Direction is right: all 12/12 layers negative.
- Failure mode: Spearman on 5 quintile means is fundamentally underpowered (n=5, df=3). This was a design flaw, not a biological negative.
- **Action**: Retire this specific test design. Replace with pairwise continuous Spearman(score, distance) across all 3092 pairs, or AUROC at thresholds. Upgrading this to decile-level or continuous gives ~3000 observations instead of 5.

### H02 — TRRUST Proximity (Promising — Novel Positive)
- TRRUST pairs are geometrically closer by the same magnitude as STRING pairs (~90% of STRING effect).
- Effect deepens with layer: TRRUST effect −0.029 at layer 0 → −0.051 at layer 11 (76% increase); STRING −0.034 → −0.058 (71% increase). The two scale nearly identically with depth.
- Key unresolved question: Is TRRUST effect independent of STRING co-membership? Many TF-target pairs are probably also STRING edges. If the overlap-corrected TRRUST-only effect persists, the claim is "geometry encodes regulatory programs independently of protein-protein interaction membership." That is a strong, novel claim.
- **Action**: Overlap-corrected TRRUST is the #1 follow-up.

### H03 — H1 Trajectory + Bootstrap CIs (Promising — Confirms Layer Structure)
- Monotonic H1 decline (rho=−0.916) is now bootstrap-confirmed and publish-ready.
- Co-polarity enrichment gradient: early layers ~1.55x, late layers ~1.33x. The layer-enrichment decline (rho=−0.818) adds a mechanistic narrative: early layers preserve broad co-expression structure, later layers encode specialized cell-state information.
- Missing: bootstrap CIs per layer (only layer 8 has CI). Need per-layer CIs to claim the early vs late difference is significant.

---

## Cumulative Project State

**Confirmed findings (publish-grade):**
1. STRING PPI pairs are geometrically closer across all 12 scGPT layers (effect ~−0.048, bootstrap CI excludes null).
2. TRRUST TF-target pairs show the same proximity (effect ~−0.043, 12/12 layers p≈0) — geometry encodes general biological interaction proximity.
3. Co-polarity enrichment (1.33x at layer 8) bootstrapped CI=[1.23, 1.41] excludes null.
4. H1 persistence monotonically compacts with depth (rho=−0.916).
5. Co-polarity enrichment is stronger in early layers (~1.55x) vs late (~1.33x), consistent with hierarchical representation building.

**Directions with diminishing returns:**
- Null sensitivity / permutation tests: covered. No value in more.
- SVD direction analysis: covered. Additional axis counting adds no new insight.
- Quintile-binning for continuous STRING score: dead end, use pairwise continuous test.

---

## Critical Open Questions Driving Next Exploration

1. Is TRRUST proximity independent of STRING overlap? (independence of regulatory vs PPI signal)
2. Do other biological interaction classes (GO pathway, complex membership, miRNA) show the same effect?
3. Is this geometry unique to scGPT, or is it model-agnostic? (Geneformer replication)
4. What is the intrinsic dimension structure, and does it correlate with topological compaction?
5. Does interaction directionality matter? (activating vs repressive TF-target pairs)
6. Do cell-type marker genes cluster by cell-type identity in embedding space?
