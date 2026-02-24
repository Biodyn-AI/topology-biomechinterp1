# Brainstormer Structured Feedback — iter_0020

## Gate status
`passed_min_research_gate: true` — two strong positives (H01, H03), one confirmed neutral (H02). Clean iteration.

---

## What was learned

### H01 — Multi-axis composite validated
The 2.18x mean-across-layers enrichment survives permutation null (z=4.88, 0/1000 exceed). Layer-8 alone is 1.33x (consistent with mean being pulled by high layers). Magnitude quintile within 3-axis stratum is monotone (Spearman r=0.900). This is a publication-ready claim pending bootstrap CIs across all 12 layers.

### H02 — Attention-SVD dissociation confirmed
SVD count AUROC=0.548 (STRING), attention AUROC=0.582 (TRRUST). Joint predictor is non-synergistic. The signals are orthogonal in biology space (PPI geometry vs TF regulation). This is a useful structural result but no new discovery. Close this line unless a better attention normalization scheme is available.

### H03 — First hard topological result
STRING pairs are closer on unit sphere at ALL 12 layers (mean d=−0.237, p<1e-72 per layer). Permutation null is flat (+0.020). The H0 mean lifetime trajectory in the PH data shows clear embedding convergence across layers (0.667 → 0.242): the cloud is geometrically collapsing as depth increases. H1 feature count is 103–133 across layers. H1 mean lifetime also decreases (0.0127 → 0.0047). These layer trajectories have not yet been analyzed.

---

## Critical observation from H03 PH data

The layer-by-layer H0 mean lifetime decrease is monotone and large (0.667 at L0 → 0.242 at L11). This is not just "embeddings are closer" — it's a progressive geometric organization. Meanwhile the STRING distance effect peaks at L7 (d=−0.254) not at the final layer (L11, d=−0.190). This divergence is scientifically meaningful: deeper layers compress the whole embedding cloud but the STRING-specific clustering signal is strongest at mid-to-deep layers. This is the single most interesting unanalyzed pattern in the iter_0020 data.

---

## Directions this feedback retires or limits

- **Joint attention-SVD predictor**: H02 confirms no synergy. Retire this exact combination. Any future ROC work should use better normalization or head-level decomposition.
- **H0 mean lifetime as a proxy for "topology"**: The trivial H0 result (all embeddings are connected at max filtration) is fully expected. The meaningful quantity is H0 co-clustering at low filtration — i.e., do specific gene communities form early connected components?
- **Further attention-STRING tests**: AUROC=0.482 (below chance). This direction is dead.

---

## What needs to happen next

Priority sequence:
1. **H1 loop membership analysis** — the H1 data is on disk. Assign genes to high-persistence H1 cycles, test GO/STRING enrichment of cycle members.
2. **Layer-stratified bootstrap CIs** for multi-axis composite enrichment — needed for publication figures.
3. **STRING confidence × distance gradient** — extend H03 by stratifying STRING pairs by confidence quintile; test if d-effect scales monotonically with confidence.
4. **CORUM complex H0 co-clustering** — at what filtration threshold do complex members merge? Compare to random gene sets.
5. **Local intrinsic dimension by community** — cheap TwoNN or PCA-decay test.
