# iter_0009 Next Steps

## Top priorities for iter_0010

### 1. SV2 top-pole biological integration (HIGH PRIORITY)
- iter_0008 found SV2 top-pole = IL-4 immune axis (GO:0032753 IL-4 production regulation, OR=Inf, p=9.5e-6, emp_p=0.000).
- Identify specific gene identities in top-52 SV2 pole.
- Test IL-4 pathway markers (IL4R, STAT6, CD23/FCER2E) presence.
- Run additional immune-signaling GO terms against SV2 top pole to fully characterize the axis.

### 2. H02 scan extension: extracellular vesicle on SV2 across all layers (HIGH PRIORITY)
- H02 ran SV1 only. H01 showed EV is on SV2. Run the same 12-layer scan for SV2/SV3 projections.
- Specifically: does EV enrichment on SV2 bottom pole appear at all layers, or only at late layers?
- Expected result: present early, strengthens late.

### 3. SV1 axis: full gene list + TRRUST TF regulatory test (NOVEL METHOD)
- SV1 top-52 genes are secretory pathway (signal peptide, ER lumen, extracellular).
- Test if known TF regulators of secretory pathway (CREB3, XBP1, ATF6, MIST1) are in the top SV1 pole.
- Test TRRUST TF→target structure: do top-SV1 targets cluster by TF (assortativity test)?

### 4. Geodesic vs Euclidean test for regulatory proximity (FAMILY 3, NOT YET TESTED)
- Retired: none in this family.
- Test: for TRRUST TF→target pairs, compare L2 Euclidean distance vs shortest-path in kNN graph.
- Hypothesis: regulatory partners are geodesically closer than expected from Euclidean alone.

### 5. Null sensitivity: annotation density confounder check
- H02 showed secreted and ER_lumen are persistently significant. Possible confounder: well-annotated genes may systematically cluster in the SV1 axis.
- Run: correlation of GO annotation count per gene vs SV1 projection value.
- If annotation density → SV1 position, results are confounded.

## Retirement status
- module_structure (standalone TF co-clustering): 2+ negative/inconclusive → RETIRED unless novel rescue
- graph_topology transitivity: 1 negative → monitor
