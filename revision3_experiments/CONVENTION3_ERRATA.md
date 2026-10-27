# CONVENTION3 errata (no specification changed)

Recorded after gate A2; CONVENTION3.md itself stays frozen (hash in FROZEN.txt).

1. CONVENTION3.md §1, `contexts.py` and `prep_cache.py` describe every context as "256 sampled cells".
   kidney_B (seed 67) and the three kidney checkpoint contexts (seed 87) used **192** cells. Every
   specification is defined per context in terms of its own extraction, so no analysis changes.
2. CONVENTION3.md §1 lists kidney_A and kidney_B as two extractions of TS kidney. Gate A2 established
   more precisely that kidney_B's file is the first 8,000 rows of kidney_A's file (one donor, TSP2):
   the two are extractions of overlapping cells with different HVG axes.
3. **Vertex-subsampling interval scale (CONVENTION3 §6; `stats3.vertex_subsample_ci`).** The convention specifies
   vertex subsampling without replacement at m = n/2 with √(m/n) rescaling. √(m/n) is the subsampling rate for
   m ≪ n; it ignores the finite-population correction, which at m = n/2 is not negligible. The variance of a
   statistic on an m-of-n subsample without replacement exceeds that of the full-sample statistic by the
   delete-d jackknife factor (n − m)/m, so the correct rescaling of deviations is √(m/(n − m)) (= 1 at m = n/2),
   not √(m/n) (= 0.71). Found by the E2 verifier: the √(m/n) intervals imply SEs 0.71× the delete-d and delete-1
   jackknife SEs in all five contexts (lung 0.0180 vs 0.0255 vs 0.0255). **All network-statistic intervals in the
   manuscript use the delete-d scaling** (equivalently, the pre-registered half-widths × √2); the pre-registered
   values are kept in the CSVs. No decision rule depends on these intervals except the replication criterion,
   which holds under either scaling by > 10 SE.
