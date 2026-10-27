#!/bin/bash
# Hash-seed spread of the submitted (round-2) residualisation: PYTHONHASHSEED 0..20 (point estimates only).
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
cd "$(dirname "$0")"
for s in $(seq 0 20); do
  [ -f continuity_hashseed_spread/hashseed_$s.csv ] && continue
  PYTHONHASHSEED=$s /Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python e5_continuity.py spread
done
echo SPREAD ALL DONE
