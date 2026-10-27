#!/bin/bash
# Lock-based follow-up worker (several may run): for each context whose main pipeline is DONE, run the extra
# sensitivity rows and the power simulation. A context is claimed by creating work/<ctx>/followup.lock.
cd "$(dirname "$0")"
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
CTXS="${1:-lung lung_external immune immune_seed43 immune_seed44 immune_TSP14 immune_TSP21 immune_TSP25 immune_Bcell immune_CD4T immune_CD8T kidney_ckpt_whole_human kidney_ckpt_kidney kidney_ckpt_brain kidney_A kidney_B}"
while true; do
  left=0
  for c in $CTXS; do
    [ -f "work/$c/power.json" ] && [ -f "work/$c/sens_extra.json" ] && continue
    left=1
    [ -f "work/$c/DONE" ] || continue
    mkdir "work/$c/followup.lock" 2>/dev/null || continue
    echo "[$(date +%H:%M:%S)] followup $c"
    $PY e4a_extra.py "$c" && $PY e4a_power.py "$c" || echo "FAILED $c"
  done
  [ $left -eq 0 ] && break
  # stop when every remaining context is claimed by another worker
  unclaimed=0
  for c in $CTXS; do [ -f "work/$c/power.json" ] || [ -d "work/$c/followup.lock" ] || unclaimed=1; done
  [ $unclaimed -eq 0 ] && break
  sleep 60
done
echo "[$(date +%H:%M:%S)] worker exits"
