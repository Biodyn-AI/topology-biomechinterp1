#!/bin/bash
# Runs the power simulation and the extra sensitivity rows for each context as soon as its main pipeline
# has finished (work/<ctx>/DONE). Resumable: both scripts skip finished contexts.
cd "$(dirname "$0")"
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
CTXS="kidney_A immune immune_seed43 immune_seed44 immune_TSP14 immune_TSP21 immune_TSP25 immune_Bcell immune_CD4T immune_CD8T kidney_B lung lung_external kidney_ckpt_whole_human kidney_ckpt_kidney kidney_ckpt_brain"
$PY e4a_power.py legacy || echo "FAILED legacy"
remaining="$CTXS"
while [ -n "$remaining" ]; do
  next=""
  for c in $remaining; do
    if [ -f "work/$c/DONE" ]; then
      echo "[$(date +%H:%M:%S)] followup $c"
      $PY e4a_extra.py "$c" && $PY e4a_power.py "$c" || echo "FAILED $c"
    else
      next="$next $c"
    fi
  done
  remaining=$(echo $next)
  [ -n "$remaining" ] && sleep 60
done
echo "[$(date +%H:%M:%S)] followup all done"
