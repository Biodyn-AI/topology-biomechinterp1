#!/bin/bash
# Resume-time addition: T2 seed 1 for lung, lung_external, kidney_B (restores the full 2-seed T2 design
# if machine load allows). Started after queues C and D exit; two processes share batches via lock files.
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 NUMBA_NUM_THREADS=2
export E4B_NJOBS=3
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
cd "$(dirname "$0")"
while kill -0 39099 2>/dev/null || kill -0 63602 2>/dev/null; do sleep 30; done
J=""
for c in lung lung_external kidney_B; do J="$J $c:T2:genie3:1 $c:T2:grnboost2:1"; done
$PY e4b_grn.py --jobs $J > grn_queueE1.log 2>&1 &
sleep 5
$PY e4b_grn.py --jobs $J > grn_queueE2.log 2>&1 &
wait
