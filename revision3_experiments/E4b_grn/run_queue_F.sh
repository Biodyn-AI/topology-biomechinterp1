#!/bin/bash
# Post-reboot restart of the remaining optional T2 seed-1 jobs (stale locks cleared before launch).
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 NUMBA_NUM_THREADS=2
export E4B_NJOBS=3
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
cd "$(dirname "$0")"
J="lung_external:T2:grnboost2:1 kidney_B:T2:genie3:1 kidney_B:T2:grnboost2:1"
$PY e4b_grn.py --jobs $J > grn_queueF1.log 2>&1 &
sleep 5
$PY e4b_grn.py --jobs $J > grn_queueF2.log 2>&1 &
wait
echo "queue F done $(date)" >> grn_queueF_done.log
