#!/bin/bash
# Job queues (see SUMMARY.md, "Compute plan").
#   A: T1 (5 seeds) for lung, lung_external, kidney_B; n_jobs=3
#   B: T1 (5 seeds) + T2 (2 seeds) for immune, kidney_A; n_jobs=2
#   C: T2 (seed 0 only -- fallback, see SUMMARY) for lung, lung_external, kidney_B; n_jobs=3; started after B ends
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
cd "$(dirname "$0")"
Q=$1
J=""
if [ "$Q" = "A" ]; then
  for c in lung lung_external kidney_B; do for s in 0 1 2 3 4; do J="$J $c:T1:genie3:$s $c:T1:grnboost2:$s"; done; done
  export E4B_NJOBS=3
elif [ "$Q" = "B" ]; then
  for c in immune kidney_A; do for s in 0 1 2 3 4; do J="$J $c:T1:genie3:$s $c:T1:grnboost2:$s"; done
    for s in 0 1; do J="$J $c:T2:genie3:$s $c:T2:grnboost2:$s"; done; done
  export E4B_NJOBS=2
elif [ "$Q" = "C" ]; then
  for c in lung lung_external kidney_B; do J="$J $c:T2:genie3:0 $c:T2:grnboost2:0"; done
  export E4B_NJOBS=3
else
  # D: the same T2 jobs in the same priority order (primary lung first), started when A ends; batches are
  #    claimed with lock files, so C and D never fit the same batch and either may assemble a job
  for c in lung lung_external kidney_B; do J="$J $c:T2:genie3:0 $c:T2:grnboost2:0"; done
  export E4B_NJOBS=3
fi
$PY e4b_grn.py --jobs $J
