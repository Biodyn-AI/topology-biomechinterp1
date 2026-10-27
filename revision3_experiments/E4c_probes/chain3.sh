#!/bin/bash
# Tuned M2/M3 (added by the resuming agent); starts after this context's chain2 has finished.
ctx=$1
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
P=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
cd "$(dirname "$0")"
while [ ! -f logs/chain_${ctx}.done ]; do sleep 60; done
$P e4c_nonlin_tuned.py $ctx mlp >> logs/mlp_tuned_${ctx}.out 2>&1
$P e4c_nonlin_tuned.py $ctx hgb >> logs/hgb_tuned_${ctx}.out 2>&1
echo "chain3 $ctx finished $(date)" >> logs/chain3_${ctx}.done
