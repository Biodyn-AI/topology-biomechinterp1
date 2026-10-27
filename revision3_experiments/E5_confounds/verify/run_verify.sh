#!/bin/bash
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 NUMBA_NUM_THREADS=2
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
cd "$(dirname "$0")"
if [ "$1" = A ]; then
  [ -f verify_clogit_kA_L11.json ] || $PY verify_clogit.py kA_L11 kidney_A 11 500 20 > log_kA_L11.log 2>&1
else
  for spec in "imm_L2 immune 2 20 10" "lung_L0 lung 0 3 2" "lungext_L1 lung_external 1 0 0" "kB_L0 kidney_B 0 0 0" "lung_L11 lung 11 0 0"; do
    set -- $spec; [ -f verify_clogit_$1.json ] || $PY verify_clogit.py $spec > log_$1.log 2>&1
  done
fi
echo ALL DONE
