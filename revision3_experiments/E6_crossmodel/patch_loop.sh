#!/bin/bash
# re-run the 1-D patch on newly finished jobs every 3 min while e6_run workers are alive, then once more
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 NUMBA_NUM_THREADS=2
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
while pgrep -f "e6_run.py" >/dev/null; do $PY e6_patch1d.py >> patch1d.log 2>&1; sleep 180; done
$PY e6_patch1d.py >> patch1d.log 2>&1
echo "patch loop done $(date)" >> patch1d.log
