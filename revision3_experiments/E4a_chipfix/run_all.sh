#!/bin/bash
# Background runner: 4 workers, BLAS threads 2 each. Usage: bash run_all.sh
cd "$(dirname "$0")"
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
SMALL="immune immune_seed43 immune_seed44 immune_TSP14 immune_TSP21 immune_TSP25 immune_Bcell immune_CD4T immune_CD8T kidney_A"
( $PY e4a_chipfix.py sens lung && $PY e4a_chipfix.py repro lung ) > w1.log 2>&1 &
( $PY e4a_chipfix.py sens lung_external && $PY e4a_chipfix.py repro lung_external ) > w2.log 2>&1 &
( $PY e4a_chipfix.py sens kidney_B kidney_ckpt_whole_human && $PY e4a_chipfix.py repro kidney_B ) > w3.log 2>&1 &
( $PY e4a_chipfix.py sens $SMALL && $PY e4a_chipfix.py repro $SMALL && $PY e4a_chipfix.py sens kidney_ckpt_kidney kidney_ckpt_brain ) > w4.log 2>&1 &
wait
echo ALLDONE > run_all.done
