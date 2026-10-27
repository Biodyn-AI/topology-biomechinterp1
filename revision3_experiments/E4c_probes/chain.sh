#!/bin/bash
# Sequential E4c job chain for one context (keeps this experiment at <= 3 concurrent processes).
ctx=$1
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
P=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
cd "$(dirname "$0")"
while pgrep -f "e4c_m1.py $ctx --phase true" >/dev/null; do sleep 30; done
$P e4c_m1.py $ctx --phase true      >> logs/m1_${ctx}_true.out 2>&1
$P e4c_m1.py $ctx --phase null      >> logs/m1_${ctx}_null.out 2>&1
$P e4c_m1.py $ctx --phase ext_node  >> logs/m1_${ctx}_ext.out 2>&1
if [ "$ctx" = "kidney_B" ]; then $P e4c_m1.py immune --phase true >> logs/m1_immune_true.out 2>&1; fi
$P e4c_nonlin.py $ctx mlp --phase true >> logs/mlp_${ctx}.out 2>&1
$P e4c_nonlin.py $ctx hgb --phase true >> logs/hgb_${ctx}.out 2>&1
$P e4c_nonlin.py $ctx mlp --phase null >> logs/mlp_${ctx}.out 2>&1
$P e4c_nonlin.py $ctx hgb --phase null --features concat >> logs/hgb_${ctx}.out 2>&1
echo "chain $ctx finished $(date)" >> logs/chain_${ctx}.done
