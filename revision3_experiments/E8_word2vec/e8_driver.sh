#!/bin/bash
# E8 orchestration after the three production training processes (keeps <= 3 python jobs at a time):
#   when processes 2 and 3 finish -> start calibration training (e8_train.py 4 and 5);
#   when process 1 finishes       -> evaluate lung_external, kidney_B, kidney_A;
#   when calibration finishes     -> evaluate lung, immune; then merge.
cd "$(dirname "$0")"
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
running() { pgrep -f "e8_train.py $1\$" >/dev/null; }
while running 2 || running 3; do sleep 30; done
echo "$(date +%T) starting calibration training"
nohup $PY e8_train.py 4 > e8_train_p4.log 2>&1 &
nohup $PY e8_train.py 5 > e8_train_p5.log 2>&1 &
sleep 10
while running 1; do sleep 30; done
for c in lung_external kidney_B kidney_A; do
  [ -f eval/$c.json ] || $PY e8_eval.py $c >> e8_eval.log 2>&1
done
while running 4 || running 5; do sleep 30; done
for c in lung immune; do
  [ -f eval/$c.json ] || $PY e8_eval.py $c >> e8_eval.log 2>&1
done
$PY e8_merge.py >> e8_merge.log 2>&1
echo "$(date +%T) DRIVER DONE"
