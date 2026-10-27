#!/bin/zsh
# usage: run_job.sh <script> <ctx> [task]   -- skips a job whose log already ends in DONE
export TMPDIR="/Volumes/Crucial X6/MacBook/biomechinterp/.tmp_r2" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
cd "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/revision3_experiments/E3_node"
PY=/Users/ihorkendiukhov/anaconda3/envs/subproject40-topology/bin/python
name=$(basename $1 .py | sed 's/e3_//')_$2${3:+_$3}
if grep -q "DONE" logs/$name.log 2>/dev/null; then echo "skip $name"; exit 0; fi
$PY -B $1 $2 $3 > logs/$name.log 2>&1
echo "$(date +%H:%M:%S) finished $name rc=$?"
