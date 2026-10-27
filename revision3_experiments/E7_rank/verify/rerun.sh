#!/bin/bash
source "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop/revision3_experiments/E7_rank/env.sh"
cd "$E7"
$PY e7_run.py curves --ctx immune_seed43 > verify/rerun_curves_immune_seed43.log 2>&1
( $PY e7_run.py eval --ctx lung > verify/rerun_eval_lung.log 2>&1; $PY e7_run.py eval --ctx kidney_B > verify/rerun_eval_kidney_B.log 2>&1 ) &
( $PY e7_run.py eval --ctx lung_external > verify/rerun_eval_lung_external.log 2>&1; $PY e7_run.py eval --ctx immune > verify/rerun_eval_immune.log 2>&1; $PY e7_run.py eval --ctx kidney_A > verify/rerun_eval_kidney_A.log 2>&1 ) &
$PY e7_static.py > verify/rerun_static.log 2>&1
wait
echo ALLDONE > verify/rerun.done
