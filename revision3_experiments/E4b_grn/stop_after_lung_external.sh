#!/bin/bash
# Resume-time compute reduction (see SUMMARY.txt): T2 for kidney_B is dropped. Queues C and D are stopped as
# soon as the last lung_external T2 job is assembled, before they spend hours on kidney_B T2.
cd "$(dirname "$0")"
until [ -f partial/lung_external_T2_grnboost2_s0.npz ]; do sleep 20; done
for p in 39099 63602; do
  kids=$(pgrep -P $p)
  kill $p 2>/dev/null; sleep 2; [ -n "$kids" ] && kill $kids 2>/dev/null
done
echo "$(date '+%H:%M:%S') stopped queues C (39099) and D (63602) after lung_external T2; kidney_B T2 dropped" >> stop_after_lung_external.log
