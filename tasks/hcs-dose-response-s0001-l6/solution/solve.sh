#!/bin/bash
# Oracle (level 4): classical pipeline without Cellpose. SOLUTION=2 selects the independent second
# implementation: Cellpose with settings validated on the annotated fields, run as stages.
set -euo pipefail
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  H="python3 /solution/hcs_l4_cellpose.py"
  $H prepare
  # counting procedure: Cellpose settings are tried on the annotated fields (SOP section 2) ...
  $H try --model nuclei --diameter 17 --flow 0.4
  $H try --model nuclei --diameter 12 --flow 0.4
  $H try --model nuclei --diameter 12 --flow 0
  # ... and the one that meets the criterion on every field with the smallest worst-case error is used
  $H choose
  $H segment
  $H measure
  $H analyse
else
  python3 /solution/hcs_l4_voting.py
fi
