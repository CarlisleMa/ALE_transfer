#!/bin/bash
# Oracle (level 7): the level-4 classical pipeline; before the plate statistics the control wells are
# checked under both orientations of the plate and the layout is applied in the one in which they behave.
# SOLUTION=2 selects the independent second implementation: Cellpose with settings validated on the
# annotated fields, orientation from the direction of the dose series, run as stages.
set -euo pipefail
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  H="python3 /solution/hcs_l7_cellpose.py"
  $H prepare
  # counting procedure: Cellpose settings are tried on the annotated fields (SOP section 2) ...
  $H try --model nuclei --diameter 17 --flow 0.4
  $H try --model nuclei --diameter 12 --flow 0.4
  $H try --model nuclei --diameter 12 --flow 0
  # ... and the one that meets the criterion on every field with the smallest worst-case error is used
  $H choose
  $H segment
  $H measure
  # which way round does the layout fit the images?
  $H orient
  $H analyse
else
  python3 /solution/hcs_l7_voting.py
fi
