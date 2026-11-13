#!/bin/bash
# Oracle of level 4. SOLUTION=2 selects an independent second implementation.
set -euo pipefail
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt_l4.py
else
  python3 /solution/flow_pipeline_l4.py
fi
