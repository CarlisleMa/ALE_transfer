#!/bin/bash
# Oracle. SOLUTION=2 selects an independent second implementation (see Harbor "multiple solutions").
set -euo pipefail
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt.py
else
  python3 /solution/flow_pipeline.py
fi
