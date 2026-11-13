#!/bin/bash
# Oracle. SOLUTION=2 selects an independent second implementation (see Harbor "multiple solutions").
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt.py
else
  python3 /solution/crispr_pipeline.py
fi
