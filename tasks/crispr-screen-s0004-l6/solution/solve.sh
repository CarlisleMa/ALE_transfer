#!/bin/bash
# Oracle, level 4. SOLUTION=2 selects an independent second implementation (Harbor "multiple solutions").
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_l4_alt.py
else
  python3 /solution/crispr_l4_pipeline.py
fi
