#!/bin/bash
# Oracle of level 7. SOLUTION=2 selects an independent second implementation.
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt_l7.py
else
  python3 /solution/amplicon_pipeline_l7.py
fi
