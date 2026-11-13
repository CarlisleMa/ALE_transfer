#!/bin/bash
# Oracle. SOLUTION=2 selects an independent second implementation (pure Python, no Picard/bcftools).
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt.py
else
  python3 /solution/trio_pipeline.py
fi
