#!/bin/bash
# Oracle. SOLUTION=2 selects an independent second implementation (different quantifier, different
# statistics, different resolution of the sample mix-up).
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt.py
else
  python3 /solution/rnaseq_pipeline.py
fi
