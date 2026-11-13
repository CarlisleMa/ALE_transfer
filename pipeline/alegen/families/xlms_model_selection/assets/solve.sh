#!/bin/bash
# Oracle. SOLUTION=2 selects an independent second implementation (no Biopython, different mapping method).
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt.py
else
  python3 /solution/xl_pipeline.py
fi
