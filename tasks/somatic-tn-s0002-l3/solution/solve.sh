#!/bin/bash
# Oracle. SOLUTION=2 selects an independent second implementation (different aligner, different
# counting engine, different pairing statistic).
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt.py
else
  python3 /solution/somatic_pipeline.py
fi
