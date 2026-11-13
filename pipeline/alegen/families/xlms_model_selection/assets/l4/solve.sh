#!/bin/bash
# Oracle for level 4. SOLUTION=2 selects an independent second implementation (no Biopython, no
# scipy; own PDB reader, block matching instead of alignment, bucketed lattice search in numpy).
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt.py
else
  python3 /solution/xl4_pipeline.py
fi
