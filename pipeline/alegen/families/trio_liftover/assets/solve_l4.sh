#!/bin/bash
# Oracle, level 4. SOLUTION=2 selects an independent second implementation (pure Python: own chain
# lift-over, own re-expression at assembly differences, own codon logic; no Picard, no bcftools).
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt_l4.py
else
  python3 /solution/trio_pipeline_l4.py
fi
