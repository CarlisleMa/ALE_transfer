#!/bin/bash
# Level-4 oracle. SOLUTION=2 selects the independent second solution (own k-mer quantifier, voom-style
# statistics, exclusion instead of relabelling, audit without re-running salmon).
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt_l4.py
else
  python3 /solution/takeover_pipeline.py
fi
