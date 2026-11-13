#!/bin/bash
# Oracle. SOLUTION=2 selects the independent second tool chain (no pdb4amber / tleap).
set -euo pipefail
cd /task
export OPENMM_CPU_THREADS="${OPENMM_CPU_THREADS:-2}"
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/solve_alt.py
else
  python3 /solution/build_tleap.py
fi
