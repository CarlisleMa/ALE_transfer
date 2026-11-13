#!/bin/bash
# Oracle of level 4. IGV is the only reader of the alignments, so both solutions go through IGV:
#   default      IGV's command-line route: igvtools count + the batch port (review_cli.py)
#   SOLUTION=2   independent second solution: the IGV desktop application operated with the
#                `desktop` command - alignment export from the track menu, session and images
#                from the File menu (igv_gui.py + sop_review.py)
set -euo pipefail
cd /task
if [ "${SOLUTION:-1}" = "2" ]; then
  python3 /solution/igv_gui.py
else
  python3 /solution/review_cli.py
fi
