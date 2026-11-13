"""Wait (read-only) until a text file has at least N lines matching a regex, or a time limit passes.

  python3 wait_re.py <file> <regex> <min_matches> <max_seconds>
"""
import re
import sys
import time
from pathlib import Path

f, rx, need, limit = Path(sys.argv[1]), re.compile(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
t0 = time.time()
hits = []
while time.time() - t0 < limit:
    hits = [l for l in (f.read_text(errors="replace").splitlines() if f.exists() else []) if rx.search(l)]
    if len(hits) >= need:
        break
    time.sleep(20)
print(f"{len(hits)} matching line(s) (wanted {need}) at", time.strftime("%H:%M:%S", time.gmtime()))
print("\n".join(l[:200] for l in hits[-12:]))
