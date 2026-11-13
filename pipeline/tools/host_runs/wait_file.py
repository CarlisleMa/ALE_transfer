"""Wait (read-only) until a text file contains a substring, or a time limit passes; then print the file's tail.

  python3 wait_file.py <file> <substring> <max_seconds> [tail_lines]
"""
import sys
import time
from pathlib import Path

f, needle, limit = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
tail = int(sys.argv[4]) if len(sys.argv) > 4 else 40
t0 = time.time()
found = False
while time.time() - t0 < limit:
    if f.exists() and needle in f.read_text(errors="replace"):
        found = True
        break
    time.sleep(15)
print(("found" if found else "NOT found") + f" {needle!r} in {f.name} at", time.strftime("%H:%M:%S", time.gmtime()))
if f.exists():
    print("\n".join(f.read_text(errors="replace").splitlines()[-tail:]))
