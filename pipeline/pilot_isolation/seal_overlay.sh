#!/bin/bash
# Non-destructive sealing for a pilot (used when the plaintext may not be deleted):
#   1. reset the access time of every file and directory under each <dir> to 2026-01-01
#      (ext4 relatime: the next read of a file, or listing of a directory, sets it to "now");
#   2. hide each <dir> under an empty tmpfs mount, so that nothing on the host can reach it by path.
# unseal_overlay.sh removes the mounts and lists every entry whose access time changed.
#   seal_overlay.sh <state_dir> <dir> [<dir> ...]
set -euo pipefail
V=$1; shift
mkdir -p "$V"
[ -e "$V/SEALED" ] && { echo "already sealed: $V"; exit 1; }
: > "$V/dirs.txt"
for D in "$@"; do D=$(readlink -f "$D"); [ -d "$D" ] || { echo "not a directory: $D"; exit 1; }; echo "$D" >> "$V/dirs.txt"; done
while read -r D; do
  find "$D" -type f -exec touch -a -d '2026-01-01 00:00:00' {} +
  # every file must now carry the old access time (this traversal reads directories, so directories are reset last)
  N=$(find "$D" -type f -newerat '2026-01-02 00:00:00' | wc -l)
  [ "$N" = 0 ] || { echo "access-time reset failed for $N files under $D"; exit 1; }
  find "$D" -depth -type d -exec touch -a -d '2026-01-01 00:00:00' {} +
  echo "$D files=$(find "$D" -type f | wc -l)" >> "$V/counts.tmp"
  find "$D" -depth -type d -exec touch -a -d '2026-01-01 00:00:00' {} +     # the count above listed directories again
  mount -t tmpfs -o size=64m,mode=755 tmpfs "$D"
done < "$V/dirs.txt"
mv "$V/counts.tmp" "$V/counts.sealed"
date -u +%FT%TZ > "$V/SEALED"
echo "sealed $(wc -l < "$V/dirs.txt") directories at $(cat "$V/SEALED")"; cat "$V/counts.sealed"
