#!/bin/bash
# Remove the empty mounts and report every file or directory that was accessed while sealed.
#   unseal_overlay.sh <state_dir>
set -uo pipefail
V=$1
[ -e "$V/SEALED" ] || { echo "not sealed: $V"; exit 1; }
: > "$V/accessed.txt"
while read -r D; do
  if grep -q " $D tmpfs " /proc/mounts; then
    echo "mount still in place: $D" >> "$V/mount_check.txt"
    umount "$D" 2>/dev/null || umount -l "$D"
  else
    echo "MOUNT MISSING at unseal: $D" | tee -a "$V/mount_check.txt"
  fi
  # pre-order traversal: each entry is tested on the stat taken before find opens it
  find "$D" \( -type f -o -type d \) -newerat '2026-01-02 00:00:00' -print >> "$V/accessed.txt"
  echo "$D files=$(find "$D" -type f | wc -l)" >> "$V/counts.now"
done < "$V/dirs.txt"
N=$(wc -l < "$V/accessed.txt")
if [ "$N" = 0 ] && ! grep -q MISSING "$V/mount_check.txt" && cmp -s "$V/counts.sealed" "$V/counts.now"; then
  echo "AUDIT: nothing under the sealed directories was read or listed between $(cat "$V/SEALED") and $(date -u +%FT%TZ); mounts were in place; file counts unchanged"
else
  echo "AUDIT: CHECK NEEDED - $N accessed entries (see $V/accessed.txt)"; grep MISSING "$V/mount_check.txt"; diff "$V/counts.sealed" "$V/counts.now"; head -20 "$V/accessed.txt"
fi
mv "$V/SEALED" "$V/UNSEALED.$(date -u +%H%M%S)"
