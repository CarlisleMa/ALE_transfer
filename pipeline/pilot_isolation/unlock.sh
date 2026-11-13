#!/bin/bash
# Report whether the sealed archive or key was read since lock.sh, then restore the plaintext.
#   unlock.sh <vault_dir> <key_file>
set -euo pipefail
V=$1; KEY=$2
[ -e "$V/LOCKED" ] || { echo "no locked vault at $V"; exit 1; }
stat -c '%n atime=%X mtime=%Y' "$V"/part.*.enc "$KEY" > "$V/atimes.now"
if cmp -s "$V/atimes.sealed" "$V/atimes.now"; then
  echo "AUDIT: archive and key were NOT read between lock ($(cat "$V/LOCKED")) and now ($(date -u +%FT%TZ))"
else
  echo "AUDIT: ACCESS DETECTED"; diff "$V/atimes.sealed" "$V/atimes.now" || true
fi
i=0
while read -r D; do
  i=$((i+1)); P=$(dirname "$D"); B=$(basename "$D")
  [ -e "$D" ] && { echo "refusing to overwrite existing $D"; exit 1; }
  mkdir -p "$P"
  openssl enc -d -aes-256-cbc -pbkdf2 -pass "file:$KEY" -in "$V/part.$i.enc" | tar -C "$P" -xf -
  (cd "$P" && find "$B" -type f -print0 | sort -z | xargs -0 -r sha256sum) | cmp -s - "$V/manifest.$i" \
    && echo "restored and verified: $D" || { echo "RESTORE MISMATCH for $D"; exit 1; }
done < "$V/dirs.txt"
mv "$V/LOCKED" "$V/UNLOCKED.$(date -u +%H%M%S)"
rm -f "$V"/part.*.enc "$KEY"
