#!/bin/bash
# Seal answer material before a pilot: encrypt the given directories, verify the archive by a full
# decrypt-and-compare, and only then remove the plaintext. Reading the archive or the key afterwards
# changes their access time (ext4 relatime), which unlock.sh reports.
#   lock.sh <vault_dir> <key_file> <dir> [<dir> ...]
set -euo pipefail
V=$1; KEY=$2; shift 2
mkdir -p "$V"
[ -e "$V/LOCKED" ] && { echo "vault at $V is already locked"; exit 1; }
head -c 48 /dev/urandom | base64 > "$KEY"; chmod 600 "$KEY"
i=0
: > "$V/dirs.txt"
for D in "$@"; do
  D=$(readlink -f "$D"); [ -d "$D" ] || { echo "not a directory: $D"; exit 1; }
  i=$((i+1)); P=$(dirname "$D"); B=$(basename "$D")
  (cd "$P" && find "$B" -type f -print0 | sort -z | xargs -0 -r sha256sum) > "$V/manifest.$i"
  tar -C "$P" -cf - "$B" | openssl enc -aes-256-cbc -pbkdf2 -pass "file:$KEY" -out "$V/part.$i.enc"
  T=$(mktemp -d "$V/verify.XXXXXX")
  openssl enc -d -aes-256-cbc -pbkdf2 -pass "file:$KEY" -in "$V/part.$i.enc" | tar -C "$T" -xf -
  (cd "$T" && find "$B" -type f -print0 | sort -z | xargs -0 -r sha256sum) | cmp -s - "$V/manifest.$i" \
    || { echo "VERIFY FAILED for $D - plaintext left untouched"; rm -rf "$T"; exit 1; }
  rm -rf "$T"
  echo "$D" >> "$V/dirs.txt"
done
# every archive verified: now remove the plaintext
while read -r D; do rm -rf "$D"; done < "$V/dirs.txt"
touch -d '2026-01-01 00:00:00' "$V"/part.*.enc "$KEY"
stat -c '%n atime=%X mtime=%Y' "$V"/part.*.enc "$KEY" > "$V/atimes.sealed"
date -u +%FT%TZ > "$V/LOCKED"
echo "locked $(wc -l < "$V/dirs.txt") directories at $(cat "$V/LOCKED")"
