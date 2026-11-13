#!/bin/bash
# Copy a pilot container's deliverables and its filesystem diff out (no answer material involved).
#   collect.sh <task_id>
set -uo pipefail
T=$1; C=pilot-$T; O=/opt/alegen_pilot/out/$T
rm -rf "$O"; mkdir -p "$O"
docker cp "$C:/task/output" "$O/output" >/dev/null 2>&1 || mkdir -p "$O/output"
docker diff "$C" 2>/dev/null | grep -v -E '__pycache__|\.pyc$' > "$O/docker_diff.txt"
docker exec "$C" bash -lc 'date -u +%FT%TZ; ls -la --time-style=+%H:%M:%S /task/output 2>/dev/null' > "$O/output_listing.txt" 2>&1
echo "$T: $(find "$O/output" -type f | wc -l) output files, $(wc -l < "$O/docker_diff.txt") changed paths"
