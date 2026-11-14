#!/bin/bash
# Quantification of all libraries of the sample sheet with salmon.
# usage: bash 01_quantify.sh <work_dir>
# (PROJECT_INPUT overrides the input location, default /task/input)
set -euo pipefail
IN=${PROJECT_INPUT:-/task/input}
WORK=${1:?usage: 01_quantify.sh <work_dir>}
mkdir -p "$WORK/quant"
salmon index -t "$IN/reference/transcripts.fa.gz" -i "$WORK/index" -k 31 -p 2 > "$WORK/index.log" 2>&1
tail -n +2 "$IN/samples.tsv" | cut -f1,6 | while IFS=$'\t' read -r sid fq; do
  salmon quant -i "$WORK/index" -l U -r "$IN/fastq/$fq" --validateMappings -p 2 \
    -o "$WORK/quant/$sid" > "$WORK/quant/$sid.log" 2>&1 < /dev/null
done
echo "quantified $(ls "$WORK/quant" | grep -vc '\.log$') libraries into $WORK/quant"
