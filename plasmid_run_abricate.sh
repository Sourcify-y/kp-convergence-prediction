#!/usr/bin/env bash
# run abricate (plasmidfinder + card) on every genome id listed in one chunk file
# usage: bash plasmid_run_abricate.sh <chunk_name>
set -u
CHUNK_NAME="$1"
CHUNK="${PHASE10_DIR:-./phase10}/assignment/abricate_chunks/$CHUNK_NAME"
OUT_DIR="${PHASE10_DIR:-./phase10}/abricate_out"
INDEX_FILE="${PHASE10_DIR:-./phase10}/genome_path_index.tsv"

mkdir -p "$OUT_DIR/plasmidfinder" "$OUT_DIR/card" "$OUT_DIR/_logs"

while IFS= read -r id; do
  [ -z "$id" ] && continue

  # skip genomes that already have both outputs
  if [ -f "$OUT_DIR/plasmidfinder/${id}.tab" ] && [ -f "$OUT_DIR/card/${id}.tab" ]; then
    continue
  fi

  genome_path=$(awk -F'\t' -v id="$id" '$1==id{print $2; exit}' "$INDEX_FILE")
  if [ -z "$genome_path" ]; then
    echo "ERROR: no path found for $id" >> "$OUT_DIR/_logs/errors.log"
    continue
  fi

  abricate --db plasmidfinder --quiet "$genome_path" > "$OUT_DIR/plasmidfinder/${id}.tab" 2>> "$OUT_DIR/_logs/errors.log"
  abricate --db card --quiet "$genome_path" > "$OUT_DIR/card/${id}.tab" 2>> "$OUT_DIR/_logs/errors.log"
done < "$CHUNK"