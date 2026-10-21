#!/usr/bin/env bash
# run mob_recon on every genome id listed in one chunk file
# usage: bash plasmid_run_mob_suite.sh <chunk_name>
set -u

MOB_ROOT="${PHASE10_DIR:-./phase10}/mob_out"
ASSIGN_ROOT="${PHASE10_DIR:-./phase10}/assignment/mac_chunks"
INDEX_FILE="${PHASE10_DIR:-./phase10}/assignment/genome_index.tsv"

CHUNK_NAME="$1"
CHUNK="$ASSIGN_ROOT/$CHUNK_NAME"
[ -f "$CHUNK" ] || { echo "ERROR: chunk not found: $CHUNK"; exit 1; }

mkdir -p "$MOB_ROOT/_logs"
STATUS_FILE="$MOB_ROOT/_logs/${CHUNK_NAME}.status.tsv"

while IFS= read -r id; do
    [ -z "$id" ] && continue

    genome_path=$(awk -F'\t' -v id="$id" '$1==id {print $2; exit}' "$INDEX_FILE")
    if [ -z "$genome_path" ]; then
        echo "ERROR: genome not found: $id"
        continue
    fi

    # contig_report.txt is the last thing written, so its presence means it finished
    if [ -f "$MOB_ROOT/$id/contig_report.txt" ]; then
        echo "SKIP $id - already complete"
        continue
    fi

    echo "===== START $id $(date +%T) ====="
    start=$(date +%s)

    mob_recon -i "$genome_path" -o "$MOB_ROOT/$id" -n 1 --force \
        > "$MOB_ROOT/_logs/${id}.log" 2>&1
    exit_code=$?

    end=$(date +%s)
    printf '%s\t%s\t%s\t%s\n' "$id" "$exit_code" "$((end-start))" "$(date -u +%FT%TZ)" >> "$STATUS_FILE"
    echo "===== FINISH $id exit=$exit_code ($((end-start))s) ====="
done < "$CHUNK"