#!/usr/bin/env bash
set -euo pipefail

# ---- user config ----
LIST="/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/pf2_list_faa.txt"
OUT_BASE="/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset"
OUT_PREFIX="result"          # produces: result1, result2, ...
BATCH_SIZE=30

PF2_BIN="pathogenfinder2"

# ---- checks ----
[[ -s "$LIST" ]] || { echo "[err] List not found or empty: $LIST"; exit 1; }
mkdir -p "$OUT_BASE"

# ---- temp workspace ----
CHUNK_DIR="$(mktemp -d)"
trap 'rm -rf "$CHUNK_DIR"' EXIT

# Normalize the list: drop blank lines, strip CRLF
CLEAN_LIST="$CHUNK_DIR/clean.txt"
awk 'NF' "$LIST" | sed 's/\r$//' > "$CLEAN_LIST"
TOTAL=$(wc -l < "$CLEAN_LIST")
[[ "$TOTAL" -gt 0 ]] || { echo "[err] No valid entries after cleaning $LIST"; exit 1; }

echo "[info] Total inputs: $TOTAL → batching by $BATCH_SIZE"

# Split into BATCH_SIZE-line chunks: pf2_0000, pf2_0001, ...
split -d -a 4 -l "$BATCH_SIZE" "$CLEAN_LIST" "$CHUNK_DIR/pf2_"

# Count chunks
NUM_CHUNKS=$(ls -1 "$CHUNK_DIR"/pf2_* 2>/dev/null | wc -l | tr -d ' ')
echo "[info] Chunks created: $NUM_CHUNKS"

# ---- resume logic: find highest existing resultN, delete it, resume from N ----
shopt -s nullglob
max_done=0
for d in "$OUT_BASE"/"$OUT_PREFIX"*; do
  [[ -d "$d" ]] || continue
  base="$(basename "$d")"
  if [[ "$base" == ${OUT_PREFIX}[0-9]* ]]; then
    n="${base#${OUT_PREFIX}}"
    [[ "$n" =~ ^[0-9]+$ ]] && { (( n > max_done )) && max_done="$n"; }
  fi
done
shopt -u nullglob

start_i=1
if (( max_done > 0 )); then
  echo "[resume] Found last result: ${OUT_PREFIX}${max_done}"
  if [[ -d "$OUT_BASE/${OUT_PREFIX}${max_done}" ]]; then
    echo "[resume] Deleting (assumed incomplete): $OUT_BASE/${OUT_PREFIX}${max_done}"
    rm -rf "$OUT_BASE/${OUT_PREFIX}${max_done}"
  fi
  start_i="$max_done"
  echo "[resume] Will resume from batch index: $start_i"
else
  echo "[resume] No existing results. Starting from batch 1."
fi

(( start_i > NUM_CHUNKS )) && { echo "[warn] Resume index > chunks. Using last chunk."; start_i="$NUM_CHUNKS"; }

# ---- run batches ----
i=1
for chunk in "$CHUNK_DIR"/pf2_*; do
  [[ -s "$chunk" ]] || { i=$((i+1)); continue; }

  if (( i < start_i )); then
    echo "[skip] Batch $i already completed → ${OUT_PREFIX}${i}"
    i=$((i+1))
    continue
  fi

  LINES=$(wc -l < "$chunk")
  OUT_DIR="${OUT_BASE}/${OUT_PREFIX}${i}"
  echo "[info] Batch $i/$NUM_CHUNKS: $LINES items → $OUT_DIR"

  "$PF2_BIN" predict \
    --multipleFiles "$chunk" \
    -f proteome \
    -o "$OUT_DIR"

  i=$((i+1))
done

echo "[done] Finished batches $start_i..$((i-1)) under $OUT_BASE (prefix '${OUT_PREFIX}')."
