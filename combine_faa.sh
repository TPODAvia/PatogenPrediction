#!/usr/bin/env bash
set -euo pipefail

# write absolute paths into this list
LIST="/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/pf2_list_faa.txt"
mkdir -p /home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset

# find *.faa under .../prokka_outputs/<ID>/<ID>.faa
find "/home/rover2/2856_genomes/prokka_outputs" \
  -mindepth 2 -maxdepth 2 -type f -iname "*.faa" -size +0 \
  -print0 | sort -z | tr '\0' '\n' > "$LIST"

# quick sanity check
wc -l "$LIST"
head -n 4 "$LIST"
