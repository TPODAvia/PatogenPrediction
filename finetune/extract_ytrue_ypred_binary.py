#!/usr/bin/env python3
import sys, os, re, csv

# Usage:
#   python extract_ytrue_ypred.py /path/to/val_predictions.tsv [out.tsv]
#
# Example:
#   python extract_ytrue_ypred.py \
#     /home/rover2/HW1_popgen/kursov/PatogenPrediction/splits_threeclass_result/trainPF2/val_predictions.tsv

in_path = sys.argv[1]
out_path = (sys.argv[2] if len(sys.argv) > 2
            else os.path.join(os.path.dirname(in_path), "val_ytrue_ypred.tsv"))

# regex to capture the final two floats on a line (supports 1.000000, 0.623047, 1e-3, etc.)
tail_two_floats = re.compile(
    r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*$'
)

with open(in_path, "r", encoding="utf-8", errors="ignore") as fin, \
     open(out_path, "w", newline="", encoding="utf-8") as fout:
    w = csv.writer(fout, delimiter="\t")
    w.writerow(["y_true", "y_pred"])  # only two columns as requested

    # skip header if present
    first = fin.readline()
    # process the rest
    for line in fin:
        if not line.strip():
            continue
        m = tail_two_floats.search(line)
        if not m:
            # could be a wrapped/continued line; skip quietly
            continue
        y_true, y_pred = m.group(1), m.group(2)
        w.writerow([y_true, y_pred])

print(f"Saved: {out_path}")
