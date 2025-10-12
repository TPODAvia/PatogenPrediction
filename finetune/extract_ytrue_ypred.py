#!/usr/bin/env python3
import sys, os, csv, numpy as np

# Usage:
#   python extract_ytrue_ypred.py /path/to/val_predictions.tsv [out.tsv]
# Default out: val_probs.tsv (includes y_true/y_pred + per-class probs)

in_path = sys.argv[1]
out_path = (sys.argv[2] if len(sys.argv) > 2
            else os.path.join(os.path.dirname(in_path), "val_probs.tsv"))

def is_float_token(tok: str) -> bool:
    try:
        float(tok); return True
    except Exception:
        return False

def try_extract(tokens):
    """
    tokens: list of tab-split fields for a (possibly wrapped) record.
    A complete record ends with 2*C floats: y_true_c0..c{C-1} then y_pred_c0..c{C-1}.
    Returns (y_true_idx, y_pred_idx, y_true_vec, y_pred_vec) or None if not ready.
    """
    k = 0
    for t in reversed(tokens):
        if is_float_token(t): k += 1
        else: break
    if k < 4 or (k % 2) != 0:
        return None
    C = k // 2
    y_true_vals = list(map(float, tokens[-2*C:-C]))
    y_pred_vals = list(map(float, tokens[-C:]))
    y_true_idx = int(np.argmax(y_true_vals))
    y_pred_idx = int(np.argmax(y_pred_vals))
    return y_true_idx, y_pred_idx, y_true_vals, y_pred_vals

with open(in_path, "r", encoding="utf-8", errors="ignore") as fin, \
     open(out_path, "w", newline="", encoding="utf-8") as fout:

    # writer
    w = csv.writer(fout, delimiter="\t")

    # skip header if present
    first_line = fin.readline()
    buffer = ""
    wrote_header = False

    for line in fin:
        if not line:
            continue
        buffer += line.rstrip("\n")
        tokens = buffer.split("\t")

        res = try_extract(tokens)
        if res is None:
            buffer += "\n"
            continue

        y_true_idx, y_pred_idx, y_true_vec, y_pred_vec = res
        C = len(y_pred_vec)

        if not wrote_header:
            header = (["y_true", "y_pred"]
                      + [f"y_true_c{i}" for i in range(C)]
                      + [f"y_pred_c{i}" for i in range(C)])
            w.writerow(header)
            wrote_header = True

        row = ([y_true_idx, y_pred_idx]
               + [f"{v:.6f}" for v in y_true_vec]
               + [f"{v:.6f}" for v in y_pred_vec])
        w.writerow(row)

        buffer = ""

print(f"Saved: {out_path}")
