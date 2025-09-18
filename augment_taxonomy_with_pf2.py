#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
augment_taxonomy_with_pf2_by_seed.py

Strict SEED_ID-based augmentation:
- Finds: /<base_root>/result<i>/predictionPF2/results_*.faa/predictions.tsv
- Extracts SEED_ID from "results_*.faa"
- Reads PF2 columns from each predictions.tsv
- Left-joins onto taxonomy by SEED_ID only
- Guarantees same number of rows as input taxonomy

Usage:
python3 augment_taxonomy_with_pf2.py \
  --base_root "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset" \
  --taxonomy_csv "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/Ref_taxonomies_2856.csv" \
  --out_csv "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/Ref_taxonomies_2856_mod.csv" \
  --unmatched_out "/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/unmatched_by_seed.csv"
"""

import os
import re
import glob
import sys
import argparse
import warnings
from typing import Dict, List, Optional

import pandas as pd

# PF2 columns we will import (missing ones will be created as empty)
PF2_COLS = [
    "Prediction Mean",
    "Prediction STD",
    "Phenotype",
    "Binary Prediction Mean",
]

# ------------------------- helpers -------------------------

RESULTS_DIR_NAME_RE = re.compile(r"results_(.+?)\.faa$", re.IGNORECASE)
RESULTS_DIR_FALLBACK_RE = re.compile(r"results_(.+)$", re.IGNORECASE)

def extract_seed_id_from_results_dir(results_dir_name: str) -> Optional[str]:
    """
    From a directory name like 'results_195099.5.faa' -> '195099.5'
    Fallback handles 'results_195099.5' if '.faa' missing.
    """
    m = RESULTS_DIR_NAME_RE.search(results_dir_name)
    if m:
        return m.group(1)
    m = RESULTS_DIR_FALLBACK_RE.search(results_dir_name)
    if m:
        return m.group(1)
    return None


NUMERIC_SEED_RE = re.compile(r'^\s*\d+(?:\.\d+)?\s*$')

def canon_seed(x) -> str:
    """
    Canonicalize SEED_ID so '1404.50' == '1404.5' and '33038.800' == '33038.8'.
    - strip spaces
    - if numeric-looking, strip trailing zeros after decimal and any trailing dot
    """
    if x is None:
        return ""
    s = str(x).strip()
    if NUMERIC_SEED_RE.match(s):
        if '.' in s:
            s = s.rstrip('0').rstrip('.')
    return s

def norm_seed(x) -> str:
    """
    Normalize SEED_ID to a comparable string:
    - strip spaces
    - keep as original text (no float formatting), but ensure it's str
    """
    if x is None:
        return ""
    s = str(x).strip()
    return s

def read_predictions_tsv(tsv_path: str) -> pd.DataFrame:
    """
    Read PF2 predictions.tsv (commented with '#').
    Expect one row per input.
    """
    return pd.read_csv(tsv_path, sep="\t", comment="#", engine="python")

def collect_pf2_by_seed(base_root: str) -> pd.DataFrame:
    """
    Scan all result{i}/predictionPF2 trees, collect predictions into a
    dataframe with columns: SEED_ID + PF2_COLS
    """
    pattern = os.path.join(base_root, "result*", "predictionPF2", "results_*", "predictions.tsv")
    tsv_paths = sorted(glob.glob(pattern))
    if not tsv_paths:
        warnings.warn(f"No predictions.tsv found under {base_root} (pattern: {pattern})")

    records: List[Dict] = []

    for tsv in tsv_paths:
        # results dir name is parent of predictions.tsv
        results_dir = os.path.basename(os.path.dirname(tsv))
        seed_id = extract_seed_id_from_results_dir(results_dir)
        if not seed_id:
            warnings.warn(f"Could not extract SEED_ID from folder '{results_dir}', skipping.")
            continue

        try:
            df = read_predictions_tsv(tsv)
        except Exception as e:
            warnings.warn(f"Failed reading {tsv}: {e}")
            continue

        if df.empty:
            warnings.warn(f"Empty predictions.tsv: {tsv}")
            continue

        # PF2 writes one row per input; take the first
        row = df.iloc[0].to_dict()

        rec = {"SEED_ID": canon_seed(seed_id)}
        # pick required PF2 columns (fill missing as None)
        for c in PF2_COLS:
            rec[c] = row.get(c, None)
        records.append(rec)

    if not records:
        return pd.DataFrame(columns=["SEED_ID"] + PF2_COLS)

    pf2 = pd.DataFrame.from_records(records)

    # Normalize SEED_ID and drop duplicates (keep first seen)
    pf2["SEED_ID"] = pf2["SEED_ID"].map(canon_seed)
    pf2 = pf2.drop_duplicates(subset=["SEED_ID"], keep="first").reset_index(drop=True)

    return pf2

# ------------------------- main -------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Augment taxonomy by SEED_ID using PF2 predictions from all result{i}/predictionPF2 trees."
    )
    ap.add_argument("--base_root", required=True,
                    help="Parent folder containing result{i}/predictionPF2 trees.")
    ap.add_argument("--taxonomy_csv", required=True,
                    help="Original taxonomy CSV (must include 'SEED_ID').")
    ap.add_argument("--out_csv", required=True,
                    help="Output CSV path.")
    ap.add_argument("--unmatched_out", required=False,
                    help="Optional CSV path to write ALL unmatched rows.")
    args = ap.parse_args()

    # ---- Load taxonomy and normalize SEED_ID ----
    try:
        tax = pd.read_csv(args.taxonomy_csv, sep=None, engine="python")
    except Exception as e:
        sys.exit(f"Failed to read taxonomy CSV: {args.taxonomy_csv}\n{e}")

    if "SEED_ID" not in tax.columns:
        sys.exit("Input taxonomy is missing required column 'SEED_ID'.")

    in_rows = len(tax)
    tax = pd.read_csv(args.taxonomy_csv, sep=None, engine="python", dtype={"SEED_ID": str})
    tax["SEED_ID"] = tax["SEED_ID"].map(canon_seed)

    # ---- Collect PF2 predictions grouped by SEED_ID ----
    pf2 = collect_pf2_by_seed(args.base_root)

    # Ensure PF2 columns exist (even if empty)
    for c in PF2_COLS:
        if c not in pf2.columns:
            pf2[c] = None

    # Deduplicate by SEED_ID (keep first)
    if not pf2.empty:
        pf2 = pf2.drop_duplicates(subset=["SEED_ID"], keep="first").reset_index(drop=True)

    # ---- Left-join strictly by SEED_ID ----
    merged = tax.merge(
        pf2[["SEED_ID"] + PF2_COLS] if not pf2.empty else pd.DataFrame(columns=["SEED_ID"] + PF2_COLS),
        on="SEED_ID",
        how="left",
        copy=False,
        validate="many_to_one"  # taxonomy may have unique SEED_ID; PF2 should be <=1 per SEED_ID
    )

    out_rows = len(merged)
    if out_rows != in_rows:
        raise RuntimeError(f"Row count changed after merge: in={in_rows}, out={out_rows}")

    # ---- Report matches & list ALL unmatched rows ----
    matched_mask = merged[PF2_COLS].notna().any(axis=1)
    unmatched_mask = ~matched_mask

    matched_count = int(matched_mask.sum())
    unmatched_count = int(unmatched_mask.sum())
    pct = (matched_count / in_rows * 100.0) if in_rows else 0.0

    print(f"[OK] Rows preserved: {in_rows} → {out_rows}")
    print(f"[INFO] Matched by SEED_ID: {matched_count} / {in_rows} ({pct:.1f}%)")
    print(f"[INFO] Unmatched rows (blank PF2 columns): {unmatched_count}")

    # Pretty-print all unmatched rows (compact view)
    unmatched = merged.loc[unmatched_mask].copy()
    print("\n[UNMATCHED ROWS BY SEED_ID]")
    if unmatched.empty:
        print("None")
    else:
        show_cols = [c for c in [
            "SEED_ID",
            "Organism (genomeID)",
            "Organism",
            "TaxID"
        ] if c in unmatched.columns]
        if not show_cols:
            show_cols = ["SEED_ID"]

        with pd.option_context(
            "display.max_rows", None,
            "display.max_colwidth", 120,
            "display.width", 200
        ):
            print(unmatched[show_cols].to_string(index=False))

        # Optionally save ALL columns of unmatched rows
        if args.unmatched_out:
            try:
                unmatched.to_csv(args.unmatched_out, index=False)
                print(f"\n[Saved] Unmatched rows (all columns) → {args.unmatched_out}")
            except Exception as e:
                print(f"[WARN] Failed to write unmatched_out CSV: {e}", file=sys.stderr)

    # ---- Write final augmented CSV ----
    try:
        merged.to_csv(args.out_csv, index=False)
        print(f"\n[Wrote] {args.out_csv}")
    except Exception as e:
        sys.exit(f"Failed to write output CSV: {args.out_csv}\n{e}")

if __name__ == "__main__":
    main()
