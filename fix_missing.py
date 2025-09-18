#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from pathlib import Path
import pandas as pd

UNMATCHED_CSV = Path("/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/unmatched_by_seed.csv")
INPUT_LIST    = Path("/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/pf2_list_fna.txt")
OUTPUT_LIST   = Path("/home/rover2/HW1_popgen/kursov/PatogenPrediction/dataset/pf2_list_fna_fix.txt")

# numeric tokens like 1404, 1404.5, 33038.800, 6666666.205260, etc.
NUMERIC_SEED_RE = re.compile(r"\d+(?:\.\d+)?")

def canon_seed(s: str) -> str:
    """Canonicalize numeric SEED_ID: strip trailing zeros after decimal and any trailing dot."""
    s = str(s).strip().strip('"').strip("'")
    if '.' in s:
        s = s.rstrip('0').rstrip('.')
    return s

def load_unmatched_seeds_pandas(csv_path: Path) -> set[str]:
    """
    Robustly load SEED_ID values from CSV:
    - auto-detect delimiter (sep=None, engine='python')
    - normalize column names by stripping spaces/quotes
    - canonicalize numeric forms (1404.50 => 1404.5)
    """
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing file: {csv_path}")

    df = pd.read_csv(csv_path, sep=None, engine="python", dtype=str)
    # normalize column names
    df.columns = [c.strip().strip('"').strip("'") for c in df.columns]

    # find a column that equals SEED_ID (case-insensitive, trimmed)
    col_map = {c.lower(): c for c in df.columns}
    if "seed_id" not in col_map:
        raise RuntimeError(f"SEED_ID column not found in {csv_path} (have: {', '.join(df.columns)})")
    col = col_map["seed_id"]

    seeds = set()
    for v in df[col].dropna().astype(str):
        v = v.strip()
        # some rows may include extra text: extract first numeric token
        m = NUMERIC_SEED_RE.search(v)
        if m:
            seeds.add(canon_seed(m.group(0)))
        else:
            # if the cell itself is a clean numeric string
            if NUMERIC_SEED_RE.fullmatch(v):
                seeds.add(canon_seed(v))
    return seeds

def line_has_target_seed(line: str, targets: set[str]) -> bool:
    """Return True if any numeric token in the line canonicalizes to a target SEED_ID."""
    for m in NUMERIC_SEED_RE.finditer(line):
        if canon_seed(m.group(0)) in targets:
            return True
    return False

def main():
    targets = load_unmatched_seeds_pandas(UNMATCHED_CSV)
    if not INPUT_LIST.exists():
        raise FileNotFoundError(f"Missing file: {INPUT_LIST}")

    seen = set()
    kept = []

    with INPUT_LIST.open("r", encoding="utf-8", errors="replace") as fin:
        for line in fin:
            if line_has_target_seed(line, targets):
                if line not in seen:
                    kept.append(line)
                    seen.add(line)

    with OUTPUT_LIST.open("w", encoding="utf-8") as fout:
        for ln in kept:
            fout.write(ln)

    print(f"[OK] Loaded SEED_IDs: {len(targets)}")
    # show a few examples
    ex = list(sorted(targets))[:10]
    print(f"[EXAMPLES] {ex}")
    print(f"[OK] Lines copied: {len(kept)}")
    print(f"[WROTE] {OUTPUT_LIST}")

if __name__ == "__main__":
    main()
