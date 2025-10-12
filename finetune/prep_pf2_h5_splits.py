#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Builds train/val pairs for PF2 from nested ProteinEmbeddings.h5 layout like:
  <H5_ROOT>/result*/predictionPF2/preprocess_<genome_id>.faa/ProteinEmbeddings.h5

Outputs:
  train_pairs.tsv, val_pairs.tsv, label_map.json, all_h5_list.txt,
  missing_genomes.tsv, found_h5.csv, duplicates_h5.csv

Label mapping:
  0.0 -> 0, 0.5 -> 1, 1.0 -> 2
"""

import argparse, os, re, json, time, stat
from collections import defaultdict
import pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit

PREPROCESS_RE = re.compile(r"^preprocess_(?P<gid>.+)\.faa$")

def _result_num(path: str) -> int:
    # try to extract result number from ".../result123/..."
    parts = path.split(os.sep)
    for p in parts:
        if p.startswith("result"):
            num = p[len("result"):]
            return int(num) if num.isdigit() else -1
    return -1

def build_h5_index(h5_root: str):
    """
    Returns:
      best: dict[genome_id] = path
      all_rows: list of dict rows for CSV (genome_id, resultN, mtime, path)
      dups: list of genome_ids that had >1 candidate
    """
    rows = []
    for root, dirs, files in os.walk(h5_root):
        if "ProteinEmbeddings.h5" in files:
            # Expect parent folder name like preprocess_<gid>.faa
            parent = os.path.basename(root)
            m = PREPROCESS_RE.match(parent)
            if not m:
                continue
            gid = m.group("gid")
            path = os.path.join(root, "ProteinEmbeddings.h5")
            try:
                st = os.stat(path)
                mtime = int(st.st_mtime)
            except OSError:
                mtime = 0
            rows.append({
                "genome_id": gid,
                "resultN": _result_num(path),
                "mtime": mtime,
                "path": path
            })

    if not rows:
        raise SystemExit(f"[err] No ProteinEmbeddings.h5 found under {h5_root}")

    df = pd.DataFrame(rows)
    # choose best per genome: max resultN, then max mtime
    df["rank_key"] = list(zip(df["resultN"], df["mtime"]))
    df_sorted = df.sort_values(["genome_id", "resultN", "mtime"], ascending=[True, False, False])

    best = {}
    dups = []
    grouped = df_sorted.groupby("genome_id")
    for gid, g in grouped:
        if len(g) > 1:
            dups.append(gid)
        best_path = g.iloc[0]["path"]
        best[gid] = best_path

    return best, df_sorted.drop(columns=["rank_key"]), dups

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True,
                    help="CSV with at least columns: genome_id, Phenotype")
    ap.add_argument("--h5-root", required=True,
                    help="Root folder that contains result*/predictionPF2/.../ProteinEmbeddings.h5")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--test-size", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[info] Scanning H5 root: {args.h5_root}")
    h5_index, found_df, dups = build_h5_index(args.h5_root)
    found_csv = os.path.join(args.outdir, "found_h5.csv")
    found_df.to_csv(found_csv, index=False)
    if dups:
        dup_csv = os.path.join(args.outdir, "duplicates_h5.csv")
        (found_df[found_df["genome_id"].isin(dups)]
         .sort_values(["genome_id", "resultN", "mtime"], ascending=[True, False, False])
        ).to_csv(dup_csv, index=False)
        print(f"[warn] Duplicates for {len(dups)} genome_ids → {dup_csv}")

    print(f"[ok] Indexed {len(h5_index)} distinct genome_ids with H5.")

    df = pd.read_csv(args.csv)
    if not {"genome_id", "Phenotype"}.issubset(df.columns):
        raise SystemExit("[err] Input CSV must have columns: genome_id, Phenotype")

    df = df[["genome_id", "Phenotype"]].copy()
    # normalize types/strings
    df["genome_id"] = df["genome_id"].astype(str).str.strip()
    # map phenotype to class id
    def to_float(x):
        try:
            return float(str(x).strip())
        except Exception:
            return None
    df["Phenotype_f"] = df["Phenotype"].apply(to_float)

    # lab_map = {0.0: 0, 0.5: 1, 1.0: 2}
    lab_map = {0.0: "No Pathogenic", 0.5: "Unknown", 1.0: "Pathogenic"}
    df["PathoPhenotype"] = df["Phenotype_f"].map(lab_map)

    bad = df[df["PathoPhenotype"].isna()]
    if len(bad):
        bad_path = os.path.join(args.outdir, "bad_phenotype_rows.csv")
        bad.to_csv(bad_path, index=False)
        print(f"[warn] Dropping {len(bad)} rows with unsupported Phenotype (saved {bad_path})")
        df = df[~df.index.isin(bad.index)].copy()

    # attach h5 path
    df["File_Embedding"] = df["genome_id"].map(h5_index)
    missing = df[df["File_Embedding"].isna()]
    if len(missing):
        miss_path = os.path.join(args.outdir, "missing_genomes.tsv")
        missing[["genome_id", "Phenotype"]].to_csv(miss_path, sep="\t", index=False)
        print(f"[warn] {len(missing)} genome_ids missing H5 → {miss_path}")

    df = df.dropna(subset=["File_Embedding", "PathoPhenotype"]).copy()
    # df["PathoPhenotype"] = df["PathoPhenotype"].astype(int)

    # save an all-h5 list (useful for predict -f embeddings)
    all_list = os.path.join(args.outdir, "all_h5_list.txt")
    with open(all_list, "w") as f:
        for p in sorted(df["File_Embedding"].unique()):
            f.write(p + "\n")

    # stratified split
    X = df["File_Embedding"].values
    y = df["PathoPhenotype"].values
    sss = StratifiedShuffleSplit(n_splits=1, test_size=args.test_size, random_state=args.seed)
    tr_idx, va_idx = next(sss.split(X, y))
    tr = df.iloc[tr_idx][["File_Embedding", "PathoPhenotype"]]
    va = df.iloc[va_idx][["File_Embedding", "PathoPhenotype"]]

    tr_path = os.path.join(args.outdir, "train_pairs.tsv")
    va_path = os.path.join(args.outdir, "val_pairs.tsv")
    tr.to_csv(tr_path, sep="\t", index=False)
    va.to_csv(va_path, sep="\t", index=False)

    # label map json
    with open(os.path.join(args.outdir, "label_map.json"), "w") as f:
        # json.dump({"0.0": 0, "0.5": 1, "1.0": 2}, f, indent=2)
        json.dump({"0.0": "No Pathogenic",  0.5: "Unknown", "1.0": "Pathogenic"}, f, indent=2)

    # quick report
    report = {
        "rows_in_csv": int(len(df) + len(missing) + len(bad)),
        "usable_rows": int(len(df)),
        "missing_h5": int(len(missing)),
        "bad_phenotype_rows": int(len(bad)),
        "train_rows": int(len(tr)),
        "val_rows": int(len(va)),
        "class_counts_full": df["PathoPhenotype"].value_counts().to_dict(),
        "class_counts_train": tr["PathoPhenotype"].value_counts().to_dict(),
        "class_counts_val": va["PathoPhenotype"].value_counts().to_dict(),
    }
    rep_path = os.path.join(args.outdir, "prep_report.json")
    with open(rep_path, "w") as f:
        json.dump(report, f, indent=2)

    print("[done] Wrote:")
    print(f"  • {tr_path}")
    print(f"  • {va_path}")
    print(f"  • {all_list}")
    print(f"  • {found_csv}")
    print(f"  • {rep_path}")
    print("  • label_map.json")
    if len(missing): print(f"  • missing_genomes.tsv")

if __name__ == "__main__":
    main()