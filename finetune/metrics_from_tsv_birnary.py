#!/usr/bin/env python3
import argparse, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve, auc, f1_score, precision_score, recall_score,
    accuracy_score, matthews_corrcoef, confusion_matrix
)

def load_df(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", engine="python")
    # Normalize headers
    lower = {c: c.strip().lower() for c in df.columns}
    df = df.rename(columns=lower)
    # Try to find y_true / y_pred
    cand_true = [c for c in df.columns if "y_true" in c or c in {"true","label","target"}]
    cand_pred = [c for c in df.columns if "y_pred" in c or c in {"pred","score","prob"} or "score" in c or "prob" in c]
    if not cand_true or not cand_pred:
        # fallback: last two columns are y_true, y_pred
        if df.shape[1] < 2:
            raise ValueError("Could not locate y_true/y_pred columns.")
        sub = df.iloc[:, -2:].copy()
        sub.columns = ["y_true","y_pred"]
    else:
        sub = df[[cand_true[0], cand_pred[0]]].copy()
        sub.columns = ["y_true","y_pred"]

    sub["y_true"] = pd.to_numeric(sub["y_true"], errors="coerce")
    sub["y_pred"] = pd.to_numeric(sub["y_pred"], errors="coerce")
    sub = sub.dropna()
    # binarize true labels and clip scores
    sub["y_true"] = (sub["y_true"] > 0.5).astype(int)
    sub["y_pred"] = sub["y_pred"].clip(0.0, 1.0)
    return sub

def metrics_at(y_true, y_score, thr: float):
    y_hat = (y_score >= thr).astype(int)
    f1 = f1_score(y_true, y_hat, zero_division=0)
    prec = precision_score(y_true, y_hat, zero_division=0)
    rec = recall_score(y_true, y_hat, zero_division=0)
    acc = accuracy_score(y_true, y_hat)
    mcc = matthews_corrcoef(y_true, y_hat) if len(np.unique(y_true)) == 2 else np.nan
    tn, fp, fn, tp = confusion_matrix(y_true, y_hat, labels=[0,1]).ravel()
    return dict(threshold=thr, f1=f1, precision=prec, recall=rec,
                accuracy=acc, mcc=mcc, tp=int(tp), fp=int(fp), tn=int(tn), fn=int(fn))

def sweep_thresholds(y_true, y_score):
    uniq = np.unique(y_score)
    grid = np.linspace(0, 1, 501)
    thr_all = np.unique(np.concatenate([uniq, grid, [0.5]]))
    rows = [metrics_at(y_true, y_score, float(t)) for t in thr_all]
    return pd.DataFrame(rows).sort_values("threshold").reset_index(drop=True)

def plot_roc(y_true, y_score, out_png: Path) -> float:
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    plt.figure(figsize=(6,5), dpi=150)
    plt.plot(fpr, tpr, linewidth=2)
    plt.plot([0,1],[0,1], linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve (AUC = {roc_auc:.4f})")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()
    return float(roc_auc)

def main():
    ap = argparse.ArgumentParser(description="Compute F1, Precision, MCC, and ROC/AUC from val_ytrue_ypred.tsv")
    ap.add_argument("tsv", help="Path to val_ytrue_ypred.tsv")
    ap.add_argument("--outdir", default=None, help="Output directory (default: alongside input)")
    args = ap.parse_args()

    outdir = Path(args.outdir) if args.outdir else Path(args.tsv).parent
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_df(args.tsv)
    y_true = df["y_true"].to_numpy(dtype=int)
    y_score = df["y_pred"].to_numpy(dtype=float)

    # ROC & AUC
    auc_val = plot_roc(y_true, y_score, outdir / "roc_auc.png")

    # Threshold sweep
    sweep = sweep_thresholds(y_true, y_score)
    sweep.to_csv(outdir / "threshold_metrics.tsv", sep="\t", index=False)

    # Metrics at 0.5 and best F1
    m05 = metrics_at(y_true, y_score, 0.5)
    best_idx = int(sweep["f1"].idxmax())
    mbest = sweep.iloc[best_idx].to_dict()

    summary = pd.DataFrame([
        {"which":"auc_only", "AUC": auc_val},
        {"which":"thr=0.5", **{k:v for k,v in m05.items() if k!="threshold"}},
        {"which":"best_f1", **{k:v for k,v in mbest.items() if k!="threshold"}},
    ])
    summary.to_csv(outdir / "metrics_summary.tsv", sep="\t", index=False)

    # Console summary
    print(f"\nRows: {len(df)}")
    print(f"AUC: {auc_val:.6f}")
    print("\n@ threshold = 0.5")
    print(f"  F1: {m05['f1']:.6f} | Precision: {m05['precision']:.6f} | Recall: {m05['recall']:.6f} | "
          f"Accuracy: {m05['accuracy']:.6f} | MCC: {m05['mcc']:.6f}")
    print(f"  Confusion: TP={m05['tp']} FP={m05['fp']} TN={m05['tn']} FN={m05['fn']}")
    print("\n@ best F1")
    print(f"  Thr*: {mbest['threshold']:.4f} | F1: {mbest['f1']:.6f} | Precision: {mbest['precision']:.6f} | "
          f"Recall: {mbest['recall']:.6f} | Accuracy: {mbest['accuracy']:.6f} | MCC: {mbest['mcc']:.6f}")
    print(f"  Confusion: TP={int(mbest['tp'])} FP={int(mbest['fp'])} TN={int(mbest['tn'])} FN={int(mbest['fn'])}")
    print(f"\nSaved files in: {outdir}")
    print(" - roc_auc.png")
    print(" - metrics_summary.tsv")
    print(" - threshold_metrics.tsv")

if __name__ == "__main__":
    main()
