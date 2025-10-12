#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import (
    precision_recall_fscore_support,
    accuracy_score,
    matthews_corrcoef,
    confusion_matrix,
    roc_curve,
    auc,
)

def _lower_cols(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={c: c.strip().lower() for c in df.columns})

def _sorted_cols(df, prefix):
    cols = [c for c in df.columns if c.startswith(prefix)]
    def keyfn(s):
        try:
            return int(s[len(prefix):])
        except Exception:
            return 10**9
    return sorted(cols, key=keyfn)

def load_any(path: str):
    """
    Returns:
      y_true_idx: (N,) int array of class indices 0..C-1
      y_pred_idx: (N,) int array of predicted class indices
      y_proba:    (N,C) float array of per-class probabilities or None if absent
      labels_sorted: list [0..C-1]
    Accepts:
      - minimal 2-col TSV: y_true, y_pred (indices), or
      - wide TSV: ... y_true_c0..y_true_c{C-1} y_pred_c0..y_pred_c{C-1}
    """
    df = pd.read_csv(path, sep="\t", engine="python")
    df = _lower_cols(df)

    pred_cols = _sorted_cols(df, "y_pred_c")
    true_cols = _sorted_cols(df, "y_true_c")

    if pred_cols:
        y_proba = df[pred_cols].to_numpy(dtype=float)
        C = y_proba.shape[1]
        if true_cols and len(true_cols) == C:
            y_true_idx = df[true_cols].to_numpy(dtype=float).argmax(axis=1)
        elif "y_true" in df.columns:
            y_true_idx = pd.to_numeric(df["y_true"], errors="coerce").astype("Int64").dropna().astype(int).to_numpy()
            if len(y_true_idx) != len(y_proba):
                raise ValueError("y_true length != probabilities length")
        else:
            raise ValueError("Could not find y_true from one-hot or index column.")
        y_pred_idx = y_proba.argmax(axis=1)
        labels_sorted = list(range(C))
        return y_true_idx, y_pred_idx, y_proba, labels_sorted

    # fallback: minimal two columns
    cand_true = [c for c in df.columns if c in {"y_true","true","label","target"} or "y_true" in c]
    cand_pred = [c for c in df.columns if c in {"y_pred","pred"} or "y_pred" in c]
    if not cand_true or not cand_pred:
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
    y_true_idx = sub["y_true"].astype(int).to_numpy()
    y_pred_idx = sub["y_pred"].astype(int).to_numpy()
    labels_sorted = sorted(set(y_true_idx) | set(y_pred_idx))
    return y_true_idx, y_pred_idx, None, labels_sorted

def plot_confusion(cm: np.ndarray, labels, out_png: Path, title: str, normalize: bool = False):
    plt.figure(figsize=(6,5), dpi=150)
    im = plt.imshow(cm, interpolation="nearest", aspect="auto")
    plt.title(title)
    plt.colorbar(im, fraction=0.046, pad=0.04)
    tick = np.arange(len(labels))
    plt.xticks(tick, labels, rotation=45, ha="right")
    plt.yticks(tick, labels)
    fmt = ".2f" if normalize else "d"
    thresh = cm.max()/2.0 if cm.size else 0.0
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            val = cm[i, j]
            plt.text(j, i, format(val, fmt),
                     ha="center", va="center",
                     color="white" if val > thresh else "black")
    plt.ylabel("True label")
    plt.xlabel("Predicted label")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def plot_multiclass_roc_perclass_only(y_true_idx, y_proba, class_names, out_png: Path):
    """
    One-vs-rest ROC: plot ONLY per-class curves (no micro/macro averages).
    Returns dict of per-class AUCs {class_idx: auc}.
    """
    from sklearn.preprocessing import label_binarize
    n_classes = len(class_names)
    y_true_bin = label_binarize(y_true_idx, classes=list(range(n_classes)))  # (N,C)

    aucs = {}
    plt.figure(figsize=(7,6), dpi=150)
    for c in range(n_classes):
        fpr, tpr, _ = roc_curve(y_true_bin[:, c], y_proba[:, c])
        auc_c = auc(fpr, tpr)
        aucs[c] = auc_c
        plt.plot(fpr, tpr, linewidth=1.8, label=f"{class_names[c]} (AUC = {auc_c:.3f})")

    plt.plot([0,1],[0,1], linestyle=":")
    plt.xlim([0.0, 1.0]); plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate"); plt.ylabel("True Positive Rate")
    plt.title("Multiclass ROC (OvR, per-class)")
    plt.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()
    return aucs

def plot_per_class_rocs(y_true_idx, y_proba, class_names, out_dir: Path):
    """Optional: one PNG per class."""
    from sklearn.preprocessing import label_binarize
    n_classes = len(class_names)
    y_true_bin = label_binarize(y_true_idx, classes=list(range(n_classes)))
    aucs = {}
    for c in range(n_classes):
        fpr, tpr, _ = roc_curve(y_true_bin[:, c], y_proba[:, c])
        aucs[c] = auc(fpr, tpr)
        plt.figure(figsize=(6,5), dpi=150)
        plt.plot(fpr, tpr, linewidth=2, label=f"AUC = {aucs[c]:.3f}")
        plt.plot([0,1],[0,1], linestyle="--")
        plt.xlim([0.0, 1.0]); plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate"); plt.ylabel("True Positive Rate")
        plt.title(f"ROC: {class_names[c]} (OvR)")
        plt.legend(loc="lower right")
        plt.tight_layout()
        plt.savefig(out_dir / f"roc_auc_{c}_{class_names[c].replace(' ','_')}.png")
        plt.close()
    return aucs

def main():
    ap = argparse.ArgumentParser(description="Multiclass metrics + per-class ROC/AUC (no micro/macro)")
    ap.add_argument("tsv", help="Path to val_probs.tsv (or wide val_predictions.tsv)")
    ap.add_argument("--outdir", default=None, help="Output directory (default: alongside input)")
    ap.add_argument("--class-names", default=None,
                    help="Comma-separated class names (e.g. 'No Pathogenic,Pathogenic,Unknown').")
    ap.add_argument("--save-individual-rocs", action="store_true",
                    help="Also save one ROC PNG per class.")
    args = ap.parse_args()

    outdir = Path(args.outdir) if args.outdir else Path(args.tsv).parent
    outdir.mkdir(parents=True, exist_ok=True)

    y_true, y_pred, y_proba, labels_sorted = load_any(args.tsv)

    # Class names
    if args.class_names:
        names = [s.strip() for s in args.class_names.split(",")]
        if len(names) != len(labels_sorted):
            raise ValueError(f"--class-names length ({len(names)}) != #labels ({len(labels_sorted)}).")
        class_names = names
    else:
        class_names = ["No Pathogenic", "Pathogenic", "Unknown"] if len(labels_sorted) == 3 else [str(i) for i in labels_sorted]

    # Confusion matrices
    cm = confusion_matrix(y_true, y_pred, labels=labels_sorted)
    with np.errstate(all="ignore"):
        row_sums = cm.sum(axis=1, keepdims=True)
        cm_norm = np.divide(cm, row_sums, where=row_sums!=0)

    # Overall metrics
    acc = accuracy_score(y_true, y_pred)
    mcc = matthews_corrcoef(y_true, y_pred) if len(labels_sorted) > 1 else np.nan

    prec_micro, rec_micro, f1_micro, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=labels_sorted, average="micro", zero_division=0
    )
    prec_macro, rec_macro, f1_macro, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=labels_sorted, average="macro", zero_division=0
    )
    prec_weighted, rec_weighted, f1_weighted, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=labels_sorted, average="weighted", zero_division=0
    )
    prec_c, rec_c, f1_c, sup_c = precision_recall_fscore_support(
        y_true, y_pred, labels=labels_sorted, average=None, zero_division=0
    )

    # Save per-class metrics
    per_class = pd.DataFrame({
        "class_id": labels_sorted,
        "class_name": class_names,
        "support": sup_c.astype(int),
        "precision": prec_c,
        "recall": rec_c,
        "f1": f1_c,
    }).sort_values("class_id")
    per_class.to_csv(outdir / "per_class_metrics.tsv", sep="\t", index=False)

    # Summary TSV (keeps micro/macro/weighted F1 values; remove here too if you don't want them)
    summary = pd.DataFrame([{
        "accuracy": acc,
        "mcc": mcc,
        "f1_micro": f1_micro,
        "precision_micro": prec_micro,
        "recall_micro": rec_micro,
        "f1_macro": f1_macro,
        "precision_macro": prec_macro,
        "recall_macro": rec_macro,
        "f1_weighted": f1_weighted,
        "precision_weighted": prec_weighted,
        "recall_weighted": rec_weighted,
        "n_samples": int(len(y_true)),
        "n_classes": int(len(labels_sorted)),
    }])
    summary.to_csv(outdir / "metrics_summary.tsv", sep="\t", index=False)

    # Confusion matrices to disk
    cm_df = pd.DataFrame(cm, index=class_names, columns=class_names)
    cm_df.to_csv(outdir / "confusion_matrix.tsv", sep="\t")
    cmn_df = pd.DataFrame(cm_norm, index=class_names, columns=class_names)
    cmn_df.to_csv(outdir / "confusion_matrix_normalized.tsv", sep="\t", float_format="%.6f")
    plot_confusion(cm, class_names, outdir / "confusion_matrix.png", "Confusion matrix (counts)", normalize=False)
    plot_confusion(cm_norm, class_names, outdir / "confusion_matrix_normalized.png", "Confusion matrix (row-normalized)", normalize=True)

    # ROC/AUC (per-class only)
    if y_proba is not None:
        aucs = plot_multiclass_roc_perclass_only(
            y_true_idx=y_true, y_proba=y_proba, class_names=class_names,
            out_png=outdir / "roc_auc_multiclass.png"
        )
        if args.save_individual_rocs:
            plot_per_class_rocs(y_true, y_proba, class_names, outdir)

        # Save ONLY per-class AUCs (no micro/macro rows)
        auc_rows = [{"which": f"class_{i}", "name": class_names[i], "auc": float(aucs[i])}
                    for i in range(len(class_names))]
        pd.DataFrame(auc_rows).to_csv(outdir / "roc_auc_values.tsv", sep="\t", index=False)
        print(f"\nROC/AUC saved → {outdir/'roc_auc_multiclass.png'}  (per-class only)")
        if args.save_individual_rocs:
            print("Per-class ROC PNGs saved (one file per class).")
    else:
        print("\nNo per-class probabilities found; skipping ROC/AUC.")
        print("Tip: run on the wide TSV with columns like y_pred_c0..y_pred_cK (your val_probs.tsv).")

    # Console summary
    print(f"\nRows: {len(y_true)}  |  Classes: {len(labels_sorted)} ({', '.join(class_names)})")
    print(f"Accuracy: {acc:.6f}  |  MCC: {mcc:.6f}")
    print(f"F1 micro: {f1_micro:.6f}  |  macro: {f1_macro:.6f}  |  weighted: {f1_weighted:.6f}")
    print("Saved files in:", outdir)
    print(" - metrics_summary.tsv")
    print(" - per_class_metrics.tsv")
    print(" - confusion_matrix.tsv / confusion_matrix.png")
    print(" - confusion_matrix_normalized.tsv / confusion_matrix_normalized.png")
    if y_proba is not None:
        print(" - roc_auc_values.tsv  (per-class only)")
        print(" - roc_auc_multiclass.png")
        if args.save_individual_rocs:
            print(" - roc_auc_<id>_<name>.png (per-class)")
if __name__ == "__main__":
    main()
