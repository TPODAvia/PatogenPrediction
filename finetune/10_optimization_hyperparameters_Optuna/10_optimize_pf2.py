#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json, os, subprocess, copy, shutil, csv
from pathlib import Path

import optuna
from optuna.samplers import TPESampler
from optuna.pruners import NopPruner
from optuna.exceptions import TrialPruned

# ---------- ПУТИ ----------
BASE_CONFIG = Path("/home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/train.json")
ROOT = Path("/home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/10_optimization_hyperparameters_Optuna"); ROOT.mkdir(parents=True, exist_ok=True)
SUFFIX = "_tuned"
STUDY_NAME = "pf2_convnext_addatt_mcc_auc" + SUFFIX
N_TRIALS = 20

# Пост-обработка из README
EXTRACT_SCRIPT = Path("/home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/extract_ytrue_ypred.py")
METRICS_SCRIPT = Path("/home/rover2/HW1_popgen/kursov/PatogenPrediction/finetune/metrics_from_tsv.py")

SUBDIR = "trainPF2"
VAL_PREDS = "val_predictions.tsv"
VAL_PROBS = "val_probs.tsv"

# ---------- УТИЛИТЫ ----------
def load_base():
    with open(BASE_CONFIG, "r", encoding="utf-8") as f:
        return json.load(f)

def trial_dirs(i: int):
    tdir = ROOT / f"trial_{i:03d}{SUFFIX}"
    if tdir.exists():
        shutil.rmtree(tdir, ignore_errors=True)
    tdir.mkdir(parents=True, exist_ok=True)
    rdir = tdir / f"results{SUFFIX}"   # ВАЖНО: PF2 сам создаст; мы не создаём заранее
    pdir = rdir / SUBDIR               # появится после train
    return tdir, rdir, pdir

def write_config(cfg: dict, tdir: Path):
    cfg_path = tdir / f"config{SUFFIX}.json"
    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=2)
    return cfg_path

def run(cmd: str, env=None):
    print(f"[RUN] {cmd}")
    proc = subprocess.run(cmd, shell=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {cmd}")

def train_with_cli(cfg_path: Path, outdir: Path):
    if outdir.exists():  # PF2 не любит существующие папки вывода
        shutil.rmtree(outdir, ignore_errors=True)
    run(f'pathogenfinder2 train -c "{cfg_path}" -o "{outdir}"')

def postprocess(pf2_dir: Path):
    preds_tsv = pf2_dir / VAL_PREDS
    if preds_tsv.exists():
        run(f'python3 "{EXTRACT_SCRIPT}" "{preds_tsv}"')
    probs_tsv = pf2_dir / VAL_PROBS
    if probs_tsv.exists():
        run(f'python3 "{METRICS_SCRIPT}" "{probs_tsv}" --outdir "{pf2_dir}"')

def try_read_metrics_json(pf2_dir: Path):
    for p in pf2_dir.glob("*.json"):
        try:
            with open(p, "r", encoding="utf-8") as f:
                data = json.load(f)
            keys = {k.lower(): k for k in data.keys()}
            mcc = data.get(keys.get("mcc")) if "mcc" in keys else None
            auc = data.get(keys.get("roc_auc")) if "roc_auc" in keys else None
            if mcc is not None and auc is not None:
                return float(mcc), float(auc)
        except Exception:
            pass
    return None

def try_read_metrics_summary(pf2_dir: Path):
    """Надёжно читаем MCC и ROC AUC из metrics_summary.tsv, чтоб не падать из-за отсутствия y_true."""
    summ = pf2_dir / "metrics_summary.tsv"
    if not summ.exists():
        return None
    with open(summ, "r", encoding="utf-8") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        return None
    header = [h.strip().lower() for h in rows[0]]
    # Вариант wide: берём последнюю строку
    if len(header) > 1 and len(rows) > 1:
        def idx(names):
            for n in names:
                if n in header:
                    return header.index(n)
            return None
        mcc_i = idx(["mcc","matthews","matthews_corrcoef"])
        auc_i = idx(["roc auc","roc_auc","auroc","roc_auc_macro","roc-auc"])
        if mcc_i is not None and auc_i is not None:
            last = rows[-1]
            try:
                return float(last[mcc_i]), float(last[auc_i])
            except:
                pass
    # Вариант key/value: две колонки
    mcc = auc = None
    for r in rows:
        if len(r) < 2:
            continue
        k, v = r[0].strip().lower(), r[1].strip()
        if "mcc" in k:
            try: mcc = float(v)
            except: pass
        if "roc" in k and "auc" in k:
            try: auc = float(v)
            except: pass
    if mcc is not None and auc is not None:
        return mcc, auc
    return None

def compute_metrics_from_files(pf2_dir: Path):
    from sklearn.metrics import matthews_corrcoef, roc_auc_score
    import numpy as np, csv as _csv
    y_true = None
    y_true_txt = pf2_dir / "y_true.txt"
    if y_true_txt.exists():
        y_true = [int(x.strip()) for x in y_true_txt.read_text().splitlines() if x.strip()!=""]
    else:
        preds_tsv = pf2_dir / VAL_PREDS
        if preds_tsv.exists():
            with open(preds_tsv, "r", encoding="utf-8") as f:
                rows = list(_csv.DictReader(f, delimiter="\t"))
            for c in ["y_true","true","label","target","class_true"]:
                if rows and c in rows[0]:
                    y_true = [int(r[c]) for r in rows]; break
    if y_true is None:
        raise FileNotFoundError("Не удалось получить y_true.")
    probs_tsv = pf2_dir / VAL_PROBS
    if not probs_tsv.exists():
        raise FileNotFoundError(f"Нет {VAL_PROBS}")
    with open(probs_tsv, "r", encoding="utf-8") as f:
        reader = _csv.reader(f, delimiter="\t")
        header = next(reader); rows = list(reader)
    header_l = [h.lower() for h in header]
    idxs = []
    for names in (["p0","p1","p2"],["prob_0","prob_1","prob_2"],["class0","class1","class2"],["prob0","prob1","prob2"]):
        try:
            idxs = [header_l.index(n) for n in names]; break
        except ValueError: idxs=[]
    if not idxs: idxs = list(range(len(header)-3, len(header)))
    import numpy as np
    probas = np.array([[float(r[i]) for i in idxs] for r in rows], dtype=float)
    y_true = np.array(y_true[:len(probas)], dtype=int)
    y_pred = probas.argmax(axis=1)
    mcc = matthews_corrcoef(y_true, y_pred)
    auc = roc_auc_score(y_true, probas, multi_class="ovr", average="macro")
    return float(mcc), float(auc)

def get_metrics(pf2_dir: Path):
    m = try_read_metrics_json(pf2_dir)
    if m is not None: return m
    m = try_read_metrics_summary(pf2_dir)   # важная вставка
    if m is not None: return m
    return compute_metrics_from_files(pf2_dir)

def copy_with_suffix_all(pf2_dir: Path):
    for p in pf2_dir.rglob("*"):
        if p.is_file():
            tuned = p.with_name(p.stem + SUFFIX + p.suffix)
            if tuned.name != p.name and not tuned.exists():
                try: shutil.copy2(p, tuned)
                except Exception: pass

# ---------- ПРОСТРАНСТВО ПОИСКА ----------
def force_valid_pf2(cfg: dict):
    # Строки, которые ожидает PF2
    cfg.setdefault("Model Parameters", {}).update({
        "Loss Function": "crossentropy",
        "Input Dimensions": 1024,
    })
    net = cfg["Model Parameters"].setdefault("Network Structure", {})
    net["Length Information"] = "concat1"     # фиксируем чтобы не ловить shape mismatch
    net["Length Dimensions"] = 30
    cfg.setdefault("Train Parameters", {}).update({
        "Loss Function": "CrossEntropyLoss",
    })
    return cfg

def suggest(trial, base_cfg):
    cfg = force_valid_pf2(copy.deepcopy(base_cfg))
    mp = cfg["Model Parameters"]; net = mp["Network Structure"]; tp = cfg["Train Parameters"]

    # Не трогаем Length Information / Length Dimensions
    net["Num Blocks"] = trial.suggest_int("num_blocks", 2, 4)
    net["Block Dimensions"] = trial.suggest_categorical("convnext_dim", [64, 128, 256, 512])
    net["Attention Dimensions"] = trial.suggest_categorical("attn_dim", [64, 128, 256, 512])

    mp["Sequence Dropout"] = trial.suggest_categorical("seq_dropout", [0.2, 0.3, 0.4])
    mp["Attention Dropout"] = trial.suggest_categorical("attn_dropout", [0.3, 0.4, 0.5, 0.6])
    mp["Stochastic Depth Prob"] = trial.suggest_categorical("stoch_depth", [0.1, 0.2, 0.3, 0.4])
    mp["Stochastic Depth Prob Att"] = trial.suggest_categorical("stoch_depth_att", [False, True])

    tp["Learning Rate"] = trial.suggest_float("lr", 5e-5, 5e-3, log=True)
    tp["Optimizer"] = trial.suggest_categorical("optimizer", ["NAdam"])  # убрали RAdam
    tp["Weight Decay"] = trial.suggest_float("weight_decay", 1e-6, 1e-3, log=True)
    tp["Warm Up"] = trial.suggest_int("warmup", 0, 10)
    tp["Stochastic Depth Prob"] = mp["Stochastic Depth Prob"]
    tp["Batch Size"] = trial.suggest_categorical("batch_size", [8, 16])
    tp["Epochs"] = 10  # по твоему требованию
    return cfg

# ---------- ЦЕЛЕВАЯ ФУНКЦИЯ ----------
def objective(trial):
    base = load_base()
    cfg = suggest(trial, base)
    tdir, rdir, pf2_dir = trial_dirs(trial.number)
    cfg["Train Parameters"]["Results dir"] = str(rdir)
    cfg_path = write_config(cfg, tdir)

    try:
        train_with_cli(cfg_path, rdir)
        postprocess(pf2_dir)
        mcc, auc = get_metrics(pf2_dir)
        copy_with_suffix_all(pf2_dir)  # делаем _tuned-копии всех артефактов
    except Exception as e:
        trial.set_user_attr("fail_reason", str(e))
        raise TrialPruned(str(e))

    trial.set_user_attr("trial_dir", str(tdir))
    trial.set_user_attr("results_dir", str(rdir))
    trial.set_user_attr("pf2_dir", str(pf2_dir))
    trial.set_user_attr("mcc", mcc); trial.set_user_attr("roc_auc", auc)
    return mcc, auc

# ---------- ЗАПУСК ----------
def main():
    sampler = TPESampler(seed=42, n_startup_trials=5)
    study = optuna.create_study(
        study_name=STUDY_NAME,
        directions=["maximize","maximize"],  # MCC, ROC AUC
        sampler=sampler,
        pruner=NopPruner()
    )
    study.optimize(objective, n_trials=N_TRIALS, show_progress_bar=True)

    # скаляризация для выбора "главного" победителя
    def score(t): return 0.7*t.values[0] + 0.3*t.values[1]
    valid_best = [t for t in study.best_trials if t.values is not None]
    if valid_best:
        best = max(valid_best, key=score)
    else:
        best = None

    # сохраняем агрегаты
    with open(ROOT / f"optuna_study{SUFFIX}.json", "w", encoding="utf-8") as f:
        json.dump([
            {
                "number": t.number,
                "state": str(t.state),
                "values": t.values,
                "mcc": (t.values[0] if t.values else None),
                "roc_auc": (t.values[1] if t.values else None),
                "params": t.params,
                "trial_dir": t.user_attrs.get("trial_dir"),
                "results_dir": t.user_attrs.get("results_dir"),
                "pf2_dir": t.user_attrs.get("pf2_dir"),
                "fail_reason": t.user_attrs.get("fail_reason")
            } for t in study.trials
        ], f, ensure_ascii=False, indent=2)

    if best is not None:
        with open(ROOT / f"optuna_best{SUFFIX}.json", "w", encoding="utf-8") as f:
            json.dump({
                "best_trial": best.number,
                "mcc": best.values[0],
                "roc_auc": best.values[1],
                "trial_dir": best.user_attrs.get("trial_dir"),
                "results_dir": best.user_attrs.get("results_dir"),
                "pf2_dir": best.user_attrs.get("pf2_dir"),
                "params": best.params
            }, f, ensure_ascii=False, indent=2)

    # доп. «харвестер»: собираем сводку в один TSV по всем успешным трейалам
    OUT = ROOT / f"optuna_harvest{SUFFIX}.tsv"
    fields = ["trial","mcc","roc_auc","num_blocks","convnext_dim","attn_dim",
              "seq_dropout","attn_dropout","stoch_depth","stoch_depth_att",
              "lr","optimizer","weight_decay","warmup","batch_size","epochs"]
    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for t in study.trials:
            pf2_dir = t.user_attrs.get("pf2_dir")
            cfg_path = (ROOT / f"trial_{t.number:03d}{SUFFIX}" / f"config{SUFFIX}.json")
            if not pf2_dir or not cfg_path.exists():
                continue
            pf2_dir = Path(pf2_dir)
            m = try_read_metrics_summary(pf2_dir) or (t.values if t.values else None)
            if not m:
                continue
            mcc, auc = (m if isinstance(m, tuple) else (m[0], m[1]))
            try:
                with open(cfg_path, "r", encoding="utf-8") as cf:
                    conf = json.load(cf)
                mp = conf.get("Model Parameters", {})
                tp = conf.get("Train Parameters", {})
                net = mp.get("Network Structure", {})
                row = {
                    "trial": f"trial_{t.number:03d}{SUFFIX}",
                    "mcc": mcc, "roc_auc": auc,
                    "num_blocks": net.get("Num Blocks"),
                    "convnext_dim": net.get("Block Dimensions"),
                    "attn_dim": net.get("Attention Dimensions"),
                    "seq_dropout": mp.get("Sequence Dropout"),
                    "attn_dropout": mp.get("Attention Dropout"),
                    "stoch_depth": mp.get("Stochastic Depth Prob"),
                    "stoch_depth_att": mp.get("Stochastic Depth Prob Att"),
                    "lr": tp.get("Learning Rate"),
                    "optimizer": tp.get("Optimizer"),
                    "weight_decay": tp.get("Weight Decay"),
                    "warmup": tp.get("Warm Up"),
                    "batch_size": tp.get("Batch Size"),
                    "epochs": tp.get("Epochs"),
                }
                w.writerow(row)
            except Exception:
                continue

    print("Saved study & harvest under:", ROOT)

if __name__ == "__main__":
    main()
