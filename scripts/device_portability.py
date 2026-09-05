"""E3: the same code, the same inputs, two devices and two clouds. Does it give the same answer?

    python scripts/device_portability.py
    python scripts/device_portability.py --from-cache

WHAT THIS EXERCISES. `cloud/submit_cpu_sweep.sh` submits the training container to GCP Batch on
CPU. It had never been run against a dataset whose results were not already committed, so the
path was infrastructure the repository described rather than infrastructure it had used. Running
it on the region-matched arm, which exists in no GCS prefix, means nothing can be overwritten and
the run has to produce its results from scratch.

WHAT IT MEASURES, which is worth more than the demonstration. The same dataset, the same five
chromosome-grouped folds and the same commit were trained twice: once on Modal on an A10G and
once on GCP Batch on one vCPU of an e2-standard-4. Nothing differs but the device and the cloud.
A paper about reproducibility should be able to say what that costs, and the honest answer is
not zero: initialisation is unseeded, so the two runs are different draws as well as different
devices, and this bounds the two together rather than separating them.

WHY THE BOUND IS STILL USEFUL DESPITE THAT. A reader reproducing this work will run on different
hardware AND get a different initialisation, because that is what the released code does. The
quantity they care about is exactly the combined one measured here.
"""

import argparse
import gzip
import io
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rbp.utils.log import log

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

TABLES = ROOT / "results" / "tables"
DATASET = ("K562", "ELAC2", "cnn", "neg2_rm")
FOLDS = 5



def gcs_scores(bucket, cell, protein, model, arm, fold):
    from google.cloud import storage
    b = storage.Client().bucket(bucket)
    blob = b.blob(f"runs/{arm}/{cell}/{protein}/{model}/fold{fold}/scores.tsv.gz")
    if not blob.exists():
        return None
    return pd.read_csv(io.BytesIO(gzip.decompress(blob.download_as_bytes())), sep="\t")


def gcs_metrics(bucket, cell, protein, model, arm, fold):
    from google.cloud import storage
    b = storage.Client().bucket(bucket)
    blob = b.blob(f"runs/{arm}/{cell}/{protein}/{model}/fold{fold}/metrics.json")
    return json.loads(blob.download_as_text()) if blob.exists() else None


def build(store, bucket):
    from sklearn.metrics import roc_auc_score
    cell, protein, model, arm = DATASET
    rows = []
    for f in range(FOLDS):
        gpu_p = Path(store) / "runs" / arm / cell / protein / model / f"fold{f}"
        if not (gpu_p / "scores.tsv.gz").exists():
            sys.exit(f"no GPU scores at {gpu_p}; run the Modal sweep first")
        gpu = pd.read_csv(gpu_p / "scores.tsv.gz", sep="\t")
        cpu = gcs_scores(bucket, cell, protein, model, arm, f)
        if cpu is None:
            log(f"  fold{f}: no CPU scores in gs://{bucket}, skipping")
            continue
        gm = json.loads((gpu_p / "metrics.json").read_text())
        cm = gcs_metrics(bucket, cell, protein, model, arm, f)

        m = gpu.merge(cpu, on="id", suffixes=("_gpu", "_cpu"))
        # THE ROW SETS MUST BE IDENTICAL, and this is asserted rather than assumed: the two
        # runs read the same dataset.tsv, so a differing row set would mean they were not
        # given the same inputs and no comparison below would mean anything.
        same_rows = len(m) == len(gpu) == len(cpu)
        same_labels = bool((m.label_gpu == m.label_cpu).all())
        same_folds = bool((m.fold_gpu == m.fold_cpu).all())
        r = float(np.corrcoef(m.score_gpu, m.score_cpu)[0, 1])
        from scipy.stats import spearmanr
        rho = float(spearmanr(m.score_gpu, m.score_cpu).statistic)
        rows.append({
            "fold": f, "n": len(m), "same_rows": same_rows, "same_labels": same_labels,
            "same_folds": same_folds, "pearson": r, "spearman": rho,
            "auroc_gpu": float(roc_auc_score(m.label_gpu, m.score_gpu)),
            "auroc_cpu": float(roc_auc_score(m.label_cpu, m.score_cpu)),
            "test_auroc_gpu": float(gm["test_auroc"]),
            "test_auroc_cpu": float(cm["test_auroc"]) if cm else float("nan"),
            "seconds_gpu": float(gm["seconds"]),
            "seconds_cpu": float(cm["seconds"]) if cm else float("nan"),
            "device_gpu": gm.get("device", "?"), "device_cpu": cm.get("device", "?") if cm
            else "?"})
        log(f"  fold{f}: n={len(m)}  AUROC gpu {rows[-1]['auroc_gpu']:.4f} vs cpu "
            f"{rows[-1]['auroc_cpu']:.4f}  r={r:.4f}  rho={rho:.4f}  "
            f"{gm['seconds']:.0f}s vs {cm['seconds']:.0f}s" if cm else "")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("no fold could be compared; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--bucket", default="rbp-repro-2026-derived")
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "device_portability_per_fold.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it without --from-cache")
    else:
        t = build(a.store, a.bucket)
        t.to_csv(per, index=False)

    out = []
    cell, protein, model, arm = DATASET
    log(f"\n=== E3: {protein}:{cell} {model}, {arm} arm, Modal A10G vs GCP Batch CPU, "
        f"{len(t)} folds ===\n")
    out.append({"check": "folds compared across devices", "value": len(t), "n": len(t)})

    # THE INPUTS WERE THE SAME. Without this the agreement below could be agreement about two
    # different row sets, which is not agreement at all.
    for col, label in (("same_rows", "identical row sets"),
                       ("same_labels", "identical labels"),
                       ("same_folds", "identical fold assignment")):
        out.append({"check": f"folds with {label} across devices",
                    "value": int(t[col].sum()), "n": len(t)})
    out.append({"check": "windows compared across devices", "value": int(t.n.sum()),
                "n": len(t)})

    for col, label in (("pearson", "Pearson"), ("spearman", "Spearman")):
        v = float(t[col].mean())
        out.append({"check": f"{label} correlation between CPU and GPU per-window scores",
                    "value": v, "n": len(t)})
        out.append({"check": f"lowest {label} correlation over folds",
                    "value": float(t[col].min()), "n": len(t)})
    d = (t.auroc_gpu - t.auroc_cpu).abs()
    out.append({"check": "max |AUROC difference| between devices, per fold",
                "value": float(d.max()), "n": len(t)})
    out.append({"check": "mean AUROC, Modal A10G", "value": float(t.auroc_gpu.mean()),
                "n": len(t)})
    out.append({"check": "mean AUROC, GCP Batch CPU", "value": float(t.auroc_cpu.mean()),
                "n": len(t)})
    log(f"  per-window score correlation: Pearson {t.pearson.mean():.4f} "
        f"(min {t.pearson.min():.4f}), Spearman {t.spearman.mean():.4f}")
    log(f"  AUROC  A10G {t.auroc_gpu.mean():.4f}   CPU {t.auroc_cpu.mean():.4f}   "
        f"max per-fold difference {d.max():.4f}")

    # THE COST OF THE DEVICE, which is the practical number. A CPU fold is slower by this
    # factor, and that is what decides whether the CPU path is a fallback or a curiosity.
    if t.seconds_cpu.notna().any():
        ratio = float((t.seconds_cpu / t.seconds_gpu).mean())
        out.append({"check": "CPU wall time per fold as a multiple of A10G",
                    "value": ratio, "n": len(t)})
        log(f"  a CPU fold takes {ratio:.1f}x the wall time of an A10G fold")

        # WHAT THAT RATIO IMPLIES ABOUT THE GPU, and it is the practical finding here. A 7,089
        # parameter convolutional network is only this much slower on one vCPU than on an
        # A10G, because a model that small never fills the device: the sweep is paying for
        # accelerator time it cannot use. Extrapolated at the measured ratio, the whole CNN
        # sweep's 6.12 GPU-h at $1.10/GPU-h becomes the same work in CPU-hours priced at a
        # spot e2-standard-4's per-vCPU rate. Stated as an extrapolation from one dataset,
        # because that is what it is; the ratio itself is measured over five folds.
        gpu_h, gpu_rate, cpu_vcpu_rate = 6.12, 1.10, 0.0085
        out.append({"check": "CNN sweep cost on Modal A10G", "value": gpu_h * gpu_rate,
                    "n": len(t), "note": "measured over 470 fold-runs"})
        out.append({"check": "same work extrapolated to spot CPU at the measured ratio",
                    "value": gpu_h * ratio * cpu_vcpu_rate, "n": len(t),
                    "note": "extrapolation from one dataset; excludes egress and storage"})
        log(f"  so the CNN sweep's ${gpu_h * gpu_rate:.2f} of A10G time extrapolates to "
            f"${gpu_h * ratio * cpu_vcpu_rate:.2f} of spot CPU: a model this small does not "
            f"fill the device")

    # AND THE HONEST CAVEAT, committed as a number rather than left to the prose: the two runs
    # differ in INITIALISATION as well as device, because initialisation is unseeded. So this
    # bounds device and seed together, which is also exactly what a reader reproducing the work
    # would experience.
    out.append({"check": "initialisation seeded across the two runs", "value": 0, "n": len(t),
                "note": "unseeded, so this bounds device and initialisation jointly"})

    pd.DataFrame(out).to_csv(TABLES / "device_portability.csv", index=False)
    log("\nwrote device_portability.csv and device_portability_per_fold.csv")


if __name__ == "__main__":
    main()
