"""Estimate the GPU cost of the training sweep before spending anything.

The budget is a hard constraint ($300 of GCP trial credit), and the fold count k is
the single biggest lever on it, so the estimate has to exist before k is chosen.

Method: count the sequences the sweep will actually push through each model, convert
to FLOPs, divide by realised accelerator throughput. FLOPs per sequence uses the
standard 6 * params * tokens for a training step (2 forward, 4 backward). LoRA models
still backpropagate through the frozen stack, so they get the same count -- only the
optimiser state shrinks, not the maths. That makes this an over-estimate for LoRA,
which is the direction an estimate should err.

The two soft numbers are MFU (fraction of peak FLOPS actually realised) and the
effective epoch count under early stopping. Both are swept, so the output is a range
with the assumptions visible rather than a single figure to be taken on faith.

    python scripts/cost_model.py --folds 5
    python scripts/cost_model.py --folds 5 --census config/panel_census.tsv,config/panel_census_HepG2.tsv
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd  # noqa: E402

from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

# Realised throughput, not spec sheet. Peak dense bf16 TFLOPS x assumed MFU.
# 101-token sequences are short, so attention is cheap but the pipeline is
# launch-bound; 25-40% MFU is the honest band for this shape of workload.
GPUS = {
    # name          peak TFLOPS   $/hr on-demand   $/hr spot
    "T4":         (65.0,  0.35, 0.11),
    "L4":         (121.0, 0.71, 0.22),
    "A100-40GB":  (312.0, 3.67, 1.10),
}
# g2/n1 host cost rides along with the accelerator; add a flat vCPU+RAM charge
HOST_PER_HR = 0.15


def sweep_sequences(census_paths, min_pairs, folds):
    """Total training-sequence-visits for one architecture, one epoch, all folds.

    Each fold trains on (k-1)/k of a protein's rows, and there are k folds, so one
    epoch over the full CV protocol sees (k-1) times the dataset -- not k times.
    """
    frames = []
    for p in census_paths:
        d = pd.read_csv(p, sep="\t")
        d["source"] = Path(p).stem
        frames.append(d)
    d = pd.concat(frames, ignore_index=True)

    passing = d[d.windows >= min_pairs].copy()
    rows = passing.windows.sum() * 2          # positives + matched negatives
    return passing, rows, rows * (folds - 1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--census", default="config/panel_census.tsv")
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--min-pairs", type=int, default=400,
                   help="inclusion threshold on out-of-fold scored pairs")
    p.add_argument("--epochs", type=int, default=None, help="override train.epochs")
    a = p.parse_args()

    cfg = cfgmod.load(a.config)
    paths = [ROOT / c for c in a.census.split(",")]
    epochs = a.epochs or cfg.train["epochs"]

    passing, rows, per_epoch = sweep_sequences(paths, a.min_pairs, a.folds)
    tokens = cfg.windows["size"]

    print(f"panel: {len(passing)} datasets clear {a.min_pairs} pairs")
    print(f"  rows (pos+neg):        {rows:,}")
    print(f"  folds:                 {a.folds}  -> {per_epoch:,} sequence-visits/epoch")
    print(f"  epochs (cap):          {epochs}, early stopping patience "
          f"{cfg.train['patience']}\n")

    # Per-architecture FLOPs for the whole sweep at the epoch cap.
    arch = []
    for name, m in cfg["models"].items():
        pm = m.get("params_m")
        if pm is None:
            pm = 0.007      # the CNN, 7,089 params
        flops = 6 * pm * 1e6 * tokens * per_epoch * epochs
        arch.append((name, m.get("label", name), pm, flops))
    total = sum(f for *_, f in arch)

    print(f"{'model':16} {'params(M)':>10} {'PFLOPs':>10} {'share':>7}")
    for name, label, pm, f in sorted(arch, key=lambda r: -r[3]):
        print(f"{label:16} {pm:10.2f} {f/1e15:10.1f} {100*f/total:6.1f}%")
    print(f"{'TOTAL':16} {'':10} {total/1e15:10.1f}\n")

    print(f"cost at the {epochs}-epoch cap (early stopping typically halves this):\n")
    print(f"{'gpu':12} {'MFU':>5} {'gpu-hrs':>9} {'on-demand':>11} {'spot':>9}")
    for g, (tf, od, sp) in GPUS.items():
        for mfu in (0.25, 0.40):
            hrs = total / (tf * 1e12 * mfu) / 3600
            print(f"{g:12} {mfu:5.0%} {hrs:9.1f} "
                  f"${hrs*(od+HOST_PER_HR):10.0f} ${hrs*(sp+HOST_PER_HR):8.0f}")

    print("\nnot included: storage (~$5), egress, pipeline orchestration overhead,")
    print("failed/preempted runs. Budget 1.5x the spot figure for those.")


if __name__ == "__main__":
    main()
