"""Stage 5: train one (cell, protein, model, fold) run.

    python scripts/train.py --protein QKI --cell K562 --model cnn --fold 0
    python scripts/train.py --protein QKI --cell K562 --model cnn --limit 400   # smoke

Array-friendly: with --all the grid is enumerated in a fixed order and SLURM_ARRAY_TASK_ID
or --index picks one, so the same command works locally, on Slurm, and as a cloud job.

WHAT DEFINES A RUN. Four things, and every one of them must appear in the output path or
runs overwrite each other: the negative arm, the cell line, the protein, the model, and
the fold. The panel, the arm and the pair threshold all come from config through
rbp.utils.panel.study -- this script used to read a stale 17-protein config/panel_final.tsv
and load from the GC directory in one hardcoded cell line, which would have swept 85 wrong
runs and called it the study.

Runs are idempotent. If metrics.json exists the task returns immediately, so re-submitting
after a preemption only redoes what was lost.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rbp.models import registry  # noqa: E402
from rbp.train import data as tdata  # noqa: E402
from rbp.train import trainer  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def grid(cfg, models, arm, folds):
    """(cell, protein, model, fold), in a fixed order.

    Sorted by dataset key so the enumeration is stable across machines and reruns: an
    array job's task index is only meaningful if index 7 means the same run everywhere.
    require_files=False for the same reason: a worker missing one dataset must not
    renumber every task after it.
    """
    paths, _ = panelmod.study(cfg, arm, require_files=False)
    out = []
    for key in sorted(paths):
        protein, cell = key.split(":")
        for m in models:
            for f in folds:
                out.append((cell, protein, m, f))
    return out


def subsample(loaders, limit):
    """Shrink every split, for a fast smoke test on real data."""
    import torch
    from torch.utils.data import DataLoader, Subset
    out = {}
    for k, dl in loaders.items():
        n = min(limit, len(dl.dataset))
        idx = torch.randperm(len(dl.dataset))[:n].tolist()
        out[k] = DataLoader(Subset(dl.dataset, idx), batch_size=dl.batch_size,
                            shuffle=(k == "train"), collate_fn=tdata.collate)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--protein")
    p.add_argument("--cell")
    p.add_argument("--model")
    p.add_argument("--all", action="store_true", help="use the grid + --index")
    p.add_argument("--index", type=int, default=None)
    p.add_argument("--models", default=None)
    p.add_argument("--arm", default=None, choices=sorted(panelmod.ARMS),
                   help="negative arm. Default: negatives.primary_arm from config")
    p.add_argument("--config", default=None)
    p.add_argument("--device", default=None)
    p.add_argument("--epochs", type=int, default=None)
    p.add_argument("--batch-size", type=int, default=None)
    p.add_argument("--limit", type=int, default=None, help="rows per split (smoke test)")
    p.add_argument("--outroot", default=str(ROOT / "results/runs"))
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--fold", type=int, default=None,
                   help="cross-validation fold to hold out. Omit for the legacy "
                        "fixed train/val/test split, which is NOT the primary protocol")
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    arm = panelmod.arm_of(cfg, a.arm)
    k = cfg.cv["k"]

    all_models = registry.names(cfg)
    models = [m.strip() for m in a.models.split(",")] if a.models else list(all_models)
    for m in models:
        if m not in all_models:
            sys.exit(f"unknown model {m!r}; config defines {list(all_models)}")

    if a.all or a.index is not None or "SLURM_ARRAY_TASK_ID" in os.environ:
        folds = [a.fold] if a.fold is not None else list(range(k))
        g = grid(cfg, models, arm, folds)
        i = a.index if a.index is not None else int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
        if i >= len(g):
            sys.exit(f"index {i} out of range: grid has {len(g)} tasks")
        cell, protein, model, fold = g[i]
        print(f"grid task {i}/{len(g)-1}: {protein} {cell} x {model} fold{fold}")
    else:
        if not (a.protein and a.model and a.cell):
            sys.exit("give --protein, --cell and --model, or --all with --index")
        cell, protein, model, fold = a.cell, a.protein, a.model, a.fold

    tcfg = cfg["train"]
    epochs = a.epochs or tcfg["epochs"]
    batch = a.batch_size or tcfg["batch_size"]

    dataset = panelmod.data_dir(cell, arm) / protein / "dataset.tsv"
    if not dataset.exists():
        sys.exit(f"missing {dataset}")

    device = registry.device_of(a.device)
    t0 = time.time()
    handle = registry.build(model, cfg)
    sizes = handle.sizes()
    print(f"{protein} {cell} {arm} fold{fold} x {handle.label}  mode={handle.mode}  "
          f"params {sizes['params_total']:,} (trainable {sizes['params_trainable']:,}, "
          f"{sizes['trainable_frac']:.2%})  device={device}  "
          f"built in {time.time()-t0:.1f}s", flush=True)

    dl = tdata.loaders(dataset, batch, seed=cfg.seed, fold=fold, k=k)
    # Row ids are only meaningful for the real loaders. --limit takes a random subset, so
    # the ids from a second pass over the file would not match, and the alignment check in
    # the trainer would (correctly) reject them.
    ids = None if a.limit else tdata.test_ids(dataset, fold=fold, k=k)
    if a.limit:
        dl = subsample(dl, a.limit)
    for split in ("train", "val", "test"):
        print(f"  {split:5} {tdata.class_balance(dl[split])}")

    outdir = Path(a.outroot) / arm / cell / protein / model
    if fold is not None:
        outdir = outdir / f"fold{fold}"
    metrics = trainer.train(handle, dl, outdir, epochs=epochs, patience=tcfg["patience"],
                            seed=cfg.seed, device=device,
                            weight_decay=float(tcfg["weight_decay"]),
                            resume=not a.no_resume, test_ids=ids)
    print(f"\n{protein} {cell} {handle.label}: val {metrics['val_auroc']:.4f}  "
          f"test AUROC {metrics['test_auroc']:.4f}  AUPRC {metrics['test_auprc']:.4f}  "
          f"({metrics['seconds']:.0f}s)")
    print(json.dumps({k_: v for k_, v in metrics.items() if k_ != "history"}, indent=2))


if __name__ == "__main__":
    main()
