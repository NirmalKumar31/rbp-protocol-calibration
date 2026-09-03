"""Training loop built to survive being killed.

Two constraints shape this file. Spot and low-priority cloud instances are reclaimed
without warning, and the HPC gpu partition caps a job at eight hours. So every epoch
writes a checkpoint holding model, optimiser, epoch and best-so-far, and a restart
resumes from it instead of starting over.

The run is also idempotent: if the final metrics file already exists, the run returns
immediately. Re-submitting a partially finished sweep therefore costs nothing for the
work already done, which matters when a preemption takes out one task in fifty.

The test split is scored exactly once, at the end, from the best validation checkpoint.

WHAT THE LOCAL FILESYSTEM DOES NOT SURVIVE. A crash leaves the disk intact and the resume
above works. A *preemption* destroys the VM, so the next attempt starts on a machine where
none of these files exist. `on_epoch` is the seam for that: the caller is handed the
checkpoint after every atomic write and can mirror it somewhere durable. This file stays
filesystem-only and therefore testable without a cloud.
"""

import json
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import average_precision_score, roc_auc_score

CHECKPOINT = "checkpoint.pt"
BEST = "best.pt"


def set_seed(seed):
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def atomic_save(obj, path):
    """Write, then rename. A kill mid-write leaves the old file, never a truncated one.

    A half-written checkpoint is worse than no checkpoint: torch.load may well succeed on
    it and training would silently continue from corrupt weights.
    """
    tmp = Path(f"{path}.tmp")
    torch.save(obj, tmp)
    tmp.replace(path)


def labels_of(loader):
    ds = loader.dataset
    if hasattr(ds, "labels"):
        return np.asarray(ds.labels, dtype=np.float32)
    return np.asarray([ds[i][1] for i in range(len(ds))], dtype=np.float32)


def check_splits(loaders):
    """Refuse to start a run that cannot produce a meaningful number.

    A split with one class makes AUROC undefined. sklearn returns nan, `nan > best` is
    False for every epoch, so `best` never updates, early stopping fires at `patience`,
    and the run reports val_auroc -1.0 and a test score from whatever weights the last
    epoch happened to leave behind. Every part of that is silent. Checked up front, before
    a single batch, so the failure costs nothing and says what is wrong.
    """
    for name in ("train", "val", "test"):
        y = labels_of(loaders[name])
        if len(y) == 0:
            raise ValueError(f"{name} split is empty")
        if len(np.unique(y)) < 2:
            raise ValueError(
                f"{name} split has one class only ({int(y.sum())} positive of {len(y)}); "
                "AUROC is undefined, so this run cannot produce a usable number")


@torch.no_grad()
def evaluate(handle, loader, device):
    handle.model.eval()
    probs, ys = [], []
    for seqs, y in loader:
        logits = handle.forward(handle.batch(seqs, device))
        probs.append(torch.sigmoid(logits.float()).cpu().numpy())
        ys.append(y.numpy())
    p = np.concatenate(probs) if probs else np.array([])
    y = np.concatenate(ys) if ys else np.array([])
    if len(y) == 0 or len(np.unique(y)) < 2:
        return {"auroc": float("nan"), "auprc": float("nan"), "n": int(len(y))}, p, y
    return ({"auroc": float(roc_auc_score(y, p)),
             "auprc": float(average_precision_score(y, p)),
             "n": int(len(y))}, p, y)


class Checkpoint:
    """Resumable state on disk. Written atomically so a kill mid-write cannot corrupt it.

    The best-so-far weights live in a SEPARATE file rather than inside this blob. They
    change only when validation improves, while this blob changes every epoch, and for a
    full fine-tune the weights are a third of the payload -- carrying them in both places
    made every epoch's write half again as large for no benefit.
    """

    def __init__(self, path):
        self.path = Path(path)

    def exists(self):
        return self.path.exists()

    def save(self, handle, optimizer, epoch, best, best_epoch, history, elapsed):
        """`best_epoch` is stored rather than derived.

        It was originally recovered by matching `best` against the val scores in
        history, but history rounds them, so the comparison found nothing and every
        resumed run crashed. Persist state you need; do not reconstruct it.
        """
        atomic_save({"model": handle.state(), "optimizer": optimizer.state_dict(),
                     "epoch": epoch, "best": best, "best_epoch": best_epoch,
                     "history": history, "elapsed": elapsed}, self.path)

    def load(self, map_location="cpu"):
        return torch.load(self.path, map_location=map_location, weights_only=False)


def train(handle, loaders, outdir, *, epochs=12, patience=4, seed=7,
          device=None, weight_decay=1e-2, resume=True, log=print,
          on_epoch=None, test_ids=None):
    """Fine-tune one model on one protein. Returns the metrics dict.

    `on_epoch(epoch, outdir)` runs after each epoch's checkpoint is safely on disk. Used
    by the cloud driver to mirror the checkpoint to object storage.

    `test_ids` is the (ids, labels) pair from rbp.train.data.test_ids. Supplying it writes
    row ids alongside the test predictions AND asserts they line up with what the loader
    actually yielded. Without ids a score vector cannot be pooled across folds or aligned
    to the composition arm, which is the whole point of scoring out-of-fold.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    final = outdir / "metrics.json"
    if final.exists():
        log(f"  [skip] metrics.json already in {outdir}")
        return json.loads(final.read_text())

    check_splits(loaders)
    set_seed(seed)
    device = device or torch.device("cpu")
    handle.to(device)
    opt = torch.optim.AdamW(handle.param_groups(), weight_decay=weight_decay)
    lossfn = nn.BCEWithLogitsLoss()
    ckpt = Checkpoint(outdir / CHECKPOINT)
    best_path = outdir / BEST

    start_epoch, best, best_epoch, history, prior = 1, -1.0, 0, [], 0.0
    if resume and ckpt.exists():
        blob = ckpt.load()
        handle.load(blob["model"])
        try:
            opt.load_state_dict(blob["optimizer"])
        except ValueError:
            log("  optimiser state incompatible, continuing with a fresh optimiser")
        start_epoch = blob["epoch"] + 1
        best = blob["best"]
        best_epoch = blob.get("best_epoch", blob["epoch"])
        history = blob["history"]
        # Time already spent, so `seconds` is the run's total and not the time since the
        # most recent resume. A preempted run otherwise under-reports its own cost.
        prior = blob.get("elapsed", 0.0)
        if "best_state" in blob and not best_path.exists():   # legacy checkpoint layout
            atomic_save(blob["best_state"], best_path)
        log(f"  resumed from epoch {blob['epoch']} "
            f"(best val AUROC {best:.4f} at epoch {best_epoch}, {prior:.0f}s already spent)")

    t0 = time.time()

    def elapsed():
        return prior + time.time() - t0

    for epoch in range(start_epoch, epochs + 1):
        handle.model.train()
        total, nb = 0.0, 0
        for seqs, y in loaders["train"]:
            opt.zero_grad()
            logits = handle.forward(handle.batch(seqs, device))
            loss = lossfn(logits.float(), y.to(device))
            loss.backward()
            opt.step()
            total += float(loss.detach())
            nb += 1
        val, _, _ = evaluate(handle, loaders["val"], device)
        row = {"epoch": epoch, "train_loss": round(total / max(nb, 1), 5),
               "val_auroc": round(val["auroc"], 5), "seconds": round(elapsed(), 1)}
        history.append(row)
        log(f"  epoch {epoch:2d}  loss {row['train_loss']:.4f}  "
            f"val AUROC {row['val_auroc']:.4f}")

        if val["auroc"] > best:
            best, best_epoch = val["auroc"], epoch
            atomic_save({k: v.detach().cpu().clone()
                         for k, v in handle.state().items()}, best_path)
        # best.pt is written BEFORE the checkpoint that refers to it. A kill between the
        # two leaves a best.pt from a later epoch than the checkpoint claims, which costs
        # a repeated epoch; the reverse would leave the checkpoint pointing at weights
        # that were never written.
        ckpt.save(handle, opt, epoch, best, best_epoch, history, elapsed())
        if on_epoch:
            on_epoch(epoch, outdir)
        if epoch - best_epoch >= patience:
            log(f"  early stop: no improvement in {patience} epochs")
            break

    if best_path.exists():
        handle.load(torch.load(best_path, map_location="cpu", weights_only=False))
    test, probs, ys = evaluate(handle, loaders["test"], device)

    ids = None
    if test_ids is not None:
        ids, id_labels = test_ids
        id_labels = np.asarray(id_labels, dtype=np.float32)
        # The ids come from a second pass over the same file, so this asserts the two
        # passes agree on both membership and ORDER. If they ever did not, every pooled
        # out-of-fold score would be attached to the wrong row and nothing downstream
        # would notice.
        if len(ids) != len(ys) or not np.array_equal(id_labels, ys):
            raise ValueError(
                f"test ids do not line up with the test loader: {len(ids)} ids against "
                f"{len(ys)} scored rows. Row ids and scores must be in the same order.")

    metrics = {
        "model": handle.name, "label": handle.label, "mode": handle.mode,
        "best_epoch": best_epoch, "val_auroc": round(best, 5),
        "test_auroc": round(test["auroc"], 5), "test_auprc": round(test["auprc"], 5),
        "n_test": test["n"], "epochs_run": history[-1]["epoch"] if history else 0,
        "seconds": round(elapsed(), 1), "seed": seed,
        "device": str(device), **handle.sizes(), "history": history,
    }
    out = {"prob": probs, "label": ys}
    if ids is not None:
        out["id"] = np.asarray(ids)
    np.savez_compressed(outdir / "test_predictions.npz", **out)
    final.write_text(json.dumps(metrics, indent=2))
    (outdir / CHECKPOINT).unlink(missing_ok=True)   # finished, so stop carrying it
    return metrics
