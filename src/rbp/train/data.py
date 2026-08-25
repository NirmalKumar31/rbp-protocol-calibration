"""Dataset plumbing.

Sequences stay as strings until the last moment, because each backbone tokenises
differently. The model handle owns that conversion, so the loader is trivial and the
same batches reach every model.
"""

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset


class WindowDataset(Dataset):
    def __init__(self, seqs, labels):
        self.seqs = list(seqs)
        self.labels = np.asarray(labels, dtype=np.float32)

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, i):
        return self.seqs[i], self.labels[i]


def collate(batch):
    seqs = [b[0] for b in batch]
    y = torch.tensor([b[1] for b in batch], dtype=torch.float32)
    return seqs, y


def split_frame(df, split, fold=None, k=5):
    """Rows belonging to one role, from either protocol.

    TWO PROTOCOLS LIVE IN THIS FILE AND ONLY ONE IS PRIMARY.

    `fold=None` uses the `split` column: a single fixed train/val/test partition. That is
    the older protocol, kept for the shuffle-arm sensitivity analysis.

    `fold=i` uses the `fold` column and derives the roles: fold i is test, fold i+1 is
    validation, the rest train. Run for i in 0..k-1 and every row is scored exactly once
    by a model that never saw it -- which is what "pooled out-of-fold" in params.yaml
    means, and what the composition control it is compared against already does.

    The trainer read only the `split` column until 2026-08-23, so the model arm and the
    composition arm would have been measured under different protocols and could not have
    been compared. Every dataset has carried a `fold` column since preprocessing; it was
    simply never wired up.
    """
    if fold is None:
        return df[df.split == split]
    from ..data.splits import split_of_fold
    roles = df.fold.map(lambda f: split_of_fold(int(f), fold, k))
    return df[roles == split]


def load_split(path, split, fold=None, k=5):
    df = pd.read_csv(path, sep="\t")
    d = split_frame(df, split, fold, k)
    return WindowDataset(d.seq_rna.values, d.label.values)


def loaders(path, batch_size, seed=7, num_workers=0, fold=None, k=5):
    """train/val/test loaders. Only train is shuffled, with a seeded generator.

    With `fold` set the three roles come from the cross-validation fold map instead of the
    frozen split column; see split_frame.
    """
    g = torch.Generator().manual_seed(seed)
    out = {}
    for split in ("train", "val", "test"):
        ds = load_split(path, split, fold, k)
        out[split] = DataLoader(
            ds, batch_size=batch_size, shuffle=(split == "train"),
            collate_fn=collate, num_workers=num_workers,
            generator=g if split == "train" else None)
    return out


def test_ids(path, fold=None, k=5):
    """Row ids of the held-out set, in loader order.

    The sweep must emit an out-of-fold score for every row, and a bare score vector is
    useless without knowing which rows it belongs to. Returned separately rather than
    carried through the Dataset so the collate function stays the same for every model.
    """
    df = pd.read_csv(path, sep="\t")
    d = split_frame(df, "test", fold, k)
    return d.id.tolist(), d.label.to_numpy()


def class_balance(loader):
    """Works on a Subset too, so smoke tests report the balance they actually saw."""
    ds = loader.dataset
    if hasattr(ds, "labels"):
        y = ds.labels
    else:                                    # torch Subset
        y = np.asarray([ds[i][1] for i in range(len(ds))], dtype=np.float32)
    return {"n": int(len(y)), "positives": int(y.sum()),
            "frac_positive": round(float(y.mean()), 4)}
