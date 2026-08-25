"""DeepBind-style convolutional baseline.

Two convolutional blocks act as learned motif detectors, then a small MLP classifies.
The second pooling step is global: it takes each filter's strongest response anywhere in
the window, which is what makes the model position invariant. The EDA showed the
discriminative signal sits ~15 nt off centre and varies per protein, so position
invariance is a requirement here rather than a convenience.

Four learnable layers, ~7K parameters. Small on purpose: it is the honest floor that
the pretrained models have to beat.
"""

import torch
import torch.nn as nn

BASES = "ACGU"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}


def one_hot(seq):
    """(4, L) float tensor. Unknown characters become all-zero columns."""
    x = torch.zeros(len(BASES), len(seq))
    for i, ch in enumerate(seq):
        j = BASE_INDEX.get(ch)
        if j is not None:
            x[j, i] = 1.0
    return x


def one_hot_batch(seqs):
    return torch.stack([one_hot(s) for s in seqs])


class DeepBindCNN(nn.Module):
    def __init__(self, channels=(16, 32), kernels=(12, 8), hidden=64,
                 pool=4, dropout=0.5):
        super().__init__()
        c1, c2 = channels
        k1, k2 = kernels
        self.features = nn.Sequential(
            nn.Conv1d(len(BASES), c1, k1, padding=k1 // 2),
            nn.ReLU(),
            nn.MaxPool1d(pool),
            nn.Conv1d(c1, c2, k2, padding=k2 // 2),
            nn.ReLU(),
            nn.AdaptiveMaxPool1d(1),          # strongest response anywhere in the window
        )
        self.head = nn.Sequential(
            nn.Flatten(),
            nn.Linear(c2, hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, 1),
        )

    def forward(self, x):
        return self.head(self.features(x)).squeeze(-1)

    @property
    def n_params(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


def first_layer_filters(model):
    """Conv1 weights as (n_filters, 4, width), for reading learned motifs later."""
    return model.features[0].weight.detach().cpu().numpy()
