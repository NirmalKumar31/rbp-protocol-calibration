"""How much do different RBPs bind the same places?

This decides the training design. If a 101-nt window built from protein A's peak
usually also covers protein B's peak, then one shared multi-task model can be given
many labels per window (cross-labelling) and is worth the extra machinery. If windows
are mostly private to one protein, per-protein models with their own matched negatives
stay the better design.

Overlap is measured at the window level, not the peak level, because the window is
what a model actually sees.
"""

import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def load_peaks(path):
    """chrom -> (starts, ends, midpoints) as sorted arrays. Reads narrowPeak.gz."""
    by_chrom = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            by_chrom.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    out = {}
    for c, iv in by_chrom.items():
        iv.sort()
        a = np.array(iv, dtype=np.int64)
        out[c] = (a[:, 0], a[:, 1], (a[:, 0] + a[:, 1]) // 2)
    return out


def merge(starts, ends):
    """Collapse overlapping intervals into a sorted, disjoint set."""
    if len(starts) == 0:
        return starts, ends
    order = np.argsort(starts)
    s, e = starts[order], ends[order]
    ms, me = [s[0]], [e[0]]
    for i in range(1, len(s)):
        if s[i] <= me[-1]:
            me[-1] = max(me[-1], e[i])
        else:
            ms.append(s[i])
            me.append(e[i])
    return np.array(ms), np.array(me)


def windows_of(peaks, size):
    """chrom -> (start, end) of the size-nt window centred on each peak midpoint."""
    half = size // 2
    return {c: (mid - half, mid + half + 1) for c, (_, _, mid) in peaks.items()}


def covered(query, target):
    """Boolean mask: does each query interval overlap any target interval?

    query/target are per-chromosome dicts of (starts, ends). Target intervals are
    merged first, so a single binary search per query interval suffices.
    """
    merged = {c: merge(s, e) for c, (s, e) in target.items()}
    out = {}
    for c, (qs, qe) in query.items():
        if c not in merged or len(merged[c][0]) == 0:
            out[c] = np.zeros(len(qs), dtype=bool)
            continue
        ts, te = merged[c]
        idx = np.searchsorted(ts, qe, side="right") - 1
        ok = idx >= 0
        hit = np.zeros(len(qs), dtype=bool)
        hit[ok] = te[idx[ok]] > qs[ok]
        out[c] = hit
    return out


def matrix(peak_paths, window=101):
    """Asymmetric overlap matrix: row A, column B = fraction of A's windows hitting B.

    Also returns per-protein 'shared' stats: how much of each protein's binding is
    also covered by at least one other protein in the panel.
    """
    peaks = {p: load_peaks(f) for p, f in peak_paths.items()}
    wins = {p: windows_of(pk, window) for p, pk in peaks.items()}
    names = sorted(peaks)

    m = pd.DataFrame(0.0, index=names, columns=names)
    n_windows = {}
    any_other = {}

    for a in names:
        qa = {c: (s, e) for c, (s, e) in wins[a].items()}
        total = sum(len(s) for s, _ in qa.values())
        n_windows[a] = total
        shared_any = {c: np.zeros(len(s), dtype=bool) for c, (s, _) in qa.items()}
        for b in names:
            tb = {c: (s, e) for c, (s, e, _) in
                  ((c, (s, e, mid)) for c, (s, e, mid) in peaks[b].items())}
            hits = covered(qa, tb)
            n_hit = sum(int(h.sum()) for h in hits.values())
            m.loc[a, b] = n_hit / total if total else np.nan
            if b != a:
                for c, h in hits.items():
                    shared_any[c] |= h
        any_other[a] = sum(int(h.sum()) for h in shared_any.values()) / total if total else np.nan

    stats = pd.DataFrame({
        "n_windows": pd.Series(n_windows),
        "frac_shared_with_any_other": pd.Series(any_other),
    })
    stats["frac_private"] = 1 - stats.frac_shared_with_any_other
    return m, stats


def peak_paths_from(dirpath):
    """protein -> peak file, from data/raw/peaks/<PROT>.<ACC>.bed.gz."""
    out = {}
    for f in sorted(Path(dirpath).glob("*.bed.gz")):
        out[f.name.split(".")[0]] = f
    return out
