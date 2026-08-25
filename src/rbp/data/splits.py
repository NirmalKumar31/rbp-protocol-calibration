"""Chromosome-level assignment: a fixed three-way split, or k cross-validation folds.

Splitting by whole chromosome, with one fixed assignment reused for every protein, is
what stops a locus appearing on both sides of the split. Nearby genomic positions share
sequence, so a random row split leaks and inflates scores.

This prevents *locus* leakage. It does not prevent similar sequences on different
chromosomes landing in different splits, which is what the leakage audit measures
later.

Cross-validation is the primary protocol. A single held-out split spends 20% of the
data on test and reports a ranking that rests on one arbitrary chromosome partition,
which is unanswerable to "would the ranking hold under a different split?". Grouped
k-fold scores every window exactly once out-of-fold, so the AUROC uses all pairs, each
fold model trains on more data, and the partition question is answered structurally.
The three-way helpers are kept because the split optimiser and its tests predate the
change and the fold code is a strict generalisation of them.
"""


def split_of(chrom, cfg_split):
    test = set(cfg_split["test"])
    val = set(cfg_split["val"])
    if chrom in test:
        return "test"
    if chrom in val:
        return "val"
    return "train"


def assign(rows, cfg_split):
    for r in rows:
        r["split"] = split_of(r["chrom"], cfg_split)
    return rows


def proportions(rows):
    n = len(rows) or 1
    out = {}
    for s in ("train", "val", "test"):
        out[s] = sum(1 for r in rows if r["split"] == s) / n
    return out


def peak_counts(peak_paths, chroms):
    """(protein x chromosome) peak-count matrix, for choosing the assignment."""
    import numpy as np

    from .windows import read_peaks
    names = sorted(peak_paths)
    idx = {c: i for i, c in enumerate(chroms)}
    m = np.zeros((len(names), len(chroms)), dtype=np.int64)
    for r, p in enumerate(names):
        for chrom, _, _, _ in read_peaks(peak_paths[p]):
            j = idx.get(chrom)
            if j is not None:
                m[r, j] += 1
    return names, m


def assignment_loss(counts, assign, target):
    """Sum over proteins of squared deviation from the target split proportions.

    Every protein contributes equally, so a small protein cannot be sacrificed to
    tidy up a large one.
    """
    import numpy as np
    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    loss = 0.0
    worst = 0.0
    for k, t in enumerate(target):
        prop = counts[:, assign == k].sum(axis=1, keepdims=True) / totals
        d = prop - t
        loss += float((d ** 2).sum())
        worst = max(worst, float(np.abs(d).max()))
    return loss, worst


def optimize_assignment(counts, target=(0.64, 0.16, 0.20), min_per_split=3,
                        restarts=40, iters=4000, seed=7):
    """Pick one chromosome -> group assignment that fits every protein well.

    A fixed chromosome split is what prevents locus leakage, but which chromosomes go
    where was originally chosen blind. Peaks are distributed very unevenly (chr19 is
    gene dense), so a blind choice leaves some proteins at 13% test and others at 24%.

    This searches for an assignment minimising deviation from the target proportions
    across all proteins at once. It optimises on peak COUNTS only, never on labels,
    model output or performance, so it introduces no leakage - it is stratification,
    the same principle as a stratified train/test split.

    The number of groups is `len(target)`, so the same search produces a three-way
    split or k equal cross-validation folds.
    """
    import numpy as np
    rng = np.random.default_rng(seed)
    n_chrom = counts.shape[1]
    n_groups = len(target)
    if n_chrom < n_groups * min_per_split:
        raise ValueError(f"{n_chrom} chromosomes cannot fill {n_groups} groups with "
                         f"at least {min_per_split} each")
    best = None

    for _ in range(restarts):
        assign = rng.integers(0, n_groups, n_chrom)
        for k in range(n_groups):               # guarantee no group is starved
            if (assign == k).sum() < min_per_split:
                assign[rng.choice(n_chrom, min_per_split, replace=False)] = k
        cur, _ = assignment_loss(counts, assign, target)
        for _ in range(iters):
            i = int(rng.integers(0, n_chrom))
            old = assign[i]
            new = int(rng.integers(0, n_groups))
            if new == old:
                continue
            assign[i] = new
            if min((assign == k).sum() for k in range(n_groups)) < min_per_split:
                assign[i] = old
                continue
            cand, _ = assignment_loss(counts, assign, target)
            if cand <= cur:
                cur = cand
            else:
                assign[i] = old
        loss, worst = assignment_loss(counts, assign, target)
        if best is None or loss < best[0]:
            best = (loss, worst, assign.copy())
    return best


def optimize_folds(counts, k=5, min_per_fold=3, restarts=40, iters=4000, seed=7):
    """Balanced chromosome -> fold assignment for grouped k-fold CV.

    Equal target mass per fold, since every fold serves as the test set exactly once
    and unequal folds would make the per-fold AUROCs incomparable.
    """
    return optimize_assignment(counts, target=(1.0 / k,) * k,
                               min_per_split=min_per_fold, restarts=restarts,
                               iters=iters, seed=seed)


def fold_roles(fold, k):
    """Which folds are test, val and train when `fold` is the held-out one.

    Validation is the next fold cyclically. Taking it from inside the k folds keeps
    every chromosome in exactly one role per iteration and costs no extra data: over
    the k iterations each fold is test once, val once, and train k-2 times.
    """
    if not 0 <= fold < k:
        raise ValueError(f"fold {fold} out of range for k={k}")
    if k < 3:
        raise ValueError("k must be at least 3 to leave a validation fold")
    val = (fold + 1) % k
    return fold, val, [f for f in range(k) if f not in (fold, val)]


def assign_folds(rows, fold_map):
    """Tag every row with its chromosome's fold index.

    The fold index is stored once in the dataset and the train/val/test roles are
    derived per iteration, so one dataset file serves all k runs. Materialising k
    copies would multiply storage and let the copies drift apart.
    """
    for r in rows:
        r["fold"] = fold_map[r["chrom"]]
    return rows


def split_of_fold(row_fold, fold, k):
    test, val, _ = fold_roles(fold, k)
    if row_fold == test:
        return "test"
    if row_fold == val:
        return "val"
    return "train"


def fold_proportions(rows, k):
    n = len(rows) or 1
    return {f: sum(1 for r in rows if r["fold"] == f) / n for f in range(k)}


def check_folds(rows, k):
    """Structural checks on a fold assignment. Returns a list of problems.

    Catches the two failures that matter: a chromosome landing in two folds (leakage),
    and an empty fold (an unusable iteration). Both are silent at training time and
    only show up as an implausibly good or missing result.
    """
    problems = []
    seen = {}
    for r in rows:
        c, f = r["chrom"], r["fold"]
        if c in seen and seen[c] != f:
            problems.append(f"{c} in folds {seen[c]} and {f}")
        seen[c] = f
    for f in range(k):
        if not any(r["fold"] == f for r in rows):
            problems.append(f"fold {f} is empty")
    return sorted(set(problems))


def check_disjoint(rows):
    """No chromosome may appear in more than one split. Returns offenders."""
    seen = {}
    bad = {}
    for r in rows:
        c, s = r["chrom"], r["split"]
        if c in seen and seen[c] != s:
            bad.setdefault(c, set()).update({seen[c], s})
        seen[c] = s
    return bad
