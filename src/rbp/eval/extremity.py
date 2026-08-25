"""Compositional extremity: how unusual a protein's binding sites are, before any model.

THE PROBLEM THIS SOLVES. Our negatives are matched on GC content, which constrains G+C as
a single total and leaves the other fifteen degrees of freedom of dinucleotide composition
free. So the obvious objection to any composition result is: "your finding is an artefact
of weak matching -- match on dinucleotides and it goes away."

The standard fix in this literature is to build negatives by dinucleotide-preserving
SHUFFLE. We tested that and abandoned it: TARDBP's motif survives the shuffle because the
motif IS a dinucleotide repeat, and shuffled sequence is detectable per se at AUROC ~0.69,
giving the arm a floor that cannot be subtracted away.

So instead of trying to build better negatives, this measures how compositionally unusual
a protein's binding sites are relative to the sequence around them:

    extremity = median L1(positive, its matched negative)
              - median L1(negative, a random other negative)

READ THE NEXT PARAGRAPH BEFORE USING THIS. The second term was originally described as
"the floor" -- the best any real-sequence matching procedure could achieve -- and that
identification is WRONG. A random pairing of windows is not the best achievable match; a
targeted nearest-neighbour search does far better. Measured on 21 datasets: random pairing
gives median L1 0.540, GC matching gives 0.500, and a directed dinucleotide search gives
**0.220**, beating the supposed "floor" on 20 of 21 datasets. Every conclusion that treated
the random-pair distance as a floor was therefore void, including a claim that 81% of the
panel was already optimally matched and that dinucleotide matching could not improve on GC
matching. It can, and it does: composition-only AUROC falls from 0.793 to 0.609 under
dinucleotide-matched negatives.

WHAT EXTREMITY IS STILL GOOD FOR. It remains a useful RELATIVE measure of how far a
protein's binding sites sit from ordinary sequence in its own neighbourhood, and the
empirical relationships are unaffected because they were measured, not derived. Across 187
datasets it correlates +0.697 with a model's eventual AUROC, +0.755 with what a
composition-only baseline achieves, and -0.372 with how much the model adds beyond
composition, and it replicates across independent cell lines at +0.857 -- more stably than
any model-derived quantity here. The categorical split (extreme vs not) separates cleanly in
both cell lines: extreme datasets keep composition-only AUROC at 0.73 after dinucleotide
matching against 0.60 for the rest.

It does NOT predict how much an individual dataset resists matching. An earlier version of
this docstring claimed +0.243 for that, measured on 21 datasets; the full 99-dataset K562
panel gives -0.119, i.e. flipped sign and weak either way. Withdrawn.

WHAT IT IS NOT. It is not a bound on achievable match quality, and it is not a substitute
for actually building matched negatives. Establishing how much of a model's performance is
composition requires the matched arm; extremity only tells you where to expect trouble.
"""

import numpy as np

ALPHABET = "ACGT"
DINUCS = [a + b for a in ALPHABET for b in ALPHABET]
_INDEX = {d: i for i, d in enumerate(DINUCS)}


def dinuc_freq(seq):
    """16-vector of dinucleotide frequencies.

    Pairs containing an unknown base are skipped rather than mapped to a real
    dinucleotide, so an N cannot inflate a count.
    """
    v = np.zeros(16, dtype=np.float64)
    s = seq.upper()
    n = 0
    for i in range(len(s) - 1):
        j = _INDEX.get(s[i:i + 2])
        if j is not None:
            v[j] += 1.0
            n += 1
    return v / n if n else v


def dinuc_matrix(seqs):
    """(n_seq x 16) frequency matrix."""
    return np.vstack([dinuc_freq(s) for s in seqs]) if len(seqs) else np.zeros((0, 16))


def l1(a, b):
    """Row-wise L1 distance between two frequency matrices.

    L1 on frequency vectors runs 0 to 2. A value of 0.16 means the average dinucleotide
    frequency differs by one percentage point.
    """
    return np.abs(np.asarray(a) - np.asarray(b)).sum(axis=1)


def random_pair_distance(neg_matrix, seed=0):
    """Typical L1 between two UNRELATED real windows from the same pool.

    A reference scale for "how different are two arbitrary windows here", not a bound on
    achievable matching. It was originally named `floor_distance` and treated as the best
    any matcher could do, which is wrong by a factor of more than two: a targeted search
    reaches 0.220 where random pairing gives 0.540. The old name is kept as an alias so
    existing callers do not break, but the concept is a reference, not a floor.
    """
    n = len(neg_matrix)
    if n < 2:
        return np.nan
    perm = np.random.default_rng(seed).permutation(n)
    # a fixed point maps a window to itself and contributes a spurious zero
    same = perm == np.arange(n)
    if same.any():
        perm[same] = (perm[same] + 1) % n
    return float(np.median(l1(neg_matrix, neg_matrix[perm])))


def extremity(pos_seqs, neg_seqs, seed=0):
    """Compositional extremity of one dataset, plus the parts it is built from."""
    P, N = dinuc_matrix(list(pos_seqs)), dinuc_matrix(list(neg_seqs))
    n = min(len(P), len(N))
    if n < 2:
        return {"n": n, "l1_pos_neg": np.nan, "l1_floor": np.nan,
                "extremity": np.nan, "matchable": None}
    pn = float(np.median(l1(P[:n], N[:n])))
    ref = random_pair_distance(N[:n], seed=seed)
    # `l1_floor` and `matchable` keep their old names so downstream tables stay readable,
    # but note what they now mean: `l1_floor` is the random-pair reference, not a bound,
    # and `matchable` means "no more different from its negatives than two arbitrary
    # windows are from each other" -- NOT "cannot be matched better".
    return {"n": n, "l1_pos_neg": round(pn, 5), "l1_floor": round(ref, 5),
            "extremity": round(pn - ref, 5), "matchable": bool(pn - ref <= 0)}


def from_dataset(df, seq_col="seq_dna", label_col="label", seed=0):
    """Extremity from a prepared dataset frame.

    Uses `seq_dna` rather than `seq_rna` by default. The RNA column is strand-corrected,
    so a minus-strand window has been reverse-complemented and its dinucleotide counts are
    those of the complement. Mixing strands in one composition comparison would blur the
    measurement; the reference-strand DNA is the consistent choice.
    """
    pos = df.loc[df[label_col] == 1, seq_col].tolist()
    neg = df.loc[df[label_col] == 0, seq_col].tolist()
    return extremity(pos, neg, seed=seed)


def panel(datasets, **kw):
    """Extremity for every dataset. `datasets` maps name -> dataframe."""
    import pandas as pd
    rows = []
    for name, df in datasets.items():
        rows.append({"dataset": name, **from_dataset(df, **kw)})
    return pd.DataFrame(rows)


# Retained so existing callers and tables keep working after the rename.
floor_distance = random_pair_distance
