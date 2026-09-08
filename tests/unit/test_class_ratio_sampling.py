"""The class-ratio subsample must keep the SAME positives in every arm, deterministically.

Why this exists. The first implementation drew with an RNG seeded on (dataset, arm, ratio) and
selected by ROW INDEX. For the ratios below one, which are reached by cutting positives, each
arm therefore kept a different random subset of positives, so the arm comparison varied both the
negative protocol and the positive realisation when the whole point is to vary only the first.
Correcting it moved the reported spans by up to 0.35 fold units.

The second implementation ranked on `chrom:start`, which is not a window: two records can share
a start and differ in end or strand, and 2 of 40 datasets sampled contain such a collision.

So the current rule is: rank on the full identity `chrom:start:end:strand`; require that
positives are unique on it, which is true in every committed arm; and give negatives a canonical
occurrence index, because duplicate negatives are real (7 in the dinucleotide arm, 300 in the
bias-aware arm, whose negatives come from a finite pool of other proteins' sites).

These tests pin all three properties on synthetic frames.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

cr = pytest.importorskip("class_ratio")

WIN = 101


def _frame(rows):
    """rows: list of (chrom, start, strand, label, fold)."""
    return pd.DataFrame({
        "chrom": [r[0] for r in rows], "start": [r[1] for r in rows],
        "end": [r[1] + WIN for r in rows], "strand": [r[2] for r in rows],
        "label": [r[3] for r in rows], "fold": [r[4] for r in rows]})


# 100 per fold, so the sparsest ratio (1:4) still leaves 25 positives per fold, above the
# MIN_PER_FOLD floor of 20. At 60 it returned None and the test proved nothing.
def _balanced(n_per_fold=100, folds=5, strand="+"):
    rows = []
    for f in range(folds):
        for i in range(n_per_fold):
            rows.append(("chr1", f * 10_000_000 + i * 200, strand, 1, f))
            rows.append(("chr2", f * 10_000_000 + i * 200, strand, 0, f))
    return _frame(rows)


def test_the_identity_is_the_full_window_not_just_the_start():
    d = _frame([("chr1", 100, "+", 1, 0), ("chr1", 100, "-", 1, 0)])
    ids = cr.window_ids(d)
    assert ids[0] != ids[1], (
        "two windows at the same start on opposite strands got one identity; chrom:start is "
        "coarser than a window and collides in real data")


def test_a_missing_identity_column_is_an_error():
    d = _frame([("chr1", 100, "+", 1, 0)]).drop(columns=["strand"])
    with pytest.raises(SystemExit):
        cr.window_ids(d)


def test_duplicate_positives_are_rejected_rather_than_ranked_ambiguously():
    """Positives are unique on the full identity in every committed arm. Assert it."""
    ids = np.array(["chr1:1:102:+", "chr1:1:102:+", "chr2:5:106:-"])
    with pytest.raises(SystemExit) as e:
        cr._rank(ids, "DS:CELL", "1:2", "pos", True)
    assert "duplicated positive" in str(e.value)


def test_duplicate_negatives_are_allowed_and_ordered_deterministically():
    """Duplicate negatives are real: 300 rows in the bias-aware arm. They must not crash, and
    the order must not depend on how the rows arrived."""
    ids = np.array(["chrA:1:102:+", "chrA:1:102:+", "chrB:2:103:-", "chrC:3:104:+"])
    a = cr._rank(ids, "DS:CELL", "1:2", "neg|dn", False)
    b = cr._rank(ids[::-1], "DS:CELL", "1:2", "neg|dn", False)
    # the SELECTED IDENTITIES must match, which is what matters; identical rows are
    # interchangeable, so the physical index may differ
    assert sorted(ids[a][:3]) == sorted(ids[::-1][b][:3])


def test_both_arms_keep_the_same_positives_when_positives_cut():
    """The confound the correction removes, stated directly."""
    shared = _balanced()
    arm_a = shared.copy()
    arm_b = shared.copy()
    # different negatives in each arm, the same positives
    arm_b.loc[arm_b.label == 0, "start"] = arm_b.loc[arm_b.label == 0, "start"] + 7
    arm_b.loc[arm_b.label == 0, "end"] = arm_b.loc[arm_b.label == 0, "start"] + WIN

    kept = {}
    for tag, d in (("a", arm_a), ("b", arm_b)):
        ids = cr.window_ids(d)
        m = cr.subsample(ids, d.label.values, d.fold.values, 0.25, "DS:CELL", tag, "1:4")
        assert m is not None
        kept[tag] = set(ids[m & (d.label.values == 1)])
    assert kept["a"] == kept["b"], (
        "the two arms retained different positives at 1:4, which is the confound this ranking "
        "exists to remove: the comparison would vary the positive sample as well as the "
        "negative protocol")


def test_shuffling_the_rows_does_not_change_which_windows_are_kept():
    d = _balanced()
    ids = cr.window_ids(d)
    m1 = cr.subsample(ids, d.label.values, d.fold.values, 0.5, "DS:CELL", "dn", "1:2")
    sh = d.sample(frac=1.0, random_state=11).reset_index(drop=True)
    ids2 = cr.window_ids(sh)
    m2 = cr.subsample(ids2, sh.label.values, sh.fold.values, 0.5, "DS:CELL", "dn", "1:2")
    assert set(ids[m1]) == set(ids2[m2]), "the subsample depends on row order"


def test_the_ratio_is_actually_reached_within_every_fold():
    d = _balanced()
    ids = cr.window_ids(d)
    m = cr.subsample(ids, d.label.values, d.fold.values, 0.5, "DS:CELL", "dn", "1:2")
    kept = d[m]
    for f in kept.fold.unique():
        s = kept[kept.fold == f]
        n1, n0 = int((s.label == 1).sum()), int((s.label == 0).sum())
        assert abs(n1 * 2 - n0) <= 1, f"fold {f} is {n1}:{n0}, not 1:2"


def test_a_fold_too_small_returns_none_rather_than_a_degenerate_split():
    d = _frame([("chr1", i * 200, "+", i % 2, 0) for i in range(10)])
    assert cr.subsample(cr.window_ids(d), d.label.values, d.fold.values,
                        0.25, "DS:CELL", "dn", "1:4") is None
