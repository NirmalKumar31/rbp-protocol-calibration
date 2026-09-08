"""The external chromosome map must be one map, and the strand metric must use strand.

Why this exists. Two defects shipped in `scripts/external_sensitivity.py` and passed 1141
assertions, because every gate asserted a weaker property than the manuscript claimed.

  * `chrom_folds()` was called once per arm. It balances fold sizes by the counts it is handed;
    the arms share positives but differ in negatives, so the maps diverged. Measured across all
    135 datasets, they were identical in ZERO of them. The verifier asserted that each arm was
    chromosome-blocked, which was true of each arm separately, and never that the two maps
    agreed. So a fold-design difference was mixed into a comparison whose entire purpose is to
    hold the folds fixed and vary only the negative construction.
  * `horlacher_arm.windows()` read the BED6 strand to reverse-complement the sequence and then
    did not return it, so the caller filled every strand with "+" and published a
    strand-agnostic number under a same-strand label.

These tests are the ones that would have failed. They run on synthetic frames rather than the
real deposit, which is 135 datasets and tens of minutes, and which a reader does not have.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

es = pytest.importorskip("external_sensitivity")


def _arm(chroms, starts, labels, strands=None):
    n = len(chroms)
    return pd.DataFrame({
        "chrom": chroms, "start": starts, "label": labels,
        "strand": strands if strands is not None else ["+"] * n,
        "seq_rna": ["ACGU" * 5] * n, "fold": [i % 5 for i in range(n)]})


def _two_arms(shared_pos, neg_a, neg_b):
    """Two arms with an IDENTICAL positive set and different negatives, as the deposit has."""
    def build(neg):
        ch = [c for c, _ in shared_pos] + [c for c, _ in neg]
        st = [s for _, s in shared_pos] + [s for _, s in neg]
        lb = [1] * len(shared_pos) + [0] * len(neg)
        return _arm(ch, st, lb)
    return {"n1": build(neg_a), "n2": build(neg_b)}


# --- the map --------------------------------------------------------------------------------

def test_the_two_arms_get_a_byte_identical_chromosome_map():
    """The defect: 0 of 135 real datasets had identical maps."""
    pos = [(f"chr{1 + i % 6}", i * 1000) for i in range(60)]
    # Deliberately lopsided negatives, which is what made the per-arm maps diverge.
    neg_a = [("chr1", 900000 + i) for i in range(50)]
    neg_b = [("chr6", 900000 + i) for i in range(50)]
    arms = _two_arms(pos, neg_a, neg_b)
    cmap = es.shared_map(arms)
    assert cmap is not None
    for a in arms.values():
        assert set(a.chrom.unique()) <= set(cmap), "a chromosome has no fold"
    # one map object, so identity is structural; assert it anyway against each arm's view
    views = [{c: cmap[c] for c in a.chrom.unique()} for a in arms.values()]
    shared = set(views[0]) & set(views[1])
    assert shared, "the arms share no chromosome, so the test proves nothing"
    assert all(views[0][c] == views[1][c] for c in shared)


def test_the_map_does_not_depend_on_the_negatives_at_all():
    """Arm-invariance, stated directly: change only the negatives, get the same map."""
    pos = [(f"chr{1 + i % 6}", i * 1000) for i in range(60)]
    a = es.shared_map(_two_arms(pos, [("chr1", 1)], [("chr1", 2)]))
    b = es.shared_map(_two_arms(pos, [("chr2", 5)] * 40, [("chr3", 7)] * 90))
    common = set(a) & set(b)
    assert {c: a[c] for c in common} == {c: b[c] for c in common}


def test_the_map_does_not_depend_on_row_order():
    pos = [(f"chr{1 + i % 6}", i * 1000) for i in range(60)]
    arms = _two_arms(pos, [("chr1", 1)], [("chr1", 2)])
    straight = es.shared_map(arms)
    shuffled = es.shared_map({k: v.sample(frac=1.0, random_state=7).reset_index(drop=True)
                              for k, v in arms.items()})
    assert straight == shuffled, "the assignment depends on row order, so it is not a function"


def test_a_chromosome_seen_only_among_negatives_still_gets_a_fold():
    pos = [(f"chr{1 + i % 5}", i * 1000) for i in range(50)]
    arms = _two_arms(pos, [("chrX", 10)], [("chrY", 20)])
    cmap = es.shared_map(arms)
    assert "chrX" in cmap and "chrY" in cmap, "a negative-only chromosome was dropped"
    assert cmap["chrX"] == es.shared_map(_two_arms(pos, [("chrY", 20)], [("chrX", 10)]))["chrX"], \
        "the negative-only assignment depends on which arm is processed first"


def test_too_few_chromosomes_returns_none_rather_than_a_broken_partition():
    pos = [("chr1", i) for i in range(40)]
    assert es.shared_map(_two_arms(pos, [("chr1", 99)], [("chr1", 98)])) is None


def test_a_mismatched_positive_set_is_a_hard_failure():
    """The map's arm-invariance rests on the positives being identical. Verify, do not assume."""
    arms = _two_arms([(f"chr{1 + i % 6}", i) for i in range(60)], [("chr1", 1)], [("chr1", 2)])
    arms["n2"] = arms["n2"].drop(index=0).reset_index(drop=True)      # remove one positive
    with pytest.raises(SystemExit):
        es.shared_map(arms)


def test_every_chromosome_lands_in_exactly_one_fold():
    pos = [(f"chr{1 + i % 8}", i * 1000) for i in range(80)]
    arms = _two_arms(pos, [("chr2", 5)], [("chr7", 9)])
    cmap = es.shared_map(arms)
    for a in arms.values():
        folds = np.array([cmap[c] for c in a.chrom.values])
        assert es.one_fold_per_chrom(a.chrom.values, folds)


# --- the strand metric ------------------------------------------------------------------------

def test_the_leakage_metric_actually_uses_strand():
    """A synthetic case whose answer DIFFERS when strand is ignored, which is the whole point.

    Two windows 100 bp apart, on opposite strands, in different folds. Ignoring strand makes
    them neighbours that cross a fold. Respecting strand makes them not neighbours at all.
    """
    chrom = ["chr1", "chr1"]
    start = [1000, 1100]
    folds = [0, 1]
    have_s, cross_s = es.leakage(chrom, start, ["+", "-"], folds)
    have_i, cross_i = es.leakage(chrom, start, ["+", "+"], folds)
    assert (have_i, cross_i) == (2, 2), "the strand-blind case should see a crossing pair"
    assert have_s == 0, (
        "with opposite strands these are not same-strand neighbours; a metric returning the "
        "same answer either way is not using strand, which is exactly the shipped defect")


def test_same_strand_neighbours_in_one_fold_are_counted_but_not_crossing():
    have, cross = es.leakage(["chr1", "chr1"], [1000, 1100], ["+", "+"], [3, 3])
    assert (have, cross) == (2, 0)


def test_windows_further_than_a_kilobase_are_not_neighbours():
    have, cross = es.leakage(["chr1", "chr1"], [1000, 5000], ["+", "+"], [0, 1])
    assert (have, cross) == (0, 0)


def test_the_deposit_reader_returns_strand():
    """horlacher_arm.windows() read strand, used it, and did not return it."""
    ha = pytest.importorskip("horlacher_arm")
    import inspect
    src = inspect.getsource(ha.windows)
    assert '"strand": strand' in src, (
        "windows() no longer returns strand. external_sensitivity.py then has no strand to "
        "use, and a same-strand statistic without strand is a different statistic")
