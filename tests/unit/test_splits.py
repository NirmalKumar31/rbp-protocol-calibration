"""Tests for chromosome-level splitting and grouped k-fold folds.

The properties here are the ones whose violation is silent. A chromosome in two folds
inflates every score and looks like success; an unbalanced fold makes per-fold AUROCs
incomparable without any error being raised.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp.data import splits  # noqa: E402

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def counts_matrix(n_proteins=12, seed=3):
    """Peak counts that are uneven across chromosomes, as the real data is."""
    rng = np.random.default_rng(seed)
    scale = rng.integers(200, 9000, n_proteins)[:, None]
    shape = rng.dirichlet(np.linspace(0.6, 2.0, len(CHROMS)), n_proteins)
    return (scale * shape).astype(np.int64)


class TestFoldRoles:
    def test_test_fold_is_the_requested_one(self):
        assert splits.fold_roles(2, 5)[0] == 2

    def test_val_is_the_next_fold_cyclically(self):
        assert splits.fold_roles(4, 5)[1] == 0

    def test_train_excludes_test_and_val(self):
        test, val, train = splits.fold_roles(1, 5)
        assert test not in train and val not in train

    def test_all_folds_accounted_for_exactly_once(self):
        test, val, train = splits.fold_roles(0, 5)
        assert sorted([test, val] + train) == list(range(5))

    @pytest.mark.parametrize("k", [3, 4, 5, 10])
    def test_every_fold_is_test_exactly_once_over_k_iterations(self, k):
        tests = [splits.fold_roles(f, k)[0] for f in range(k)]
        assert sorted(tests) == list(range(k))

    @pytest.mark.parametrize("k", [3, 4, 5, 10])
    def test_every_fold_is_val_exactly_once_over_k_iterations(self, k):
        vals = [splits.fold_roles(f, k)[1] for f in range(k)]
        assert sorted(vals) == list(range(k))

    def test_rejects_out_of_range_fold(self):
        with pytest.raises(ValueError):
            splits.fold_roles(5, 5)

    def test_rejects_k_below_three(self):
        """k=2 leaves no fold for validation, so early stopping would peek at test."""
        with pytest.raises(ValueError):
            splits.fold_roles(0, 2)


class TestAssignFolds:
    def test_tags_every_row(self):
        fmap = {"chr1": 0, "chr2": 1}
        rows = [{"chrom": "chr1"}, {"chrom": "chr2"}, {"chrom": "chr1"}]
        splits.assign_folds(rows, fmap)
        assert [r["fold"] for r in rows] == [0, 1, 0]

    def test_same_chromosome_always_same_fold(self):
        fmap = {c: i % 5 for i, c in enumerate(CHROMS)}
        rows = [{"chrom": c} for c in CHROMS * 4]
        splits.assign_folds(rows, fmap)
        per = {}
        for r in rows:
            per.setdefault(r["chrom"], set()).add(r["fold"])
        assert all(len(v) == 1 for v in per.values())

    def test_missing_chromosome_raises(self):
        with pytest.raises(KeyError):
            splits.assign_folds([{"chrom": "chrZ"}], {"chr1": 0})


class TestSplitOfFold:
    def test_roles_are_consistent_with_fold_roles(self):
        k = 5
        for fold in range(k):
            test, val, train = splits.fold_roles(fold, k)
            assert splits.split_of_fold(test, fold, k) == "test"
            assert splits.split_of_fold(val, fold, k) == "val"
            for t in train:
                assert splits.split_of_fold(t, fold, k) == "train"

    def test_every_row_is_scored_as_test_exactly_once_across_folds(self):
        """The defining property of k-fold: all data is out-of-fold exactly once."""
        k = 5
        for row_fold in range(k):
            n = sum(1 for fold in range(k)
                    if splits.split_of_fold(row_fold, fold, k) == "test")
            assert n == 1


class TestCheckFolds:
    def test_clean_assignment_has_no_problems(self):
        rows = [{"chrom": c, "fold": i % 5} for i, c in enumerate(CHROMS)]
        assert splits.check_folds(rows, 5) == []

    def test_detects_chromosome_in_two_folds(self):
        rows = [{"chrom": "chr1", "fold": 0}, {"chrom": "chr1", "fold": 3}]
        assert any("chr1" in p for p in splits.check_folds(rows, 5))

    def test_detects_empty_fold(self):
        rows = [{"chrom": "chr1", "fold": 0}, {"chrom": "chr2", "fold": 1}]
        probs = splits.check_folds(rows, 5)
        assert sum("empty" in p for p in probs) == 3

    def test_reports_each_problem_once(self):
        rows = [{"chrom": "chr1", "fold": 0}] * 3 + [{"chrom": "chr1", "fold": 1}] * 3
        assert len([p for p in splits.check_folds(rows, 2 + 1) if "chr1" in p]) == 1


@pytest.fixture(scope="module")
def solved():
    counts = counts_matrix()
    loss, worst, assign = splits.optimize_folds(counts, k=5, restarts=8, iters=800)
    return counts, loss, worst, assign


class TestOptimizeFolds:
    def test_returns_one_fold_per_chromosome(self, solved):
        counts, _, _, assign = solved
        assert len(assign) == counts.shape[1]

    def test_uses_all_k_folds(self, solved):
        _, _, _, assign = solved
        assert set(assign.tolist()) == set(range(5))

    def test_respects_min_per_fold(self, solved):
        _, _, _, assign = solved
        assert min((assign == f).sum() for f in range(5)) >= 3

    def test_beats_a_naive_round_robin(self, solved):
        counts, loss, _, _ = solved
        naive = np.array([i % 5 for i in range(counts.shape[1])])
        naive_loss, _ = splits.assignment_loss(counts, naive, (0.2,) * 5)
        assert loss < naive_loss

    def test_worst_protein_deviation_is_small(self, solved):
        """Every protein, not just the average, must land near 20% per fold."""
        _, _, worst, _ = solved
        assert worst < 0.10

    def test_is_deterministic_for_a_fixed_seed(self):
        counts = counts_matrix()
        a = splits.optimize_folds(counts, k=5, restarts=4, iters=400, seed=11)[2]
        b = splits.optimize_folds(counts, k=5, restarts=4, iters=400, seed=11)[2]
        assert (a == b).all()

    def test_different_seeds_are_allowed_to_differ(self):
        counts = counts_matrix()
        a = splits.optimize_folds(counts, k=5, restarts=2, iters=200, seed=1)[0]
        b = splits.optimize_folds(counts, k=5, restarts=2, iters=200, seed=2)[0]
        assert a >= 0 and b >= 0

    @pytest.mark.parametrize("k", [3, 4, 5])
    def test_works_for_several_k(self, k):
        counts = counts_matrix()
        _, _, assign = splits.optimize_folds(counts, k=k, restarts=4, iters=400)
        assert set(assign.tolist()) == set(range(k))

    def test_refuses_impossible_request(self):
        """23 chromosomes cannot give 10 folds of at least 3."""
        with pytest.raises(ValueError):
            splits.optimize_folds(counts_matrix(), k=10, min_per_fold=3)


class TestThreeWayStillWorks:
    """The fixed split is retained; generalising the optimiser must not break it."""

    def test_default_target_gives_three_groups(self):
        counts = counts_matrix()
        _, _, assign = splits.optimize_assignment(counts, restarts=4, iters=400)
        assert set(assign.tolist()) == {0, 1, 2}

    def test_split_of_is_unchanged(self):
        cfg = {"test": ["chr1"], "val": ["chr9"]}
        assert splits.split_of("chr1", cfg) == "test"
        assert splits.split_of("chr9", cfg) == "val"
        assert splits.split_of("chr7", cfg) == "train"
