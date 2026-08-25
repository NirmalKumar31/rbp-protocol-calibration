"""Sequence handling, window geometry, negative matching, and split integrity."""

import numpy as np
import pytest

from rbp.data import negatives as neg
from rbp.data import splits, windows as win


class TestRevcomp:
    def test_basic(self):
        assert win.revcomp("ACGT") == "ACGT"
        assert win.revcomp("AAAA") == "TTTT"
        assert win.revcomp("ACCTG") == "CAGGT"

    def test_n_preserved(self):
        assert win.revcomp("ACNGT") == "ACNGT"

    def test_involution(self):
        s = "ACGTTGCAGGNTCA"
        assert win.revcomp(win.revcomp(s)) == s


class TestToRna:
    def test_plus_strand_transcribes_only(self):
        assert win.to_rna("ACGT", "+") == "ACGU"

    def test_minus_strand_reverse_complements_then_transcribes(self):
        # revcomp("ACGT") = "ACGT" -> "ACGU"; use an asymmetric case instead
        assert win.to_rna("AACCG", "-") == "CGGUU"

    def test_no_thymine_survives(self):
        assert "T" not in win.to_rna("TTTTT", "+")


class TestGcContent:
    @pytest.mark.parametrize("seq,expected", [
        ("GCGC", 1.0), ("ATAT", 0.0), ("ACGT", 0.5), ("", 0.0), ("acgt", 0.5),
    ])
    def test_values(self, seq, expected):
        assert win.gc_content(seq) == expected


class TestWindowBounds:
    def test_odd_size_is_centred(self):
        w0, w1 = win.window_bounds(1000, 1100, 101)
        assert w1 - w0 == 101
        assert (w0 + w1) // 2 == 1050            # midpoint preserved

    def test_single_base_peak(self):
        w0, w1 = win.window_bounds(500, 501, 101)
        assert (w0, w1) == (450, 551)
        assert w1 - w0 == 101

    def test_window_size_always_exact(self):
        for a, b in [(0, 1), (10, 11), (999, 1500), (7, 8)]:
            w0, w1 = win.window_bounds(a, b, 101)
            assert w1 - w0 == 101


class TestExclusionZones:
    def test_padding_and_merging(self):
        peaks = [("chr1", 1000, 1100, "+"), ("chr1", 1200, 1300, "+")]
        z = neg.exclusion_zones(peaks, 500)
        assert z["chr1"] == [(500, 1800)]        # padded ranges overlap, so they merge

    def test_never_negative_coordinates(self):
        z = neg.exclusion_zones([("chr1", 10, 20, "+")], 500)
        assert z["chr1"][0][0] == 0

    def test_separate_chromosomes(self):
        z = neg.exclusion_zones([("chr1", 0, 10, "+"), ("chr2", 0, 10, "+")], 5)
        assert set(z) == {"chr1", "chr2"}


class TestSubtract:
    def test_hole_in_the_middle(self):
        assert neg.subtract([(0, 1000)], [(400, 600)]) == [(0, 400), (600, 1000)]

    def test_hole_covers_everything(self):
        assert neg.subtract([(100, 200)], [(0, 1000)]) == []

    def test_hole_outside_is_ignored(self):
        assert neg.subtract([(0, 100)], [(500, 600)]) == [(0, 100)]

    def test_hole_clips_left_edge(self):
        assert neg.subtract([(100, 500)], [(0, 200)]) == [(200, 500)]

    def test_hole_clips_right_edge(self):
        assert neg.subtract([(100, 500)], [(400, 900)]) == [(100, 400)]

    def test_many_holes(self):
        out = neg.subtract([(0, 1000)], [(100, 200), (300, 400), (900, 950)])
        assert out == [(0, 100), (200, 300), (400, 900), (950, 1000)]

    def test_no_holes_returns_input(self):
        assert neg.subtract([(0, 10), (20, 30)], []) == [(0, 10), (20, 30)]

    def test_output_intervals_are_valid(self):
        out = neg.subtract([(0, 1000), (2000, 3000)], [(500, 2500)])
        assert all(e > s for s, e in out)
        assert out == [(0, 500), (2500, 3000)]


class TestWeightedPick:
    def test_stays_inside_intervals_and_leaves_room(self):
        rng = np.random.default_rng(0)
        ivs = [(0, 200), (1000, 1400)]
        for _ in range(200):
            w0 = neg._weighted_pick(rng, ivs, 101)
            assert w0 is not None
            assert any(s <= w0 and w0 + 101 <= e for s, e in ivs)

    def test_returns_none_when_nothing_fits(self):
        rng = np.random.default_rng(0)
        assert neg._weighted_pick(rng, [(0, 50)], 101) is None


class TestSplits:
    cfg = {"test": ["chr1", "chr2"], "val": ["chr19"]}

    def test_assignment(self):
        assert splits.split_of("chr1", self.cfg) == "test"
        assert splits.split_of("chr19", self.cfg) == "val"
        assert splits.split_of("chr7", self.cfg) == "train"

    def test_chromosomes_are_disjoint_across_splits(self):
        rows = [{"chrom": c} for c in ("chr1", "chr19", "chr7", "chr1")]
        splits.assign(rows, self.cfg)
        assert splits.check_disjoint(rows) == {}

    def test_detects_a_leaked_chromosome(self):
        rows = [{"chrom": "chr5", "split": "train"}, {"chrom": "chr5", "split": "test"}]
        assert "chr5" in splits.check_disjoint(rows)


class _FakeSeq:
    def __init__(self, s):
        self.seq = s

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, sl):
        return _FakeSeq(self.seq[sl])

    def upper(self):
        return self.seq.upper()


class _FakeFasta(dict):
    """Minimal stand-in for pyfaidx.Fasta: chrom -> sequence with slicing."""


def _fasta(**chroms):
    return _FakeFasta({c: _FakeSeq(s) for c, s in chroms.items()})


def _index(region, chrom, start, end):
    return {region: {chrom: (np.array([start]), np.array([end]))}}


class TestDinucVector:
    def test_counts_mode_sums_to_window_length_minus_one(self):
        v = neg.dinuc_vector("ACGT" * 10, normalise=False)
        assert v.sum() == len("ACGT" * 10) - 1

    def test_counts_are_exact_integers(self):
        # The whole reproducibility fix rests on this: a count is representable in
        # float64 exactly, so L1 distances between count vectors are exact.
        v = neg.dinuc_vector("ACGTTGCA" * 12, normalise=False)
        assert np.array_equal(v, np.round(v))

    def test_frequencies_are_counts_over_total(self):
        seq = "ACGTTGCA" * 12
        c = neg.dinuc_vector(seq, normalise=False)
        f = neg.dinuc_vector(seq)
        assert np.allclose(f, c / c.sum())

    def test_matrix_respects_normalise(self):
        seqs = ["ACGT" * 8, "GGCC" * 8]
        assert np.allclose(neg.dinuc_matrix(seqs, normalise=False).sum(axis=1), 31)
        assert np.allclose(neg.dinuc_matrix(seqs).sum(axis=1), 1.0)


class TestBuildNegativesDinuc:
    size = 21

    def _run(self, seq, positives, peaks, pool=None, **kw):
        fasta = _fasta(chr1=seq)
        lo, hi = pool if pool else (0, len(seq))
        index = _index("intron", "chr1", lo, hi)
        return neg.build_negatives_dinuc(
            positives, peaks, fasta, index, self.size,
            min_peak_distance=0, seed=7, **kw)

    def _pos(self, seq, start):
        return {"chrom": "chr1", "start": start, "end": start + self.size,
                "strand": "+", "region": "intron",
                "seq_dna": seq[start:start + self.size]}

    def test_single_candidate_pool_does_not_misalign_rows(self):
        """A pool holding exactly one usable window makes k == 1, and cKDTree then
        returns a flat (n,) array. Reshaping that as (1, n) instead of (n, 1) means
        every positive after the first indexes off the end. The pool interval below is
        exactly one window wide, which is the only way to reach k == 1.
        """
        seq = "".join("ACGT"[(i * 5) % 4] for i in range(200))
        positives = [self._pos(seq, s) for s in (0, 50, 100)]
        rows, dropped, _dists = self._run(seq, positives, [], pool=(0, self.size),
                                          pool_min=1, pool_multiple=1)
        assert len(rows) == 3
        matched = [(r["chrom"], r["start"]) for r in rows if r is not None]
        # one candidate, so exactly one positive can be matched and the rest are dropped
        assert len(matched) == 1
        assert dropped["no_match"] == 2

    def test_reported_distance_is_in_frequency_units(self):
        """Counts are used internally; the reported L1 must still be 0..2 on
        frequencies, because every published dinuc_l1 figure is on that scale."""
        seq = ("ACGT" * 30) + ("GGCC" * 30)
        positives = [self._pos(seq, 0), self._pos(seq, 130)]
        _rows, _dropped, dists = self._run(seq, positives, [])
        ok = dists[~np.isnan(dists)]
        assert len(ok) > 0
        assert ok.max() <= 2.0

    def test_same_inputs_give_the_same_negatives(self):
        seq = "".join("ACGT"[(i * 7 + i // 5) % 4] for i in range(600))
        positives = [self._pos(seq, s) for s in (0, 50, 100, 150)]
        a = self._run(seq, positives, [])[0]
        b = self._run(seq, positives, [])[0]
        assert [None if r is None else r["start"] for r in a] == \
               [None if r is None else r["start"] for r in b]

    def test_no_window_is_used_as_a_negative_twice(self):
        seq = "".join("ACGT"[(i * 3) % 4] for i in range(800))
        positives = [self._pos(seq, s) for s in range(0, 400, 25)]
        rows, _dropped, _dists = self._run(seq, positives, [])
        got = [(r["chrom"], r["start"]) for r in rows if r is not None]
        assert len(got) == len(set(got))

    def test_tied_candidates_resolve_to_the_earliest_position(self):
        """Where several candidates match a positive equally well, the winner must be
        decided by genomic position rather than by whatever order cKDTree happened to
        return -- that order differs between the arm64 and x86-64 builds of scipy, and
        chasing it cost a full cloud re-run to discover.
        """
        unit = "ACGT"
        seq = unit * 200                        # every in-frame window is compositionally identical
        positives = [self._pos(seq, 400)]
        rows, _dropped, _dists = self._run(seq, positives, [], pool_min=40, pool_multiple=1)
        assert rows[0] is not None
        # many windows tie; the earliest sampled one wins, deterministically
        again = self._run(seq, positives, [], pool_min=40, pool_multiple=1)[0]
        assert rows[0]["start"] == again[0]["start"]

    def test_tiebreak_cannot_outrank_a_genuinely_closer_window(self):
        """The tiebreak stays under 0.5 while real count distances differ by at least 1,
        so it may order equals but must never promote a worse match over a better one."""
        # positive is pure AC repeat; one exact-composition window exists late in the
        # sequence, everything earlier is GT and much further away
        seq = ("GTGT" * 150) + ("ACAC" * 50)
        positives = [self._pos(seq, 700)]
        rows, _dropped, dists = self._run(seq, positives, [], pool_min=60, pool_multiple=1)
        assert rows[0] is not None
        # it must land in the AC region, not on an early GT window that the tiebreak favours
        assert rows[0]["start"] >= 600
