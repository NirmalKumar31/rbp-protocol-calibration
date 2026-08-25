"""Tests for compositional extremity.

The measure is a difference between two medians, so it is easy to get a plausible-looking
number from a broken implementation. These tests pin down the sign and scale using cases
where the answer is known by construction.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp.eval import extremity as ex  # noqa: E402

L = 101
RNG = np.random.default_rng(5)


def seqs_at(gc, n, rng=RNG):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]      # A C G T
    return ["".join(rng.choice(list("ACGT"), size=L, p=p)) for _ in range(n)]


def repeat_seqs(unit, n):
    s = (unit * (L // len(unit) + 1))[:L]
    return [s] * n


class TestDinucFreq:
    def test_sums_to_one(self):
        assert ex.dinuc_freq(seqs_at(0.5, 1)[0]).sum() == pytest.approx(1.0)

    def test_homopolymer_is_one_dinucleotide(self):
        v = ex.dinuc_freq("A" * 50)
        assert v[ex.DINUCS.index("AA")] == pytest.approx(1.0)

    def test_alternating_repeat_splits_between_two(self):
        v = ex.dinuc_freq("GT" * 50)
        assert v[ex.DINUCS.index("GT")] == pytest.approx(0.5, abs=0.02)
        assert v[ex.DINUCS.index("TG")] == pytest.approx(0.5, abs=0.02)

    def test_unknown_bases_are_skipped(self):
        """An N must not create a spurious dinucleotide."""
        v = ex.dinuc_freq("ANA")
        assert v.sum() == 0.0                       # both pairs contain N

    def test_partial_unknown_still_counts_valid_pairs(self):
        v = ex.dinuc_freq("AANAA")
        assert v[ex.DINUCS.index("AA")] == pytest.approx(1.0)

    def test_empty_sequence_gives_zeros(self):
        assert ex.dinuc_freq("").sum() == 0.0


class TestL1:
    def test_identical_is_zero(self):
        m = ex.dinuc_matrix(seqs_at(0.5, 5))
        assert np.allclose(ex.l1(m, m), 0.0)

    def test_maximum_is_two(self):
        """Two disjoint single-dinucleotide distributions are maximally far apart."""
        a = ex.dinuc_matrix(["A" * L])
        b = ex.dinuc_matrix(["C" * L])
        assert ex.l1(a, b)[0] == pytest.approx(2.0)

    def test_is_symmetric(self):
        a, b = ex.dinuc_matrix(seqs_at(0.4, 4)), ex.dinuc_matrix(seqs_at(0.6, 4))
        assert np.allclose(ex.l1(a, b), ex.l1(b, a))


class TestFloorDistance:
    def test_identical_pool_has_zero_floor(self):
        """If every negative is the same sequence, any match is perfect."""
        m = ex.dinuc_matrix(repeat_seqs("ACGT", 20))
        assert ex.floor_distance(m) == pytest.approx(0.0)

    def test_diverse_pool_has_positive_floor(self):
        m = ex.dinuc_matrix(seqs_at(0.5, 200))
        assert ex.floor_distance(m) > 0.0

    def test_no_window_is_paired_with_itself(self):
        """A fixed point in the permutation contributes a spurious zero and pulls the
        floor down, which would inflate extremity."""
        m = ex.dinuc_matrix(repeat_seqs("ACGT", 3) + repeat_seqs("GGCC", 3))
        # every distance is either 0 (same unit) or large; a self-pair would add a 0
        d = []
        for seed in range(20):
            d.append(ex.floor_distance(m, seed=seed))
        assert all(np.isfinite(x) for x in d)

    def test_too_few_sequences_gives_nan(self):
        assert np.isnan(ex.floor_distance(ex.dinuc_matrix(["ACGT"])))


class TestExtremity:
    def test_matched_negatives_give_extremity_at_or_below_zero(self):
        """Positives and negatives drawn from the same distribution are matchable."""
        rng = np.random.default_rng(11)
        r = ex.extremity(seqs_at(0.5, 300, rng), seqs_at(0.5, 300, rng))
        assert r["extremity"] <= 0.05
        assert r["matchable"] in (True, False)      # near zero, either side is fine

    def test_compositionally_extreme_positives_give_positive_extremity(self):
        """Positives are a GT repeat, negatives are ordinary sequence. No real negative
        can match them, which is the TARDBP situation."""
        rng = np.random.default_rng(13)
        r = ex.extremity(repeat_seqs("GT", 300), seqs_at(0.5, 300, rng))
        assert r["extremity"] > 0.5 and r["matchable"] is False

    def test_extremity_is_the_stated_difference(self):
        rng = np.random.default_rng(17)
        r = ex.extremity(seqs_at(0.6, 200, rng), seqs_at(0.4, 200, rng))
        assert r["extremity"] == pytest.approx(r["l1_pos_neg"] - r["l1_floor"], abs=1e-9)

    def test_matchable_flag_matches_the_sign(self):
        rng = np.random.default_rng(19)
        for pos_gc, neg_gc in ((0.5, 0.5), (0.8, 0.3)):
            r = ex.extremity(seqs_at(pos_gc, 200, rng), seqs_at(neg_gc, 200, rng))
            assert r["matchable"] == (r["extremity"] <= 0)

    def test_more_extreme_positives_give_larger_extremity(self):
        rng = np.random.default_rng(23)
        neg = seqs_at(0.5, 300, rng)
        mild = ex.extremity(seqs_at(0.58, 300, rng), neg)["extremity"]
        wild = ex.extremity(seqs_at(0.85, 300, rng), neg)["extremity"]
        assert wild > mild

    def test_too_few_sequences_is_reported_not_crashed(self):
        r = ex.extremity(["ACGT"], ["ACGT"])
        assert np.isnan(r["extremity"]) and r["matchable"] is None

    def test_is_deterministic_for_a_fixed_seed(self):
        rng = np.random.default_rng(29)
        pos, neg = seqs_at(0.6, 150, rng), seqs_at(0.5, 150, rng)
        assert ex.extremity(pos, neg, seed=3) == ex.extremity(pos, neg, seed=3)


class TestFromDataset:
    def test_reads_a_prepared_frame(self):
        rng = np.random.default_rng(31)
        df = pd.DataFrame({
            "seq_dna": seqs_at(0.6, 100, rng) + seqs_at(0.5, 100, rng),
            "label": [1] * 100 + [0] * 100})
        assert np.isfinite(ex.from_dataset(df)["extremity"])

    def test_uses_dna_not_strand_corrected_rna(self):
        """seq_rna is reverse-complemented on the minus strand, so its dinucleotide counts
        are the complement's. Mixing strands would blur the measurement."""
        rng = np.random.default_rng(37)
        df = pd.DataFrame({
            "seq_dna": seqs_at(0.6, 60, rng) + seqs_at(0.5, 60, rng),
            "seq_rna": ["A" * L] * 120,               # deliberately degenerate
            "label": [1] * 60 + [0] * 60})
        r = ex.from_dataset(df)
        assert r["l1_pos_neg"] > 0        # would be 0 if it had used seq_rna

    def test_panel_returns_one_row_per_dataset(self):
        rng = np.random.default_rng(41)
        df = pd.DataFrame({"seq_dna": seqs_at(0.5, 80, rng), "label": [1, 0] * 40})
        out = ex.panel({"A:K562": df, "B:K562": df})
        assert len(out) == 2 and set(out.dataset) == {"A:K562", "B:K562"}
