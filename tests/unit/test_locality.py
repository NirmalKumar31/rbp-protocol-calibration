"""Tests for the repeat measures and the data-driven locality control.

The locality control decides whether we trust a 0.10 AUROC drop as confound removal, so a
version of it that reports a large effect for a composition-only model would quietly
validate the wrong conclusion. The constructed cases below pin that down.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp.eval import locality as loc  # noqa: E402
from rbp.eval import repeats as rp  # noqa: E402

L = 101
RNG = np.random.default_rng(3)


def rand_seqs(n, rng=RNG, gc=0.5):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return ["".join(rng.choice(list("ACGU"), size=L, p=p)) for _ in range(n)]


def implant(seqs, motif, rng=RNG):
    out = []
    for s in seqs:
        i = int(rng.integers(0, L - len(motif)))
        out.append(s[:i] + motif + s[i + len(motif):])
    return out


class TestMinimalPeriod:
    @pytest.mark.parametrize("s,want", [
        ("UUUUU", 1), ("UGUGU", 2), ("UGUUGU", 3), ("A", 1), ("AC", 2),
        # GCAUG has G at index 0 and 4, so it IS periodic at p=4. My first version of
        # this test asserted 5 and was simply wrong about the string.
        ("GCAUG", 4), ("GCAUC", 5),
    ])
    def test_known_values(self, s, want):
        assert rp.minimal_period(s) == want

    def test_period_never_exceeds_length(self):
        assert rp.minimal_period("ACGUACGUA") <= 9


class TestRepeatScore:
    def test_homopolymer_is_one(self):
        assert rp.repeat_score("UUUUU") == pytest.approx(1.0)

    def test_dinucleotide_repeat_is_three_quarters_at_k5(self):
        assert rp.repeat_score("UGUGU") == pytest.approx(0.75)

    def test_non_repeat_is_zero(self):
        assert rp.repeat_score("GCAUC") == pytest.approx(0.0)

    def test_accidental_periodicity_is_reported(self):
        """GCAUG repeats G at distance 4, so it scores 0.25 rather than 0. The measure is
        structural and does not know which repeats are biologically meaningful."""
        assert rp.repeat_score("GCAUG") == pytest.approx(0.25)

    def test_single_base_is_zero(self):
        assert rp.repeat_score("A") == 0.0

    def test_bounded_zero_to_one(self):
        for s in ("AAAA", "ACAC", "ACGU", "ACGUA"):
            assert 0.0 <= rp.repeat_score(s) <= 1.0


class TestLowComplexity:
    def test_homopolymer_is_fully_covered(self):
        assert rp.low_complexity("UUUUU") == pytest.approx(1.0)

    def test_alternating_has_no_runs(self):
        assert rp.low_complexity("UGUGU") == pytest.approx(0.0)

    def test_catches_a_run_the_period_measure_misses(self):
        """AUUUA is not periodic but is obviously low-complexity."""
        assert rp.repeat_score("AUUUA") < 0.3 and rp.low_complexity("AUUUA") > 0.5

    def test_runs_shorter_than_three_do_not_count(self):
        assert rp.low_complexity("AACGU") == pytest.approx(0.0)


class TestKmerEnrichment:
    def test_implanted_motif_is_the_top_kmer(self):
        rng = np.random.default_rng(11)
        pos = implant(rand_seqs(300, rng), "GCAUG", rng)
        km, enr = loc.top_kmer(pos, rand_seqs(300, rng), k=5)
        assert km == "GCAUG" and enr > 1.0

    def test_no_signal_gives_small_enrichment(self):
        rng = np.random.default_rng(13)
        _, enr = loc.top_kmer(rand_seqs(300, rng), rand_seqs(300, rng), k=5)
        assert enr < 1.5

    def test_enrichment_is_symmetric_under_swap(self):
        rng = np.random.default_rng(17)
        a, b = rand_seqs(200, rng), implant(rand_seqs(200, rng), "ACUAA", rng)
        e1 = rp.kmer_enrichment(a, b, k=5)
        e2 = rp.kmer_enrichment(b, a, k=5)
        assert np.allclose(e1, -e2, atol=1e-9)


class TestBuildPairs:
    def test_both_mutants_change_exactly_one_base(self):
        rng = np.random.default_rng(19)
        seqs = implant(rand_seqs(50, rng), "GCAUG", rng)
        for p in loc.build_pairs(seqs, "GCAUG", rng=rng):
            assert sum(a != b for a, b in zip(p["ref"], p["disruptive"])) == 1
            assert sum(a != b for a, b in zip(p["ref"], p["neutral"])) == 1

    def test_both_mutants_make_the_same_substitution(self):
        """Otherwise the comparison confounds position with which base changed."""
        rng = np.random.default_rng(23)
        seqs = implant(rand_seqs(50, rng), "GCAUG", rng)
        for p in loc.build_pairs(seqs, "GCAUG", rng=rng):
            d = [(a, b) for a, b in zip(p["ref"], p["disruptive"]) if a != b][0]
            n = [(a, b) for a, b in zip(p["ref"], p["neutral"]) if a != b][0]
            assert d == n

    def test_disruptive_mutation_lands_inside_the_kmer(self):
        rng = np.random.default_rng(29)
        seqs = implant(rand_seqs(40, rng), "GCAUG", rng)
        for p in loc.build_pairs(seqs, "GCAUG", rng=rng):
            pos = [i for i, (a, b) in enumerate(zip(p["ref"], p["disruptive"])) if a != b][0]
            start = p["ref"].find("GCAUG")
            assert start <= pos < start + 5

    def test_neutral_mutation_is_far_from_every_occurrence(self):
        rng = np.random.default_rng(31)
        seqs = implant(rand_seqs(40, rng), "GCAUG", rng)
        for p in loc.build_pairs(seqs, "GCAUG", min_distance=25, rng=rng):
            pos = [i for i, (a, b) in enumerate(zip(p["ref"], p["neutral"])) if a != b][0]
            for h in loc._occurrences(p["ref"], "GCAUG"):
                assert abs(pos - (h + 2)) >= 25

    def test_windows_without_the_kmer_are_skipped(self):
        rng = np.random.default_rng(37)
        assert loc.build_pairs(rand_seqs(20, rng), "GGGGGGGG", rng=rng) == []


class TestLocality:
    def _frame(self, pos, neg, n_folds=5, rng=RNG):
        seqs = pos + neg
        y = [1] * len(pos) + [0] * len(neg)
        folds = list(rng.integers(0, n_folds, len(seqs)))
        return pd.DataFrame({"seq_rna": seqs, "label": y, "fold": folds})

    def test_detects_a_real_local_motif(self):
        rng = np.random.default_rng(41)
        df = self._frame(implant(rand_seqs(400, rng), "GCAUG", rng),
                         rand_seqs(400, rng), rng=rng)
        r = loc.locality(df, k=4, kmer_k=5)
        assert r is not None and r["cohens_d"] > 0.5

    def test_composition_only_signal_ALSO_scores_high_which_is_the_known_flaw(self):
        """THE important test, and it documents a FAILURE rather than a guarantee.

        Positives differ from negatives only in global composition -- there is no local
        feature anywhere. A valid locality probe should return a small effect. This one
        returns about 1.8.

        The reason is structural and not fixable here: a bag-of-k-mers model IS local by
        construction. It represents "composition" as weights on k-mers, so mutating the
        most enriched k-mer always moves the score more than mutating elsewhere, whether
        the underlying signal is a motif or global composition. And the disruptive site is
        SELECTED as the most enriched k-mer, i.e. exactly where the model's weight is
        concentrated -- so the comparison is biased by construction.

        Consequence: absolute locality cannot be read as "this model uses a motif", and the
        187-dataset stratification built on it in scripts/correction_quality.py is not
        interpretable. Answering that question needs a model that could in principle be
        NON-local -- a CNN or a transformer -- which is one of the concrete reasons to run
        the trained models rather than assume their results.
        """
        rng = np.random.default_rng(43)
        df = self._frame(rand_seqs(400, rng, gc=0.68), rand_seqs(400, rng, gc=0.32),
                         rng=rng)
        r = loc.locality(df, k=4, kmer_k=5)
        if r is None:
            pytest.skip("no usable pairs")
        assert r["cohens_d"] > 0.5, (
            "if this ever fails, the flaw documented here has gone away and the "
            "stratification may become usable -- re-derive it before trusting it")

    def test_returns_none_when_nothing_is_usable(self):
        rng = np.random.default_rng(47)
        df = self._frame(rand_seqs(12, rng), rand_seqs(12, rng), rng=rng)
        assert loc.locality(df, k=4, kmer_k=5) is None

    def test_reports_the_kmer_it_used(self):
        rng = np.random.default_rng(53)
        df = self._frame(implant(rand_seqs(400, rng), "ACUAA", rng),
                         rand_seqs(400, rng), rng=rng)
        r = loc.locality(df, k=4, kmer_k=5)
        assert r["kmers"] and all(len(km) == 5 for km in r["kmers"])

    def test_is_deterministic_for_a_fixed_seed(self):
        rng = np.random.default_rng(59)
        df = self._frame(implant(rand_seqs(300, rng), "GCAUG", rng),
                         rand_seqs(300, rng), rng=rng)
        a = loc.locality(df, k=4, kmer_k=5, seed=5)
        b = loc.locality(df, k=4, kmer_k=5, seed=5)
        assert a["cohens_d"] == pytest.approx(b["cohens_d"])

    def test_disruptive_change_exceeds_neutral_for_a_real_motif(self):
        rng = np.random.default_rng(61)
        df = self._frame(implant(rand_seqs(400, rng), "GCAUG", rng),
                         rand_seqs(400, rng), rng=rng)
        r = loc.locality(df, k=4, kmer_k=5)
        assert r["delta_disruptive"] > r["delta_neutral"]
