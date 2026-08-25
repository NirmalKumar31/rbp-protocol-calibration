"""Constructed controls for the ISM locality probe.

THE POINT OF THIS FILE. The previous probe (locality.py) validated at r = +0.96 against a
literature positive control on nine proteins and was still invalid: a constructed
pure-global-composition signal, with no local feature anywhere, scored Cohen's d about 1.8
where a valid probe must score ~0. Agreement with ground truth on a small sample proved
nothing, because both were driven by the same third thing -- how strong the protein's signal
was at all.

So this probe is judged only on cases where the answer is known BY CONSTRUCTION:

    pure global composition, no motif  ->  sensitivity profile must be FLAT   (low Gini)
    implanted motif                    ->  profile must be SPIKY             (high Gini)
    and the gap between them must be large enough to separate real datasets

The old probe fails the first of those. The test that pins it lives in test_locality.py and
asserts the failure so it cannot be forgotten. Here we require the opposite.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp.eval import locality_ism as loc  # noqa: E402

BASES = "ACGU"


def rand_seqs(n, rng, length=60, gc=0.5):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]      # A C G U
    return ["".join(rng.choice(list(BASES), length, p=p)) for _ in range(n)]


def implant(seqs, motif, rng):
    out = []
    for s in seqs:
        i = int(rng.integers(0, len(s) - len(motif)))
        out.append(s[:i] + motif + s[i + len(motif):])
    return out


# --- score functions with KNOWN locality -------------------------------------------------

def gc_score(seqs):
    """Pure global composition. Every position contributes equally; no local feature."""
    return np.array([(s.count("G") + s.count("C")) / len(s) for s in seqs])


def motif_score(motif):
    """Pure local feature: presence of one short motif, position-independent."""
    def f(seqs):
        return np.array([float(motif in s) for s in seqs])
    return f


def motif_count_score(motif):
    """Local, but graded -- so the probe is not only tested on a step function."""
    def f(seqs):
        return np.array([float(s.count(motif)) for s in seqs])
    return f


class TestGini:
    def test_flat_is_zero(self):
        assert loc.gini(np.ones(50)) == pytest.approx(0.0, abs=0.02)

    def test_single_spike_is_near_one(self):
        x = np.zeros(50)
        x[7] = 1.0
        assert loc.gini(x) > 0.95

    def test_all_zero_is_zero_not_nan(self):
        """A model that ignores the sequence must not produce a NaN Gini."""
        assert loc.gini(np.zeros(30)) == 0.0

    def test_scale_free(self):
        """Locality is a question about SHAPE, so scaling every value must not move it."""
        rng = np.random.default_rng(0)
        x = rng.random(40)
        assert loc.gini(x) == pytest.approx(loc.gini(x * 1000.0))


class TestIsmProfile:
    def test_length_matches_sequence(self):
        p = loc.ism_profile(gc_score, "ACGUACGUAC")
        assert p.shape == (10,)

    def test_motif_positions_dominate(self):
        seq = "AAAAAAAAAAGCAUGAAAAAAAAAA"
        p = loc.ism_profile(motif_score("GCAUG"), seq)
        inside = p[10:15].sum()
        outside = p.sum() - inside
        assert inside > 0 and outside == pytest.approx(0.0, abs=1e-9)

    def test_composition_profile_is_flat(self):
        """The control that matters. Under a GC model every position moves the score by a
        comparable amount, because the mean over alternative bases removes base identity."""
        rng = np.random.default_rng(1)
        seq = rand_seqs(1, rng, length=40)[0]
        p = loc.ism_profile(gc_score, seq)
        assert p.std() / p.mean() < 0.30

    def test_a_sequence_independent_model_gives_zeros(self):
        p = loc.ism_profile(lambda s: np.zeros(len(s)), "ACGUACGU")
        assert p.sum() == 0.0


class TestTheControlThatKilledTheOldProbe:
    """A pure global-composition signal MUST look non-local. This is the whole test."""

    def test_composition_only_signal_is_FLAT(self):
        rng = np.random.default_rng(7)
        r = loc.locality(gc_score, rand_seqs(12, rng, length=50), max_windows=12)
        assert r is not None
        assert r["gini"] < 0.25, (
            f"a pure composition signal scored gini={r['gini']:.3f}; the probe is measuring "
            "something other than locality and must not be used")

    def test_motif_signal_is_SPIKY(self):
        rng = np.random.default_rng(11)
        seqs = implant(rand_seqs(12, rng, length=50), "GCAUG", rng)
        r = loc.locality(motif_score("GCAUG"), seqs, max_windows=12)
        assert r["gini"] > 0.80

    def test_the_two_are_far_apart(self):
        """A probe that separates them by a hair cannot classify real datasets."""
        rng = np.random.default_rng(13)
        flat = loc.locality(gc_score, rand_seqs(12, rng, length=50), max_windows=12)
        seqs = implant(rand_seqs(12, rng, length=50), "GCAUG", rng)
        spiky = loc.locality(motif_score("GCAUG"), seqs, max_windows=12)
        assert spiky["gini"] - flat["gini"] > 0.55

    def test_graded_local_signal_also_registers(self):
        rng = np.random.default_rng(17)
        seqs = implant(implant(rand_seqs(12, rng, length=60), "GCAUG", rng), "GCAUG", rng)
        r = loc.locality(motif_count_score("GCAUG"), seqs, max_windows=12)
        assert r["gini"] > 0.60


class TestMixedSignal:
    def test_motif_plus_composition_lands_between(self):
        """Real models are mixtures, so the probe has to be graded, not binary."""
        rng = np.random.default_rng(19)
        seqs = implant(rand_seqs(12, rng, length=50), "GCAUG", rng)

        def mixed(s):
            return 0.5 * gc_score(s) + 0.5 * motif_score("GCAUG")(s)

        mid = loc.locality(mixed, seqs, max_windows=12)["gini"]
        lo = loc.locality(gc_score, seqs, max_windows=12)["gini"]
        hi = loc.locality(motif_score("GCAUG"), seqs, max_windows=12)["gini"]
        assert lo < mid < hi


class TestReporting:
    def test_top10_frac_is_a_fraction(self):
        rng = np.random.default_rng(23)
        r = loc.locality(gc_score, rand_seqs(8, rng, length=40), max_windows=8)
        assert 0.0 <= r["top10_frac"] <= 1.0

    def test_top10_frac_is_high_for_a_motif(self):
        rng = np.random.default_rng(29)
        seqs = implant(rand_seqs(8, rng, length=50), "GCAUG", rng)
        r = loc.locality(motif_score("GCAUG"), seqs, max_windows=8)
        assert r["top10_frac"] > 0.7

    def test_deterministic_for_a_fixed_seed(self):
        rng = np.random.default_rng(31)
        seqs = rand_seqs(20, rng, length=40)
        a = loc.locality(gc_score, seqs, max_windows=8, seed=3)
        b = loc.locality(gc_score, seqs, max_windows=8, seed=3)
        assert a["gini"] == pytest.approx(b["gini"])

    def test_none_when_there_is_nothing_to_score(self):
        assert loc.locality(gc_score, [], max_windows=8) is None

    def test_none_when_the_model_ignores_sequence(self):
        rng = np.random.default_rng(37)
        r = loc.locality(lambda s: np.zeros(len(s)), rand_seqs(5, rng, length=30),
                         max_windows=5)
        assert r is None
