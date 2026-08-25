"""The conservation-controlled test, FDR correction, and the positive control.

These test the statistical machinery that the project's central claim rests on, so the
emphasis is on: does it recover an effect that is really there, does it stay quiet when
nothing is there, and does it survive the degenerate cases (separation, one class, few
events) that our restricted variant sets actually produce.
"""

import numpy as np
import pandas as pd
import pytest

from rbp.variants import conservation as cons
from rbp.variants import phylop
from rbp.variants import positive_control as pc


class TestStandardise:
    def test_mean_zero_sd_one(self):
        z = cons._z([1, 2, 3, 4, 5])
        assert z.mean() == pytest.approx(0.0)
        assert z.std() == pytest.approx(1.0)

    def test_constant_input_does_not_divide_by_zero(self):
        # a model producing identical delta everywhere must not blow up the fit
        assert np.all(cons._z([3, 3, 3]) == 0)


class TestFirth:
    def test_stays_finite_under_separation(self):
        """Maximum likelihood diverges here; Firth is the reason we use it."""
        X = np.array([[-2.], [-1.], [-0.5], [0.5], [1.], [2.]])
        y = np.array([0, 0, 0, 1, 1, 1])
        beta, se = cons.firth_fit(X, y)
        assert np.all(np.isfinite(beta))
        assert beta[1] > 0                      # correct direction
        assert np.all(se > 0)

    def test_recovers_a_known_coefficient(self):
        rng = np.random.default_rng(0)
        n = 4000
        x = rng.standard_normal(n)
        y = rng.binomial(1, 1 / (1 + np.exp(-(0.0 + 1.5 * x))))
        beta, _ = cons.firth_fit(x.reshape(-1, 1), y)
        assert beta[1] == pytest.approx(1.5, abs=0.2)

    def test_intercept_is_first(self):
        rng = np.random.default_rng(1)
        n = 2000
        x = rng.standard_normal(n)
        y = rng.binomial(1, 1 / (1 + np.exp(-(-2.0 + 0.0 * x))))   # rare events
        beta, _ = cons.firth_fit(x.reshape(-1, 1), y)
        assert beta[0] < -1.0                   # intercept reflects low prevalence
        assert abs(beta[1]) < 0.3               # no slope was injected


class TestBenjaminiHochberg:
    def test_monotone_and_no_smaller_than_p(self):
        p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212])
        q = cons.benjamini_hochberg(p)
        assert np.all(np.diff(q) >= -1e-12)     # non-decreasing
        assert np.all(q >= p - 1e-12)           # adjustment never lowers a p-value
        assert np.all(q <= 1.0)

    def test_smallest_p_scales_by_m(self):
        p = np.array([0.01, 0.5, 0.5, 0.5, 0.5])
        q = cons.benjamini_hochberg(p)
        assert q[0] == pytest.approx(0.05)      # 0.01 * 5 / 1

    def test_order_is_preserved(self):
        p = np.array([0.5, 0.01, 0.2])
        q = cons.benjamini_hochberg(p)
        assert np.argmin(q) == 1                # the smallest p keeps its position

    def test_nan_passes_through(self):
        q = cons.benjamini_hochberg([0.01, np.nan, 0.02])
        assert np.isnan(q[1])
        assert np.all(np.isfinite(q[[0, 2]]))

    def test_all_nan(self):
        assert np.all(np.isnan(cons.benjamini_hochberg([np.nan, np.nan])))


class TestFitDeltaCoef:
    def _data(self, n=800, beta=1.2, gamma=1.0, r=0.0, seed=0):
        """delta correlated with conservation at r; label driven by both.

        delta is generated NON-NEGATIVE, because the real delta is a magnitude and the
        fit standardises |delta|. Drawing a signed normal and injecting the effect on
        the signed value would break the relationship the fit is looking for -- the
        mistake that was also present in the power simulation.
        """
        rng = np.random.default_rng(seed)
        cons_x = rng.standard_normal(n)
        latent = r * cons_x + np.sqrt(max(0.0, 1 - r * r)) * rng.standard_normal(n)
        delta = np.abs(latent)
        logit = -0.5 + beta * cons._z(delta) + gamma * cons._z(cons_x)
        y = rng.binomial(1, 1 / (1 + np.exp(-logit)))
        return delta, y, cons_x

    def test_recovers_a_real_effect(self):
        d, y, c = self._data(beta=1.5, r=0.0)
        f = cons.fit_delta_coef(d, y, c, n_boot=300, method="firth")
        assert f.survives
        assert f.ci_low > 0 and f.coef > 0

    def test_stays_quiet_when_there_is_no_effect(self):
        d, y, c = self._data(beta=0.0, r=0.0, seed=3)
        f = cons.fit_delta_coef(d, y, c, n_boot=300, method="firth")
        assert not f.survives                   # interval must include zero

    def test_control_removes_a_purely_confounded_signal(self):
        """The whole point of the project, as a unit test.

        A latent "functional importance" drives everything: conserved positions are
        functional, models react at functional positions, and functional positions are
        more often pathogenic. delta has no independent effect on the label. Here
        conservation is the CLEANER measurement of that latent, so controlling for it
        should remove delta's apparent signal.
        """
        rng = np.random.default_rng(7)
        n = 2000
        f = rng.standard_normal(n)
        c = f + 0.10 * rng.standard_normal(n)                     # precise proxy for f
        delta = np.log1p(np.exp(f)) + 0.60 * rng.random(n)        # noisier proxy
        y = rng.binomial(1, 1 / (1 + np.exp(-(-0.3 + 1.8 * f))))  # only f matters

        alone = cons.fit_delta_coef(delta, y, None, n_boot=300, method="firth")
        ctrl = cons.fit_delta_coef(delta, y, c, n_boot=300, method="firth")
        assert alone.survives                    # naive test sees a signal
        assert abs(ctrl.coef) < abs(alone.coef)  # controlling attenuates it
        assert not ctrl.survives                 # and it does not survive

    def test_residual_confounding_when_delta_measures_the_latent_better(self):
        """A KNOWN LIMITATION, asserted so it is never forgotten.

        Adjusting for a mismeasured confounder does not fully remove confounding. If
        delta happens to track the latent "functional importance" more precisely than
        phyloP does, delta keeps a significant coefficient despite having no independent
        causal effect -- a false positive our test cannot rule out.

        This is why a surviving signal is reported as a lead rather than a result, and
        why the delta-conservation correlation must be reported alongside it: the
        proteins where a model survives AND is most conservation-entangled are exactly
        where this failure mode is most likely.
        """
        rng = np.random.default_rng(7)
        n = 2000
        f = rng.standard_normal(n)
        c = f + 0.60 * rng.standard_normal(n)                     # noisy proxy
        delta = np.log1p(np.exp(f)) + 0.10 * rng.random(n)        # precise proxy
        y = rng.binomial(1, 1 / (1 + np.exp(-(-0.3 + 1.8 * f))))  # still only f matters

        ctrl = cons.fit_delta_coef(delta, y, c, n_boot=300, method="firth")
        assert ctrl.survives, (
            "delta has no causal effect here, yet it survives, because conservation is "
            "the noisier measurement. This is residual confounding, not a real signal.")

    def test_bootstrap_p_is_floored(self):
        d, y, c = self._data(beta=3.0)
        f = cons.fit_delta_coef(d, y, c, n_boot=100, method="firth")
        assert f.p_boot >= 1.0 / f.n_boot_ok    # cannot resolve below one resample

    def test_reports_when_bootstrap_mostly_failed(self):
        y = np.array([1] * 3 + [0] * 3)
        d = np.array([1.0, 1.1, 1.2, 0.1, 0.2, 0.3])
        c = np.zeros(6)
        f = cons.fit_delta_coef(d, y, c, n_boot=20, method="firth")
        # tiny n: interval may be unavailable, but it must not raise
        assert f.n_boot_ok <= 20


class TestTestGroup:
    def test_refuses_insufficient_data(self):
        out = cons.test_group(np.arange(5.0), np.array([1, 0, 1, 0, 1]),
                              np.zeros(5), n_boot=20)
        assert out["note"] == "insufficient data"
        assert "controlled_coef" not in out

    def test_refuses_single_class(self):
        out = cons.test_group(np.arange(50.0), np.ones(50, dtype=int),
                              np.zeros(50), n_boot=20)
        assert out["note"] == "insufficient data"

    def test_drops_nan_conservation_and_counts_it(self):
        rng = np.random.default_rng(0)
        n = 200
        d = rng.random(n)
        y = rng.integers(0, 2, n)
        c = rng.standard_normal(n)
        c[:20] = np.nan
        out = cons.test_group(d, y, c, n_boot=100)
        assert out["n_dropped_nan"] == 20
        assert out["n"] == n - 20

    def test_reports_conservation_benchmark(self):
        rng = np.random.default_rng(2)
        n = 400
        c = rng.standard_normal(n)
        y = rng.binomial(1, 1 / (1 + np.exp(-2.0 * c)))
        out = cons.test_group(rng.random(n), y, c, n_boot=100)
        # conservation alone should be strongly predictive here, and we must say so
        assert out["conservation_auroc"] > 0.8


class TestRunAndFdr:
    def test_adds_q_values_and_fdr_flags(self):
        rng = np.random.default_rng(0)
        rows = []
        for prot in ("A", "B", "C"):
            n = 300
            c = rng.standard_normal(n)
            d = rng.random(n)
            y = rng.binomial(1, 1 / (1 + np.exp(-(0.5 * c))))
            rows.append(pd.DataFrame({"protein": prot, "m1": d, "m2": rng.random(n),
                                      "label": y, "conservation": c}))
        df = pd.concat(rows, ignore_index=True)
        res = cons.run(df, ["m1", "m2"], group_col="protein", n_boot=120)
        assert len(res) == 6
        for col in ("controlled_q", "controlled_survives_fdr",
                    "alone_q", "alone_survives_fdr"):
            assert col in res.columns
        assert res.controlled_q.dropna().between(0, 1).all()

    def test_fdr_is_stricter_than_the_raw_interval(self):
        # q-values can only make a call harder, never easier
        res = pd.DataFrame({"controlled_coef": [0.5, 0.4, 0.3],
                            "controlled_p_wald": [0.04, 0.045, 0.048],
                            "alone_coef": [0.5, 0.4, 0.3],
                            "alone_p_wald": [0.04, 0.045, 0.048],
                            "note": ["", "", ""]})
        out = cons.add_fdr(res, alpha=0.05)
        assert (out.controlled_q >= out.controlled_p_wald).all()

    def test_pooled_mode(self):
        rng = np.random.default_rng(1)
        n = 400
        df = pd.DataFrame({"m1": rng.random(n), "label": rng.integers(0, 2, n),
                           "conservation": rng.standard_normal(n)})
        res = cons.run(df, ["m1"], group_col=None, n_boot=100)
        assert len(res) == 1 and res.iloc[0]["group"] == "ALL"


class TestPower:
    def test_larger_n_detects_smaller_effects(self):
        small = cons.min_detectable_effect(120, 0.3, n_sim=25, n_boot=80, seed=0)
        large = cons.min_detectable_effect(1200, 0.3, n_sim=25, n_boot=80, seed=0)
        assert large <= small

    def test_returns_inf_when_hopeless(self):
        mde = cons.min_detectable_effect(25, 0.1, effects=[0.1, 0.2],
                                         n_sim=15, n_boot=60, seed=0)
        assert mde == float("inf")


class TestGenomicPos:
    """Index-to-genome mapping. Getting the minus strand wrong silently shifts every
    conservation lookup to the wrong end of the window."""

    def test_plus_strand_counts_forward(self):
        assert pc.genomic_pos(1000, 1101, "+", 0) == 1001
        assert pc.genomic_pos(1000, 1101, "+", 50) == 1051

    def test_minus_strand_counts_backward(self):
        # the window sequence is reverse-complemented, so index 0 is the highest coordinate
        assert pc.genomic_pos(1000, 1101, "-", 0) == 1101
        assert pc.genomic_pos(1000, 1101, "-", 50) == 1051

    def test_centre_agrees_on_both_strands(self):
        assert (pc.genomic_pos(1000, 1101, "+", 50)
                == pc.genomic_pos(1000, 1101, "-", 50))

    def test_zero_based_option(self):
        assert pc.genomic_pos(1000, 1101, "+", 0, one_based=False) == 1000


class TestFindAll:
    def test_finds_every_occurrence(self):
        assert pc.find_all("AUGCAUGCAUG", "CAUG") == [3, 7]

    def test_finds_overlapping(self):
        assert pc.find_all("AAAA", "AA") == [0, 1, 2]

    def test_absent(self):
        assert pc.find_all("ACGU", "GGGG") == []


class TestBuildPairs:
    def _df(self, seqs):
        return pd.DataFrame([
            {"id": f"w{i}", "chrom": "chr1", "start": 1000 + i * 200,
             "end": 1101 + i * 200, "strand": "+", "seq_rna": s}
            for i, s in enumerate(seqs)])

    def test_makes_one_disruptive_and_one_neutral(self):
        seq = "GCAUG" + "A" * 96
        pairs, dropped = pc.build_pairs(self._df([seq]), "GCAUG", 2, "C",
                                        min_distance=25)
        assert len(pairs) == 2
        assert set(pairs.kind) == {"disruptive", "neutral"}
        assert list(pairs.label) == [1, 0]

    def test_both_mutants_use_the_same_substitution(self):
        seq = "GCAUG" + "A" * 96
        pairs, _ = pc.build_pairs(self._df([seq]), "GCAUG", 2, "C", min_distance=25)
        assert pairs.ref_base.nunique() == 1     # same base changed
        assert pairs.alt_base.nunique() == 1     # to the same replacement

    def test_neutral_site_is_far_from_every_motif_hit(self):
        seq = "GCAUG" + "A" * 40 + "GCAUG" + "A" * 51
        pairs, _ = pc.build_pairs(self._df([seq]), "GCAUG", 2, "C", min_distance=25)
        neutral = pairs[pairs.kind == "neutral"].iloc[0]
        for h in pc.find_all(seq, "GCAUG"):
            assert abs(neutral.mut_index - h) >= 25

    def test_alt_seq_differs_by_exactly_one_base(self):
        seq = "GCAUG" + "A" * 96
        pairs, _ = pc.build_pairs(self._df([seq]), "GCAUG", 2, "C", min_distance=25)
        for r in pairs.itertuples():
            diffs = sum(1 for a, b in zip(r.ref_seq, r.alt_seq) if a != b)
            assert diffs == 1
            assert len(r.ref_seq) == len(r.alt_seq)

    def test_counts_why_windows_were_dropped(self):
        no_motif = "A" * 101
        crowded = "GCAUG" * 20 + "A"          # motif everywhere: no neutral site
        _, dropped = pc.build_pairs(self._df([no_motif, crowded]), "GCAUG", 2, "C",
                                    min_distance=25)
        assert dropped["no_motif"] == 1
        assert dropped["no_neutral_site"] == 1

    def test_skips_when_base_is_already_the_replacement(self):
        seq = "GCCUG" + "A" * 96              # index 2 is already C
        _, dropped = pc.build_pairs(self._df([seq]), "GCCUG", 2, "C", min_distance=25)
        assert dropped["already_replacement"] == 1


class TestScoreAndEffect:
    def test_reference_scored_once_per_unique_window(self):
        seq = "GCAUG" + "A" * 96
        pairs, _ = pc.build_pairs(
            pd.DataFrame([{"id": "w", "chrom": "chr1", "start": 0, "end": 101,
                           "strand": "+", "seq_rna": seq}]), "GCAUG", 2, "C")
        calls = []

        def fake_score(seqs, batch_size=64):
            calls.append(len(seqs))
            return np.linspace(0.1, 0.9, len(seqs))

        out = pc.score_pairs(pairs, fake_score)
        assert calls == [1, 2]                  # 1 unique ref, 2 alts
        assert (out.delta >= 0).all()

    def test_effect_size_direction(self):
        scored = pd.DataFrame({
            "kind": ["disruptive"] * 4 + ["neutral"] * 4,
            "delta": [0.5, 0.6, 0.55, 0.45, 0.02, 0.03, 0.01, 0.04]})
        e = pc.effect_size(scored)
        assert e["cohens_d"] > 2
        assert e["delta_disruptive_mean"] > e["delta_neutral_mean"]


class TestPhylopHelpers:
    def test_chrom_normalisation_is_idempotent(self):
        assert phylop._norm("1") == "chr1"
        assert phylop._norm("chr1") == "chr1"

    def test_cache_roundtrip_including_nan(self, tmp_path):
        p = tmp_path / "cache.tsv"
        vals = {("chr1", 100): 1.5, ("chr2", 200): float("nan")}
        phylop.save_cache(str(p), vals)
        back = phylop.load_cache(str(p))
        assert back[("chr1", 100)] == 1.5
        assert np.isnan(back[("chr2", 200)])

    def test_missing_cache_is_empty_not_an_error(self, tmp_path):
        assert phylop.load_cache(str(tmp_path / "nope.tsv")) == {}

    def test_coverage_reports_missing(self):
        df = pd.DataFrame({"conservation": [1.0, np.nan, 2.0, np.nan]})
        c = phylop.coverage(df)
        assert c == {"n": 4, "missing": 2, "frac_missing": 0.5}
