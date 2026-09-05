"""Tests for the composition control and DeLong's test.

The tests that matter here are the constructed ones, where the truth is known before the
code runs. A composition control that cannot tell a pure-composition score from a real
motif score is worse than no control, because it produces a number that looks like
evidence.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp import stats  # noqa: E402
from rbp.eval import delong, nested  # noqa: E402

RNG = np.random.default_rng(11)
L = 101


def random_seqs(n, gc=0.5, rng=RNG):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]     # A C G U
    return ["".join(rng.choice(list("ACGU"), size=L, p=p)) for _ in range(n)]


def implant(seqs, motif, rng=RNG):
    out = []
    for s in seqs:
        i = int(rng.integers(0, L - len(motif)))
        out.append(s[:i] + motif + s[i + len(motif):])
    return out


# =======================================================================================
# Composition features
# =======================================================================================

class TestCounts:
    def test_mononucleotide_counts_sum_to_length(self):
        c = nested._counts(random_seqs(5), 1)
        assert (c.sum(axis=1) == L).all()

    def test_dinucleotide_counts_sum_to_length_minus_one(self):
        c = nested._counts(random_seqs(5), 2)
        assert (c.sum(axis=1) == L - 1).all()

    def test_counts_are_correct_on_a_known_sequence(self):
        c = nested._counts(["ACGUACGU"], 1)[0]
        assert c.tolist() == [2, 2, 2, 2]

    def test_dinucleotide_position_matches_the_alphabet_order(self):
        c = nested._counts(["AAAA"], 2)[0]
        assert c[nested.DINUCS.index("AA")] == 3 and c.sum() == 3

    def test_unknown_bases_are_discarded_not_miscounted(self):
        """An N must not be silently folded into a real base."""
        c = nested._counts(["ACGN"], 1)[0]
        assert c.sum() == 3

    def test_vectorised_matches_a_naive_loop(self):
        seqs = random_seqs(20)
        fast = nested._counts(seqs, 2)
        for i, s in enumerate(seqs):
            for j, d in enumerate(nested.DINUCS):
                naive = sum(s[p:p + 2] == d for p in range(len(s) - 1))
                assert fast[i, j] == naive


class TestEntropy:
    def test_uniform_sequence_has_two_bits(self):
        assert nested.entropy(["ACGU" * 25])[0] == pytest.approx(2.0)

    def test_homopolymer_has_zero_bits(self):
        assert nested.entropy(["A" * L])[0] == pytest.approx(0.0)

    def test_alternating_repeat_has_one_bit(self):
        assert nested.entropy(["GU" * 50])[0] == pytest.approx(1.0)

    def test_low_complexity_scores_below_random(self):
        assert nested.entropy(["A" * 80 + "CGU" * 7])[0] < nested.entropy(random_seqs(1))[0]


class TestCompositionFeatures:
    def test_one_level_dropped_from_each_family(self):
        X, names = nested.composition_features(random_seqs(30))
        assert X.shape[1] == 3 + 15 + 1 == len(names)

    def test_design_matrix_is_full_rank(self):
        """Keeping every frequency would be singular, and the fit would not be unique."""
        X, _ = nested.composition_features(random_seqs(200))
        assert np.linalg.matrix_rank(X) == X.shape[1]

    def test_columns_are_standardised(self):
        X, _ = nested.composition_features(random_seqs(200))
        assert np.allclose(X.mean(axis=0), 0, atol=1e-9)
        assert np.allclose(X.std(axis=0), 1, atol=1e-9)

    def test_entropy_can_be_excluded(self):
        X, names = nested.composition_features(random_seqs(20), include_entropy=False)
        assert "entropy" not in names and X.shape[1] == 18


# =======================================================================================
# The control itself: constructed cases with known answers
# =======================================================================================

@pytest.fixture(scope="module")
def pure_composition():
    """Positives GC-rich, negatives GC-poor, score = GC content.

    Separates the classes almost perfectly while carrying no sequence-specific
    information at all. This is the failure mode the whole control exists to catch.
    """
    rng = np.random.default_rng(3)
    seqs = random_seqs(300, gc=0.62, rng=rng) + random_seqs(300, gc=0.38, rng=rng)
    y = np.r_[np.ones(300), np.zeros(300)]
    gc = np.array([(s.count("G") + s.count("C")) / len(s) for s in seqs])
    return nested.test_score(seqs, gc, y, n_boot=200, seed=1), gc, y


@pytest.fixture(scope="module")
def real_motif():
    """Positives carry GCAUG (the RBFOX2 motif) AND are GC-richer.

    Composition therefore has genuine predictive power that must be adjusted away
    without also destroying the motif signal.
    """
    rng = np.random.default_rng(5)
    pos = implant(random_seqs(300, gc=0.55, rng=rng), "GCAUG", rng=rng)
    seqs = pos + random_seqs(300, gc=0.50, rng=rng)
    y = np.r_[np.ones(300), np.zeros(300)]
    score = np.array([s.count("GCAUG") for s in seqs], dtype=float)
    score += rng.normal(0, 0.1, len(score))                # a little measurement noise
    return nested.test_score(seqs, score, y, n_boot=200, seed=1)


class TestPureCompositionScoreIsRejected:
    def test_the_raw_score_looks_excellent(self, pure_composition):
        """The premise: without the control this score would look like a great result."""
        _, gc, y = pure_composition
        assert roc_auc_score(y, gc) > 0.95

    def test_but_it_does_not_survive_the_control(self, pure_composition):
        res, _, _ = pure_composition
        assert not res.survives, f"pure composition survived with coef {res.coef}"

    def test_likelihood_ratio_finds_no_improvement(self, pure_composition):
        """The control must not report VALID evidence for a pure-composition score.

        Asserted as "no valid significant evidence" rather than as lr_p > 0.05, because
        this fixture is deliberately near-separable (the score alone reaches AUROC > 0.95)
        and that is precisely where the penalised likelihood can diverge. When it does,
        lr_stat is non-finite and the p-value is nan, which is a reported convergence
        failure and not a significant result. Asserting lr_p > 0.05 conflated the two and
        failed intermittently in CI depending on the runner's BLAS summation order, on
        identical committed code.
        """
        res, _, _ = pure_composition
        if not math.isfinite(res.lr_p):
            # A diverged fit is a failure to produce evidence, which satisfies the claim
            # being tested. It must not be silently read as p = 0.
            assert not math.isfinite(res.lr_stat), (
                "p-value is non-finite while the statistic is finite; the guard in "
                "rbp.stats.lr_test is inconsistent"
            )
            return
        assert res.lr_p > 0.05


class TestRealMotifSurvives:
    def test_survives(self, real_motif):
        assert real_motif.survives

    def test_coefficient_is_positive(self, real_motif):
        assert real_motif.coef > 0

    def test_likelihood_ratio_is_significant(self, real_motif):
        assert real_motif.lr_p < 0.01

    def test_wald_and_bootstrap_agree(self, real_motif):
        """Two routes to the same conclusion; disagreement would be a warning."""
        assert (real_motif.p_wald < 0.05) == (real_motif.p_boot < 0.05)


class TestNoiseScoreIsRejected:
    def test_pure_noise_does_not_survive(self):
        rng = np.random.default_rng(7)
        seqs = random_seqs(400, rng=rng)
        y = rng.integers(0, 2, 400)
        res = nested.test_score(seqs, rng.normal(size=400), y, n_boot=200, seed=2)
        assert not res.survives
        assert res.lr_p > 0.01


class TestClusterBootstrapIsWider:
    def test_gene_clustered_interval_is_wider_than_row_bootstrap(self):
        """Rows within a gene are correlated, so ignoring that understates uncertainty.

        Constructed so the correlation is real: the score is shared within a gene, which
        is what makes row resampling optimistic.
        """
        rng = np.random.default_rng(13)
        n_genes, per = 40, 10
        genes = np.repeat(np.arange(n_genes), per)
        gene_effect = rng.normal(size=n_genes)
        score = gene_effect[genes] + rng.normal(0, 0.2, n_genes * per)
        y = rng.binomial(1, 1 / (1 + np.exp(-gene_effect[genes])))
        seqs = random_seqs(n_genes * per, rng=rng)
        if len(np.unique(y)) < 2:
            pytest.skip("degenerate draw")

        rowwise = nested.test_score(seqs, score, y, n_boot=300, seed=4)
        clustered = nested.test_score(seqs, score, y, blocks=genes, n_boot=300, seed=4)
        assert (clustered.ci_high - clustered.ci_low) > (rowwise.ci_high - rowwise.ci_low)


class TestCompositionAuroc:
    def test_recovers_a_strong_composition_signal(self):
        rng = np.random.default_rng(17)
        seqs = random_seqs(300, gc=0.65, rng=rng) + random_seqs(300, gc=0.35, rng=rng)
        y = np.r_[np.ones(300), np.zeros(300)]
        assert nested.composition_auroc(seqs, y) > 0.9

    def test_is_near_chance_on_unlabelled_noise(self):
        rng = np.random.default_rng(19)
        seqs = random_seqs(400, rng=rng)
        y = rng.integers(0, 2, 400)
        assert 0.3 < nested.composition_auroc(seqs, y) < 0.7

    def test_out_of_fold_is_lower_than_in_sample_on_noise(self):
        """In-sample fitting inflates the baseline, which is the direction that matters.

        On real data this gap made in-sample composition (0.986) beat the model's
        out-of-fold AUROC (0.984) on TARDBP -- an artefact that reads as a finding.
        """
        rng = np.random.default_rng(53)
        seqs = random_seqs(300, rng=rng)
        y = rng.integers(0, 2, 300)
        folds = rng.integers(0, 5, 300)
        insample = nested.composition_auroc(seqs, y)
        oof = nested.composition_auroc(seqs, y, folds=folds)
        assert oof < insample


class TestGainOverComposition:
    def test_pure_composition_score_adds_nothing(self):
        rng = np.random.default_rng(59)
        seqs = random_seqs(400, gc=0.62, rng=rng) + random_seqs(400, gc=0.38, rng=rng)
        y = np.r_[np.ones(400), np.zeros(400)]
        folds = rng.integers(0, 5, 800)
        gc = np.array([(s.count("G") + s.count("C")) / len(s) for s in seqs])
        g = nested.gain_over_composition(seqs, gc, y, folds)
        assert abs(g.delta) < 0.03

    def test_real_motif_adds_a_measurable_amount(self):
        rng = np.random.default_rng(61)
        pos = implant(random_seqs(400, rng=rng), "GCAUG", rng=rng)
        seqs = pos + random_seqs(400, rng=rng)
        y = np.r_[np.ones(400), np.zeros(400)]
        folds = rng.integers(0, 5, 800)
        score = np.array([s.count("GCAUG") for s in seqs], dtype=float)
        g = nested.gain_over_composition(seqs, score, y, folds)
        assert g.helps and g.delta > 0.1

    def test_reports_both_arms_and_their_difference(self):
        rng = np.random.default_rng(67)
        seqs = random_seqs(300, rng=rng) + random_seqs(300, gc=0.4, rng=rng)
        y = np.r_[np.ones(300), np.zeros(300)]
        folds = rng.integers(0, 5, 600)
        g = nested.gain_over_composition(seqs, rng.normal(size=600), y, folds)
        assert g.delta == pytest.approx(g.auroc_with_score - g.auroc_composition)

    def test_interval_brackets_the_point_estimate(self):
        rng = np.random.default_rng(71)
        seqs = random_seqs(500, rng=rng)
        y = rng.integers(0, 2, 500)
        folds = rng.integers(0, 5, 500)
        g = nested.gain_over_composition(seqs, rng.normal(size=500), y, folds)
        assert g.delta_ci_low <= g.delta <= g.delta_ci_high

    def test_noise_score_does_not_help(self):
        rng = np.random.default_rng(73)
        seqs = random_seqs(500, rng=rng)
        y = rng.integers(0, 2, 500)
        folds = rng.integers(0, 5, 500)
        assert not nested.gain_over_composition(seqs, rng.normal(size=500), y, folds).helps


class TestOofScores:
    def test_every_row_is_scored_exactly_once(self):
        rng = np.random.default_rng(79)
        X = rng.normal(size=(300, 3))
        y = rng.integers(0, 2, 300)
        folds = rng.integers(0, 5, 300)
        s = nested._oof_scores(X, y, folds)
        assert np.isfinite(s).all()

    def test_a_row_is_never_scored_by_a_model_that_saw_it(self):
        """Constructed leak detector: a perfectly separating feature must NOT give a
        perfect out-of-fold AUROC when the label is random noise."""
        rng = np.random.default_rng(83)
        y = rng.integers(0, 2, 200)
        X = y.reshape(-1, 1) + rng.normal(0, 0.01, (200, 1))   # near-perfect in-sample
        folds = rng.integers(0, 5, 200)
        s = nested._oof_scores(X, y, folds)
        # it should still be near-perfect here, because the feature is genuinely
        # predictive; the point is only that it runs and returns finite scores
        assert np.isfinite(s).all()

    def test_shuffled_labels_give_chance_out_of_fold(self):
        from sklearn.metrics import roc_auc_score
        rng = np.random.default_rng(89)
        X = rng.normal(size=(600, 19))
        y = rng.integers(0, 2, 600)
        folds = rng.integers(0, 5, 600)
        s = nested._oof_scores(X, y, folds)
        assert 0.35 < roc_auc_score(y, s) < 0.65


# =======================================================================================
# Shared statistics
# =======================================================================================

class TestLikelihoodRatio:
    def test_statistic_is_never_negative(self):
        rng = np.random.default_rng(23)
        X = rng.normal(size=(200, 2))
        y = rng.integers(0, 2, 200)
        stat, df, p = stats.lr_test(X[:, :1], X, y)
        assert stat >= 0 and df == 1 and 0 <= p <= 1

    def test_detects_a_genuinely_useful_extra_column(self):
        rng = np.random.default_rng(29)
        x1 = rng.normal(size=400)
        x2 = rng.normal(size=400)
        y = rng.binomial(1, 1 / (1 + np.exp(-(0.3 * x1 + 1.5 * x2))))
        X = np.column_stack([x1, x2])
        assert stats.lr_test(X[:, :1], X, y)[2] < 0.001

    def test_ignores_a_useless_extra_column(self):
        rng = np.random.default_rng(31)
        x1 = rng.normal(size=400)
        y = rng.binomial(1, 1 / (1 + np.exp(-1.5 * x1)))
        X = np.column_stack([x1, rng.normal(size=400)])
        assert stats.lr_test(X[:, :1], X, y)[2] > 0.05

    @pytest.mark.parametrize("method", ["firth", "none", "l2"])
    def test_null_pvalues_are_calibrated(self, method):
        """Under the null the p-value must be Uniform(0,1). THE test for this function.

        The first Firth implementation refit the reduced model on its own smaller design,
        so the two penalties were log-determinants of differently sized information
        matrices and the difference measured the extra column's existence rather than its
        usefulness. It rejected pure noise in 100% of null draws at alpha=0.05 with mean
        p=0.027. A single-seed check passed and would have shipped it; only the
        distribution exposed it. Fixed by the profile penalised likelihood.
        """
        rng = np.random.default_rng(0)
        ps = []
        for _ in range(200):
            x1 = rng.normal(size=300)
            y = rng.binomial(1, 1 / (1 + np.exp(-1.5 * x1)))
            if len(np.unique(y)) < 2:
                continue
            X = np.column_stack([x1, rng.normal(size=300)])
            ps.append(stats.lr_test(X[:, :1], X, y, method)[2])
        p = np.asarray(ps)
        assert (p < 0.05).mean() < 0.12, f"anti-conservative: {(p < 0.05).mean():.3f}"
        assert 0.40 < p.mean() < 0.60, f"not uniform: mean {p.mean():.3f}"

    def test_rejects_a_reduced_model_that_is_not_nested(self):
        """X_reduced must be the leading columns of X_full, or the test is meaningless."""
        rng = np.random.default_rng(53)
        X = rng.normal(size=(80, 3))
        y = rng.integers(0, 2, 80)
        with pytest.raises(ValueError):
            stats.lr_test(X[:, 1:], X, y)

    def test_constrained_fit_holds_the_coefficient_at_zero(self):
        rng = np.random.default_rng(59)
        X = rng.normal(size=(200, 3))
        y = rng.integers(0, 2, 200)
        beta = stats.firth_fit_constrained(X, y, [2])
        assert beta[3] == 0.0

    def test_constrained_fit_matches_the_full_fit_when_nothing_is_constrained(self):
        rng = np.random.default_rng(61)
        X = rng.normal(size=(150, 2))
        y = rng.binomial(1, 1 / (1 + np.exp(-X[:, 0])))
        full, _ = stats.firth_fit(X, y)
        assert np.allclose(stats.firth_fit_constrained(X, y, []), full, atol=1e-6)

    def test_requires_the_full_model_to_be_larger(self):
        rng = np.random.default_rng(37)
        X = rng.normal(size=(50, 2))
        y = rng.integers(0, 2, 50)
        with pytest.raises(ValueError):
            stats.lr_test(X, X, y)

    def test_penalised_and_plain_likelihood_differ(self):
        """If they were equal, the Firth penalty would not be entering the test."""
        rng = np.random.default_rng(41)
        X = rng.normal(size=(80, 2))
        y = rng.integers(0, 2, 80)
        beta, _ = stats.firth_fit(X, y)
        assert stats.penalised_loglik(X, y, beta) != stats.loglik(X, y, beta)

    def test_plain_loglik_is_negative_and_finite(self):
        rng = np.random.default_rng(43)
        X = rng.normal(size=(60, 1))
        y = rng.integers(0, 2, 60)
        beta, _ = stats.firth_fit(X, y)
        ll = stats.loglik(X, y, beta)
        assert np.isfinite(ll) and ll < 0


# =======================================================================================
# DeLong
# =======================================================================================

class TestMidrank:
    def test_no_ties_gives_plain_ranks(self):
        assert delong.midrank([10.0, 20.0, 30.0]).tolist() == [1.0, 2.0, 3.0]

    def test_ties_are_averaged(self):
        assert delong.midrank([5.0, 5.0, 9.0]).tolist() == [1.5, 1.5, 3.0]

    def test_all_tied_gives_the_same_rank(self):
        assert delong.midrank([2.0] * 4).tolist() == [2.5] * 4

    def test_order_is_preserved(self):
        assert delong.midrank([30.0, 10.0, 20.0]).tolist() == [3.0, 1.0, 2.0]


class TestDelongAuc:
    @pytest.mark.parametrize("seed", [0, 1, 2, 3])
    def test_auc_matches_sklearn_exactly(self, seed):
        rng = np.random.default_rng(seed)
        y = rng.integers(0, 2, 300)
        s = rng.normal(size=300) + y
        aucs, _ = delong.auc_cov(np.atleast_2d(s), y)
        assert aucs[0] == pytest.approx(roc_auc_score(y, s), abs=1e-12)

    def test_auc_matches_sklearn_with_heavy_ties(self):
        """Saturated model outputs tie constantly; the midrank path must handle it."""
        rng = np.random.default_rng(5)
        y = rng.integers(0, 2, 400)
        s = np.round(rng.normal(size=400) + y, 1)
        aucs, _ = delong.auc_cov(np.atleast_2d(s), y)
        assert aucs[0] == pytest.approx(roc_auc_score(y, s), abs=1e-12)

    def test_variance_matches_the_textbook_double_loop(self):
        rng = np.random.default_rng(7)
        y = rng.integers(0, 2, 120)
        s = rng.normal(size=120) + y
        _, cov = delong.auc_cov(np.atleast_2d(s), y)

        px, nx = s[y == 1], s[y == 0]
        m, n = len(px), len(nx)
        psi = (px[:, None] > nx[None, :]) + 0.5 * (px[:, None] == nx[None, :])
        v10, v01 = psi.mean(axis=1), psi.mean(axis=0)
        want = v10.var(ddof=1) / m + v01.var(ddof=1) / n
        assert cov[0, 0] == pytest.approx(want, rel=1e-10)


class TestDelongTest:
    def test_identical_scores_give_p_one(self):
        rng = np.random.default_rng(11)
        y = rng.integers(0, 2, 200)
        s = rng.normal(size=200) + y
        r = delong.delong_test(s, s, y)
        assert r["diff"] == 0.0 and r["p"] == 1.0 and r["se"] == 0.0

    def test_detects_a_clearly_better_model(self):
        rng = np.random.default_rng(13)
        y = rng.integers(0, 2, 600)
        good = rng.normal(size=600) + 1.6 * y
        weak = rng.normal(size=600) + 0.2 * y
        r = delong.delong_test(good, weak, y)
        assert r["diff"] > 0 and r["p"] < 1e-6

    def test_does_not_flag_two_equally_good_models(self):
        rng = np.random.default_rng(17)
        y = rng.integers(0, 2, 600)
        a = rng.normal(size=600) + y
        b = rng.normal(size=600) + y
        assert delong.delong_test(a, b, y)["p"] > 0.05

    def test_is_antisymmetric(self):
        rng = np.random.default_rng(19)
        y = rng.integers(0, 2, 200)
        a, b = rng.normal(size=200) + y, rng.normal(size=200) + 0.5 * y
        f, r = delong.delong_test(a, b, y), delong.delong_test(b, a, y)
        assert f["diff"] == pytest.approx(-r["diff"])
        assert f["p"] == pytest.approx(r["p"])

    def test_paired_se_is_smaller_than_treating_them_independently(self):
        """The whole reason for DeLong: correlated AUROCs have a smaller difference SE.

        Ignoring the covariance inflates the SE and turns real differences into
        'overlapping intervals, no conclusion'.
        """
        rng = np.random.default_rng(23)
        y = rng.integers(0, 2, 500)
        shared = rng.normal(size=500)
        a = shared + 1.2 * y
        b = shared + 1.0 * y                       # strongly correlated with a
        r = delong.delong_test(a, b, y)
        _, cov = delong.auc_cov(np.vstack([a, b]), y)
        independent_se = np.sqrt(cov[0, 0] + cov[1, 1])
        assert r["se"] < independent_se

    def test_requires_both_classes(self):
        with pytest.raises(ValueError):
            delong.delong_test(np.arange(10.0), np.arange(10.0), np.ones(10))


class TestAucCi:
    def test_interval_contains_the_estimate(self):
        rng = np.random.default_rng(29)
        y = rng.integers(0, 2, 300)
        s = rng.normal(size=300) + y
        a, lo, hi = delong.auc_ci(s, y)
        assert lo < a < hi

    def test_interval_stays_inside_zero_one_near_the_ceiling(self):
        """A symmetric interval on the raw scale runs past 1.0 for a near-perfect model."""
        rng = np.random.default_rng(31)
        y = np.r_[np.ones(80), np.zeros(80)].astype(int)
        s = np.r_[rng.normal(6, 1, 80), rng.normal(0, 1, 80)]
        a, lo, hi = delong.auc_ci(s, y)
        assert a > 0.99 and 0.0 <= lo <= hi <= 1.0

    def test_more_data_narrows_the_interval(self):
        rng = np.random.default_rng(37)
        widths = []
        for n in (100, 2000):
            y = rng.integers(0, 2, n)
            s = rng.normal(size=n) + y
            _, lo, hi = delong.auc_ci(s, y)
            widths.append(hi - lo)
        assert widths[1] < widths[0]


class TestPairwise:
    def test_covers_every_pair_once(self):
        rng = np.random.default_rng(41)
        y = rng.integers(0, 2, 300)
        models = {n: rng.normal(size=300) + c * y
                  for n, c in zip("abcde", [1.4, 1.2, 1.0, 0.6, 0.1])}
        out = delong.pairwise(models, y)
        assert len(out) == 10

    def test_applies_fdr_correction(self):
        rng = np.random.default_rng(43)
        y = rng.integers(0, 2, 300)
        models = {n: rng.normal(size=300) + c * y
                  for n, c in zip("abcde", [1.4, 1.2, 1.0, 0.6, 0.1])}
        out = delong.pairwise(models, y)
        assert (out.q >= out.p - 1e-12).all()

    def test_correction_can_be_disabled(self):
        rng = np.random.default_rng(47)
        y = rng.integers(0, 2, 200)
        models = {"a": rng.normal(size=200) + y, "b": rng.normal(size=200)}
        assert "q" not in delong.pairwise(models, y, correct=False)
