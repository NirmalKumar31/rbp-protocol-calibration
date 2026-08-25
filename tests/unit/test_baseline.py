"""Tests for the k-mer baseline under cross-validation folds.

The property worth guarding is the fold discipline. A model that scores rows it trained
on produces an inflated AUROC and a variant delta that partly memorises its own training
data, and neither failure raises anything.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp.eval import baseline  # noqa: E402

L = 60
K = 5


def seqs_with_motif(n, motif, rng):
    out = []
    for _ in range(n):
        s = list(rng.choice(list("ACGU"), size=L))
        i = int(rng.integers(0, L - len(motif)))
        s[i:i + len(motif)] = list(motif)
        out.append("".join(s))
    return out


def random_seqs(n, rng):
    return ["".join(rng.choice(list("ACGU"), size=L)) for _ in range(n)]


@pytest.fixture(scope="module")
def toy():
    rng = np.random.default_rng(7)
    pos = seqs_with_motif(250, "GCAUG", rng)
    neg = random_seqs(250, rng)
    seqs = pos + neg
    y = np.r_[np.ones(250), np.zeros(250)].astype(int)
    folds = np.tile(np.arange(K), 100)
    return seqs, y, folds


class TestKmerMatrix:
    def test_shape_is_sequences_by_vocabulary(self):
        rng = np.random.default_rng(1)
        X, vec = baseline.kmer_matrix(random_seqs(10, rng), 3)
        assert X.shape[0] == 10 and X.shape[1] == len(vec.get_feature_names_out())

    def test_counts_sum_to_the_number_of_kmers(self):
        rng = np.random.default_rng(2)
        X, _ = baseline.kmer_matrix(random_seqs(10, rng), 4)
        assert (np.asarray(X.sum(axis=1)).ravel() == L - 3).all()

    def test_vectoriser_is_reusable_on_new_sequences(self):
        """Variant scoring must transform with the SAME vocabulary or the coefficients
        refer to the wrong columns -- silently, since the shapes can still match."""
        rng = np.random.default_rng(3)
        X, vec = baseline.kmer_matrix(random_seqs(20, rng), 3)
        assert vec.transform(random_seqs(5, rng)).shape[1] == X.shape[1]


class TestFoldModels:
    def test_one_model_per_fold(self, toy):
        seqs, y, folds = toy
        models, _ = baseline.fit_fold_models(seqs, y, folds, k=K)
        assert set(models) == set(range(K))

    def test_a_fold_model_never_saw_its_own_fold(self, toy):
        """Constructed leak detector: poison one fold's labels. If the model for that
        fold trained on it, its predictions there would follow the poison."""
        seqs, y, folds = toy
        poisoned = y.copy()
        poisoned[folds == 0] = 1 - poisoned[folds == 0]
        models, vec = baseline.fit_fold_models(seqs, poisoned, folds, k=K)
        X = vec.transform(seqs)
        sel = folds == 0
        s = models[0].decision_function(X[sel])
        # scored against the POISONED labels it never saw, it should be no better than
        # chance and in fact anti-correlated, because the real signal is inverted there
        assert roc_auc_score(poisoned[sel], s) < 0.5

    def test_single_class_fold_is_skipped_not_crashed(self):
        rng = np.random.default_rng(5)
        seqs = random_seqs(100, rng)
        y = np.zeros(100, dtype=int)
        y[:50] = 1
        folds = np.zeros(100, dtype=int)
        folds[y == 1] = 1                     # fold 0 is all-negative, fold 1 all-positive
        models, _ = baseline.fit_fold_models(seqs, y, folds, k=2)
        assert models == {}


class TestOofScores:
    def test_every_row_scored_exactly_once(self, toy):
        seqs, y, folds = toy
        s, _, _ = baseline.oof_scores(seqs, y, folds, k=K)
        assert np.isfinite(s).all()

    def test_recovers_a_real_motif(self, toy):
        seqs, y, folds = toy
        s, _, _ = baseline.oof_scores(seqs, y, folds, k=K)
        assert roc_auc_score(y, s) > 0.9

    def test_is_chance_on_shuffled_labels(self):
        """The out-of-fold AUROC must collapse when the labels carry no information.

        An in-sample fit on 1024 5-mer features and 500 rows would score near-perfectly
        here, so this test is what distinguishes out-of-fold from in-sample.
        """
        rng = np.random.default_rng(11)
        seqs = random_seqs(500, rng)
        y = rng.integers(0, 2, 500)
        folds = np.tile(np.arange(K), 100)
        s, _, _ = baseline.oof_scores(seqs, y, folds, k=K)
        assert 0.35 < roc_auc_score(y, s) < 0.65

    def test_out_of_fold_is_lower_than_in_sample(self, toy):
        from sklearn.linear_model import LogisticRegression
        seqs, y, folds = toy
        s, _, _ = baseline.oof_scores(seqs, y, folds, k=K)
        X, _ = baseline.kmer_matrix(seqs, K)
        ins = LogisticRegression(max_iter=3000).fit(X, y).decision_function(X)
        assert roc_auc_score(y, s) <= roc_auc_score(y, ins)


class TestVariantDelta:
    def test_identical_sequences_give_zero_delta(self, toy):
        seqs, y, folds = toy
        models, vec = baseline.fit_fold_models(seqs, y, folds, k=K)
        d = baseline.variant_delta(models, vec, seqs[:20], seqs[:20], folds[:20])
        assert np.allclose(d, 0.0)

    def test_destroying_the_motif_gives_a_positive_delta(self, toy):
        """score(reference) - score(alternate) must be positive when the alternate
        allele removes the motif."""
        seqs, y, folds = toy
        models, vec = baseline.fit_fold_models(seqs, y, folds, k=K)
        refs = [s for s in seqs[:120] if "GCAUG" in s]
        alts = [s.replace("GCAUG", "AAAAA", 1) for s in refs]
        d = baseline.variant_delta(models, vec, refs, alts,
                                   np.zeros(len(refs), dtype=int))
        assert np.nanmean(d) > 0

    def test_creating_the_motif_gives_a_negative_delta(self, toy):
        seqs, y, folds = toy
        models, vec = baseline.fit_fold_models(seqs, y, folds, k=K)
        refs = [s for s in seqs[250:350] if "GCAUG" not in s][:40]
        alts = [s[:20] + "GCAUG" + s[25:] for s in refs]
        d = baseline.variant_delta(models, vec, refs, alts,
                                  np.zeros(len(refs), dtype=int))
        assert np.nanmean(d) < 0

    def test_variants_with_no_fold_model_come_back_nan(self, toy):
        seqs, y, folds = toy
        models, vec = baseline.fit_fold_models(seqs, y, folds, k=K)
        d = baseline.variant_delta(models, vec, seqs[:5], seqs[:5],
                                  np.full(5, 99))          # fold 99 has no model
        assert np.isnan(d).all()

    def test_nan_fold_is_not_scored(self, toy):
        seqs, y, folds = toy
        models, vec = baseline.fit_fold_models(seqs, y, folds, k=K)
        f = np.array([0.0, np.nan, 1.0, np.nan, 2.0])
        d = baseline.variant_delta(models, vec, seqs[:5], seqs[:5], f)
        assert np.isnan(d[[1, 3]]).all() and np.isfinite(d[[0, 2, 4]]).all()

    def test_each_variant_uses_its_own_fold_model(self, toy):
        """Two variants with identical sequences but different folds should generally get
        different deltas, because different models score them."""
        seqs, y, folds = toy
        models, vec = baseline.fit_fold_models(seqs, y, folds, k=K)
        refs = [seqs[0]] * K
        alts = [seqs[0].replace("GCAUG", "AAAAA", 1)] * K
        d = baseline.variant_delta(models, vec, refs, alts, np.arange(K))
        assert len(set(np.round(d, 9))) > 1


class TestEvaluate:
    def test_reports_auroc_with_an_interval(self, toy):
        seqs, y, folds = toy
        df = pd.DataFrame({"seq_rna": seqs, "label": y, "fold": folds})
        r = baseline.evaluate(df, k=K)
        assert r["ci_low"] < r["auroc"] < r["ci_high"]

    def test_n_counts_only_scored_rows(self, toy):
        seqs, y, folds = toy
        df = pd.DataFrame({"seq_rna": seqs, "label": y, "fold": folds})
        assert baseline.evaluate(df, k=K)["n"] == len(df)

    def test_sweep_covers_every_k(self, toy):
        seqs, y, folds = toy
        df = pd.DataFrame({"seq_rna": seqs, "label": y, "fold": folds})
        out = baseline.sweep_k(df, ks=(3, 4))
        assert set(out) == {3, 4} and all("auroc" in v for v in out.values())
