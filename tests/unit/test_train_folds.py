"""Fold-based splitting for the training path.

An earlier revision of the trainer read only the frozen `split` column, so a GPU sweep would
have been measured under a different protocol from the composition control it is compared
against. These tests pin the cross-validation path that params.yaml calls primary.

No torch import here: split_frame is pure pandas, and the CPU image has no torch.
"""

import pandas as pd
import pytest

from rbp.data.splits import split_of_fold
from rbp.train.data import split_frame

K = 5


def frame(n_per_fold=6):
    """Rows tagged with both protocols, so the two can be told apart."""
    rows = []
    for f in range(K):
        for j in range(n_per_fold):
            rows.append({"id": f"r{f}_{j}", "label": j % 2, "fold": f,
                         "split": ("test" if f == 0 else "val" if f == 1 else "train"),
                         "seq_rna": "ACGU" * 25})
    return pd.DataFrame(rows)


class TestLegacySplit:
    def test_fold_none_uses_the_split_column(self):
        df = frame()
        assert set(split_frame(df, "test").fold) == {0}
        assert set(split_frame(df, "val").fold) == {1}
        assert set(split_frame(df, "train").fold) == {2, 3, 4}

    def test_legacy_and_fold_paths_disagree(self):
        """If these ever matched, the test below would prove nothing."""
        df = frame()
        assert set(split_frame(df, "test").id) != set(split_frame(df, "test", fold=3).id)


class TestFoldSplit:
    @pytest.mark.parametrize("fold", range(K))
    def test_roles_follow_split_of_fold(self, fold):
        df = frame()
        test, val, train = (split_frame(df, s, fold=fold, k=K) for s in
                            ("test", "val", "train"))
        assert set(test.fold) == {fold}
        assert set(val.fold) == {(fold + 1) % K}
        assert set(train.fold) == {f for f in range(K) if f not in (fold, (fold + 1) % K)}

    @pytest.mark.parametrize("fold", range(K))
    def test_the_three_roles_partition_the_data(self, fold):
        df = frame()
        got = pd.concat([split_frame(df, s, fold=fold, k=K)
                         for s in ("train", "val", "test")])
        assert len(got) == len(df)
        assert set(got.id) == set(df.id)

    def test_every_row_is_held_out_exactly_once_across_all_folds(self):
        """The property the whole protocol rests on. If a row were held out twice, or
        never, the pooled out-of-fold AUROC would double-count it or silently drop it."""
        df = frame()
        seen = []
        for fold in range(K):
            seen += split_frame(df, "test", fold=fold, k=K).id.tolist()
        assert len(seen) == len(df)
        assert set(seen) == set(df.id)

    def test_a_row_is_never_in_train_and_test_at_once(self):
        df = frame()
        for fold in range(K):
            tr = set(split_frame(df, "train", fold=fold, k=K).id)
            te = set(split_frame(df, "test", fold=fold, k=K).id)
            va = set(split_frame(df, "val", fold=fold, k=K).id)
            assert not (tr & te) and not (tr & va) and not (te & va)

    def test_validation_is_the_next_fold_cyclically(self):
        assert split_of_fold(0, K - 1, K) == "val"     # wraps around

    def test_empty_when_a_fold_has_no_rows(self):
        df = frame()
        df = df[df.fold != 2]
        assert len(split_frame(df, "test", fold=2, k=K)) == 0
