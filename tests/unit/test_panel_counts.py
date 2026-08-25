"""Pin the four panel sizes and their nesting, so 95/187/189 can never drift again.

WHY THIS FILE EXISTS. Four legitimately different panels, each with one cause, got referred
to by their sizes rather than their roles, and the result read like an inconsistency every
time it came up. docs/PANELS.md is the prose; this is the executable version. If a number
here changes, either the data changed or something is wrong, and either way somebody has to
look.

These are integration-flavoured: they read the committed result tables rather than mocking,
because the thing being asserted IS the shape of the committed results.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
T = ROOT / "results" / "tables"

FULL, MATCHED, DEEP, VARIANT = 189, 187, 95, 94


def _load(name):
    p = T / name
    if not p.exists():
        pytest.skip(f"{name} not built in this checkout")
    return pd.read_csv(p)


def _sets():
    full = set(_load("rehearsal_binding_dinuc.csv").dataset)
    matched = set(_load("rehearsal_binding_gc.csv").dataset)
    deep = set(_load("matched95_four_models.csv").dataset)
    v = _load("variant_scores_splicebert.csv")
    variant = set(v.protein + ":" + v.cell)
    return full, matched, deep, variant


class TestSizes:
    def test_the_four_counts(self):
        full, matched, deep, variant = _sets()
        assert (len(full), len(matched), len(deep), len(variant)) == \
            (FULL, MATCHED, DEEP, VARIANT)

    def test_no_duplicate_datasets_anywhere(self):
        for name in ("rehearsal_binding_dinuc.csv", "rehearsal_binding_gc.csv",
                     "matched95_four_models.csv", "locality_ism.csv"):
            d = _load(name)
            assert not d.dataset.duplicated().any(), f"{name} has duplicate datasets"


class TestNesting:
    """Every panel is a subset of FULL. A stray dataset would mean a panel was built from
    a different source than the one its name implies."""

    def test_matched_subset_of_full(self):
        full, matched, _, _ = _sets()
        assert matched <= full

    def test_deep_subset_of_full(self):
        full, _, deep, _ = _sets()
        assert deep <= full

    def test_variant_subset_of_deep(self):
        _, _, deep, variant = _sets()
        assert variant <= deep


class TestTheDifferencesHaveTheDocumentedCause:
    def test_the_two_datasets_missing_from_the_gc_arm(self):
        """They fall below min_pairs under GC matching, and only under GC matching."""
        full, matched, _, _ = _sets()
        assert full - matched == {"DDX51:K562", "NCBP2:K562"}

    def test_those_two_clear_the_floor_under_dinucleotide_matching(self):
        """The stricter control KEPT datasets the looser one dropped, which is counter-
        intuitive enough that it gets asserted rather than remembered."""
        d = _load("rehearsal_binding_dinuc.csv").set_index("dataset")
        for ds in ("DDX51:K562", "NCBP2:K562"):
            assert d.loc[ds, "pairs"] >= 400

    def test_the_dataset_missing_from_the_variant_panel(self):
        _, _, deep, variant = _sets()
        assert deep - variant == {"NCBP2:K562"}

    def test_deep_is_a_systematic_sample_not_a_size_threshold(self):
        """A size-thresholded subset would confound the panel with dataset size, which
        correlates with AUROC at r = +0.53 to +0.67. Systematic sampling by pair rank does
        not, so DEEP must span essentially the same size range as FULL."""
        full_d = _load("rehearsal_binding_dinuc.csv")
        deep_d = _load("matched95_four_models.csv")
        assert deep_d.pairs.min() <= full_d.pairs.quantile(0.05)
        assert deep_d.pairs.max() >= full_d.pairs.quantile(0.95)


class TestClaimsUseTheRightPanel:
    def test_cost_of_matching_uses_matched(self):
        assert len(_load("cost_of_matching.csv")) == MATCHED

    def test_locality_uses_deep(self):
        assert len(_load("locality_ism.csv")) == DEEP

    def test_variant_arms_are_like_for_like(self):
        """Both score models must be reported on the SAME variants, or the comparison is
        model confounded with panel."""
        r = _load("variant_results_splicebert.csv")
        assert r.n.nunique() == 1
        assert r.n_pathogenic.nunique() == 1
