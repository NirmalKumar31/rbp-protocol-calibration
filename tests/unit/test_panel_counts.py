"""Pin this pipeline's panel sizes and their nesting, so 94/95/189 can never drift again.

WHY THIS FILE EXISTS. Legitimately different panels, each with one cause, got referred to by
their sizes rather than their roles, and the result read like an inconsistency every time it
came up. docs/PANELS.md is the prose; this is the executable version. If a number here
changes, either the data changed or something is wrong, and either way somebody has to look.

THESE ARE THIS PIPELINE'S NUMBERS, NOT THE EARLIER STUDY'S. The file previously asserted a
189-dataset panel losing two datasets to GC matching, which described a different and now
discarded build. Here the study panel is 95 datasets selected from the candidate pool, and
189 is a COUNT OF TASKS rather than of datasets: 95 in the dinucleotide arm plus 94 in the GC
arm, one task per dataset per arm. Conflating the two is what made the old assertions look
like facts about this run.

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

# STUDY   the selected panel, and the dinucleotide arm covers all of it
# MATCHED the GC arm, one smaller because NCBP2:K562 drops below the min_pairs floor
# VARIANT the ClinVar arm, which needs a dataset to have variants near its peaks
STUDY, MATCHED, VARIANT = 95, 94, 94

# What stage 7 submits: one task per dataset per arm. Not a panel size.
STAGE7_TASKS = STUDY + MATCHED

# The floor a dataset must clear in an arm to be analysable in it.
MIN_PAIRS = 400

# The one dataset that clears the floor under dinucleotide matching and not under GC.
THE_FALLER = "NCBP2:K562"


def _load(name):
    p = T / name
    if not p.exists():
        pytest.skip(f"{name} not built in this checkout")
    return pd.read_csv(p)


# cloud_analysis.four_models() writes matched_four_models.csv; the earlier study called the
# same table matched95_four_models.csv. Both are accepted so these tests exercise whichever
# one a checkout actually has, instead of skipping and certifying nothing.
def _load_deep():
    for name in ("matched_four_models.csv", "matched95_four_models.csv"):
        if (T / name).exists():
            return pd.read_csv(T / name)
    pytest.skip("no four-model table built in this checkout")


def _arms():
    """The two rehearsal arms, which are the only panels always present locally."""
    study = set(_load("rehearsal_binding_dinuc.csv").dataset)
    matched = set(_load("rehearsal_binding_gc.csv").dataset)
    return study, matched


def _sets():
    study, matched = _arms()
    deep = set(_load_deep().dataset)
    v = _load("variant_scores_splicebert.csv")
    variant = set(v.protein + ":" + v.cell)
    return study, matched, deep, variant


class TestSizes:
    def test_the_two_rehearsal_arms(self):
        study, matched = _arms()
        assert (len(study), len(matched)) == (STUDY, MATCHED)

    def test_stage7_task_count_is_the_sum_of_the_arms(self):
        """189 is a task count. Asserting it here is what stops it being read as a panel."""
        study, matched = _arms()
        assert len(study) + len(matched) == STAGE7_TASKS

    def test_the_deep_and_variant_panels(self):
        _, _, deep, variant = _sets()
        assert (len(deep), len(variant)) == (STUDY, VARIANT)

    def test_no_duplicate_datasets_anywhere(self):
        for name in ("rehearsal_binding_dinuc.csv", "rehearsal_binding_gc.csv",
                     "locality_ism.csv"):
            d = _load(name)
            assert not d.dataset.duplicated().any(), f"{name} has duplicate datasets"
        deep = _load_deep()
        assert not deep.dataset.duplicated().any(), "four-model table has duplicate datasets"


class TestNesting:
    """Every panel is a subset of the study panel. A stray dataset would mean a panel was
    built from a different source than the one its name implies."""

    def test_matched_subset_of_study(self):
        study, matched = _arms()
        assert matched <= study

    def test_deep_equals_study(self):
        """The deep panel IS the study panel; all four models run on all 95."""
        study, _, deep, _ = _sets()
        assert deep == study

    def test_variant_subset_of_study(self):
        study, _, _, variant = _sets()
        assert variant <= study


class TestTheDifferenceHasTheDocumentedCause:
    def test_exactly_one_dataset_is_missing_from_the_gc_arm(self):
        study, matched = _arms()
        assert study - matched == {THE_FALLER}

    def test_it_clears_the_floor_under_dinucleotide_matching(self):
        """The stricter control KEPT a dataset the looser one dropped, which is counter-
        intuitive enough that it gets asserted rather than remembered. NCBP2:K562 matches
        406 pairs under dinucleotide control and 384 under GC, so it clears 400 in exactly
        one arm. This is also why R1's n is 94 and not 95."""
        d = _load("rehearsal_binding_dinuc.csv").set_index("dataset")
        assert d.loc[THE_FALLER, "pairs"] >= MIN_PAIRS

    def test_nothing_appears_in_the_gc_arm_that_is_absent_from_the_dinuc_arm(self):
        """The arms are built from one candidate list, so GC cannot contain a surprise."""
        study, matched = _arms()
        assert not matched - study

    def test_the_dataset_missing_from_the_variant_panel(self):
        study, _, _, variant = _sets()
        assert study - variant == {THE_FALLER}

    def test_study_is_a_systematic_sample_not_a_size_threshold(self):
        """A size-thresholded subset would confound the panel with dataset size, which
        correlates with AUROC at r = +0.53 to +0.67. Systematic sampling by pair rank does
        not, so the study panel must span essentially the whole size range it was drawn
        from. scripts/select_panel.py asserts this at selection time; this is the check on
        the committed result."""
        sweep = _load("sweep_dinuc.csv")
        study_d = _load("rehearsal_binding_dinuc.csv")
        if "pairs" not in sweep.columns:
            pytest.skip("sweep_dinuc.csv carries no pairs column in this checkout")
        assert study_d.pairs.min() <= sweep.pairs.quantile(0.05)
        assert study_d.pairs.max() >= sweep.pairs.quantile(0.95)


class TestClaimsUseTheRightPanel:
    def test_cost_of_matching_uses_the_matched_arm(self):
        assert len(_load("cost_of_matching.csv")) == MATCHED

    def test_locality_uses_the_study_panel(self):
        assert len(_load("locality_ism.csv")) == STUDY

    def test_variant_arms_are_like_for_like(self):
        """Both score models must be reported on the SAME variants, or the comparison is
        model confounded with panel."""
        r = _load("variant_results_splicebert.csv")
        assert r.n.nunique() == 1
        assert r.n_pathogenic.nunique() == 1
