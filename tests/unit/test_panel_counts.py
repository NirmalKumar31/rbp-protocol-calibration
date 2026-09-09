"""Pin this pipeline's panel sizes and their nesting, so 94/95/189 can never drift again.

Why this file exists. Legitimately different panels, each with one cause, got referred to by
their sizes rather than their roles, and the result read like an inconsistency every time it
came up. docs/PANELS.md is the prose; this is the executable version. If a number here
changes, either the data changed or something is wrong, and either way somebody has to look.

These are this pipeline's numbers, not the earlier study's. The file previously asserted a
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
#
# VARIANT read 94 and was wrong, and nothing found out because the table it was checked
# against, variant_scores_splicebert.csv, is produced by a full pipeline run but is not
# committed: four tests here skipped on every run and certified nothing. docs/PANELS.md, which
# declares itself the single source for these counts, said 95. The released
# variant_assignments.csv settles it: all 95 study datasets carry ClinVar assignments,
# NCBP2:K562 among them. The doc was right and the test was wrong.
STUDY, MATCHED, VARIANT = 95, 94, 95

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
    # variant_assignments.csv IS in the release and covers the whole candidate pool, so the
    # variant panel is its intersection with the study panel. The superseded
    # variant_scores_splicebert.csv is not shipped, and reading it made every test below skip.
    v = _load("variant_assignments.csv")
    variant = set(v.protein.astype(str) + ":" + v.cell.astype(str)) & study
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

    def test_no_study_dataset_is_missing_from_the_variant_panel(self):
        """THE_FALLER drops out of the GC arm, not out of the variant arm. This asserted the
        opposite while skipping, so the mistake cost nothing until it was read."""
        study, _, _, variant = _sets()
        assert study - variant == set()

    def test_study_is_a_systematic_sample_not_a_size_threshold(self):
        """A size-thresholded subset would confound the panel with dataset size, which
        correlates with AUROC at r = +0.53 to +0.67. Systematic sampling by pair rank does
        not, so the study panel must span essentially the whole size range it was drawn
        from. scripts/select_panel.py asserts this at selection time; this is the same
        check on the committed result.

        THE POPULATION HAS TO BE ON THE OTHER SIDE OF THE COMPARISON, and for as long as
        this test has existed it was not. It read sweep_dinuc.csv, which holds the 95
        SELECTED datasets and has never carried a `pairs` column in any commit, so the
        guard below fired on every checkout that has ever existed and the assertions never
        ran once. Had the column been there the test would have compared the panel against
        its own quantiles, and min <= q05 is true of any set whatsoever, so it would have
        passed while certifying nothing. candidate_sizes.csv is the 189-dataset pool the
        panel was drawn from, it is committed, and it is exactly what select_panel.py calls
        `full` when it makes this decision.
        """
        pool = _load("candidate_sizes.csv")
        study_d = _load("rehearsal_binding_dinuc.csv")
        assert set(study_d.dataset) < set(pool.dataset), "the panel is a strict subset"
        assert study_d.pairs.min() <= pool.pairs.quantile(0.05)
        assert study_d.pairs.max() >= pool.pairs.quantile(0.95)


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
