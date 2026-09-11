"""Every headline the dashboard displays must equal the committed table it claims to come from.

Read with the standard library rather than through `rbp_dashboard.data`, so a bug in the
loader cannot agree with itself. If a table is regenerated and a displayed number does not
follow, these fail.
"""

from __future__ import annotations

import csv

import pytest
from rbp_dashboard import copy as dash_copy
from rbp_dashboard import data

ARMS = ("gc", "dn", "neg2")


def raw(name: str, label: str) -> dict:
    with (data.TABLES / name).open() as handle:
        hits = [r for r in csv.DictReader(handle) if r["check"] == label]
    assert len(hits) == 1, f"{name}: {label!r} matched {len(hits)} rows"
    return hits[0]


def raw_value(name: str, label: str) -> float:
    return float(raw(name, label)["value"])


# --------------------------------------------------------------------- the four Overview tiles

def test_primary_span_is_the_cross_fitted_four_mer_value():
    """4.84x, and it must come from the cross-fitted row, not the two-stage one.

    These two rows differ by 0.58 and sit adjacent in the same file. Quoting the wrong one is
    the single most likely way for this dashboard to misstate the paper.
    """
    label = "4-mer three-arm span, fully cross-fitted"
    assert data.value("cross_fitting.csv", label) == raw_value("cross_fitting.csv", label)
    assert round(data.value("cross_fitting.csv", label), 2) == 4.84


def test_the_two_stage_span_is_shown_as_a_separate_number():
    label = "4-mer three-arm span, as published"
    assert round(raw_value("cross_fitting.csv", label), 2) == 5.42
    assert raw_value("cross_fitting.csv", label) != raw_value(
        "cross_fitting.csv", "4-mer three-arm span, fully cross-fitted")


def test_dataset_count_is_ninety_four_and_the_panel_is_ninety_five():
    """The Overview tile says 94. The panel file has 95 rows. Both statements must stay true."""
    with (data.TABLES / "panel_summary.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 95, "panel_summary.csv is no longer the 95-dataset selected panel"
    in_all_three = [r for r in rows if r["in_both_arms"] == "True"]
    assert len(in_all_three) == 94, "the three-protocol denominator is no longer 94"
    assert "95" in dash_copy.PANEL_CAVEAT and "94" in dash_copy.PANEL_CAVEAT


def test_known_null_bias_reduction_is_at_least_ninety_five_percent():
    """The tile claims >=95%. The smallest of the three arms decides whether that is true."""
    fractions = [
        raw_value("cross_fitting.csv",
                  f"2-mer fraction of the published value removed by cross-fitting, {arm} arm")
        for arm in ARMS
    ]
    assert min(fractions) >= 0.95, f"the >=95% claim no longer holds: {fractions}"
    assert round(min(fractions) * 100) == 96


# --------------------------------------------------------------- protocol sensitivity figures

@pytest.mark.parametrize("arm", ARMS)
def test_protocol_contribution_matches_the_contrast_table(arm):
    label = f"nested contribution, {arm} arm"
    assert data.value("three_arm_contrast.csv", label) == raw_value(
        "three_arm_contrast.csv", label)


def test_the_inverse_relation_holds_in_the_panel():
    """The bias-aware arm has the highest apparent AUROC and the lowest contribution.

    This is the counter-intuitive claim the Protocol view makes in prose. If the tables ever
    stop supporting it, the prose is wrong and must change.
    """
    apparent = {a: raw_value("three_arm_contrast.csv",
                             f"apparent AUROC (composition + score), {a} arm") for a in ARMS}
    gain = {a: raw_value("three_arm_contrast.csv", f"nested contribution, {a} arm") for a in ARMS}
    assert max(apparent, key=apparent.get) == "neg2"
    assert min(gain, key=gain.get) == "neg2"
    assert min(apparent, key=apparent.get) == "dn"
    assert max(gain, key=gain.get) == "dn"


def test_every_three_arm_row_carries_an_interval():
    for arm in ARMS:
        val, low, high = data.interval("three_arm_contrast.csv", f"nested contribution, {arm} arm")
        assert low is not None and high is not None, f"{arm} lost its interval"
        assert low < val < high, f"{arm} value sits outside its own interval"


# ------------------------------------------------------------------------- model comparison

@pytest.mark.parametrize("model,expected", [("kmer", 5.42), ("cnn", 7.42), ("splicebert", 3.72)])
def test_model_spans_match_the_table(model, expected):
    label = f"{model} three-protocol span"
    assert round(raw_value("three_arm_models.csv", label), 2) == expected


def test_the_estimator_caveat_names_the_neural_models_as_exploratory():
    """The spans above are two-stage. The notice shown beside them must say so.

    Without this, a reader compares 7.42x for the CNN against the paper's primary 4.84x and
    concludes the CNN was cross-fitted. It was not.
    """
    text = dash_copy.ESTIMATOR_CAVEAT.lower()
    assert "two-stage" in text
    assert "exploratory" in text
    assert "not cross-fitted" in text or "were not cross-fitted" in text
    assert "cnn" in text and "splicebert" in text


def test_all_three_model_classes_put_their_lowest_contribution_in_the_bias_aware_arm():
    stated = raw_value("three_arm_models.csv",
                       "model classes whose lowest contribution is the bias-aware arm")
    assert stated == 3.0
    for model in ("kmer", "cnn", "splicebert"):
        gains = {a: raw_value("three_arm_models.csv", f"{model} nested contribution, {a} arm")
                 for a in ARMS}
        assert min(gains, key=gains.get) == "neg2", f"{model} no longer bottoms out at neg2"


# ------------------------------------------------------------------------------- cross-fitting

@pytest.mark.parametrize("arm", ARMS)
def test_the_two_mer_truth_is_zero_and_two_stage_reports_above_it(arm):
    """The known-null test only means anything if the two-stage value is positive."""
    published = raw_value("cross_fitting.csv", f"2-mer contribution as published, {arm} arm")
    assert published > 0, f"{arm}: two-stage no longer reports a positive value on the null"
    note = raw("estimator_floor.csv", f"order-2 noise floor, {arm} arm")["note"]
    assert "zero by construction" in note


@pytest.mark.parametrize("arm", ARMS)
def test_cross_fitting_raises_rather_than_lowers_the_four_mer(arm):
    """Stated in prose on the Cross-fitting view. Worth pinning: it is counter-intuitive."""
    published = raw_value("cross_fitting.csv", f"4-mer contribution as published, {arm} arm")
    fitted = raw_value("cross_fitting.csv", f"4-mer contribution fully cross-fitted, {arm} arm")
    assert fitted > published, f"{arm}: cross-fitting no longer raises the 4-mer estimate"


# -------------------------------------------------------------------------- external validation

def test_external_span_and_its_replication_verdict():
    label = "SPAN across their two negative-set constructions"
    val, low, high = data.interval("external_replication.csv", label)
    assert round(val, 2) == 1.69
    assert val > 1.5 and low > 1.2, "the pre-fixed replication criteria no longer pass"
    assert raw_value("external_replication.csv",
                     "protocol verdict for Claim A on a dataset-disjoint external sample") == 1.0


def test_external_panel_is_disjoint_and_sized_as_stated():
    assert raw_value("external_replication.csv", "datasets") == 135.0
    assert raw_value("external_replication.csv", "proteins") == 108.0
    assert raw_value("external_replication.csv", "dataset overlap with our study panel") == 0.0


def test_the_inverse_relation_does_not_replicate_externally():
    """The External view states this failure. The tables must still show it.

    Externally the higher baseline goes with the higher contribution, which is the same
    direction, not the opposite one seen in the panel.
    """
    base_1 = raw_value("external_replication.csv", "composition baseline, negative-1")
    base_2 = raw_value("external_replication.csv", "composition baseline, negative-2")
    gain_1 = raw_value("external_replication.csv",
                       "nested contribution, negative-1 (bias-agnostic)")
    gain_2 = raw_value("external_replication.csv", "nested contribution, negative-2 (bias-aware)")
    assert base_1 > base_2 and gain_1 > gain_2, (
        "the external arms no longer move in the same direction, so the 'did not replicate' "
        "statement on the External view is now wrong")


# ------------------------------------------------------------------------------ reproducibility

def test_provenance_classes_are_the_ones_the_view_explains():
    with (data.TABLES / "PROVENANCE.csv").open() as handle:
        statuses = {r["status"] for r in csv.DictReader(handle)}
    explained = {"raw-reproducible", "evidence-recomputable", "frozen-cache", "frozen-only",
                 "cloud-produced", "unattributed"}
    assert statuses <= explained, f"unexplained provenance class: {statuses - explained}"


def test_the_reproducibility_caveat_refuses_the_stronger_claim():
    text = dash_copy.REPRO_CAVEAT.lower()
    assert "verification" in text
    assert "not complete raw-to-result" in text or "not complete" in text
    assert "bit-for-bit" in text


def test_the_dashboard_does_not_present_itself_as_monitoring_or_cloud_connected():
    text = dash_copy.NOT_MONITORING.lower()
    assert "not production monitoring" in text
    assert "not connected to any cloud" in text


def test_the_not_only_negatives_caveat_is_stated():
    """The headline invites 'only the negatives changed'. That is not quite true and is said so."""
    text = dash_copy.NOT_ONLY_NEGATIVES.lower()
    assert "refitted" in text
    assert "retained positive subsets differ" in text or "positives" in text


@pytest.mark.parametrize("arm", ARMS)
def test_positive_floor_counts_come_from_the_floor_table_not_the_crossfitting_one(arm):
    """Two tables carry a near-identical count of different things. Cite the right one.

    `estimator_floor.csv` counts datasets where the two-stage estimator reports a positive value
    on the known-null, which is 92, 93 and 91. `cross_fitting.csv` counts datasets where
    cross-fitting LOWERS the 2-mer estimate, which is 91, 93 and 91. They agree in two arms out
    of three, so a wrong citation would survive casual checking in two cases out of three. The
    Cross-fitting view quotes the first; this pins which.
    """
    floor = raw_value("estimator_floor.csv", f"datasets with a positive floor, {arm} arm")
    assert floor > 90, f"{arm}: the 'over 90 of 94' statement no longer holds ({floor})"
    lowers = raw_value("cross_fitting.csv",
                       f"2-mer datasets where cross-fitting LOWERS the contribution, {arm} arm")
    assert floor >= lowers - 1, "the two counts have diverged; recheck which one the view quotes"
