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


# --------------------------------------------------------------- how the negatives are built

@pytest.mark.parametrize("arm,expected", [("gc", 0.95), ("dn", 0.85), ("neg2", 0.20)])
def test_match_quality_is_what_the_construction_view_claims(arm, expected):
    """The view states 95%, 85% and 20% matched within 0.05 GC. Those drive the whole story.

    The bias-aware figure being low is not a defect in that protocol, it is the protocol: it
    corrects for assay bias instead of matching composition. But it is the reason its baseline
    is highest and its measured contribution smallest, so if the number moves the explanation
    on that view is wrong.
    """
    got = raw_value("match_quality.csv", f"fraction of pairs within |dGC| 0.05, {arm} arm")
    assert abs(got - expected) < 0.02, f"{arm} matched {got:.3f}, the view says about {expected}"


def test_the_bias_aware_arm_is_the_least_composition_matched():
    """The ordering is the mechanism. It must survive a regeneration or the prose is wrong."""
    fracs = {a: raw_value("match_quality.csv", f"fraction of pairs within |dGC| 0.05, {a} arm")
             for a in ARMS}
    assert min(fracs, key=fracs.get) == "neg2"
    assert max(fracs, key=fracs.get) == "gc"


def test_dinucleotide_matching_is_strictly_the_harder_constraint():
    """The view claims dn achieves a tighter GC gap AND a better dinucleotide match than gc.

    That is a strong claim, and it is what makes the arm a fair comparison rather than a looser
    one dressed up. Both halves are checked.
    """
    gc_gap_dn = raw_value("match_quality.csv", "gc_gap_median, dn arm")
    gc_gap_gc = raw_value("match_quality.csv", "gc_gap_median, gc arm")
    assert gc_gap_dn < gc_gap_gc, "dinucleotide matching no longer achieves a tighter GC gap"
    improvement = raw_value("match_quality.csv",
                            "dinucleotide L1 improvement factor, dn vs gc arm")
    assert improvement > 2.0, f"the dinucleotide improvement fell to {improvement:.2f}"


def test_stricter_matching_lowers_the_raw_score_on_most_datasets():
    """The cost chart's annotation counts datasets below zero. Recomputed here independently."""
    import csv as _csv

    with (data.TABLES / "cost_of_matching.csv").open() as handle:
        rows = list(_csv.DictReader(handle))
    lower = sum(1 for r in rows if float(r["cost"]) < 0)
    assert lower > len(rows) * 0.75, (
        f"only {lower} of {len(rows)} datasets score lower under stricter matching; the "
        "construction view says almost all of them do")


def test_the_result_does_not_depend_on_one_random_draw():
    """Five seeds. The view says they land inside one standard error; check that they do."""
    seeds = [7, 11, 23, 42, 101]
    vals = [raw_value("negative_draws.csv", f"panel-mean contribution, seed {s}") for s in seeds]
    se = raw_value("negative_draws.csv",
                   "between-protein standard error, the published draw")
    mean = sum(vals) / len(vals)
    assert max(abs(v - mean) for v in vals) < se, (
        "a redraw now falls outside one standard error of the others")
    ratio = raw_value("negative_draws.csv", "ratio of combined to published interval width")
    assert ratio < 1.10, f"draw uncertainty now widens the interval by {(ratio - 1) * 100:.1f}%"


# ---------------------------------------------------------------------- why the study exists

def test_no_surveyed_method_reports_a_composition_baseline():
    """Zero of seven. That number is the entire justification for the paper.

    If a regenerated survey ever finds one, the motivating claim on the opening view stops being
    true and must be rewritten rather than quietly left standing.
    """
    surveyed = raw_value("negative_set_survey.csv", "methods and benchmarks surveyed")
    reporting = raw_value("negative_set_survey.csv",
                          "surveyed sources reporting a composition-only baseline")
    assert surveyed == 7.0
    assert reporting == 0.0, f"{reporting:.0f} surveyed sources now report a baseline"

    import csv as _csv
    with (data.TABLES / "negative_set_survey_per_method.csv").open() as handle:
        rows = list(_csv.DictReader(handle))
    assert len(rows) == 7
    assert all(r["composition_baseline"] == "False" for r in rows)
    assert all(r["url"].startswith("http") for r in rows), "a source lost its citation"
    assert all(len(r["quote"]) > 40 for r in rows), "a source lost its supporting quote"


def test_the_variant_ladder_still_shows_the_controls_beating_the_model():
    """The pivot story. Conservation with no model beats the model, and a wrong protein clears
    chance by a wide margin. Both are the reason the earlier line of work stopped.
    """
    import csv as _csv

    with (data.TABLES / "variant_ladder.csv").open() as handle:
        arms = {r["arm"]: float(r["auroc"]) for r in _csv.DictReader(handle)}
    assert arms["conservation"] > arms["matched"], (
        "conservation no longer beats the matched model; the opening narrative is wrong")
    assert arms["mismatched"] > 0.60, "the mismatched-protein control no longer clears chance"
    assert arms["matched"] > arms["mismatched"] > arms["kmer"] > 0.5, "the ladder order changed"
    assert round(arms["conservation"], 3) == 0.908
    assert round(arms["matched"], 3) == 0.829


# ------------------------------------------------------------------ the recommendation's test

def test_the_recommendation_improves_every_pair_but_clears_zero_on_none():
    """Both halves matter. The direction is consistent; the effect is not established.

    Stating only the first half would be advertising, and stating only the second would throw
    away a real signal. The view says both, so both are pinned.
    """
    pairs = ["gc vs dn", "gc vs neg2", "dn vs neg2"]
    assert raw_value("recommendation_works.csv",
                     "protocol pairs where rank agreement improves") == 3.0
    for pair in pairs:
        assert raw_value("recommendation_works.csv", f"rank agreement gain, {pair}") > 0, (
            f"{pair} no longer improves; the recommendation's direction has changed")

    cleared = 0
    for pair in pairs:
        _, low, _ = data.interval(
            "recommendation_works.csv",
            f"rank agreement gain, {pair}, Bonferroni over 3 pairs")
        if low is not None and low > 0:
            cleared += 1
    assert cleared == 0, (
        f"{cleared} pair(s) now clear zero after Bonferroni correction; the view claims none do")


def test_exactly_one_pair_clears_zero_before_correction():
    """The uncorrected picture, stated separately so the correction's effect is visible."""
    cleared = []
    for pair in ("gc vs dn", "gc vs neg2", "dn vs neg2"):
        _, low, _ = data.interval("recommendation_works.csv", f"rank agreement gain, {pair}")
        if low is not None and low > 0:
            cleared.append(pair)
    assert cleared == ["gc vs dn"], (
        f"uncorrected, {cleared} clear zero; the view says only the GC/dinucleotide pair does")


def test_between_dataset_spread_exceeds_between_protocol_spread():
    """The site claimed the opposite for several releases. It is not true.

    The Map view said "columns differ from each other more than rows do, which means the
    protocol matters more than the dataset". Computed from the committed per-dataset table,
    between-dataset variance is about 1.6x the between-protocol variance, so the dataset is the
    larger source of spread.

    Nothing in the paper depended on the false version: the published claim is that the PANEL
    MEAN moves 4.84-fold, which is a statement about systematic shift and not about variance
    share. Both hold at once. This pins the direction so the stronger, wrong sentence cannot
    come back.
    """
    import statistics

    frame = data.table("three_arm_per_dataset.csv")
    columns = {arm: frame[f"gain_{arm}"].dropna().tolist() for arm in ARMS}
    n = min(len(v) for v in columns.values())
    rows = [[columns[arm][i] for arm in ARMS] for i in range(n)]

    between_protocol = statistics.pvariance([statistics.mean(columns[a]) for a in ARMS])
    between_dataset = statistics.pvariance([statistics.mean(r) for r in rows])

    assert between_dataset > between_protocol, (
        "between-protocol spread now exceeds between-dataset spread; the wording on the Map "
        "view was written for the opposite and must be revisited")
    ratio = between_dataset / between_protocol
    assert 1.2 < ratio < 2.2, f"the ratio moved to {ratio:.2f}; the view states about 1.6"


def test_the_trajectory_counter_matches_the_committed_table():
    """The chart prints "59 of 94 follow the panel ordering". Recomputed here.

    A count rendered into an annotation is exactly the kind of number that survives a data
    change unnoticed, because the chart still draws.
    """
    frame = data.table("three_arm_per_dataset.csv").dropna(
        subset=[f"gain_{a}" for a in ARMS])
    follows = sum(
        1 for _, r in frame.iterrows()
        if r["gain_dn"] > r["gain_gc"] > r["gain_neg2"]
    )
    assert len(frame) == 94
    assert follows == 59, f"{follows} of 94 now follow the ordering; the guide says 59"
    assert follows < len(frame), (
        "if every dataset followed the ordering the guide's second check, which says 35 do not, "
        "would be wrong")
