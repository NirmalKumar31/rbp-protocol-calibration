"""Every arm's membership must be readable from the repository, not from a private disk.

WHY. An external review reported that the frozen panel artefact is not in the release. It was
wrong twice and right once, and all three are worth writing down because the shape recurs.

Wrong the first time: `manifest/study_panel.tsv` is a GCS object key, not a path. The reviewer
looked for it on disk, did not find it, and called the panel missing.

Wrong the second time: the 94-dataset study panel IS committed, as
`results/tables/supplementary_table_s1.csv`, whose `in_three_arm_panel` column is exactly the
inclusion flag the review asked for, alongside the ENCODE accession and experiment for each
row. Looking for one filename and concluding the artefact does not exist is the same error as
grepping for one project id and concluding no project id is hardcoded.

Right, underneath both: `config/panel_final_{cell}_{arm}.tsv` existed for `gc` and `dinuc` and
for neither bias-aware arm. Their membership came from listing directories in `../rbp-store`,
which is 2.9 GB, uncommitted, and on one laptop. The study's third and fourth arms were defined
by something no reader could see.

WHAT THIS CHECKS. That the panels exist for every arm, and that they agree with the per-window
out-of-fold scores committed under `data/evidence/` -- which are what the published AUROCs are
recomputed from, so a disagreement means one of the two is wrong. Checking a panel against the
store would prove nothing to anyone without the store, which is the problem being fixed.

NOTE THE TWO DIFFERENT PANELS, because conflating them is this project's oldest recurring bug.
`config/panel_final_*_{gc,dinuc}.tsv` are the CANDIDATE panels: every dataset that cleared the
pair floor in that arm, 187 of them. The STUDY panel is the systematic subsample of 95, of
which 94 carry all three arms. The bias-aware arms were only ever built for the study panel, so
their panel files hold 94 and the composition-matched ones hold 187. Both are correct; they
answer different questions. docs/PANELS.md is the long version.
"""

import csv
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
CELLS = ("K562", "HepG2")
S1 = ROOT / "results" / "tables" / "supplementary_table_s1.csv"

# Arms whose panel file is the STUDY panel, with the per-dataset score directories that must
# agree with it. dinuc's evidence is stored per protein PAIR rather than per dataset, so it has
# no directory listing to compare against and is checked for existence only.
STUDY_ARMS = {"neg2": "scores_neg2", "neg2_rm": None}
CANDIDATE_ARMS = ("gc", "dinuc")

EXPECTED_STUDY = 94
EXPECTED_PAIRS = 456_734


def _panel(cell, arm):
    f = ROOT / "config" / f"panel_final_{cell}_{arm}.tsv"
    if not f.exists():
        return None
    return {ln.split("\t")[0]: int(ln.split("\t")[2])
            for ln in f.read_text().strip().splitlines()[1:]}


@pytest.mark.parametrize("arm", sorted(STUDY_ARMS) + list(CANDIDATE_ARMS))
def test_every_arm_has_a_committed_panel(arm):
    for cell in CELLS:
        assert _panel(cell, arm) is not None, (
            f"config/panel_final_{cell}_{arm}.tsv is missing, so this arm's membership is "
            "defined only by whatever directories happen to exist on the author's disk")


@pytest.mark.parametrize("arm", sorted(STUDY_ARMS))
def test_the_study_panel_is_the_size_the_paper_reports(arm):
    total = sum(len(_panel(c, arm) or {}) for c in CELLS)
    assert total == EXPECTED_STUDY, (
        f"{arm}: committed panels hold {total} datasets, the paper reports {EXPECTED_STUDY}")
    pairs = sum(sum((_panel(c, arm) or {}).values()) for c in CELLS)
    assert pairs == EXPECTED_PAIRS, (
        f"{arm}: committed panels hold {pairs:,} pairs, the paper reports {EXPECTED_PAIRS:,}")


@pytest.mark.parametrize("arm", sorted(a for a in STUDY_ARMS if STUDY_ARMS[a]))
def test_the_panel_matches_the_committed_scores(arm):
    root = ROOT / "data" / "evidence" / STUDY_ARMS[arm]
    if not root.exists():
        pytest.skip(f"{root.relative_to(ROOT)} not in this checkout")
    for cell in CELLS:
        panel = _panel(cell, arm)
        if panel is None:
            pytest.skip(f"no committed panel for {cell}/{arm}")
        d = root / cell
        scored = {p.name for p in d.iterdir() if p.is_dir()} if d.exists() else set()
        assert set(panel) == scored, (
            f"{arm}/{cell}: panel and committed scores disagree.\n"
            f"  in the panel but unscored: {sorted(set(panel) - scored)}\n"
            f"  scored but not in the panel: {sorted(scored - set(panel))}")


def test_supplementary_table_s1_is_the_study_panel_of_record():
    """The artefact the review said was absent, with the fields it said were needed."""
    if not S1.exists():
        pytest.skip("supplementary_table_s1.csv not in this checkout")
    rows = list(csv.DictReader(S1.open()))
    for field in ("dataset", "protein", "cell", "accession", "experiment", "pairs",
                  "in_three_arm_panel"):
        assert field in rows[0], f"supplementary_table_s1.csv has no {field} column"
    n = sum(1 for r in rows if r["in_three_arm_panel"] == "True")
    assert n == EXPECTED_STUDY, (
        f"{n} rows flagged in_three_arm_panel, the paper reports {EXPECTED_STUDY}")


def test_s1_and_the_neg2_panel_name_the_same_datasets():
    """Two independently written artefacts, one membership. They have disagreed before."""
    if not S1.exists():
        pytest.skip("supplementary_table_s1.csv not in this checkout")
    s1 = {r["dataset"] for r in csv.DictReader(S1.open())
          if r["in_three_arm_panel"] == "True"}
    panel = {f"{p}:{c}" for c in CELLS for p in (_panel(c, "neg2") or {})}
    assert s1 == panel, (
        f"  in S1 only:    {sorted(s1 - panel)}\n  in the panel only: {sorted(panel - s1)}")
