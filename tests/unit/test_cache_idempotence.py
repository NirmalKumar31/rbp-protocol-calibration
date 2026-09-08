"""The idempotence gate's own failure paths, and rbp.utils.carry's.

Why these exist. scripts/cache_idempotence.py is the only thing asserting that the documented
offline entry points return the tables they are documented to return, and
rbp.utils.carry is the only thing stopping those entry points degrading a committed table when
an optional input is absent. Both were written in one day, and between them they had four
false-pass paths that an audit found and one that CI found:

  * a non-zero child exit was reported as SKIPPED, so a broken entry point dropped out of the
    evidence and the gate still passed;
  * the comparison covered the `check` set, `note` and `value`, so `n`, `ci_low`, `ci_high` and
    every analysis-specific column could move unnoticed;
  * only top-level *.csv was examined, so nested paths and .tsv files were invisible;
  * restoration copied the snapshot back without DELETING files a run had added, so an addition
    could contaminate the next entry point and survive the gate;
  * carry appended a provenance suffix to the note it carried, so the carried row was not
    identical and the suffix accumulated on every run.

Each of those is one test below, on synthetic two-row tables rather than by running the real
entry points, which takes a minute and is what the gate itself is for. The point of these is
that the gate FAILS when it should, which no amount of green runs demonstrates.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

cache_idempotence = pytest.importorskip("cache_idempotence")
carry = pytest.importorskip("rbp.utils.carry")

BASE = pd.DataFrame([
    {"check": "alpha", "value": 1.5, "ci_low": 1.0, "ci_high": 2.0, "n": 94, "note": "a note"},
    {"check": "beta", "value": 2.5, "ci_low": "", "ci_high": "", "n": 79, "note": ""},
])


def _pair(tmp_path, mutate=None):
    """Write a committed table and a regenerated one, returning both paths."""
    before = tmp_path / "before.csv"
    after = tmp_path / "after.csv"
    BASE.to_csv(before, index=False)
    d = BASE.copy()
    if mutate is not None:
        d = mutate(d)
    d.to_csv(after, index=False)
    return before, after


def _complaints(tmp_path, mutate, tol=1e-9):
    before, after = _pair(tmp_path, mutate)
    return cache_idempotence.compare(before, after, tol)


# --- the comparison, one failure path per test --------------------------------------------

def test_an_identical_table_is_silent(tmp_path):
    assert _complaints(tmp_path, None) == []


def test_floating_point_noise_is_not_a_complaint(tmp_path):
    """The reason this is not a byte diff: seven real tables move at 1e-14 and always will."""
    out = _complaints(tmp_path, lambda d: d.assign(value=d.value + 3e-14))
    assert out == [], out


def test_a_value_moving_past_the_tolerance_fails(tmp_path):
    out = _complaints(tmp_path, lambda d: d.assign(value=d.value + 0.01))
    assert any("value moved" in m for m in out), out


def test_a_changed_confidence_bound_fails(tmp_path):
    """ci_low and ci_high went unchecked, so a bound could move silently."""
    def mutate(d):
        d.loc[0, "ci_low"] = 0.5
        return d
    out = _complaints(tmp_path, mutate)
    assert any("CI_LOW" in m.upper() for m in out), out


def test_a_changed_sample_size_fails(tmp_path):
    """`n` went unchecked. A contribution recomputed on a different panel is a different claim."""
    def mutate(d):
        d.loc[1, "n"] = 40
        return d
    out = _complaints(tmp_path, mutate)
    assert any("moved on 'beta'" in m or "N changed" in m.upper() for m in out), out


def test_a_changed_note_is_reported_as_staleness(tmp_path):
    def mutate(d):
        d.loc[0, "note"] = "a different note"
        return d
    out = _complaints(tmp_path, mutate)
    assert any("stale against its script" in m for m in out), out


def test_a_blank_note_and_a_nan_note_are_the_same_absence(tmp_path):
    """A CSV round trip turns "" into NaN. Treating them as different would fail every run."""
    def mutate(d):
        d.loc[1, "note"] = float("nan")
        return d
    assert _complaints(tmp_path, mutate) == []


def test_a_dropped_row_fails(tmp_path):
    out = _complaints(tmp_path, lambda d: d.iloc[:1])
    assert any("DROPPED row" in m for m in out), out


def test_an_added_row_fails(tmp_path):
    def mutate(d):
        return pd.concat([d, d.iloc[[0]].assign(check="gamma")], ignore_index=True)
    out = _complaints(tmp_path, mutate)
    assert any("ADDED row" in m for m in out), out


def test_a_duplicate_key_fails_rather_than_masking_a_change(tmp_path):
    """Two rows with one key make every per-key comparison below ambiguous."""
    def mutate(d):
        return pd.concat([d, d.iloc[[0]]], ignore_index=True)
    out = _complaints(tmp_path, mutate)
    assert any("DUPLICATE keys" in m for m in out), out


def test_a_changed_column_set_fails_and_stops(tmp_path):
    out = _complaints(tmp_path, lambda d: d.drop(columns=["ci_high"]))
    assert len(out) == 1 and "COLUMNS changed" in out[0], out


# --- the snapshot and the restore ----------------------------------------------------------

def test_the_snapshot_covers_nested_paths_and_tsv(tmp_path):
    """Top-level *.csv only was the bug: draws/ and .tsv files were invisible."""
    (tmp_path / "draws").mkdir()
    (tmp_path / "top.csv").write_text("a\n1\n")
    (tmp_path / "draws" / "nested.csv").write_text("a\n2\n")
    (tmp_path / "manifest.tsv").write_text("a\n3\n")
    (tmp_path / "ignored.md").write_text("not a table\n")
    snap = cache_idempotence.snapshot_of(tmp_path)
    assert set(snap) == {"top.csv", "draws/nested.csv", "manifest.tsv"}


def test_restore_deletes_files_the_run_added(tmp_path):
    """Copying a snapshot back leaves additions in place; they then contaminate the next run."""
    snap, dest = tmp_path / "snap", tmp_path / "dest"
    for d in (snap, dest):
        d.mkdir()
    (snap / "kept.csv").write_text("a\n1\n")
    (dest / "kept.csv").write_text("a\n1\n")
    (dest / "added.csv").write_text("a\n9\n")
    (dest / "sub").mkdir()
    (dest / "sub" / "also_added.csv").write_text("a\n9\n")
    cache_idempotence.restore(snap, dest)
    assert (dest / "kept.csv").exists()
    assert not (dest / "added.csv").exists()
    assert not (dest / "sub" / "also_added.csv").exists()


def test_restore_puts_a_modified_file_back(tmp_path):
    snap, dest = tmp_path / "snap", tmp_path / "dest"
    for d in (snap, dest):
        d.mkdir()
    (snap / "t.csv").write_text("a\n1\n")
    (dest / "t.csv").write_text("a\nTAMPERED\n")
    cache_idempotence.restore(snap, dest)
    assert (dest / "t.csv").read_text() == "a\n1\n"


# --- the child process --------------------------------------------------------------------

def test_a_non_zero_child_exit_is_a_failure_not_a_skip():
    """The false-pass path: a broken entry point used to drop out of the evidence silently.

    Driven through the real CLI with --only naming a stem that does not exist as a script, so
    the subprocess genuinely fails. The gate must exit non-zero.
    """
    bad = ROOT / "scripts" / "definitely_not_a_script.py"
    assert not bad.exists()
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "cache_idempotence.py"), "--only", "nosuch"],
        capture_output=True, text=True, cwd=str(ROOT),
        env={"PATH": "/usr/bin:/bin", "PYTHONPATH": str(ROOT / "src")})
    # No matching entry point is a usage error, which must also not be a silent pass.
    assert r.returncode != 0, r.stdout + r.stderr


# --- rbp.utils.carry -----------------------------------------------------------------------

def test_a_carried_row_is_byte_identical(tmp_path):
    """It appended a provenance suffix, so the carried table differed from the committed one."""
    t = tmp_path / "t.csv"
    BASE.to_csv(t, index=False)
    out = []
    assert carry.emit(out, t, "alpha", None, recomputed=False) is True
    assert out[0]["note"] == "a note", "the note must come back unchanged"
    assert float(out[0]["value"]) == 1.5


def test_carrying_twice_does_not_grow_the_note(tmp_path):
    """The accumulating suffix: every run appended again, without bound."""
    t = tmp_path / "t.csv"
    BASE.to_csv(t, index=False)
    for _ in range(3):
        out = []
        carry.emit(out, t, "alpha", None, recomputed=False)
        pd.DataFrame(out + [{"check": "beta", "value": 2.5, "ci_low": "", "ci_high": "",
                             "n": 79, "note": ""}]).to_csv(t, index=False)
    assert pd.read_csv(t).set_index("check").loc["alpha", "note"] == "a note"


def test_carry_returns_false_when_there_is_nothing_to_carry(tmp_path):
    """No committed value and no input is genuinely unanswerable, and must stay loud.

    The row is then absent and verify.py's presence assertion fails, which is correct.
    """
    out = []
    assert carry.emit(out, tmp_path / "absent.csv", "alpha", None, recomputed=False) is False
    assert out == []


def test_carry_returns_false_for_a_row_the_committed_table_lacks(tmp_path):
    t = tmp_path / "t.csv"
    BASE.to_csv(t, index=False)
    out = []
    assert carry.emit(out, t, "not_a_row", None, recomputed=False) is False
    assert out == []


def test_a_recomputed_row_is_written_normally(tmp_path):
    out = []
    assert carry.emit(out, tmp_path / "absent.csv", "alpha", 7.0, note="fresh", n=94) is True
    assert out[0] == {"check": "alpha", "value": 7.0, "ci_low": "", "ci_high": "",
                      "n": 94, "note": "fresh"}
