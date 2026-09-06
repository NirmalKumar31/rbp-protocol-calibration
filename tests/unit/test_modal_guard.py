"""The cost guard's failure modes, mocked. There were no tests of this file at all.

Why this matters more than its size. `cloud/modal/guard.py` is the only thing between an
unattended GPU sweep and an unbounded bill, and SECURITY.md advertises it as failing closed.
Every defect found in it so far has been of one kind: a state where the guard COULD NOT SEE was
converted, silently, into a state where it saw nothing to worry about.

  * `modal app list` exited non-zero and the returncode was never checked, so an expired token
    read as "no running app" and the guard exited 0.
  * The app-name pattern matched one literal description while the arms that actually ran were
    named differently, so three of four sweeps ran unguarded.
  * `report()` multiplied ONE app's elapsed time by the container cap and called the product an
    upper bound, while each of several concurrent apps carries its own cap.
  * A live app whose `created_at` did not parse became `(id, None)`, and the upper bound summed
    `if s`, so that app contributed $0 for as long as it ran.

Only the last one is new here; the tests cover all four, because a regression in any of them
looks exactly like a healthy run. Everything is mocked: no Modal account, no network, no cost.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "cloud" / "modal"))


@pytest.fixture
def guard(monkeypatch):
    monkeypatch.setenv("GOOGLE_CLOUD_PROJECT", "ci-no-such-project")
    return pytest.importorskip("guard")


def _cli(monkeypatch, guard, payload, returncode=0, stderr=""):
    """Replace the `modal app list` subprocess with a canned response."""
    def fake(cmd, *a, **k):
        out = payload if isinstance(payload, str) else json.dumps(payload)
        return subprocess.CompletedProcess(cmd, returncode, out, stderr)
    monkeypatch.setattr(guard.subprocess, "run", fake)


def _app(app_id, created, desc="rbp-gc-sweep", state="running"):
    return {"app_id": app_id, "created_at": created, "Description": desc, "State": state}


# --- what the CLI says, and what the guard is allowed to conclude from it ----------------

def test_a_failed_cli_call_is_not_an_empty_app_list(monkeypatch, guard):
    _cli(monkeypatch, guard, "", returncode=1, stderr="token expired")
    with pytest.raises(guard.Unobservable):
        guard.live_apps()


def test_unparseable_output_is_not_an_empty_app_list(monkeypatch, guard):
    _cli(monkeypatch, guard, "<html>login required</html>")
    with pytest.raises(guard.Unobservable):
        guard.live_apps()


def test_a_genuinely_empty_list_is_allowed_to_be_empty(monkeypatch, guard):
    _cli(monkeypatch, guard, [])
    assert guard.live_apps() == []


def test_every_arms_app_name_matches_not_just_one(monkeypatch, guard):
    """The names that actually ran. This pattern once matched only 'rbp-sweep'."""
    for desc in ("rbp-sweep", "rbp-gc-sweep", "rbp-neg2-sweep", "rbp-neg2-rm-sweep"):
        _cli(monkeypatch, guard, [_app("ap-1", "2026-09-06T10:00:00+00:00", desc=desc)])
        assert len(guard.live_apps()) == 1, f"{desc} is not being guarded"


def test_a_finished_app_is_not_charged_for(monkeypatch, guard):
    _cli(monkeypatch, guard, [_app("ap-1", "2026-09-06T10:00:00+00:00", state="stopped")])
    assert guard.live_apps() == []


# --- the bound itself --------------------------------------------------------------------

def test_a_live_app_with_an_unreadable_start_time_fails_closed(monkeypatch, guard):
    """THE FAIL-OPEN. It used to be kept in the list with a None and then filtered out."""
    _cli(monkeypatch, guard, [_app("ap-1", "not-a-timestamp")])
    with pytest.raises(guard.UnparseableStart) as e:
        guard.live_apps()
    assert "ap-1" in str(e.value)


def test_one_bad_timestamp_among_good_ones_still_fails_closed(monkeypatch, guard):
    """The dangerous case: the bound looks plausible because the other apps are fine."""
    _cli(monkeypatch, guard, [_app("ap-good", "2026-09-06T10:00:00+00:00"),
                              _app("ap-bad", "")])
    with pytest.raises(guard.UnparseableStart):
        guard.live_apps()


def test_report_refuses_to_bound_a_list_containing_a_startless_app(monkeypatch, guard):
    """Second line of defence, in case a caller assembles the list itself."""
    monkeypatch.setattr(guard, "work_done", lambda m: (0, 0.0, []))
    with pytest.raises(guard.UnparseableStart):
        guard.report("cnn", 40.0, None, [("ap-1", None)])


def test_the_upper_bound_sums_over_every_live_app(monkeypatch, guard):
    """It used to use the oldest app's elapsed time alone, which two apps can exceed."""
    monkeypatch.setattr(guard, "work_done", lambda m: (0, 0.0, []))
    now = 1_000_000.0
    monkeypatch.setattr(guard.time, "time", lambda: now)
    one, _ = guard.report("cnn", 40.0, None, [("a", now - 3600)])
    two, _ = guard.report("cnn", 40.0, None, [("a", now - 3600), ("b", now - 3600)])
    assert two == pytest.approx(2 * one), (
        "two concurrent apps each hold MAX_CONTAINERS, so the bound must double")
    assert one == pytest.approx(guard.MAX_CONTAINERS * guard.RATE)


def test_the_upper_bound_is_above_the_lower_bound(monkeypatch, guard):
    monkeypatch.setattr(guard, "work_done", lambda m: (10, 1800.0, []))
    now = 1_000_000.0
    monkeypatch.setattr(guard.time, "time", lambda: now)
    upper, lower = guard.report("cnn", 40.0, None, [("a", now - 3600)])
    assert upper > lower > 0


# --- stopping ------------------------------------------------------------------------------

def test_a_failed_stop_command_is_reported_as_a_failure(monkeypatch, guard):
    monkeypatch.setattr(guard.subprocess, "run",
                        lambda cmd, *a, **k: subprocess.CompletedProcess(cmd, 1, "", "denied"))
    assert guard.stop("ap-1") is False


def test_a_stop_that_leaves_the_app_live_is_not_confirmed(monkeypatch, guard):
    """`modal app stop` exiting 0 is a request accepted, not a container released."""
    monkeypatch.setattr(guard.subprocess, "run",
                        lambda cmd, *a, **k: subprocess.CompletedProcess(cmd, 0, "ok", ""))
    monkeypatch.setattr(guard, "live_apps", lambda: [("ap-1", 1.0)])
    monkeypatch.setattr(guard.time, "sleep", lambda s: None)
    assert guard.stop("ap-1") is False


def test_a_stop_is_confirmed_only_when_the_app_is_gone(monkeypatch, guard):
    monkeypatch.setattr(guard.subprocess, "run",
                        lambda cmd, *a, **k: subprocess.CompletedProcess(cmd, 0, "ok", ""))
    monkeypatch.setattr(guard, "live_apps", lambda: [])
    monkeypatch.setattr(guard.time, "sleep", lambda s: None)
    assert guard.stop("ap-1") is True


def test_a_stop_cannot_be_confirmed_while_modal_is_unreadable(monkeypatch, guard):
    monkeypatch.setattr(guard.subprocess, "run",
                        lambda cmd, *a, **k: subprocess.CompletedProcess(cmd, 0, "ok", ""))
    def blind():
        raise guard.Unobservable("cli down")
    monkeypatch.setattr(guard, "live_apps", blind)
    monkeypatch.setattr(guard.time, "sleep", lambda s: None)
    assert guard.stop("ap-1") is False


def test_app_state_returns_the_oldest_because_it_bounds_the_spend(monkeypatch, guard):
    _cli(monkeypatch, guard, [_app("ap-new", "2026-09-06T12:00:00+00:00"),
                              _app("ap-old", "2026-09-06T09:00:00+00:00")])
    app_id, started = guard.app_state()
    assert app_id == "ap-old", "the longest-running app is the one that bounds the spend"
    assert started is not None
