"""The suite must not silently shrink.

golden.yaml carries `integrity.min_tests_passing`, which read 480 while the suite had grown
to 576, and nothing anywhere read the key. A floor that no one enforces and no one updates is
worse than no floor: it sits in the config looking like a guarantee.

WHY THIS SKIPS WITHOUT TORCH, which is not a way of dodging the check. The CPU image
deliberately ships no torch -- that is the whole reason there are two images -- so the
torch-dependent tests cannot be COLLECTED there and the container legitimately sees 514 tests
against 576 locally. Asserting a single floor across both environments fails the image build
for a reason that has nothing to do with the suite shrinking, which is exactly the kind of
gate that gets disabled and then trusted anyway.

So the floor is enforced where the whole suite is collectable -- a developer machine and CI --
and skipped where it provably is not.
"""

from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]


def _has_torch():
    try:
        import torch  # noqa: F401
        return True
    except Exception:
        return False


def test_suite_is_at_least_the_documented_size(pytestconfig):
    if not _has_torch():
        pytest.skip("no torch: this is the CPU image, which cannot collect the full suite")
    floor = yaml.safe_load((ROOT / "config" / "golden.yaml").read_text())["integrity"][
        "min_tests_passing"]
    tr = pytestconfig.pluginmanager.get_plugin("terminalreporter")
    n = getattr(tr, "_numcollected", None) if tr is not None else None
    if n is None:
        pytest.skip("cannot determine collected test count in this runner")
    assert n >= floor, (
        f"the suite collects {n} tests but config/golden.yaml requires at least {floor}. "
        f"Tests were removed, or the floor is stale. Both need saying out loud.")
