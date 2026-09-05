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

So the floor is enforced where the whole suite is collectable and skipped where it provably is
not. "Collectable" has two ways of failing, and only one of them was handled. The absence of
torch is the first. The second is a caller who passes --ignore, which is exactly what
.github/workflows/ci.yml does for the two torch-importing modules: the collected count is then
a subset by construction and comparing it to a whole-suite floor is meaningless. On a machine
that has torch AND passes those ignores, which is what scripts/ci_local.sh does in order to
reproduce CI faithfully, the old test reported the suite as having shrunk by twelve tests when
nothing had been removed.
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
    ignored = pytestconfig.getoption("ignore") or []
    if ignored:
        pytest.skip(f"--ignore given ({len(ignored)} path(s)): the collected count is a "
                    f"subset by construction and cannot be compared to a whole-suite floor")
    # THE THIRD WAY COLLECTION SHRINKS WITHOUT ANY TEST BEING REMOVED: a partial file tree.
    # Several modules parametrise over files -- the figure checks over manuscript/figures/*.pdf,
    # the hardcoded-id checks over the tracked source -- so where those directories are absent
    # the same modules collect fewer cases. The container image copies config, src, scripts and
    # tests and nothing else, and collects 680 against 699 here. Raising the floor to match a
    # full checkout therefore failed the image build, which is the exact failure this file's
    # docstring warns about, and it took scripts/check_image_tree.sh to surface it.
    if not (ROOT / "manuscript").exists():
        pytest.skip("no manuscript/ in this tree: the file-parametrised tests collect fewer "
                    "cases by construction, so a whole-suite floor does not apply")
    floor = yaml.safe_load((ROOT / "config" / "golden.yaml").read_text())["integrity"][
        "min_tests_passing"]
    tr = pytestconfig.pluginmanager.get_plugin("terminalreporter")
    n = getattr(tr, "_numcollected", None) if tr is not None else None
    if n is None:
        pytest.skip("cannot determine collected test count in this runner")
    assert n >= floor, (
        f"the suite collects {n} tests but config/golden.yaml requires at least {floor}. "
        f"Tests were removed, or the floor is stale. Both need saying out loud.")
