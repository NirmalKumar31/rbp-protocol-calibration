"""The suite must not silently shrink.

golden.yaml carries `integrity.min_tests_passing`, which for a long time read 480 while the
suite had grown to 575, and nothing anywhere read the key. A floor that no one enforces and
no one updates is worse than no floor: it appears in the config as a guarantee.

This asserts the floor against the collected suite itself, so deleting tests to make a build
green fails the build.
"""

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def test_suite_is_at_least_the_documented_size(pytestconfig):
    floor = yaml.safe_load((ROOT / "config" / "golden.yaml").read_text())["integrity"][
        "min_tests_passing"]
    collected = pytestconfig.pluginmanager.get_plugin("terminalreporter")
    n = getattr(pytestconfig, "_collected_count", None)
    if n is None and collected is not None:
        n = getattr(collected, "_numcollected", None)
    assert n is not None, "cannot determine collected test count"
    assert n >= floor, (
        f"the suite collects {n} tests but config/golden.yaml requires at least {floor}. "
        f"Tests were removed, or the floor is stale. Both need saying out loud.")
