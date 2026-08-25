"""No source file may hardcode the author's project or bucket names.

WHY. Eighteen files did. It is invisible while you only ever run in one project, and it is
a total failure the first time anyone reproduces the work elsewhere -- which is the entire
point of the pipeline being reproducible. This test is cheap and it is the only thing that
stops the habit coming back, because the failure mode is silent: the code runs, it just
runs against somebody else's account, or reads a stale bucket and does nothing.

Docs and Terraform are exempt. Docs describe a specific historical run and should name it;
Terraform declares the id as a variable with a default, which is the correct place for one.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FORBIDDEN = re.compile(r"rbp-composition-2026")

# docker/ was missing from this list, which is exactly how two hardcoded project ids
# survived in cloudbuild.cpu.yaml and cloudbuild.gpu.yaml and sent the first image build
# looking for a cache layer in the OLD project's registry. A test that does not look
# everywhere is a test that certifies the places it looked.
SEARCH = ["scripts", "src", "cloud", "docker"]
EXEMPT_SUFFIX = {".md", ".tf", ".tfvars", ".example"}
EXEMPT_NAME = {"cloud.py"}          # the resolver's own docstring explains the convention


def _files():
    for d in SEARCH:
        for p in (ROOT / d).rglob("*"):
            if not p.is_file() or p.suffix in EXEMPT_SUFFIX or p.name in EXEMPT_NAME:
                continue
            if "tfplan" in p.name or "__pycache__" in str(p) or p.suffix == ".json":
                continue
            if p.suffix in {".py", ".sh", ".yaml", ".yml"}:
                yield p


@pytest.mark.parametrize("path", sorted(_files(), key=str), ids=lambda p: str(p.name))
def test_no_hardcoded_project_id(path):
    hits = [f"{i}: {ln.strip()}"
            for i, ln in enumerate(path.read_text(errors="ignore").splitlines(), 1)
            if FORBIDDEN.search(ln)]
    assert not hits, (
        f"{path.relative_to(ROOT)} hardcodes the project id:\n  " + "\n  ".join(hits) +
        "\nResolve it through rbp.utils.cloud (python) or $PROJECT_ID (shell) instead.")


def test_the_resolver_refuses_to_guess(monkeypatch):
    """No silent default. An unconfigured environment must stop, not pick something."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from rbp.utils import cloud

    monkeypatch.delenv(cloud.ENV_PROJECT, raising=False)
    monkeypatch.setattr(cloud, "_from_config", lambda: {})
    with pytest.raises(RuntimeError, match="not configured"):
        cloud.project()


def test_env_overrides_config(monkeypatch):
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from rbp.utils import cloud

    monkeypatch.setenv(cloud.ENV_PROJECT, "someone-elses-project")
    assert cloud.project() == "someone-elses-project"
    assert cloud.derived_bucket() == "someone-elses-project-derived"
    assert cloud.raw_bucket() == "someone-elses-project-raw"
