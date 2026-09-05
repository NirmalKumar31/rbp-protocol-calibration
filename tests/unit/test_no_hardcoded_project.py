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
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
# A PATTERN, NOT A NAME. This was the literal string `rbp-composition-2026` -- the project the
# study started in -- and it stayed that way after the move to `rbp-repro-2026`. So the test
# whose docstring says it "is the only thing that stops the habit coming back" spent the whole
# second half of the project guarding a name nothing used any more, while the habit came back
# under the new one: scripts/device_portability.py carried `rbp-repro-2026-derived` as an
# argparse default, and docs/REPRODUCE.md claimed in the same breath that no such literal
# existed anywhere in the source. An external review found it; this test could not.
#
# The pattern matches this project's id shape in any generation, plus the -derived and -raw
# bucket suffixes built from it. It deliberately does NOT try to match every possible GCP
# project id: a generic matcher would fire on `ci-no-such-project` in the CI workflow and on
# every hyphenated word in a docstring, and a test that cries wolf gets an exemption added
# rather than a bug fixed.
FORBIDDEN = re.compile(r"\brbp-[a-z0-9]+-20\d\d(-derived|-raw)?\b")
# A BILLING ACCOUNT ID IS THE THING THAT ACTUALLY MATTERS, and this test did not look for it.
# The real one sat as a shell default in cloud/cost.sh and in a gcloud command in
# cloud/killswitch/main.py, in a repo whose README claims "no hardcoded project id... a test
# fails the build if a literal reappears". It was true of the project id and false of the
# credential-adjacent identifier next to it.
BILLING = re.compile(r"\b\d{6}-[0-9A-F]{6}-[0-9A-F]{6}\b")

# docker/ was missing from this list, which is exactly how two hardcoded project ids
# survived in cloudbuild.cpu.yaml and cloudbuild.gpu.yaml and sent the first image build
# looking for a cache layer in the OLD project's registry. A test that does not look
# everywhere is a test that certifies the places it looked.
SEARCH = ["scripts", "src", "cloud", "docker"]
EXEMPT_SUFFIX = {".md", ".tf", ".tfvars", ".example"}
EXEMPT_NAME = {"cloud.py"}          # the resolver's own docstring explains the convention


def _tracked():
    """Files git would actually publish.

    Scanning the working directory is the wrong scope: cloud/jobs/rendered/ holds 30 untracked
    build artifacts carrying the old project id, and they are gitignored precisely so they
    never become public. What matters is what `git ls-files` would ship.

    THE FALLBACK EXISTS BECAUSE THE CONTAINER HAS NO GIT. This ran `git ls-files` at
    COLLECTION time and raised FileNotFoundError inside the image, which aborted collection and
    failed every Cloud Build of the GPU image -- so the published image had been stale since
    this test was added, and a Batch job died on an arm the image had never heard of. Where git
    is unavailable the fallback walks the search directories instead and takes its exclusions
    FROM .gitignore rather than from a hand-written list. Two ignored trees would otherwise
    trip it: `cloud/jobs/rendered/`, which holds 31 Batch specs rendered against whichever
    project submitted them, and `cloud/terraform/terraform.tfvars`, which holds the real
    billing account id and is ignored for exactly that reason. Reading the same file git reads
    keeps the two paths from drifting; a fallback that fails the moment it is exercised is a
    trap, not a safety net. The matching is substring-crude, so it can only over-exclude, and
    the git path already covers the superset.
    """
    try:
        out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True,
                             text=True)
    except FileNotFoundError:
        # .gitignore IS ITSELF ABSENT IN THE IMAGE, which is the second way this fallback
        # failed a build. Treat it as optional: where it exists there are ignored trees to
        # exclude, and where it does not there is nothing to exclude either, because the image
        # contains only what the Dockerfile copies.
        gi = ROOT / ".gitignore"
        pats = [ln.strip().strip("/") for ln in gi.read_text().splitlines()
                if ln.strip() and not ln.startswith("#") and "*" not in ln] \
            if gi.exists() else []
        for d in SEARCH:
            for p in (ROOT / d).rglob("*"):
                s = str(p.relative_to(ROOT))
                if p.is_file() and not any(q in s for q in pats):
                    yield p
        return
    for name in out.stdout.split("\0"):
        if not name:
            continue
        p = ROOT / name
        if p.is_file() and ".terraform/" not in name:
            yield p


_TRACKED = None


def _files():
    global _TRACKED
    if _TRACKED is None:
        _TRACKED = set(_tracked())
    for d in SEARCH:
        for p in (ROOT / d).rglob("*"):
            if not p.is_file() or p.suffix in EXEMPT_SUFFIX or p.name in EXEMPT_NAME:
                continue
            # .json was exempt, which is where cloud/jobs/ingest.json kept a service-account
            # address and a bucket name. Scan tracked .json; untracked build artifacts under
            # cloud/jobs/rendered/ are gitignored and never published.
            if "tfplan" in p.name or "__pycache__" in str(p):
                continue
            if p.suffix == ".json" and p not in _TRACKED:
                continue
            if p.suffix in {".py", ".sh", ".yaml", ".yml", ".json"}:
                yield p


def _executable_lines(path):
    """Line numbers whose content is code rather than narration.

    THE RULE IS ABOUT WHAT RUNS, NOT ABOUT WHAT IS WRITTEN. Widening the pattern from one
    historical project name to the id's shape turned up three more hits, and all three were
    docstrings explaining why the old bucket returns 403 -- `src/rbp/utils/localstore.py`
    exists BECAUSE that bucket died, and a docstring that cannot say which bucket is a
    docstring that has been lobotomised to please a regex. The test already grants docs this
    exemption for exactly this reason: "Docs describe a specific historical run and should
    name it."

    A comment cannot redirect a pipeline to the wrong account. An argparse default can, and
    that is the one this test failed to catch. So comments and docstrings are narration, and
    every other string literal is code.
    """
    text = path.read_text(errors="ignore")
    if path.suffix != ".py":
        # Shell, YAML and JSON: whole-line `#` comments only. An inline `#` inside a quoted
        # string would be mis-stripped, so it is not attempted; over-strictness here is safe.
        return {i for i, ln in enumerate(text.splitlines(), 1)
                if not ln.lstrip().startswith("#")}
    import ast
    import io
    import tokenize
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return set(range(1, text.count("\n") + 2))
    docstring_spans = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = getattr(node, "body", [])
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                docstring_spans.append((body[0].lineno, body[0].end_lineno))
    narration = {i for a, b in docstring_spans for i in range(a, b + 1)}
    for tok in tokenize.generate_tokens(io.StringIO(text).readline):
        if tok.type == tokenize.COMMENT:
            narration.add(tok.start[0])
    return set(range(1, text.count("\n") + 2)) - narration


@pytest.mark.parametrize("path", sorted(_files(), key=str), ids=lambda p: str(p.name))
def test_no_hardcoded_project_id(path):
    code = _executable_lines(path)
    hits = [f"{i}: {ln.strip()}"
            for i, ln in enumerate(path.read_text(errors="ignore").splitlines(), 1)
            if i in code and FORBIDDEN.search(ln)]
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


# --- the cloudbuild files must actually be valid YAML -----------------------------------
#
# Editing a comment into a substitutions block dropped the leading '#' from one line, which
# turned a comment into a bare YAML key and made `gcloud builds submit` fail before it even
# uploaded. A yaml file that looks fine in a diff is not a yaml file that parses.

@pytest.mark.parametrize("name", ["cloudbuild.cpu.yaml", "cloudbuild.gpu.yaml"])
def test_cloudbuild_yaml_parses(name):
    yaml = pytest.importorskip("yaml")
    p = ROOT / "docker" / name
    if not p.exists():
        pytest.skip(f"{name} not present")
    d = yaml.safe_load(p.read_text())
    assert "steps" in d, f"{name} has no steps"
    assert "_IMAGE" in d.get("substitutions", {}), f"{name} lost its _IMAGE substitution"



def test_no_billing_account_id_in_tracked_files():
    """A billing account ID must never be committed, in ANY tracked file, including docs.

    Not scoped to SEARCH, because docs/ carried it too and a public repo is public in all of
    its directories. This is the check that was missing while README claimed "a test fails the
    build if a literal reappears" -- true of the project id, false of the billing account.
    """
    hits = []
    for p in _tracked():
        try:
            text = p.read_text(errors="ignore")
        except (OSError, UnicodeDecodeError):
            continue
        for m in BILLING.finditer(text):
            hits.append(f"{p.relative_to(ROOT)}: {m.group(0)}")
    assert not hits, "billing account ID in tracked files:\n  " + "\n  ".join(hits)
