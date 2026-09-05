"""The provenance manifest must not overclaim reproducibility.

WHY A TEST. The first version of scripts/provenance.py classified 31 tables as more
reproducible than they are, and it did so while its docstring claimed to close the
table-to-script link for every committed table. Nothing caught that, because `--check` compared
hashes and names and never asked whether a status was true. These are the three defects an
external audit found, each pinned to the smallest input that reproduces it.
"""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
MANIFEST = ROOT / "results" / "tables" / "PROVENANCE.csv"


@pytest.fixture(scope="module")
def prov():
    return pytest.importorskip("provenance")


def test_a_commented_invocation_is_not_an_invocation(prov, tmp_path, monkeypatch):
    """run.sh says strand_audit.py is deliberately NOT run. It came out raw-reproducible."""
    fake = tmp_path / "run.sh"
    fake.write_text(
        's13_analysis() {\n'
        '  "$PY" scripts/real_one.py || die "real"\n'
        '  # scripts/ghost.py is deliberately NOT here. It needs the window store.\n'
        '}\n'
        'STAGES=(s13_analysis)\n')
    monkeypatch.setattr(prov, "RUN", fake)
    inv = prov.invocations()
    assert "real_one" in inv, "a real invocation must still be found"
    assert "ghost" not in inv, (
        "a script named only inside a comment was read as an invocation, which is how "
        "run.sh's own statement that a script cannot run became a reproducibility guarantee")


def test_a_cache_is_not_recomputable_from_itself(prov):
    """--from-cache rebuilds a summary FROM a per-dataset table. The table is an input."""
    assert prov.classify("x_per_dataset.csv", "x", prov.EVID) == prov.CACHE
    # And where no per-dataset sibling exists, the summary re-reads itself, so it is the cache.
    assert prov.classify("no_such_table_xyz.csv", "no_such_table_xyz",
                         prov.EVID) == prov.CACHE


def test_status_ranking_takes_the_hardest_requirement(prov):
    assert prov.RANK[prov.RAW] < prov.RANK[prov.EVID] < prov.RANK[prov.CACHE] \
        < prov.RANK[prov.FROZEN]


def test_the_manifest_covers_every_committed_table():
    """Top-level .csv only was the old scope; .tsv and nested files were silently omitted."""
    if not MANIFEST.exists():
        pytest.skip("PROVENANCE.csv not in this checkout")
    import csv
    listed = {r["table"] for r in csv.DictReader(MANIFEST.open())}
    tables = ROOT / "results" / "tables"
    meta = {"PROVENANCE.csv", "manuscript_orphans.csv", "release_facts.csv",
            "verify_summary.csv"}
    on_disk = {str(p.relative_to(tables)) for pat in ("*.csv", "*.tsv", "*/*.csv", "*/*.tsv")
               for p in tables.glob(pat) if p.name not in meta and p.name != "README.md"}
    assert on_disk <= listed, f"not in the manifest: {sorted(on_disk - listed)}"


def test_no_table_is_left_unattributed_by_accident():
    """UNKNOWN means nobody can regenerate it. Only unattributed/ may carry that."""
    if not MANIFEST.exists():
        pytest.skip("PROVENANCE.csv not in this checkout")
    import csv
    bad = [r["table"] for r in csv.DictReader(MANIFEST.open())
           if r["status"] == "UNKNOWN"]
    assert not bad, f"tables with no producing script outside unattributed/: {bad}"


def test_check_mode_agrees_with_the_committed_manifest():
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "provenance.py"), "--check"],
                       cwd=ROOT, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin", "PYTHONPATH": str(ROOT / "src")})
    assert r.returncode == 0, r.stdout + r.stderr
