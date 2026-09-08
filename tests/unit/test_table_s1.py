"""The typeset Table S1 must be what its CSV regenerates.

Why this exists. The supplement used to point at supplementary_table_s1.csv and print nothing,
so Table S1 was a reference rather than a table. It is now typeset, which trades one risk for
another: 95 rows in a .tex file go stale the moment the panel changes, and nothing would say so.

build.sh regenerates it, but build.sh cannot be the gate. scripts/pdf_freshness.py copies
manuscript/ alone into a temp tree to rebuild the PDFs, so the generator under scripts/ is not
reachable there and the call is conditional. A step that is skipped in the environment that
matters is not a check. This test is the check, and it runs everywhere.
"""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
TEX = ROOT / "manuscript" / "sections" / "table_s1.tex"
CSV = ROOT / "results" / "tables" / "supplementary_table_s1.csv"
GEN = ROOT / "scripts" / "table_s1_tex.py"


def test_the_committed_table_is_what_the_csv_regenerates(tmp_path):
    """Byte-for-byte. The generator is deterministic, so anything else means one of them moved."""
    for f in (TEX, CSV, GEN):
        if not f.exists():
            pytest.skip(f"{f.name} not in this checkout")
    before = TEX.read_bytes()
    backup = tmp_path / "table_s1.tex"
    backup.write_bytes(before)
    try:
        r = subprocess.run([sys.executable, str(GEN)], capture_output=True, text=True)
        assert r.returncode == 0, f"the generator failed: {r.stderr[-400:]}"
        assert TEX.read_bytes() == before, (
            "manuscript/sections/table_s1.tex is not what scripts/table_s1_tex.py produces from "
            "results/tables/supplementary_table_s1.csv. Either the CSV changed and the table was "
            "not regenerated, or the table was edited by hand. Run the generator and commit.")
    finally:
        TEX.write_bytes(backup.read_bytes())


def test_every_csv_row_reaches_the_typeset_table():
    """A generator that silently dropped rows would still be byte-identical to itself."""
    if not (TEX.exists() and CSV.exists()):
        pytest.skip("not in this checkout")
    import csv as _csv
    rows = list(_csv.DictReader(CSV.open(newline="")))
    body = [ln for ln in TEX.read_text().splitlines() if ln.rstrip().endswith(r"\\")
            and not ln.startswith((r"\multicolumn", r"\textbf"))]
    data = [ln for ln in body if "&" in ln and r"\textbf" not in ln]
    assert len(data) == len(rows), (
        f"{len(rows)} rows in the CSV, {len(data)} typeset. The table must be complete: a reader "
        "who cannot find an accession in it has no way to know it was dropped rather than absent")


def test_the_three_arm_flag_is_rendered_not_leaked_raw():
    """`True`/`False` in a printed table is a leaked serialisation, not a value a reader wants."""
    if not TEX.exists():
        pytest.skip("not in this checkout")
    t = TEX.read_text()
    assert "& True" not in t and "& False" not in t, "raw booleans reached the typeset table"
    assert "& yes" in t and "& no" in t, "the flag no longer renders both of its states"
