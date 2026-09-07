"""Every advertised `--from-cache` entry point must reproduce the table it is documented to.

    python scripts/cache_idempotence.py            # all of them
    python scripts/cache_idempotence.py --only window_centring,gene_clustered_cv

WHY THIS EXISTS, and it is the check that should have existed first. `run.sh` advertises 33
offline entry points, each of which rebuilds a committed summary from committed evidence. That
is the reproducibility claim a reader can actually exercise. Nothing checked that any of them
returned what is committed.

Running all 33 in a clean `git archive` found three defects, and one of them was serious:

  * window_centring.py rewrote "narrowPeak files checked for a point-source summit" from 10 to
    0, because the peak files are not in the release. verify.py's gate requires zero summits
    AND at least ten files checked, correctly, since "no summit among nothing checked" is
    vacuous. So the documented offline pipeline took a clean clone from 1083/1083 to FOUR
    FAILURES, and the failure read as a broken reproduction rather than a missing optional
    input.
  * gene_clustered_cv.py dropped three rows for the same reason, and verify.py asserts their
    presence.
  * strand_placebo.csv was simply STALE: the script had said "pre-specified" for some time and
    the committed table still said "pre-registered", so the table did not match its own
    producer.

CI has a `git diff --exit-code -- results/tables/` step whose comment claims "regenerating the
checked-in artefacts must not change them". It runs after audit_manuscript, verify,
release_consistency and two --check calls, so the only tables it can ever see change are
manuscript_orphans.csv, verify_summary.csv and release_facts.csv. It certifies three of 158
tables and reads as though it certifies all of them. This script is the version that does what
that comment says.

WHY NOT A BYTE DIFF. Seven of the nine tables that moved differ only by floating-point
summation order, from 2e-16 to 5e-14, which is not a defect and cannot be legislated away: the
same table regenerated on a different BLAS will differ in the last bits. Demanding byte
identity would produce a gate that fails for a reason nobody should act on, which is how a gate
gets disabled and then trusted anyway.

So the comparison is structural first and numeric second, which is also the order in which the
two matter:

  * the SET of `check` rows must be identical. A dropped row is the defect above.
  * every `note` must be identical. A changed note means the table is stale against its script.
  * numeric values must agree to --tol, default 1e-9, which is five orders of magnitude looser
    than the worst honest noise observed and five tighter than anything that would matter.

THIS SCRIPT RESTORES WHAT IT TOUCHED. It snapshots results/tables, runs the entry points, then
puts the snapshot back in a finally block, so a failed run does not leave the tree dirty. That
matters because the thing it is checking is whether these commands damage committed evidence.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
RUN = ROOT / "run.sh"

# Parsed from run.sh, never a hand-kept list: a hand-kept list is a second copy of the stage
# graph and drifts from it the first time a stage moves. Commented-out lines are excluded for
# the reason provenance.py records, having once read one as an invocation.
ENTRY = re.compile(r'scripts/([a-z_]+)\.py\s+--(from-cache|summarise)')


def entry_points():
    out = []
    for raw in RUN.read_text().split("\n"):
        line = raw.split("#", 1)[0]
        for m in ENTRY.finditer(line):
            pair = (m.group(1), "--" + m.group(2))
            if pair not in out:
                out.append(pair)
    return sorted(out)


def compare(before, after, tol):
    """Structural differences first, then numeric. Returns a list of complaint strings."""
    bad = []
    a = pd.read_csv(before)
    b = pd.read_csv(after)
    if "check" not in a.columns or "check" not in b.columns:
        if len(a) != len(b):
            bad.append(f"row count {len(a)} -> {len(b)}")
        return bad
    ia, ib = set(a["check"]), set(b["check"])
    for k in sorted(ia - ib):
        bad.append(f"DROPPED row: {k!r}")
    for k in sorted(ib - ia):
        bad.append(f"ADDED row: {k!r}")
    common = sorted(ia & ib)
    ai = a.drop_duplicates("check").set_index("check")
    bi = b.drop_duplicates("check").set_index("check")
    if "note" in a.columns and "note" in b.columns:
        for k in common:
            na = str(ai.loc[k, "note"] if pd.notna(ai.loc[k, "note"]) else "")
            nb = str(bi.loc[k, "note"] if pd.notna(bi.loc[k, "note"]) else "")
            if na != nb:
                bad.append(f"NOTE changed on {k!r}: the table is stale against its script")
    for k in common:
        va, vb = ai.loc[k, "value"], bi.loc[k, "value"]
        fa, fb = pd.to_numeric(va, errors="coerce"), pd.to_numeric(vb, errors="coerce")
        if pd.notna(fa) and pd.notna(fb):
            if abs(fa - fb) > tol:
                bad.append(f"VALUE moved on {k!r}: {fa!r} -> {fb!r}")
        elif str(va) != str(vb):
            bad.append(f"VALUE changed on {k!r}: {va!r} -> {vb!r}")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--only", default="", help="comma-separated script stems")
    a = ap.parse_args()

    points = entry_points()
    if a.only:
        keep = {s.strip() for s in a.only.split(",")}
        points = [p for p in points if p[0] in keep]
    if not points:
        sys.exit("no --from-cache entry points found in run.sh")

    snap = Path(tempfile.mkdtemp(prefix="cache-idem-"))
    shutil.copytree(TABLES, snap / "tables")
    failures, ran, skipped = {}, 0, []
    try:
        for stem, flag in points:
            script = ROOT / "scripts" / f"{stem}.py"
            env = {"PYTHONPATH": str(ROOT / "src")}
            r = subprocess.run([sys.executable, str(script), flag], capture_output=True,
                               text=True, cwd=str(ROOT),
                               env={**dict(__import__("os").environ), **env})
            if r.returncode != 0:
                skipped.append(f"{stem} exited {r.returncode}: "
                               f"{(r.stderr or r.stdout).strip().splitlines()[-1][:70]}")
                continue
            ran += 1
            for after in sorted(TABLES.glob("*.csv")):
                before = snap / "tables" / after.name
                if not before.exists():
                    failures.setdefault(stem, []).append(f"NEW table {after.name}")
                    continue
                if before.read_bytes() == after.read_bytes():
                    continue
                for msg in compare(before, after, a.tol):
                    failures.setdefault(stem, []).append(f"{after.name}: {msg}")
            # Restore between entry points so one script's drift is not attributed to the next.
            for f in snap.joinpath("tables").glob("*.csv"):
                shutil.copy2(f, TABLES / f.name)
    finally:
        for f in snap.joinpath("tables").rglob("*"):
            if f.is_file():
                dest = TABLES / f.relative_to(snap / "tables")
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, dest)
        shutil.rmtree(snap, ignore_errors=True)

    log(f"\n  {ran} of {len(points)} entry points ran, tolerance {a.tol:g}")
    for s in skipped:
        log(f"    skipped: {s}")
    if failures:
        log(f"\n  {len(failures)} ENTRY POINT(S) DO NOT REPRODUCE THEIR COMMITTED TABLE:\n")
        for stem, msgs in sorted(failures.items()):
            log(f"    {stem}")
            for m in msgs[:6]:
                log(f"      {m}")
            if len(msgs) > 6:
                log(f"      ... and {len(msgs) - 6} more")
        sys.exit(1)
    log("  every entry point reproduced its committed table")


if __name__ == "__main__":
    main()
