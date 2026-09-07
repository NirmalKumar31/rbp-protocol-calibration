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
import os
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
    """Every column, every field, and the row keys must be unique. Returns complaint strings.

    An earlier version compared the SET of `check` rows, the `note` column and `value`, and its
    docstring called that structural and numeric comparison. It was neither complete: `n`,
    `ci_low`, `ci_high` and every analysis-specific column went unchecked, so a run could move
    a confidence bound or a sample size and pass. An audit was right to call that a false-pass
    path in a gate whose whole purpose is to have none.
    """
    bad = []
    a = pd.read_csv(before)
    b = pd.read_csv(after)

    if list(a.columns) != list(b.columns):
        bad.append(f"COLUMNS changed: {list(a.columns)} -> {list(b.columns)}")
        return bad                      # nothing below is meaningful across different schemas

    key = "check" if "check" in a.columns else None
    if key is None:
        if len(a) != len(b):
            bad.append(f"ROW COUNT {len(a)} -> {len(b)} in a table with no `check` column")
            return bad
        ia = ib = range(len(a))
        ai, bi = a, b
    else:
        for lab, frame in (("committed", a), ("regenerated", b)):
            dup = frame[key][frame[key].duplicated()].tolist()
            if dup:
                bad.append(f"DUPLICATE keys in the {lab} table: {sorted(set(dup))[:3]}")
        if bad:
            return bad
        ia, ib = set(a[key]), set(b[key])
        for k in sorted(ia - ib):
            bad.append(f"DROPPED row: {k!r}")
        for k in sorted(ib - ia):
            bad.append(f"ADDED row: {k!r}")
        ai, bi = a.set_index(key), b.set_index(key)
        ia = ib = sorted(ia & ib)

    fields = [c for c in a.columns if c != key]
    for k in ia:
        for col in fields:
            va, vb = ai.loc[k, col], bi.loc[k, col]
            fa = pd.to_numeric(va, errors="coerce")
            fb = pd.to_numeric(vb, errors="coerce")
            if pd.notna(fa) and pd.notna(fb):
                if abs(fa - fb) > tol:                       # numeric, tolerance applies
                    bad.append(f"{col} moved on {k!r}: {fa!r} -> {fb!r}")
                continue
            # Text, or one side blank. Blank and NaN are the same absence in a CSV round trip.
            sa = "" if pd.isna(va) else str(va)
            sb = "" if pd.isna(vb) else str(vb)
            if sa != sb:
                what = "NOTE" if col == "note" else col.upper()
                extra = "; the table is stale against its script" if col == "note" else ""
                bad.append(f"{what} changed on {k!r}: {sa[:40]!r} -> {sb[:40]!r}{extra}")
    return bad


def snapshot_of(root):
    """Relative path -> bytes, over every table the release ships, not just top-level CSVs."""
    out = {}
    for f in sorted(root.rglob("*")):
        if f.is_file() and f.suffix in (".csv", ".tsv"):
            out[str(f.relative_to(root))] = f.read_bytes()
    return out


def restore(snap_dir, dest):
    """Put `dest` back exactly as `snap_dir`, including REMOVING files the run added.

    Copying the snapshot back leaves any new file in place, so it could contaminate the next
    entry point and survive the gate. The comment "restores what it touched" was therefore too
    strong, which an audit caught.
    """
    keep = set()
    for f in snap_dir.rglob("*"):
        if f.is_file():
            rel = f.relative_to(snap_dir)
            keep.add(rel)
            target = dest / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists() or target.read_bytes() != f.read_bytes():
                shutil.copy2(f, target)
    for f in sorted(dest.rglob("*"), reverse=True):
        if f.is_file() and f.relative_to(dest) not in keep:
            f.unlink()


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
    before = snapshot_of(snap / "tables")
    failures, ran = {}, 0
    try:
        for stem, flag in points:
            script = ROOT / "scripts" / f"{stem}.py"
            r = subprocess.run([sys.executable, str(script), flag], capture_output=True,
                               text=True, cwd=str(ROOT),
                               env={**os.environ, "PYTHONPATH": str(ROOT / "src")})
            # A NON-ZERO EXIT IS A FAILURE, NOT A SKIP. This reported them as skipped, so a
            # broken entry point dropped out of the evidence and the gate still passed: a
            # false-pass path in the one script written to have none.
            if r.returncode != 0:
                tail = (r.stderr or r.stdout).strip().splitlines()
                failures.setdefault(stem, []).append(
                    f"EXITED {r.returncode}: {tail[-1][:80] if tail else 'no output'}")
            else:
                ran += 1
            # Compared over the whole tree, so a NEW or DELETED table is caught too, and nested
            # paths and .tsv files are covered rather than top-level *.csv only.
            after = snapshot_of(TABLES)
            for rel in sorted(set(before) | set(after)):
                if rel not in after:
                    failures.setdefault(stem, []).append(f"DELETED table {rel}")
                elif rel not in before:
                    failures.setdefault(stem, []).append(f"NEW table {rel}")
                elif before[rel] != after[rel]:
                    for msg in compare(snap / "tables" / rel, TABLES / rel, a.tol):
                        failures.setdefault(stem, []).append(f"{rel}: {msg}")
            restore(snap / "tables", TABLES)
    finally:
        restore(snap / "tables", TABLES)
        shutil.rmtree(snap, ignore_errors=True)

    log(f"\n  {ran} of {len(points)} entry points ran clean, tolerance {a.tol:g}")
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
