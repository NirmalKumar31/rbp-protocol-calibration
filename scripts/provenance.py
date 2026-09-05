"""Every committed table, what produced it, and what a clean clone can actually rebuild.

    python scripts/provenance.py            # write results/tables/PROVENANCE.csv, check it
    python scripts/provenance.py --check    # check only, for CI

WHY. An external review made a distinction the release had blurred and it was right to. There
are two different claims:

  frozen-result verification   a clone can check that committed evidence is internally
                               consistent and matches frozen expectations. This repository does
                               that well: 959 assertions, offline, no credentials.

  computational reproduction   a clone can regenerate the evidence from documented inputs. This
                               repository does that in part. `run.sh all` rebuilds the
                               dinucleotide arm end to end; the GC and bias-aware sweeps ran on
                               Modal and their per-window scores are committed rather than
                               rebuilt; a few analyses need the 2.9 GB window store, which is
                               not redistributable because it contains genomic sequence.

README said "reproduce it in full", which claimed the second for work that had done the first.
The wording is fixed, but a sentence is not a manifest. This writes the manifest: one row per
committed table, naming its producing script, the stage that runs it, and which of three
statuses it has.

  raw-reproducible       rebuilt by `run.sh all` from downloaded inputs on any machine
  evidence-recomputable  rebuilt from committed per-window scores, no store and no cloud
  frozen-only            needs the window store, a GPU, or cloud credentials

THE STATUS IS DERIVED, NOT DECLARED. It is read from how run.sh invokes the script: a
`--from-cache` invocation recomputes a summary from a committed per-dataset table and is
therefore evidence-recomputable; a bare invocation inside the default stage list is
raw-reproducible; anything run.sh explicitly excludes from the default path, or that takes
`--store`, is frozen-only. Declaring the status by hand is how the README sentence went wrong.

WHAT THIS DOES NOT DO. It records the producing command and the reproduction class. It does not
hash upstream inputs -- the ENCODE and GENCODE accessions are in
`results/tables/supplementary_table_s1.csv` and the download code records URL, size and MD5 at
fetch time, but those records were never committed for the published run, so claiming input
integrity here would be claiming something not held.
"""

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
RUN = ROOT / "run.sh"
OUT = TABLES / "PROVENANCE.csv"

RAW = "raw-reproducible"
EVID = "evidence-recomputable"
FROZEN = "frozen-only"

# Tables written by this script or by the auditors, which have no upstream stage.
META = {"PROVENANCE.csv", "manuscript_orphans.csv", "release_facts.csv", "verify_summary.csv"}


def invocations():
    """script stem -> (status, the run.sh line that decides it).

    Parsed from run.sh rather than from a hand-kept list, because a hand-kept list is a second
    copy of the stage graph and would drift from it the first time a stage moved.
    """
    text = RUN.read_text()
    out = {}
    # The default stage list. A stage not named here does not run under `run.sh all`, which is
    # what "excluded from the default path" means concretely.
    m = re.search(r"STAGES=\(([^)]*)\)", text) or re.search(r"for s in (s\d[\w\s]*)\)", text)
    default = set(re.findall(r"s\d+\w*", m.group(1))) if m else set()

    stage = None
    for line in text.splitlines():
        fn = re.match(r"(s\d+\w*)\(\)\s*\{", line)
        if fn:
            stage = fn.group(1)
        for mm in re.finditer(r'scripts/([a-z0-9_]+)\.py([^|]*)', line):
            name, args = mm.group(1), mm.group(2)
            if "--from-cache" in args:
                st = EVID
            elif "--store" in args:
                st = FROZEN
            elif stage and default and stage not in default:
                st = FROZEN
            else:
                st = RAW
            # Worst status wins: a script invoked twice is only as reproducible as its
            # hardest requirement.
            rank = {RAW: 0, EVID: 1, FROZEN: 2}
            if name not in out or rank[st] > rank[out[name][0]]:
                out[name] = (st, stage or "?", line.strip()[:100])
    return out


def owner(table, inv):
    """Which script writes this table. Matched on the stem, then on the script's own text."""
    stem = table.replace("_per_dataset", "").replace("_per_fold", "").removesuffix(".csv")
    if stem in inv:
        return stem
    for name in inv:
        f = ROOT / "scripts" / f"{name}.py"
        if f.exists() and table in f.read_text():
            return name
    for f in sorted((ROOT / "scripts").glob("*.py")):
        if table in f.read_text():
            return f.stem
    return ""


def build():
    inv = invocations()
    rows = []
    for p in sorted(TABLES.glob("*.csv")):
        if p.name in META:
            continue
        who = owner(p.name, inv)
        status, stage, cmd = inv.get(who, (FROZEN, "", "not invoked by run.sh"))
        rows.append({
            "table": p.name,
            "producing_script": f"scripts/{who}.py" if who else "UNKNOWN",
            "run_sh_stage": stage,
            "status": status if who else "UNKNOWN",
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest()[:16],
            "bytes": p.stat().st_size,
            "invocation": cmd,
        })
    return rows


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--check", action="store_true",
                   help="fail if the committed manifest is stale or any table is unattributed")
    a = a.parse_args()

    rows = build()
    unknown = [r["table"] for r in rows if r["producing_script"] == "UNKNOWN"]
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    log(f"  {len(rows)} committed tables")
    for k in (RAW, EVID, FROZEN, "UNKNOWN"):
        if counts.get(k):
            log(f"    {k:22s} {counts[k]:3d}")

    if unknown:
        log("")
        log("  NO PRODUCING SCRIPT FOUND, which means a committed table nobody can regenerate:")
        for u in unknown:
            log(f"    {u}")

    if a.check:
        if not OUT.exists():
            log(f"\n  {OUT.name} is not committed; run scripts/provenance.py")
            return 1
        old = list(csv.DictReader(OUT.open()))
        cur = {r["table"]: r["sha256"] for r in rows}
        prev = {r["table"]: r["sha256"] for r in old}
        moved = sorted(t for t in cur if prev.get(t) not in (None, cur[t]))
        added = sorted(set(cur) - set(prev))
        gone = sorted(set(prev) - set(cur))
        if moved or added or gone:
            log("\n  PROVENANCE.csv is stale. Rerun scripts/provenance.py and commit it.")
            for t in moved:
                log(f"    changed: {t}")
            for t in added:
                log(f"    added:   {t}")
            for t in gone:
                log(f"    removed: {t}")
            return 1
        log("\n  manifest matches the committed tables")
        return 0

    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    log(f"\n  wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
