"""Every committed table, what produced it, and what a clean clone can actually rebuild.

    python scripts/provenance.py            # write results/tables/PROVENANCE.csv, check it
    python scripts/provenance.py --check    # check only, for CI

WHY. An external review made a distinction the release had blurred and it was right to. There
are two different claims:

  frozen-result verification   a clone can check that committed evidence is internally
                               consistent and matches frozen expectations. This repository does
                               that well, offline and with no credentials.

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
  evidence-recomputable  a summary rebuilt from a committed per-dataset table
  frozen-cache           the committed table such a run READS. An input, not an output
  frozen-only            needs the window store, a GPU, or cloud credentials
  unattributed           no producing script; see results/tables/unattributed/README.md

THE FIRST VERSION OF THIS FILE OVERCLAIMED FOR 31 TABLES, in three ways an audit found and this
docstring should keep naming. It read commented-out lines of run.sh as invocations, so the one
script run.sh explicitly says is not runnable came out raw-reproducible. It gave a --from-cache
summary and the cache it reads the same status, which told a reader the input reconstructs
itself. And it attributed a table to any script whose text mentioned the filename, which is
whichever script READS it. Fixing all three moved 52 raw-reproducible down to 21.

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

RAW = "raw-reproducible"       # rebuilt by `run.sh all` from downloaded inputs
EVID = "evidence-recomputable"  # a summary rebuilt from a committed per-dataset table
CACHE = "frozen-cache"          # the committed table a --from-cache run READS. An input.
FROZEN = "frozen-only"          # needs the window store, a GPU, or cloud credentials
RANK = {RAW: 0, EVID: 1, CACHE: 2, FROZEN: 3}

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
    for raw_line in text.splitlines():
        fn = re.match(r"(s\d+\w*)\(\)\s*\{", raw_line)
        if fn:
            stage = fn.group(1)
        # COMMENTS ARE NOT INVOCATIONS. This scanned every line, so run.sh's comment saying
        # "scripts/strand_audit.py is deliberately NOT here, it needs --gtf and --datasets"
        # was read as a bare invocation inside a default stage and classified
        # raw-reproducible -- a false reproducibility claim, generated by the tool whose job is
        # to stop false reproducibility claims, about the one file run.sh says is not runnable.
        line = raw_line.split("#", 1)[0]
        if not line.strip():
            continue
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
            if name not in out or RANK[st] > RANK[out[name][0]]:
                out[name] = (st, stage or "?", line.strip()[:100])
    return out


def classify(table, who, status):
    """Refine a script's status down to the individual table.

    A --from-cache invocation recomputes a SUMMARY from a committed per-dataset table. The
    summary is genuinely recomputable from released evidence; the table it reads is not, it is
    a frozen input that came from the window store or a GPU sweep. Giving both the same status
    told a reader that the input reconstructs itself.

    The rule is the naming convention the repository already follows: `X_per_dataset.csv` and
    `X_per_fold.csv` are the cache for `X.csv`. Where a --from-cache script has no such file,
    it re-reads its own summary, and that summary is the cache.
    """
    if status != EVID:
        return status
    stem = f"{who}"
    if table.startswith(stem) and ("_per_dataset" in table or "_per_fold" in table):
        return CACHE
    if table == f"{stem}.csv":
        sibling = TABLES / f"{stem}_per_dataset.csv"
        sibling_f = TABLES / f"{stem}_per_fold.csv"
        return EVID if (sibling.exists() or sibling_f.exists()) else CACHE
    return status


def owner(table, inv):
    """Which script writes this table.

    ORDER MATTERS AND THE OLD ORDER WAS WRONG. This used to fall back to "any script whose text
    mentions this filename", which attributes a table to whichever file happens to READ it. That
    is how strand_audit.csv came out attributed to an invoked script and therefore
    raw-reproducible, when scripts/strand_audit.py is not invoked at all and run.sh says in as
    many words that it cannot be.

    Now: the naming convention first, since this repository follows it without exception; then a
    quoted-literal write of the exact filename; then nothing. Returning "" is a better answer
    than a plausible wrong one, because "" is reported as UNKNOWN and a wrong attribution is
    reported as a reproducibility guarantee.
    """
    stem = table.replace("_per_dataset", "").replace("_per_fold", "").removesuffix(".csv")
    if (ROOT / "scripts" / f"{stem}.py").exists():
        return stem
    if stem in inv:
        return stem
    for f in sorted((ROOT / "scripts").glob("*.py")):
        src = f.read_text()
        if f'"{table}"' in src and "to_csv" in src:
            return f.stem
    return ""


def build():
    inv = invocations()
    rows = []
    # EVERY COMMITTED TABLE, which previously meant top-level *.csv only. That omitted
    # variant_tasks.tsv and the two files under unattributed/, while the docstring claimed
    # coverage of all of them.
    files = sorted(list(TABLES.glob("*.csv")) + list(TABLES.glob("*.tsv"))
                   + list(TABLES.glob("*/*.csv")) + list(TABLES.glob("*/*.tsv")),
                   key=lambda q: str(q.relative_to(TABLES)))
    for p in files:
        if p.name in META:
            continue
        rel = str(p.relative_to(TABLES))
        if rel.startswith("unattributed/"):
            rows.append({"table": rel, "producing_script": "", "run_sh_stage": "",
                         "status": "unattributed",
                         "sha256": hashlib.sha256(p.read_bytes()).hexdigest()[:16],
                         "bytes": p.stat().st_size,
                         "invocation": "no producing script; see unattributed/README.md"})
            continue
        who = owner(p.name, inv)
        status, stage, cmd = inv.get(who, (FROZEN, "", "not invoked by run.sh"))
        status = classify(p.name, who, status)
        rows.append({
            "table": rel,
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
    for k in (RAW, EVID, CACHE, FROZEN, "unattributed", "UNKNOWN"):
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
