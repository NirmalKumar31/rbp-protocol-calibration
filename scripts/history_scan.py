"""Scan every commit on every ref for secret material, and write the result as an artefact.

    python scripts/history_scan.py            # write results/tables/history_scan.csv
    python scripts/history_scan.py --check    # rescan and fail if the finding changed

Why it is a script and not a paragraph. SECURITY.md described a scan run by hand on a date, over
"206 commits at the time of writing". By the time an external audit read it the repository was
at 210 refs-wide, the audit itself counted 208, and the file still said 206. Three numbers, three
sources, no way to tell which described the release. A count that has to be retyped is a count
that is wrong between the retyping and the next commit, and this one sits in the document that
tells a reader whether the history is safe to clone.

So the scan is code, the counts come out of git, and SECURITY.md points here instead of quoting.
`--check` rescans and fails if the set of findings has changed, which is the property that
matters: not that the number is stable, but that nothing new has appeared.

What it looks for. Private-key blocks, service-account key material, AWS access keys, GitHub and
Slack tokens, Google API keys, OAuth tokens, and GCP billing account IDs. It reports counts by
kind and the distinct commits each kind appears in. It does NOT print the matched strings: the
one live finding is documented in SECURITY.md and reprinting a billing account ID into a
committed table would be the same disclosure this exists to measure.

What it cannot do. Forks and existing clones keep whatever they were made from, so a clean scan
of this repository is a statement about this repository. SECURITY.md says so.
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

OUT = ROOT / "results" / "tables" / "history_scan.csv"

PATTERNS = {
    "private key block": re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
    "service-account key material": re.compile(r'"type"\s*:\s*"service_account"'),
    "AWS access key": re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
    "GitHub token": re.compile(r"\bgh[pousr]_[A-Za-z0-9]{36,}\b"),
    "Slack token": re.compile(r"\bxox[baprs]-[0-9A-Za-z-]{10,}\b"),
    "Google API key": re.compile(r"\bAIza[0-9A-Za-z_-]{35}\b"),
    "OAuth access token": re.compile(r"\bya29\.[0-9A-Za-z_-]{20,}\b"),
    "GCP billing account ID": re.compile(r"\b\d{6}-[0-9A-F]{6}-[0-9A-F]{6}\b"),
}


def _git(*args):
    """Run git, and treat a non-zero exit as a failure rather than as an empty result.

    The bug class this repository keeps hitting. cost.sh reported an auth failure as zero spend;
    the Modal guard read an expired token as no app running; test_no_hardcoded_project.py caught
    only FileNotFoundError from git and so stopped checking JSON inside an unpacked archive. An
    unchecked return code turns "could not look" into "found nothing", and in a secret scan that
    is the worst possible direction to be wrong in.
    """
    r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} exited {r.returncode}: "
                           f"{(r.stderr or r.stdout).strip()[:200]}")
    return r.stdout


def scan():
    n_head = int(_git("rev-list", "--count", "HEAD").strip())
    n_all = int(_git("rev-list", "--all", "--count").strip())

    # `git log --all -p` is the whole history as patches. Large, but this runs in seconds on a
    # repository of this size and the alternative is scanning blobs, which misses deletions.
    patch = _git("log", "--all", "-p", "--no-color", "--format=commit %H")
    commit, path = None, ""
    hits = {k: set() for k in PATTERNS}          # distinct commits
    places = {k: set() for k in PATTERNS}        # distinct (commit, file) pairs
    lines = {k: 0 for k in PATTERNS}             # diff lines, both sides
    for line in patch.splitlines():
        if line.startswith("commit ") and len(line) == 47:
            commit, path = line[7:], ""
            continue
        if line.startswith("+++ b/"):
            path = line[6:]
            continue
        for kind, pat in PATTERNS.items():
            if pat.search(line):
                hits[kind].add(commit)
                places[kind].add((commit, path))
                lines[kind] += 1

    rows = [{"check": "commits on HEAD", "value": n_head, "note": ""},
            {"check": "commits on all refs", "value": n_all,
             "note": "the scan covers this set"}]
    for kind in PATTERNS:
        rows.append({"check": f"commits containing a {kind}", "value": len(hits[kind]),
                     "note": "the matched strings are deliberately not printed here"})
        # Three different numbers, all true, AND SECURITY.md QUOTED ONE WITHOUT SAYING
        # WHICH. It said the billing ID "appears in 11 places in the history". That is the
        # DIFF LINE count, and it includes the removal side of the scrub commit; the ID was
        # introduced in 4 distinct commits and touches 9 commit-and-file pairs. All three are
        # recorded, each named for what it counts, so the document can be exact instead of
        # picking one and calling it "places".
        rows.append({"check": f"commit-and-file pairs containing a {kind}",
                     "value": len(places[kind]),
                     "note": "a file touched in two commits counts twice"})
        rows.append({"check": f"diff lines containing a {kind}", "value": lines[kind],
                     "note": "both sides of the diff, so a scrub commit contributes its "
                             "removal line as well as the original addition"})
    total = sum(len(v) for k, v in hits.items() if k != "GCP billing account ID")
    rows.append({"check": "commits containing credential material of any kind", "value": total,
                 "note": "everything except the billing account ID, which is an identifier "
                         "rather than a credential and is documented in SECURITY.md"})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    # An unpacked archive has no history to scan, and this crashed with a traceback there.
    # A git export, a Zenodo deposit, a tarball: the files are present and .git is not, so
    # `git rev-list` exits 128 and _git raises. run.sh gates on this script, so the documented
    # pipeline died in exactly the archival case the release is meant to be checked in --
    # which is the third time this repository has been bitten by the boundary between "cannot
    # look" and "found nothing", and the first time in code written to warn about it.
    #
    # The right answer is neither a crash nor a pass. There is no history here to be clean or
    # dirty, so say so, and fall back to reporting the committed finding, which IS in the
    # archive. If neither git nor the table is present, that is a real failure: nothing
    # establishes what the history contains.
    if not (ROOT / ".git").exists():
        if not OUT.exists():
            log("  no .git and no committed history_scan.csv: nothing establishes what this "
                "repository's history contains. Run this in a clone.")
            sys.exit(1)
        with OUT.open(newline="") as fh:
            have = {r["check"]: r["value"] for r in csv.DictReader(fh)}
        cred = have.get("commits containing credential material of any kind")
        log("  no .git in this tree, so the history cannot be rescanned here. Reporting the "
            f"committed scan: {have.get('commits on all refs', '?')} commits, "
            f"{cred} carrying credential material.")
        if cred not in ("0", 0):
            log("  THE COMMITTED SCAN RECORDS CREDENTIAL MATERIAL. Stop and rotate.")
            sys.exit(1)
        return

    rows = scan()
    for r in rows:
        log(f"  {r['check']:52s} {r['value']}")

    cred = next(r for r in rows
                if r["check"] == "commits containing credential material of any kind")
    if cred["value"]:
        log("  CREDENTIAL MATERIAL IS PRESENT IN THE HISTORY. This is not a documented finding; "
            "stop and rotate before doing anything else.")
        sys.exit(1)

    if a.check:
        if not OUT.exists():
            log("  history_scan.csv is missing")
            sys.exit(1)
        with OUT.open(newline="") as fh:
            have = {r["check"]: r["value"] for r in csv.DictReader(fh)}
        # The counts of commits are expected to move. What must not move is any finding.
        moved = [r["check"] for r in rows
                 if "commits containing" in r["check"]
                 and have.get(r["check"]) != str(r["value"])]
        if moved:
            log(f"  the scan's FINDINGS changed: {moved}. Rerun and read the diff before "
                "committing it.")
            sys.exit(1)
        log("  no finding has changed")
        return

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["check", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    log(f"  wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
