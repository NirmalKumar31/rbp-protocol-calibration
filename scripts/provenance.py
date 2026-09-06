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
CLOUD = "cloud-produced"        # written by a cloud stage, not by any local script
FROZEN = "frozen-only"          # needs the window store, a GPU, or cloud credentials
RANK = {RAW: 0, EVID: 1, CACHE: 2, FROZEN: 3, CLOUD: 4}

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

    A --from-cache invocation recomputes a SUMMARY from a committed table. The summary is
    recomputable from released evidence; the table it reads is not, it is a frozen input. Giving
    both the same status told a reader the input reconstructs itself.

    DECIDED BY WHAT THE SCRIPT ACTUALLY READS, not by the filename. The first version matched
    `X_per_dataset.csv` against the producer's stem, so an alias name defeated it:
    horlacher_arm.py writes horlacher_per_dataset.csv, the stems differ, and the frozen cache
    was labelled recomputable.
    """
    if status != EVID:
        return status
    src = ROOT / "scripts" / f"{who}.py"
    if src.exists() and table in _literals_read(src.read_text()):
        return CACHE
    return status


def _literals_written(src):
    """Filenames this source actually WRITES, from its `.to_csv(...)` calls.

    PARSED, NOT GREPPED, because the grep version was wrong. It attributed a table to the first
    alphabetically scanned script whose text contained the quoted filename and a `to_csv`
    anywhere -- so deep_contrast_per_dataset.csv, written at deep_model_contrast.py:376, came
    out owned by baseline_order_models.py, which READS it as an anchor. Ten scripts mention that
    filename. Only one writes it.

    Resolves three shapes, which is every shape this repository uses:
        df.to_csv(TABLES / "x.csv")            a literal in the call
        p = TABLES / "x.csv"; df.to_csv(p)     a name bound to a literal
        df.to_csv(out.with_suffix(".partial.csv"))   ignored; a partial is not a release table
    """
    import ast
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return set()

    def strings(node):
        """Every string constant reachable in an expression."""
        return {n.value for n in ast.walk(node)
                if isinstance(n, ast.Constant) and isinstance(n.value, str)}

    # Names bound to an expression containing a .csv/.tsv literal.
    bound = {}
    for n in ast.walk(tree):
        if isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0],
                                                                            ast.Name):
            s = {x for x in strings(n.value) if x.endswith((".csv", ".tsv"))}
            if s:
                bound.setdefault(n.targets[0].id, set()).update(s)

    out, patterns = set(), []
    for n in ast.walk(tree):
        if not (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)):
            continue
        # write_text is here because cloud_rehearsal.py writes rehearsal_binding_{arm}.csv with
        # it and an f-string; a to_csv-only parser reported that table as having no producer.
        if n.func.attr not in ("to_csv", "to_string", "write_text"):
            continue
        # For write_text the path is the RECEIVER, `(dir / name).write_text(body)`.
        targets = list(n.args[:1])
        if n.func.attr == "write_text":
            targets = [n.func.value]
        for arg in targets:
            out |= {x for x in strings(arg) if x.endswith((".csv", ".tsv"))}
            for nm in (x.id for x in ast.walk(arg) if isinstance(x, ast.Name)):
                out |= bound.get(nm, set())
            # An f-string names a FAMILY of files. Keep its literal fragments and match on
            # them, so `rehearsal_binding_{a.arm}.csv` claims rehearsal_binding_gc.csv.
            for js in (x for x in ast.walk(arg) if isinstance(x, ast.JoinedStr)):
                frag = [v.value for v in js.values
                        if isinstance(v, ast.Constant) and isinstance(v.value, str)]
                if frag and frag[-1].endswith((".csv", ".tsv")) and len(frag[0]) >= 6:
                    patterns.append((frag[0], frag[-1]))
    return {x for x in out if ".partial." not in x}, patterns


def _literals_read(src):
    """Filenames this source READS, from its `read_csv` calls. Same parse as the writer side."""
    import ast
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return set()

    def strings(node):
        return {n.value for n in ast.walk(node)
                if isinstance(n, ast.Constant) and isinstance(n.value, str)}

    bound = {}
    for n in ast.walk(tree):
        if isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0],
                                                                            ast.Name):
            s = {x for x in strings(n.value) if x.endswith((".csv", ".tsv"))}
            if s:
                bound.setdefault(n.targets[0].id, set()).update(s)
    out = set()
    for n in ast.walk(tree):
        if not (isinstance(n, ast.Call) and getattr(n.func, "attr", "") in
                ("read_csv", "read_table")):
            continue
        if not n.args:
            continue
        out |= {x for x in strings(n.args[0]) if x.endswith((".csv", ".tsv"))}
        for nm in (x.id for x in ast.walk(n.args[0]) if isinstance(x, ast.Name)):
            out |= bound.get(nm, set())
    return out


_WRITERS = None
_READERS = None


def writers():
    """filename -> [scripts that write it]. Built once."""
    global _WRITERS
    if _WRITERS is None:
        _WRITERS = {"_patterns": []}
        for f in sorted((ROOT / "scripts").glob("*.py")):
            names, pats = _literals_written(f.read_text())
            for name in names:
                _WRITERS.setdefault(name, []).append(f.stem)
            _WRITERS["_patterns"] += [(a_, b_, f.stem) for a_, b_ in pats]
    return _WRITERS


def readers():
    """filename -> [scripts that read it]. Built once."""
    global _READERS
    if _READERS is None:
        _READERS = {}
        for f in sorted((ROOT / "scripts").glob("*.py")):
            for name in _literals_read(f.read_text()):
                _READERS.setdefault(name, []).append(f.stem)
    return _READERS


def owner(table, inv):
    """Which script writes this table, established from its write calls.

    Falls back to the naming convention only when no script demonstrably writes the file,
    which covers tables produced by the cloud stages rather than by a local script. Returns ""
    rather than a guess when neither applies: UNKNOWN is a better answer than a wrong one,
    because a wrong one is reported as a reproducibility guarantee.
    """
    W = writers()
    w = W.get(table, [])
    if not w:
        w = [stem for pre, suf, stem in W["_patterns"]
             if table.startswith(pre) and table.endswith(suf)
             and len(table) >= len(pre) + len(suf)]
    if len(w) == 1:
        return w[0]
    if len(w) > 1:
        # Two scripts writing one name is a real problem; prefer the stem match, flag otherwise.
        stem = table.replace("_per_dataset", "").replace("_per_fold", "").removesuffix(".csv")
        return stem if stem in w else sorted(w)[0]
    stem = table.replace("_per_dataset", "").replace("_per_fold", "").removesuffix(".csv")
    if (ROOT / "scripts" / f"{stem}.py").exists() or stem in inv:
        return stem
    # A CLOUD STAGE IS A PRODUCER TOO. sweep_dinuc.csv and variant_tasks.tsv are written by the
    # Batch and sweep jobs and fetched from the bucket; reporting them as "nobody can regenerate
    # this" is false, and reporting them as locally reproducible would be worse.
    for f in sorted((ROOT / "scripts").glob("cloud_*.py")) + sorted((ROOT / "cloud").rglob("*.py")):
        if table in f.read_text():
            return f"CLOUD:{f.stem}"
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
        if who.startswith("CLOUD:"):
            status, stage, cmd = CLOUD, "", f"written by {who[6:]}, a cloud stage"
            who = who[6:]
        else:
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
    for k in (RAW, EVID, CACHE, FROZEN, CLOUD, "unattributed", "UNKNOWN"):
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
        # EVERY FIELD, not the name and the hash. Comparing those two let a row keep a wrong
        # producer and a wrong reproducibility class forever, because neither is a hash: an
        # audit found four rows naming a script that only READS the table, and --check passed.
        cur = {r["table"]: r for r in rows}
        prev = {r["table"]: r for r in old}
        added = sorted(set(cur) - set(prev))
        gone = sorted(set(prev) - set(cur))
        drift = []
        for name in sorted(set(cur) & set(prev)):
            for field in ("producing_script", "run_sh_stage", "status", "sha256", "bytes"):
                a_, b_ = str(prev[name].get(field, "")), str(cur[name][field])
                if a_ != b_:
                    drift.append(f"    {name}: {field} committed={a_!r} derived={b_!r}")
        if drift or added or gone:
            log("\n  PROVENANCE.csv is stale. Rerun scripts/provenance.py and commit it.")
            for x in drift:
                log(x)
            for x in added:
                log(f"    added:   {x}")
            for x in gone:
                log(f"    removed: {x}")
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
