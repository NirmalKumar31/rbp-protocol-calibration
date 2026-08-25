"""Stage 3 in the cloud: one (cell line, protein, negative arm) per Batch task.

WHY THE WORK IS SPLIT THIS WAY. Preprocessing is embarrassingly parallel across proteins
and completely serial within one, so the natural unit is a single dataset. 244 datasets x 2
negative arms = 488 tasks. The two arms are separate tasks rather than one task doing both
because their costs differ by 27x -- pairing them would hide the fast one behind the slow
one on every node.

WHY THE OUTPUT IS COMPARED BYTE FOR BYTE. Nothing downstream is trustworthy if the
container preprocesses differently from the laptop, and "differently" here would be a
handful of windows out of millions, invisible in any summary statistic. md5 of dataset.tsv
is the only check that actually catches it.

    python scripts/cloud_prep.py index       # build regions.pkl from the GTF, publish
    python scripts/cloud_prep.py manifest    # write the task list (run locally)
    python scripts/cloud_prep.py prep        # one task, indexed by BATCH_TASK_INDEX
    python scripts/cloud_prep.py finalize    # collect reports -> panel_final / excluded
"""

import argparse
import fcntl
import os
import pickle
import resource
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import prepare  # noqa: E402
from rbp.data import annotation as ann  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
WORK = Path(os.environ.get("WORK_DIR", "/tmp/rbp"))

GENOME = "GRCh38.primary_assembly.genome.fa"
GTF = "gencode.v45.annotation.gtf.gz"
INDEX = "interim/regions.pkl"
MANIFEST = "manifest/prep_tasks.tsv"
ARMS = ("gc", "dinuc")
CELLS = {"K562": "config/panel_full.tsv", "HepG2": "config/panel_full_HepG2.tsv"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def peak_rss_gb():
    """ru_maxrss is KB on Linux. Reported so the next run can be sized on measurement
    rather than on a guess about how much memory a cKDTree over 16 dimensions wants."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def stage(bucket, name, dest):
    """Download an object to the node's disk, once per node rather than once per task.

    Four tasks share a VM and all four need the same 3.1 GB genome. Without the lock they
    race and pull it four times; with it, one pulls and three wait about 25 seconds. The
    size check makes a rerun after preemption free.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    blob = bucket.get_blob(name)
    if blob is None:
        raise FileNotFoundError(f"gs://{bucket.name}/{name}")
    with open(f"{dest}.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if dest.exists() and dest.stat().st_size == blob.size:
            return dest
        tmp = Path(f"{dest}.part")
        blob.download_to_filename(str(tmp))
        os.replace(tmp, dest)          # atomic, so a killed task cannot leave a half file
    return dest


def read_tsv(text):
    lines = text.strip().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def buckets(a):
    """On a Batch VM the metadata server supplies the project, so `project` is None and
    that is correct. Off it, `gcloud config` is invisible to the client library -- the CLI
    being authenticated says nothing about the library -- so it has to be passed in."""
    from google.cloud import storage
    c = storage.Client(project=a.project)
    return c.bucket(a.raw), c.bucket(a.derived)


# --- mode: index ------------------------------------------------------------------------

def do_index(a):
    """Parse the GTF into region intervals in the cloud instead of uploading the local
    pickle. The parse is deterministic, so the resulting md5 is itself a reproducibility
    gate -- and it is the same reasoning that had the cloud build the genome .fai."""
    import hashlib

    raw, derived = buckets(a)
    out = derived.blob(INDEX)
    if out.exists() and not a.force:
        log(f"present: gs://{a.derived}/{INDEX}")
        return
    gtf = stage(raw, GTF, WORK / "data/raw" / GTF)
    log(f"parsing {GTF}")
    t0 = time.time()
    index = ann.build_index(str(gtf))
    blob = pickle.dumps(index, protocol=5)
    log(f"  built in {time.time() - t0:.0f}s, {len(blob) / 1e6:.1f} MB, "
        f"md5 {hashlib.md5(blob).hexdigest()}")
    for r, s in ann.stats(index).items():
        log(f"  {r:9} {s['intervals']:>9,} intervals  {s['bases']:>13,} bases")
    out.upload_from_string(blob, content_type="application/octet-stream")
    log(f"published gs://{a.derived}/{INDEX}  peak rss {peak_rss_gb():.2f} GB")


# --- mode: manifest ---------------------------------------------------------------------

def do_manifest(a):
    """Freeze the task list. BATCH_TASK_INDEX is only meaningful against a fixed ordering,
    so the list is written once and read by every task rather than recomputed per task
    from an API that could return a different order.

    Ordered biggest-peak-file first. Batch hands contiguous blocks of indices to nodes, so
    with alphabetical order the heaviest datasets land wherever they happen to fall and the
    job finishes when the unluckiest node does; longest-first is the standard scheduling fix.

    Peak file size is the proxy rather than the local `seconds`, because those timings turn
    out to measure page-cache state as much as work -- QKI's GC arm is recorded at 971s and
    reruns in 7s from identical inputs to a byte-identical file. File size is a property of
    the data, available for all 244, and needs no local state.
    """
    raw, derived = buckets(a)
    rows = []
    for cell, panel in CELLS.items():
        sizes = {b.name.split("/")[-1].split(".")[0]: b.size
                 for b in raw.client.list_blobs(a.raw, prefix=f"peaks/{cell}/")}
        for r in read_tsv(raw.blob(panel).download_as_text()):
            name = r["protein"]
            if name not in sizes:
                sys.exit(f"no peak file in GCS for {name} {cell}")
            rows.extend((cell, name, arm, sizes[name]) for arm in ARMS)
    rows.sort(key=lambda r: (-r[3], r[0], r[1], r[2]))
    body = "idx\tcell_line\tprotein\tarm\tpeak_bytes\n" + "".join(
        f"{i}\t{c}\t{p}\t{m}\t{n}\n" for i, (c, p, m, n) in enumerate(rows))
    derived.blob(MANIFEST).upload_from_string(body)
    log(f"{len(rows)} tasks -> gs://{a.derived}/{MANIFEST}")
    log(f"  peak files {min(r[3] for r in rows) / 1e3:.0f} KB to "
        f"{max(r[3] for r in rows) / 1e6:.1f} MB")
    for arm in ARMS:
        log(f"  {arm:6} {sum(1 for r in rows if r[2] == arm)} tasks")
    for cell in CELLS:
        log(f"  {cell:6} {sum(1 for r in rows if r[0] == cell)} tasks")


# --- mode: prep -------------------------------------------------------------------------

def do_prep(a):
    raw, derived = buckets(a)
    tasks = read_tsv(derived.blob(MANIFEST).download_as_text())
    idx = a.index if a.index is not None else int(os.environ.get("BATCH_TASK_INDEX", 0))
    # TASK_LIST maps the array's 0..n-1 onto arbitrary manifest rows. Needed for the smoke
    # test, which wants two specific datasets, and again afterwards: spot preemption leaves
    # a scattered handful of failures and rerunning all 488 to catch nine of them is silly.
    if os.environ.get("TASK_LIST"):
        picks = [int(x) for x in os.environ["TASK_LIST"].split(",")]
        idx = picks[idx]
    if idx >= len(tasks):
        sys.exit(f"index {idx} beyond manifest of {len(tasks)}")
    t = tasks[idx]
    cell, name, arm = t["cell_line"], t["protein"], t["arm"]
    prefix = f"processed/{arm}/{cell}/{name}"
    log(f"task {idx}: {name} {cell} {arm}-matched  ({int(t['peak_bytes']) / 1e3:.0f} KB of peaks)")

    # The REPORT is the completion marker, because it is uploaded last. Keying resume off
    # dataset.tsv instead would let a task preempted between the two uploads come back, see
    # the dataset, skip the work -- and leave finalize with no report, silently dropping
    # that dataset from the panel with nothing anywhere saying so.
    if derived.blob(f"{prefix}/prep_report.txt").exists() and not a.force:
        log("already present, nothing to do")
        return

    # Everything above is manifest bookkeeping and costs nothing; everything below moves
    # gigabytes. --dry-run stops here so the 488-task array can be rehearsed for free. It
    # exists because the first smoke run died on a renamed manifest column two seconds in:
    # a compile check cannot see a dict key, and nothing else exercised this path locally.
    if a.dry_run:
        log("dry run: manifest resolved, stopping before staging")
        return

    (WORK / "data/interim").mkdir(parents=True, exist_ok=True)
    cfg = cfgmod.load(a.config)
    cfg["encode"]["cell_line"] = cell

    # The pipeline addresses everything relative to the repo root, so prepare's notion of
    # the root is pointed at the node's scratch disk. The frozen fold map is read FIRST,
    # while the root still points into the image: config/ is code and already ships there.
    # Copying it onto the shared disk instead would put four concurrent tasks in a race on
    # the same copytree, where the second one sees the directory exist and reads a
    # half-written folds.tsv.
    fold_map = prepare.load_fold_map(cfg)
    prepare.ROOT = WORK

    panel = read_tsv(raw.blob(CELLS[cell]).download_as_text())
    acc = next((r["accession"] for r in panel if r["protein"] == name), None)
    if acc is None:
        sys.exit(f"{name} not in {CELLS[cell]}")

    t0 = time.time()
    stage(raw, GENOME, WORK / "data/raw" / GENOME)
    stage(raw, f"{GENOME}.fai", WORK / "data/raw" / f"{GENOME}.fai")
    stage(raw, f"peaks/{cell}/{name}.{acc}.bed.gz",
          WORK / f"data/raw/peaks/{cell}/{name}.{acc}.bed.gz")
    index_path = WORK / "data/interim/regions.pkl"
    stage(derived, INDEX, index_path)
    log(f"  staged inputs in {time.time() - t0:.0f}s")

    from pyfaidx import Fasta
    index = pickle.loads(index_path.read_bytes())
    fasta = Fasta(str(WORK / "data/raw" / GENOME))
    outdir = WORK / "out" / arm / cell

    rep = prepare.prepare_one(name, cfg, fasta, index, outdir, fold_map, arm)
    log(f"  {rep['pairs']} pairs, {rep['rows']} rows in {rep['seconds']}s")

    d = outdir / name
    for f in ("dataset.tsv", "prep_report.txt"):
        derived.blob(f"{prefix}/{f}").upload_from_filename(str(d / f))
    log(f"published gs://{a.derived}/{prefix}/  peak rss {peak_rss_gb():.2f} GB")


# --- mode: finalize ---------------------------------------------------------------------

def do_finalize(a):
    """Decide panel membership from the reports the array produced.

    prepare.py does this at the end of a whole-panel loop, which no single Batch task can
    see. Splitting it out keeps the rule identical -- min_pairs, and a broken fold
    assignment on a dataset that clears the threshold is still fatal.
    """
    _raw, derived = buckets(a)
    cfg = cfgmod.load(a.config)
    k, min_pairs = cfg.cv["k"], cfg.cv["min_pairs"]
    arm = a.arm

    reports = {}
    for blob in derived.client.list_blobs(a.derived, prefix=f"processed/{arm}/"):
        if not blob.name.endswith("prep_report.txt"):
            continue
        cell, name = blob.name.split("/")[2:4]
        r = {}
        for ln in blob.download_as_text().splitlines():
            if ":" in ln:
                key, val = ln.split(":", 1)
                r[key.strip()] = val.strip()
        reports.setdefault(cell, {})[name] = r

    for cell in sorted(reports):
        rs = reports[cell]
        keep, drop = [], []
        for name in sorted(rs):
            r = rs[name]
            pairs = int(r["pairs"])
            problems = eval(r.get("fold_problems", "[]"), {"__builtins__": {}})
            folds = eval(r["fold_proportions"], {"__builtins__": {}})
            if len(folds) != k:
                drop.append((name, pairs, f"fold count {len(folds)} != {k}"))
            elif pairs < min_pairs:
                drop.append((name, pairs, f"pairs<{min_pairs}"))
            elif problems:
                drop.append((name, pairs, "; ".join(problems)))
            else:
                keep.append((name, pairs))
        broken = [(n, p, w) for n, p, w in drop if p >= min_pairs]
        if broken:
            sys.exit(f"FATAL {cell}: clears the threshold but the fold assignment is "
                     f"broken, so the fold map and the data are out of step: {broken}")
        for fname, body in (
            (f"panel_final_{cell}_{arm}.tsv",
             "protein\tcell_line\tpairs\n"
             + "".join(f"{n}\t{cell}\t{p}\n" for n, p in keep)),
            (f"panel_excluded_{cell}_{arm}.tsv",
             "protein\tcell_line\tpairs\treason\n"
             + "".join(f"{n}\t{cell}\t{p}\t{w}\n" for n, p, w in drop)),
        ):
            derived.blob(f"panel/{arm}/{fname}").upload_from_string(body)
        log(f"{cell} {arm}: {len(keep)} of {len(rs)} clear min_pairs={min_pairs}, "
            f"{sum(p for _, p in keep):,} pairs -> panel/{arm}/")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["index", "manifest", "prep", "finalize"])
    p.add_argument("--raw", default=os.environ.get("RAW_BUCKET"))
    p.add_argument("--derived", default=os.environ.get("DERIVED_BUCKET"))
    p.add_argument("--project", default=os.environ.get("GOOGLE_CLOUD_PROJECT"))
    p.add_argument("--config", default=None)
    p.add_argument("--index", type=int, default=None, help="override BATCH_TASK_INDEX")
    p.add_argument("--arm", default="dinuc", choices=list(ARMS), help="finalize only")
    p.add_argument("--force", action="store_true")
    p.add_argument("--dry-run", action="store_true",
                   help="prep only: resolve the task and stop before staging anything")
    a = p.parse_args()
    if not (a.raw and a.derived):
        sys.exit("--raw/--derived or RAW_BUCKET/DERIVED_BUCKET required")
    {"index": do_index, "manifest": do_manifest,
     "prep": do_prep, "finalize": do_finalize}[a.mode](a)


if __name__ == "__main__":
    main()
