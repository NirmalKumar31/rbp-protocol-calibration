"""Stage 11. The ClinVar arm's inputs, built in a container instead of on a laptop.

WHAT WAS LAPTOP-ONLY AND WHY IT MATTERED. rehearsal_variants.py reads the ClinVar VCF, the
3.1 GB genome and every processed dataset straight off local disk, and writes to
results/tables/. That is fine for exploration and fatal for reproducibility: two of the
paper's four results depended on files that existed only on one machine.

THE DESIGN, and it is deliberately boring: STAGE IN, RUN THE EXISTING CODE UNCHANGED, STAGE
OUT. Inputs are pulled from GCS onto container disk at the paths the existing code already
expects, the tested functions run untouched, and the outputs are uploaded. Reimplementing
the assignment and conservation logic against a cloud filesystem would have meant a second
copy of subtle code -- strand handling, window offsets, ref-allele checks -- diverging
silently from the copy that produced the published numbers. Copying bytes is cheap; copying
logic is how two versions of a result appear.

WHY THIS STAGE NEEDS AN EXTERNAL IP. phyloP conservation is fetched from UCSC by HTTP range
request over a bigwig. Workers get Private Google Access only, which routes to Google APIs
and nothing else, so a normal worker cannot reach UCSC at all. This stage runs on a single
short-lived VM with a public IP -- the same posture as ingest -- rather than adding a Cloud
NAT that would bill per hour and hand internet access to every other worker for no reason.

    python scripts/cloud_variants.py --what assign     variants near peaks
    python scripts/cloud_variants.py --what score      k-mer delta per variant
    python scripts/cloud_variants.py --what phylop     conservation
    python scripts/cloud_variants.py --what all        all three, in order
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402

from rbp.utils import cloud as cloudcfg  # noqa: E402

# The container's working tree. The existing code takes these paths as given, so we
# materialise GCS into them rather than teaching it about buckets.
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"

# Inputs the variant stages read, and where they live in the raw bucket.
RAW_OBJECTS = [
    "clinvar.vcf.gz",
    "GRCh38.primary_assembly.genome.fa",
    "GRCh38.primary_assembly.genome.fa.fai",
]

# Outputs this stage is responsible for. The completion marker is written LAST, after every
# payload, for the same reason as everywhere else in this project: a task preempted between
# two uploads must redo its work, not be skipped by a marker that arrived early.
OUTPUTS = ["variant_assignments.csv", "variant_scores.csv", "variant_conservation.csv"]
MARKER = "variants-complete.json"

# Which sub-stage produces which file. The chain is assign -> score -> phylop, and each link
# costs real time: assign ~10 min over 189 datasets, score ~12 min over 66,010 variants,
# phylop an unknown number of HTTP round trips against UCSC.
#
# WHY THIS MAPPING EXISTS AT ALL. It did not, and every phyloP failure threw away the
# twenty-two minutes of correct assign and score work that preceded it, because run_existing()
# raises on a nonzero return code and that exception skipped stage_out() entirely. Four
# consecutive runs recomputed identical assignments and identical scores to reach the same
# failing line. The results were never wrong; they were simply never saved.
#
# rehearsal_variants.py's own docstring promises each stage is "cached to disk so a rerun is
# cheap", and on a laptop that was true -- the disk persisted. A container's disk does not, so
# the cache has to live in GCS or it does not exist.
STAGE_OUTPUT = {
    "assign": "variant_assignments.csv",
    "score": "variant_scores.csv",
    "phylop": "variant_conservation.csv",
}

# assign writes this too; it is not part of the completion contract but the analysis reads it.
EXTRA_OUTPUTS = ["variant_availability_panel.csv"]

# phyloP's own cache: (chrom, pos) -> score, one line each. Round-tripping it through GCS
# means a run that dies halfway through 27,492 remote lookups resumes at the position it
# reached rather than at zero. This is the only stage whose cost is external and rate-limited,
# so it is the one stage where partial progress is worth carrying.
PHYLOP_CACHE = ROOT / "data" / "interim" / "phylop_cache.tsv"
PHYLOP_CACHE_OBJ = "interim/phylop_cache.tsv"


def log(m):
    print(f"[cloud_variants] {m}", flush=True)


def stage_in(bucket, raw_bucket):
    """Pull raw inputs and the study panel's processed datasets onto container disk."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    for name in RAW_OBJECTS:
        dest = RAW_DIR / name
        if dest.exists():
            log(f"have {name}")
            continue
        # get_blob(), not blob() + exists().
        #
        # bucket.blob(name) builds a LOCAL reference and fetches nothing; exists() issues a
        # HEAD and returns a bool without populating metadata. So b.size stayed None and the
        # progress line -- b.size / 1e9 -- raised TypeError before a single byte downloaded.
        # A log line took the stage down. get_blob() does one GET and returns a hydrated blob,
        # or None if it is absent, which also collapses the existence check into the same call.
        b = raw_bucket.get_blob(name)
        if b is None:
            raise SystemExit(f"missing gs://{raw_bucket.name}/{name} -- stage 3 incomplete")
        log(f"downloading {name} ({(b.size or 0) / 1e9:.2f} GB)")
        b.download_to_filename(str(dest))

    panel = pd.read_csv(
        __import__("io").StringIO(
            bucket.blob("manifest/study_panel.tsv").download_as_text()), sep="\t")
    log(f"study panel: {len(panel)} datasets")

    # THE ENCODE PEAK FILES, which this stage needs and I forgot.
    #
    # rehearsal_variants.py --what assign reads peaks to find where each protein binds before
    # it can decide which ClinVar variants sit near a site. RAW_OBJECTS listed the genome, its
    # index and the VCF -- the three files the phyloP half needs -- and nothing else, so the
    # assign half died with "no peak file for AATF in /app/data/raw/peaks/K562" after
    # successfully staging 95 datasets. A stage-in list is a contract with every function the
    # stage calls, not just the one you had in mind while writing it.
    # ALL the peaks, not just the study panel's.
    #
    # My first attempt staged peaks only for the 95 datasets in the study panel, which failed
    # on ADAT1: rehearsal_variants.py --what assign walks the FULL candidate panel, not the
    # study subset, so it needs peaks for proteins the study never scores. Two ways to fix
    # that -- teach the assign stage about the study panel, or stage everything it might read.
    # The whole design of this file is to run the existing code UNCHANGED, so it stages
    # everything. The peaks are a few hundred MB against a 3.1 GB genome already being pulled;
    # narrowing them saved nothing and cost two failed runs.
    peaks = 0
    for b in raw_bucket.client.list_blobs(raw_bucket.name, prefix="peaks/"):
        dest = RAW_DIR / b.name
        if dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        b.download_to_filename(str(dest))
        peaks += 1
    log(f"staged {peaks} peak files")

    # EVERY processed dataset in the dinuc arm, not just the study panel's 95.
    #
    # Third time this exact mistake bit, so it is worth stating as a rule rather than a fix.
    # stage_score walks variant_assignments.csv, which covers every dataset that has a variant
    # near a peak -- not the study subset -- because a variant's nearest binding site can
    # belong to a protein the study never scores. Staging only the panel's datasets fails on
    # the first protein outside it, exactly as staging only the panel's peaks failed on ADAT1.
    #
    # 462 MiB against a 3.15 GB genome already being pulled. Every attempt to be selective
    # about staging has cost a failed run and saved nothing measurable.
    staged = 0
    for b in bucket.client.list_blobs(bucket.name, prefix="processed/dinuc/"):
        if not b.name.endswith("dataset.tsv"):
            continue
        # processed/dinuc/{cell}/{protein}/dataset.tsv -> data/processed/{cell}/{protein}/
        parts = b.name.split("/")
        dest = PROC_DIR / parts[2] / parts[3] / "dataset.tsv"
        if dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        b.download_to_filename(str(dest))
        staged += 1
    log(f"staged {sum(1 for _ in PROC_DIR.rglob('dataset.tsv'))} processed datasets")

    # Anything a previous attempt already finished. Downloading these is what turns a rerun
    # into a resume: main() skips any sub-stage whose output is already on disk.
    resumed = []
    for name in OUTPUTS + EXTRA_OUTPUTS:
        b = bucket.get_blob(f"results/tables/{name}")
        if b is None:
            continue
        b.download_to_filename(str(TABLES / name))
        resumed.append(name)
    if resumed:
        log(f"resuming with {len(resumed)} table(s) from a previous attempt: "
            + ", ".join(resumed))

    cb = bucket.get_blob(PHYLOP_CACHE_OBJ)
    if cb is not None:
        PHYLOP_CACHE.parent.mkdir(parents=True, exist_ok=True)
        cb.download_to_filename(str(PHYLOP_CACHE))
        log(f"phyloP cache restored ({(cb.size or 0) / 1e6:.1f} MB)")

    return panel


def stage_out(bucket):
    """Upload whatever this stage produced, marker last.

    Called even when a sub-stage failed, which is the point: partial work is still correct
    work. The marker is what encodes completeness, so uploading two of three tables is safe
    and saves the next attempt from recomputing them.
    """
    sent = []
    for name in OUTPUTS + EXTRA_OUTPUTS:
        p = TABLES / name
        if not p.exists():
            log(f"WARNING expected output missing: {name}")
            continue
        bucket.blob(f"results/tables/{name}").upload_from_filename(
            str(p), content_type="text/csv")
        if name in OUTPUTS:
            sent.append(f"{name} ({p.stat().st_size / 1e6:.1f} MB)")
        log(f"uploaded {name}")

    if PHYLOP_CACHE.exists():
        bucket.blob(PHYLOP_CACHE_OBJ).upload_from_filename(
            str(PHYLOP_CACHE), content_type="text/tab-separated-values")
        log(f"uploaded phyloP cache ({PHYLOP_CACHE.stat().st_size / 1e6:.1f} MB)")
    if len(sent) == len(OUTPUTS):
        import json
        bucket.blob(MARKER).upload_from_string(
            json.dumps({"outputs": sent, "git_sha": os.environ.get("GIT_SHA", "unknown")}),
            content_type="application/json")
        log(f"wrote completion marker {MARKER}")
    else:
        log("NOT writing the completion marker: some outputs are missing, so a rerun "
            "must redo the work rather than skip it")


def run_existing(what):
    """Call rehearsal_variants.py in-place. Same code, same numbers, different filesystem."""
    cmd = [sys.executable, str(ROOT / "scripts" / "rehearsal_variants.py"), "--what", what]
    log("running: " + " ".join(cmd[1:]))
    rc = subprocess.run(cmd, cwd=str(ROOT)).returncode
    if rc != 0:
        raise SystemExit(f"rehearsal_variants.py --what {what} failed (rc={rc})")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--what", default="all",
                   choices=["assign", "score", "phylop", "all"])
    p.add_argument("--force", action="store_true")
    a = p.parse_args()

    bucket = cloudcfg.bucket()
    raw_bucket = cloudcfg.client().bucket(cloudcfg.raw_bucket())
    log(cloudcfg.describe())

    if bucket.blob(MARKER).exists() and not a.force:
        log(f"{MARKER} already present, nothing to do")
        return

    TABLES.mkdir(parents=True, exist_ok=True)
    stage_in(bucket, raw_bucket)

    # ASSIGN, SCORE, PHYLOP -- and score is not optional.
    #
    # I ran assign then phylop, and phylop died on a missing variant_scores.csv. The stages
    # form a chain: assign finds which variants sit near a binding site, score computes the
    # k-mer disruption delta for each, and phylop reads that file to know which distinct
    # positions need a conservation lookup. Skipping the middle one breaks the third.
    #
    # score also produces the k-mer arm of the R4 ladder -- the bottom rung, against which the
    # mismatched and matched SpliceBERT heads are compared. Leaving it out would have removed
    # the baseline the whole result is measured from.
    # stage_out runs whether or not the chain completes. Without the finally, a failure in
    # phylop discarded a finished assign and a finished score -- see STAGE_OUTPUT above.
    failure = None
    try:
        for what in (["assign", "score", "phylop"] if a.what == "all" else [a.what]):
            done = TABLES / STAGE_OUTPUT[what]
            if done.exists() and not a.force:
                log(f"skipping {what}: {done.name} already exists ({done.stat().st_size / 1e6:.1f} MB)")
                continue
            run_existing(what)
    except SystemExit as e:
        failure = e
    finally:
        stage_out(bucket)

    if failure is not None:
        log("chain did not complete; finished stages were saved and will be skipped on rerun")
        raise failure
    log("done")


if __name__ == "__main__":
    main()
