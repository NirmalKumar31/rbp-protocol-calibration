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
    python scripts/cloud_variants.py --what phylop     conservation
    python scripts/cloud_variants.py --what all        both, in order
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
OUTPUTS = ["variant_assignments.csv", "variant_conservation.csv"]
MARKER = "variants-complete.json"


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
    peaks = 0
    for r in panel.itertuples():
        cell = getattr(r, "cell_line", None) or r.cell
        for b in raw_bucket.client.list_blobs(raw_bucket.name,
                                              prefix=f"peaks/{cell}/{r.protein}"):
            dest = RAW_DIR / b.name
            if dest.exists():
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            b.download_to_filename(str(dest))
            peaks += 1
    log(f"staged {peaks} peak files")

    # Only the study panel's datasets. Pulling all 189 would move data the stage never reads.
    for r in panel.itertuples():
        dest = PROC_DIR / r.cell_line / r.protein / "dataset.tsv"
        if dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        b = bucket.blob(f"processed/dinuc/{r.cell_line}/{r.protein}/dataset.tsv")
        if not b.exists():
            log(f"WARNING no processed dataset for {r.protein}:{r.cell_line}, skipping")
            continue
        b.download_to_filename(str(dest))
    log(f"staged {sum(1 for _ in PROC_DIR.rglob('dataset.tsv'))} processed datasets")
    return panel


def stage_out(bucket):
    """Upload whatever this stage produced, marker last."""
    sent = []
    for name in OUTPUTS:
        p = TABLES / name
        if not p.exists():
            log(f"WARNING expected output missing: {name}")
            continue
        bucket.blob(f"results/tables/{name}").upload_from_filename(
            str(p), content_type="text/csv")
        sent.append(f"{name} ({p.stat().st_size / 1e6:.1f} MB)")
        log(f"uploaded {name}")
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
    p.add_argument("--what", default="all", choices=["assign", "phylop", "all"])
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

    for what in (["assign", "phylop"] if a.what == "all" else [a.what]):
        run_existing(what)

    stage_out(bucket)
    log("done")


if __name__ == "__main__":
    main()
