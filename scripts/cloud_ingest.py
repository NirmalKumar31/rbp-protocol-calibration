"""Stage 1+2 in the cloud: discover the panel, fetch every raw input, publish to GCS.

WHY THE CLOUD FETCHES THIS RATHER THAN THE LAPTOP UPLOADING IT. The genome alone is 3.1 GB.
A home connection uploads that in tens of minutes; a GCP VM pulls it from GENCODE in about
one. Ingress to GCS is free. And the path exercised here is the one the pipeline actually
uses, so it is tested rather than bypassed.

WHY PYTHON RATHER THAN A SHELL SCRIPT CALLING gcloud. The CPU image is built on
python:3.13-slim, which has no gcloud CLI. Installing it would add ~200 MB to an image that
374 tasks will pull. `google-cloud-storage` is already a dependency, is faster (resumable
uploads, no subprocess per file), and can be unit tested.

IDEMPOTENT BY DESIGN. Every artifact is checked for existence in GCS before downloading. A
rerun after a spot preemption fetches only what is missing, which matters because this runs
on a preemptible VM.

    python scripts/cloud_ingest.py --bucket $PROJECT_ID-raw
"""

import argparse
import hashlib
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import requests  # noqa: E402

from rbp.data import encode  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
WORK = Path(os.environ.get("WORK_DIR", "/tmp/raw"))

REFERENCE = {
    "GRCh38.primary_assembly.genome.fa.gz":
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/"
        "GRCh38.primary_assembly.genome.fa.gz",
    "gencode.v45.annotation.gtf.gz":
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/"
        "gencode.v45.annotation.gtf.gz",
    "clinvar.vcf.gz":
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def stream_to_gcs(url, blob, chunk=8 << 20, retries=6):
    """Download a URL straight into a GCS object, without holding it all in memory.

    A 3.1 GB genome would not fit comfortably in the container's memory and writing it to
    the boot disk first doubles the I/O. `upload_from_file` on a streaming response hands
    GCS the socket directly.
    """
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                r.raw.decode_content = True
                blob.upload_from_file(r.raw, content_type="application/octet-stream")
            return True
        except Exception as e:                                   # noqa: BLE001
            wait = min(60, 2 ** attempt)
            log(f"    attempt {attempt}/{retries} failed ({e}); retrying in {wait}s")
            time.sleep(wait)
    return False


def ensure(bucket, name, url):
    """Fetch `url` to `name` in the bucket unless it is already there."""
    blob = bucket.blob(name)
    if blob.exists():
        log(f"  skip (present): {name}")
        return False
    log(f"  downloading {name}")
    if not stream_to_gcs(url, blob):
        raise RuntimeError(f"failed to fetch {name} after retries")
    blob.reload()
    log(f"  published {name}  ({blob.size / 1e6:.1f} MB)")
    return True


def main():
    from google.cloud import storage

    p = argparse.ArgumentParser()
    p.add_argument("--bucket", default=os.environ.get("RAW_BUCKET"))
    p.add_argument("--config", default=None)
    p.add_argument("--skip-reference", action="store_true",
                   help="peaks and panel only; useful for a cheap smoke test")
    a = p.parse_args()
    if not a.bucket:
        sys.exit("--bucket or RAW_BUCKET required")

    cfg = cfgmod.load(a.config)
    client = storage.Client()
    bucket = client.bucket(a.bucket)
    WORK.mkdir(parents=True, exist_ok=True)

    # ---- 1. panel discovery, from the live API -----------------------------------------
    log("=== 1. panel discovery from the live ENCODE API ===")
    panels = {}
    for cell, fname in (("K562", "panel_full.tsv"), ("HepG2", "panel_full_HepG2.tsv")):
        cfg["encode"]["cell_line"] = cell
        rows, _dupes, _skipped = encode.build_panel(cfg)
        panels[cell] = rows
        cols = ["protein", "accession", "cell_line", "experiment", "n_replicates"]
        body = "\t".join(cols) + "\n" + "\n".join(
            "\t".join(r[c] for c in cols) for r in rows) + "\n"
        bucket.blob(f"config/{fname}").upload_from_string(body)
        log(f"  {cell}: {len(rows)} proteins -> config/{fname}")

    # ---- 2. reference data --------------------------------------------------------------
    if not a.skip_reference:
        log("=== 2. reference data ===")
        for name, url in REFERENCE.items():
            ensure(bucket, name, url)
    else:
        log("=== 2. reference data SKIPPED (--skip-reference) ===")

    # ---- 2b. decompress and index the genome ONCE ---------------------------------------
    # pyfaidx needs an uncompressed, indexed FASTA for random access. Doing this here means
    # it happens once; leaving it to preprocessing would mean 374 tasks each inflating
    # 845 MB to 3.1 GB and rebuilding the same index.
    if not a.skip_reference:
        log("=== 2b. decompress and index the genome ===")
        fa_name = "GRCh38.primary_assembly.genome.fa"
        if bucket.blob(f"{fa_name}.fai").exists():
            log("  skip (present): indexed genome")
        else:
            import gzip
            import shutil
            gz_local = WORK / f"{fa_name}.gz"
            fa_local = WORK / fa_name
            log("  fetching compressed genome to local disk")
            bucket.blob(f"{fa_name}.gz").download_to_filename(str(gz_local))
            log("  decompressing")
            with gzip.open(gz_local, "rb") as fh_in, open(fa_local, "wb") as fh_out:
                shutil.copyfileobj(fh_in, fh_out, length=32 << 20)
            gz_local.unlink()
            log(f"  {fa_local.stat().st_size / 1e9:.2f} GB uncompressed; indexing")
            from pyfaidx import Fasta
            fasta = Fasta(str(fa_local))
            names = list(fasta.keys())
            # The same assertion the local pipeline makes. Ensembl names chromosomes "1",
            # GENCODE and the ENCODE peaks use "chr1"; mixing them makes every sequence
            # lookup silently return nothing.
            if not names[0].startswith("chr"):
                raise SystemExit(f"FATAL: chromosome names are '{names[0]}'-style, not chr-prefixed")
            log(f"  indexed {len(names)} sequences, first is {names[0]}")
            for f in (fa_local, Path(str(fa_local) + ".fai")):
                log(f"  uploading {f.name} ({f.stat().st_size / 1e9:.2f} GB)")
                blob = bucket.blob(f.name)
                blob.chunk_size = 32 << 20        # resumable upload for a multi-GB object
                blob.upload_from_filename(str(f))
                f.unlink()

    # ---- 3. peaks -----------------------------------------------------------------------
    log("=== 3. peak files, both cell lines ===")
    api = cfg.encode["api"]
    for cell, rows in panels.items():
        new = 0
        for i, r in enumerate(rows, 1):
            name = f"peaks/{cell}/{r['protein']}.{r['accession']}.bed.gz"
            if bucket.blob(name).exists():
                continue
            url = encode.peak_url(r["accession"], api)
            resp = requests.get(url, timeout=120)
            resp.raise_for_status()
            bucket.blob(name).upload_from_string(resp.content)
            new += 1
            if new % 25 == 0:
                log(f"    {cell}: {new} downloaded ({i}/{len(rows)} checked)")
        log(f"  {cell}: {new} newly downloaded, {len(rows) - new} already present")

    # ---- 4. manifest --------------------------------------------------------------------
    # Provenance. "We used ClinVar" is unfalsifiable without this: ClinVar is re-released
    # weekly and the file behind that URL changes. The md5 pins exactly which one we used.
    log("=== 4. manifest with checksums ===")
    lines = ["path\tsize_bytes\tmd5_base64\tupdated_utc"]
    total = 0
    for blob in client.list_blobs(a.bucket):
        if blob.name == "manifest.tsv":
            continue
        total += blob.size or 0
        lines.append(f"{blob.name}\t{blob.size}\t{blob.md5_hash}\t{blob.updated:%Y-%m-%dT%H:%M:%SZ}")
    bucket.blob("manifest.tsv").upload_from_string("\n".join(lines) + "\n")
    log(f"  {len(lines) - 1} objects, {total / 1e9:.2f} GB total")

    log("=== ingest complete ===")


if __name__ == "__main__":
    main()
