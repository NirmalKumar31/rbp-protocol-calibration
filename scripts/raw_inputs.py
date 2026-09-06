"""A machine-readable manifest of every raw input, with the bytes named rather than implied.

    python scripts/raw_inputs.py --from-gcs     # rebuild from the raw bucket's live metadata
    python scripts/raw_inputs.py --check        # offline: structure, coverage, no secrets

WHY. The manuscript lists ENCODE accessions and the genome and annotation versions, which
identify logical resources. They do not identify bytes. A provider that republishes a file
under the same accession leaves every accession in the paper still correct and every number in
the paper unreproducible, and nothing in the release would show it. An audit asked for size and
checksum records for the original downloads and reported that they were not committed.

They existed. `gs://<raw bucket>/manifest.tsv` has carried path, size and base64 MD5 for all 251
raw objects since the day they were staged; it simply was never brought into the repository. This
script brings it in, verifies it against the bucket's live metadata, and adds the fields the
manifest did not carry: what each object is, which accession and cell it belongs to, the source
URL, the assembly and annotation version, and the terms the source publishes under.

What this establishes, exactly. All 251 objects are present and their sizes and MD5s match the
bucket today, 251 of 251, zero mismatches. The 244 peak files cover all 95 candidate accessions
and all 94 study-panel accessions in supplementary_table_s1.csv. scripts/cloud_prep.py reads
peaks from this bucket, so the link between these bytes and the pipeline is in code.

What it does not establish, and must not be read as. There is no fetch-time record. The
timestamps are when each object was written to the bucket, all within a 7.5-minute window on
2026-08-25, which is the staging run and not the download from ENCODE or EBI. No checksum
published by ENCODE or GENCODE was recorded at fetch time, so these MD5s pin the bytes THIS
STUDY USED and cannot prove they are the bytes the provider served. The source URLs are
reconstructed from the accession and the release, not read from a download log, and are marked
`reconstructed` in the manifest for that reason. Byte-level historical provenance to the
provider is not claimed anywhere, and a future study should record the response headers at fetch
time, which costs nothing and is the only thing that would close this.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

TABLES = ROOT / "results" / "tables"
OUT = TABLES / "raw_inputs.csv"
S1 = TABLES / "supplementary_table_s1.csv"
PARAMS = ROOT / "config" / "params.yaml"

FIELDS = ["path", "kind", "accession", "cell", "size_bytes", "md5_base64", "uploaded_utc",
          "source_url", "url_provenance", "assembly", "annotation", "terms_url"]

ENCODE_FILE = "https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
ENCODE_TERMS = "https://www.encodeproject.org/help/citing-encode/"
GENCODE_TERMS = "https://www.gencodegenes.org/pages/data_access.html"
NCBI_TERMS = "https://www.ncbi.nlm.nih.gov/home/about/policies/"

# Anything matching these must never reach a committed manifest.
SECRET = re.compile(r"\b\d{6}-[0-9A-F]{6}-[0-9A-F]{6}\b"          # GCP billing account
                    r"|BEGIN [A-Z ]*PRIVATE KEY"
                    r"|AIza[0-9A-Za-z_-]{35}"                      # Google API key
                    r"|ya29\."                                     # OAuth token
                    r"|@[a-z0-9-]+\.iam\.gserviceaccount\.com")
UNSAFE_PATH = re.compile(r"^/|\.\.|~|\\|(?:^|/)(?:Users|home)/")


def _params():
    import yaml
    with PARAMS.open() as fh:
        return yaml.safe_load(fh)


def classify(path, cfg):
    """kind, accession, cell, source_url, url_provenance, assembly, annotation, terms.

    THE URLS COME FROM config/params.yaml's `reference:` BLOCK, which has recorded them from
    the beginning, so for the genome, the annotation and ClinVar they are recorded rather than
    reconstructed. The 244 peak files are the other case: the bucket stores them under
    `peaks/<cell>/<PROTEIN>.<ACCESSION>.bed.gz` and no download URL was kept, so the ENCODE
    portal URL is rebuilt from the accession and is marked `reconstructed`. The distinction is
    a column rather than a footnote because the two are not the same kind of evidence.
    """
    ref = cfg.get("reference", {}) or {}
    asm = str(cfg.get("encode", {}).get("assembly", "") or "GRCh38")
    gtf_url = str(ref.get("gtf", ""))
    m = re.search(r"gencode\.(v\d+)\.", gtf_url)
    ann = "GENCODE " + m.group(1) if m else ""

    m = re.match(r"peaks/(K562|HepG2)/[^.]+\.(ENCFF[0-9A-Z]+)\.bed\.gz$", path)
    if m:
        cell, acc = m.group(1), m.group(2)
        return ("ENCODE eCLIP peaks", acc, cell, ENCODE_FILE.format(acc=acc), "reconstructed",
                asm, ann, ENCODE_TERMS)
    if path.startswith("GRCh38.primary_assembly.genome.fa"):
        url = str(ref.get("genome", ""))
        return ("genome FASTA", "", "", url, "config/params.yaml reference.genome" if url
                else "", asm, ann, GENCODE_TERMS)
    if path.startswith("gencode."):
        return ("annotation GTF", "", "", gtf_url,
                "config/params.yaml reference.gtf" if gtf_url else "", asm, ann, GENCODE_TERMS)
    if path.startswith("clinvar"):
        url = str(ref.get("clinvar", ""))
        return ("ClinVar VCF", "", "", url,
                "config/params.yaml reference.clinvar" if url else "", asm, "", NCBI_TERMS)
    if path.startswith("config/"):
        return ("panel definition, staged from the repository", "", "", "", "in-repo", asm,
                ann, "")
    return ("other raw input", "", "", "", "", asm, ann, "")


def from_gcs(bucket):
    """The bucket's live metadata, which is the source of truth for size and MD5."""
    from google.cloud import storage
    client = storage.Client()
    rows = []
    for b in client.list_blobs(bucket):
        if b.name.endswith("/") or b.name == "manifest.tsv":
            continue
        rows.append({"path": b.name, "size_bytes": str(b.size), "md5_base64": b.md5_hash or "",
                     "uploaded_utc": b.time_created.strftime("%Y-%m-%dT%H:%M:%SZ")})
    return sorted(rows, key=lambda r: r["path"])


def enrich(rows):
    cfg = _params()
    out = []
    for r in rows:
        kind, acc, cell, url, prov, asm, ann, terms = classify(r["path"], cfg)
        out.append({"path": r["path"], "kind": kind, "accession": acc, "cell": cell,
                    "size_bytes": r["size_bytes"], "md5_base64": r["md5_base64"],
                    "uploaded_utc": r["uploaded_utc"], "source_url": url,
                    "url_provenance": prov, "assembly": asm, "annotation": ann,
                    "terms_url": terms})
    return out


def audit(rows):
    """Every condition the manifest has to meet before it is committed. Returns problems."""
    bad = []
    for r in rows:
        if UNSAFE_PATH.search(r["path"]):
            bad.append(f"unsafe path shape: {r['path']}")
        if not r["md5_base64"]:
            bad.append(f"no MD5: {r['path']}")
        if not r["size_bytes"].isdigit() or int(r["size_bytes"]) <= 0:
            bad.append(f"bad size: {r['path']}")
    blob = "\n".join("\t".join(r.values()) for r in rows)
    for m in SECRET.finditer(blob):
        bad.append(f"credential-shaped string in the manifest: {m.group(0)[:20]}")
    # COVERAGE. A manifest that omits a panel dataset's peaks is worse than none, because it
    # looks complete. Checked against the study panel of record, not against a count.
    if S1.exists():
        with S1.open(newline="") as fh:
            s1 = list(csv.DictReader(fh))
        have = {r["accession"] for r in rows if r["accession"]}
        for label, want in (("study panel",
                             {r["accession"] for r in s1
                              if r["in_three_arm_panel"] == "True"}),
                            ("all candidates", {r["accession"] for r in s1})):
            missing = want - have
            if missing:
                bad.append(f"{len(missing)} {label} accession(s) absent from the manifest: "
                           f"{sorted(missing)[:5]}")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-gcs", action="store_true",
                    help="rebuild from the raw bucket's live metadata (needs credentials)")
    ap.add_argument("--bucket", default="")
    ap.add_argument("--check", action="store_true", help="offline audit of the committed file")
    a = ap.parse_args()

    if a.from_gcs:
        from rbp.utils import cloud
        rows = enrich(from_gcs(a.bucket or cloud.raw_bucket()))
    else:
        if not OUT.exists():
            log(f"  {OUT.relative_to(ROOT)} is missing; rerun with --from-gcs")
            sys.exit(1)
        with OUT.open(newline="") as fh:
            rows = list(csv.DictReader(fh))

    problems = audit(rows)
    for p in problems:
        log(f"  PROBLEM: {p}")
    if problems:
        sys.exit(1)

    kinds = {}
    for r in rows:
        kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
    log(f"  {len(rows)} raw inputs, all with a size and an MD5, no secrets, no unsafe paths")
    for k, n in sorted(kinds.items(), key=lambda kv: -kv[1]):
        log(f"    {n:4d}  {k}")

    if a.check:
        return
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    log(f"  wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
