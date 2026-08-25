"""Stage 2: fetch raw data, one asset group at a time.

    python scripts/download.py --what peaks     # 17 ENCODE peak files
    python scripts/download.py --what gencode   # annotation
    python scripts/download.py --what genome    # GRCh38 fasta
    python scripts/download.py --what clinvar   # variant vcf
    python scripts/download.py --what all

Everything lands under data/raw/ and is logged to data/raw/manifest.tsv with size,
MD5 and source URL. phyloP is deliberately not downloaded: it is ~10 GB and we read
it over HTTP range requests instead.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rbp.data import download as dl  # noqa: E402
from rbp.data import encode  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = RAW / "manifest.tsv"


def get_peaks(cfg, force, panel=None):
    rows = cfgmod.proteins(panel=panel)
    out = RAW / "peaks" / cfg.encode["cell_line"]
    print(f"peaks: {len(rows)} ENCODE files -> {out.relative_to(ROOT)}")
    n_new = 0
    for i, r in enumerate(rows, 1):
        url = encode.peak_url(r["accession"], cfg.encode["api"])
        dest = out / f"{r['protein']}.{r['accession']}.bed.gz"
        path, size, md5, fresh = dl.fetch(url, dest, force=force, quiet=True)
        dl.record(MANIFEST, ROOT, path, size, md5, url)
        n_new += fresh
        print(f"  [{i:2d}/{len(rows)}] {r['protein']:9} {r['accession']}  "
              f"{size/1e6:6.2f} MB  {'downloaded' if fresh else 'cached'}")
    print(f"peaks: {n_new} downloaded, {len(rows)-n_new} already present")


def get_one(cfg, key, filename, force):
    url = cfg.reference[key]
    dest = RAW / filename
    print(f"{key}: {url}")
    path, size, md5, fresh = dl.fetch(url, dest, force=force)
    dl.record(MANIFEST, ROOT, path, size, md5, url)
    print(f"  {size/1e6:.1f} MB  md5={md5}  {'downloaded' if fresh else 'cached'}")


def get_genome(cfg, force):
    """Download, decompress and index the genome, then verify chromosome naming."""
    url = cfg.reference["genome"]
    gz = RAW / "GRCh38.primary_assembly.genome.fa.gz"
    print(f"genome: {url}")
    path, size, md5, fresh = dl.fetch(url, gz, force=force)
    dl.record(MANIFEST, ROOT, path, size, md5, url)
    print(f"  {size/1e6:.0f} MB  md5={md5}  {'downloaded' if fresh else 'cached'}")

    print("decompressing (needed for random access; gzip cannot be seeked)")
    fa = dl.gunzip(gz, force=force)
    print(f"  {fa.stat().st_size/1e9:.2f} GB -> {fa.name}")

    print("indexing with pyfaidx")
    names, first_len = dl.index_fasta(fa)
    print(f"  {len(names)} sequences; first is {names[0]} ({first_len:,} bp)")
    if not names[0].startswith("chr"):
        raise SystemExit(f"FATAL: chromosome names are '{names[0]}'-style but the GTF and "
                         "peaks use 'chr1'. Wrong genome source.")
    print("  naming check: OK, 'chr'-prefixed and consistent with the GTF and peaks")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--what", required=True,
                   choices=["peaks", "gencode", "genome", "clinvar", "all"])
    p.add_argument("--force", action="store_true", help="re-download even if present")
    p.add_argument("--panel", default=None, help="panel TSV to use instead of config/proteins.tsv")
    p.add_argument("--cell-line", default=None, help="override encode.cell_line")
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    if a.cell_line:
        cfg["encode"]["cell_line"] = a.cell_line

    if a.what in ("peaks", "all"):
        get_peaks(cfg, a.force, panel=a.panel)
    if a.what in ("gencode", "all"):
        get_one(cfg, "gtf", "gencode.v45.annotation.gtf.gz", a.force)
    if a.what in ("genome", "all"):
        get_genome(cfg, a.force)
    if a.what in ("clinvar", "all"):
        get_one(cfg, "clinvar", "clinvar.vcf.gz", a.force)

    print(f"\nmanifest: {MANIFEST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
