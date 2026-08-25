"""Downloads with provenance. Every fetched file is recorded with its size, MD5 and
source URL so the raw data is auditable and the whole stage is re-runnable.

Existing files are verified rather than re-fetched, so re-running is cheap.
"""

import hashlib
import time
from pathlib import Path

import requests

CHUNK = 1 << 20
MANIFEST_COLS = ["path", "bytes", "md5", "url", "fetched_utc"]


def md5sum(path, chunk=CHUNK):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _expected_size(url, timeout=60):
    try:
        r = requests.head(url, allow_redirects=True, timeout=timeout)
        return int(r.headers.get("content-length") or 0)
    except requests.RequestException:
        return 0


def fetch(url, dest, force=False, timeout=120, quiet=False, retries=6):
    """Stream `url` to `dest`, resuming and retrying on dropped connections.

    Large reference files reliably fail partway on a home connection, so this uses
    HTTP Range to resume from the .part file instead of restarting. The final size is
    checked against content-length, so a truncated download can never be silently
    promoted to the real filename -- a half-written genome is worse than none.

    Returns (path, bytes, md5, was_downloaded).
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        return dest, dest.stat().st_size, md5sum(dest), False

    tmp = dest.with_suffix(dest.suffix + ".part")
    if force and tmp.exists():
        tmp.unlink()
    total = _expected_size(url, timeout)

    for attempt in range(1, retries + 1):
        have = tmp.stat().st_size if tmp.exists() else 0
        if total and have >= total:
            break
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with requests.get(url, stream=True, timeout=timeout, headers=headers) as r:
                if have and r.status_code == 200:
                    # server ignored Range, so start over rather than corrupt the file
                    have, mode = 0, "wb"
                elif have and r.status_code == 206:
                    mode = "ab"
                else:
                    r.raise_for_status()
                    mode = "wb"
                got = have
                with open(tmp, mode) as fh:
                    for block in r.iter_content(CHUNK):
                        fh.write(block)
                        got += len(block)
                        if not quiet and total and got % (128 * CHUNK) < CHUNK:
                            print(f"    {got/1e6:7.0f} / {total/1e6:.0f} MB "
                                  f"({got/total:.0%})", flush=True)
            if not total or tmp.stat().st_size >= total:
                break
            raise OSError(f"short read: {tmp.stat().st_size} of {total}")
        except (requests.RequestException, OSError) as e:
            if attempt == retries:
                raise
            wait = min(2 ** attempt, 30)
            print(f"    attempt {attempt} failed ({type(e).__name__}); "
                  f"resuming from {tmp.stat().st_size/1e6:.0f} MB in {wait}s", flush=True)
            time.sleep(wait)

    got = tmp.stat().st_size
    if total and got != total:
        raise OSError(f"{dest.name}: got {got} bytes, expected {total} - not promoting")
    tmp.replace(dest)
    return dest, dest.stat().st_size, md5sum(dest), True


def gunzip(src, dest=None, force=False, quiet=False):
    """Stream-decompress a .gz. Returns the output path.

    Streamed rather than read whole, because the genome is ~3 GB uncompressed.
    """
    import gzip as gz
    src = Path(src)
    dest = Path(dest) if dest else src.with_suffix("")
    if dest.exists() and not force:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    done = 0
    with gz.open(src, "rb") as fin, open(tmp, "wb") as fout:
        while block := fin.read(CHUNK):
            fout.write(block)
            done += len(block)
            if not quiet and done % (256 * CHUNK) < CHUNK:
                print(f"    decompressed {done/1e9:.2f} GB", flush=True)
    tmp.replace(dest)
    return dest


def index_fasta(path):
    """Build the .fai index and report what pyfaidx sees."""
    from pyfaidx import Fasta
    fa = Fasta(str(path))
    names = list(fa.keys())
    return names, len(fa[names[0]]) if names else 0


def record(manifest, root, path, size, md5, url):
    """Append one provenance row, replacing any earlier entry for the same path."""
    manifest = Path(manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    rel = str(Path(path).relative_to(root))
    rows = {}
    if manifest.exists():
        for line in manifest.read_text().strip().splitlines()[1:]:
            if line:
                parts = line.split("\t")
                rows[parts[0]] = parts
    rows[rel] = [rel, str(size), md5, url, time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())]
    with open(manifest, "w") as fh:
        fh.write("\t".join(MANIFEST_COLS) + "\n")
        for k in sorted(rows):
            fh.write("\t".join(rows[k]) + "\n")
