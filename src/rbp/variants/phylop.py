"""Fetch UCSC phyloP conservation scores for variant positions.

Reads the hg38 phyloP100way bigWig over HTTP range requests, so nothing is downloaded:
the full track is ~10 GB and we need a few thousand single positions from it. bigWig
supports random access, so this is the right trade.

Positions are 1-based here, matching VCF and ClinVar. Everything else in the project is
0-based, so callers convert at the boundary -- `clinvar.load` provides `pos_vcf` for
exactly this.
"""

import os
import time

import numpy as np

DEFAULT_URL = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/"
               "hg38.phyloP100way.bw")
SPAN = 4096

# UCSC IS A SHARED PUBLIC SERVER AND IT FLAKES. This is the only external dependency in the
# whole pipeline, and a single refused connection used to kill the stage two seconds in --
# taking twenty-two minutes of completed assign and score work with it, because the failure
# propagated before anything was saved.
#
# Proven transient, not a misconfiguration: a diagnostic job on the identical image digest,
# region and network posture opened the same URL and read real values minutes later. So the
# right response is to retry rather than to redesign.
#
# Two fallbacks, tried in rotation. hgdownload.cse is the older hostname for the same host,
# and plain http avoids a TLS handshake; both resolve to the same file, so a score fetched
# through either is the same number.
OPEN_ATTEMPTS = 6
BLOCK_ATTEMPTS = 3

# How much missing conservation is tolerable before the stage is a failure rather than a
# result. The reference run annotated 27,491 positions with ZERO misses, so anything above a
# fraction of a percent means the track was unreachable for part of the run, not that the
# genome lacks scores there.
#
# This threshold exists because the block-level except clause below turns a network error into
# NaN. That is correct for a genuinely unscorable position and silently corrupting for a
# dropped connection -- and conservation is the CONTROL variable in R4, so quietly losing it
# would weaken the headline result instead of failing loudly. A tolerant inner guard needs a
# strict outer check.
MAX_FRAC_MISSING = 0.01


def _candidates(url):
    out = [url]
    if "soe.ucsc.edu" in url:
        out.append(url.replace("soe.ucsc.edu", "cse.ucsc.edu"))
    if url.startswith("https://"):
        out.append("http://" + url[len("https://"):])
    return out


def _open(url):
    """Open the remote bigwig, retrying with backoff across mirror hostnames."""
    import pyBigWig

    if not pyBigWig.remote:
        raise RuntimeError(
            "pyBigWig was built WITHOUT remote support (pyBigWig.remote == 0), so it cannot "
            "open a URL at all. The PyPI manylinux wheel is built without libcurl; install "
            "libcurl4-openssl-dev and reinstall with --no-binary pyBigWig.")

    urls = _candidates(url)
    last = None
    for i in range(OPEN_ATTEMPTS):
        u = urls[i % len(urls)]
        try:
            bw = pyBigWig.open(u)
            if bw is not None:
                if i:
                    print(f"  phyloP: opened {u} on attempt {i + 1}", flush=True)
                return bw
            last = f"pyBigWig.open returned None for {u}"
        except Exception as e:                      # noqa: BLE001 - any curl error retries
            last = f"{type(e).__name__}: {e}"
        if i < OPEN_ATTEMPTS - 1:
            wait = min(30, 2 ** i)
            print(f"  phyloP: open failed ({last}); retrying in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(
        f"could not open the phyloP track after {OPEN_ATTEMPTS} attempts. Last error: {last}")


def _norm(chrom):
    return chrom if str(chrom).startswith("chr") else f"chr{chrom}"


def load_cache(path):
    """Cached (chrom, pos) -> value. Remote reads are the slow part, so we keep them."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            c, p, v = line.rstrip("\n").split("\t")
            out[(c, int(p))] = float(v) if v not in ("", "nan") else np.nan
    return out


def save_cache(path, values):
    if not path:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        for (c, p), v in sorted(values.items()):
            fh.write(f"{c}\t{p}\t{v}\n")


def fetch(positions, url=DEFAULT_URL, cache=None, progress=True):
    """Map (chrom, pos) -> phyloP score. Unmappable positions give NaN.

    Nearby positions are grouped into one remote read of SPAN bases and indexed into,
    which cuts HTTP requests by roughly two orders of magnitude. Variants cluster in
    genes, so the grouping wins often.
    """
    import pyBigWig

    want = {(_norm(c), int(p)) for c, p in positions}
    values = {k: v for k, v in load_cache(cache).items() if k in want}
    todo = sorted(want - set(values))
    if not todo:
        return values

    bw = _open(url)
    try:
        chrom_sizes = bw.chroms()
        done, i, failed_blocks = 0, 0, 0
        while i < len(todo):
            chrom = todo[i][0]
            size = chrom_sizes.get(chrom)
            if size is None:                       # chromosome absent from the track
                while i < len(todo) and todo[i][0] == chrom:
                    values[todo[i]] = np.nan
                    i += 1
                continue

            start = todo[i][1] - 1
            end = min(start + SPAN, size)
            batch = []
            while i < len(todo) and todo[i][0] == chrom and todo[i][1] - 1 < end:
                batch.append(todo[i])
                i += 1

            # A bad range must not kill the run, but a dropped connection must not silently
            # become NaN either: retry the block first, and count what still fails so the
            # caller can tell "unscorable" from "unreachable".
            block = None
            for attempt in range(BLOCK_ATTEMPTS):
                try:
                    block = bw.values(chrom, start, end)
                    break
                except (RuntimeError, OverflowError):
                    if attempt < BLOCK_ATTEMPTS - 1:
                        time.sleep(1 + attempt)
            if block is None:
                failed_blocks += 1

            for c, p in batch:
                if block is None:
                    values[(c, p)] = np.nan
                else:
                    off = (p - 1) - start
                    v = block[off] if 0 <= off < len(block) else None
                    values[(c, p)] = np.nan if v is None else float(v)

            done += len(batch)
            if progress and done % 2000 < len(batch):
                print(f"  phyloP: {done}/{len(todo)}", flush=True)
    finally:
        bw.close()
        # Save even on the way out of an exception. The cache is the only thing that makes a
        # partial run worth anything, and remote reads are the expensive part.
        save_cache(cache, values)

    if failed_blocks:
        print(f"  phyloP: {failed_blocks} block(s) unreadable after {BLOCK_ATTEMPTS} tries",
              flush=True)
    return values


def annotate(df, chrom_col="chrom", pos_col="pos_vcf", out_col="conservation",
             url=DEFAULT_URL, cache=None):
    """Add a phyloP column to a variant dataframe.

    Defaults to `pos_vcf` because the track is 1-based. Passing a 0-based column here
    shifts every score by one position, which is silent and wrong, so the default points
    at the column that is already in the right convention.
    """
    vals = fetch(zip(df[chrom_col], df[pos_col]), url=url, cache=cache)
    df = df.copy()
    df[out_col] = [vals.get((_norm(c), int(p)), np.nan)
                   for c, p in zip(df[chrom_col], df[pos_col])]
    return df


def coverage(df, col="conservation"):
    """How much of the data actually got a score. Reported rather than assumed."""
    n = len(df)
    miss = int(df[col].isna().sum())
    return {"n": n, "missing": miss,
            "frac_missing": round(miss / n, 4) if n else 0.0}


def assert_coverage(df, col="conservation", max_frac=MAX_FRAC_MISSING):
    """Refuse a conservation table that is too sparse to serve as R4's control.

    Called by the pipeline rather than by fetch(), because the library should report and the
    pipeline should decide. See MAX_FRAC_MISSING for why this is strict.
    """
    cov = coverage(df, col)
    if cov["frac_missing"] > max_frac:
        raise SystemExit(
            f"phyloP coverage too low to use: {cov['missing']:,}/{cov['n']:,} positions "
            f"missing ({cov['frac_missing']:.2%}), limit {max_frac:.2%}. The reference run "
            f"missed zero, so this means the track was unreachable for part of the run. "
            f"Rerun -- the cache makes it resume.")
    return cov
