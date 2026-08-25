"""Fetch UCSC phyloP conservation scores for variant positions.

Reads the hg38 phyloP100way bigWig over HTTP range requests, so nothing is downloaded:
the full track is ~10 GB and we need a few thousand single positions from it. bigWig
supports random access, so this is the right trade.

Positions are 1-based here, matching VCF and ClinVar. Everything else in the project is
0-based, so callers convert at the boundary -- `clinvar.load` provides `pos_vcf` for
exactly this.
"""

import os

import numpy as np

DEFAULT_URL = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/"
               "hg38.phyloP100way.bw")
SPAN = 4096


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

    bw = pyBigWig.open(url)
    try:
        chrom_sizes = bw.chroms()
        done, i = 0, 0
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

            try:
                block = bw.values(chrom, start, end)
            except (RuntimeError, OverflowError):
                block = None                       # a bad range must not kill the run

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

    save_cache(cache, values)
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
