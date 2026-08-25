"""Turn the GENCODE GTF into region intervals we can match negatives against.

Two things make this less trivial than it looks:

1. GENCODE emits a generic `UTR` feature, not 5'/3'. Which side a UTR is on has to be
   worked out by comparing it to the transcript's CDS bounds *and* the strand.
2. Introns are not in the file at all. They are the gaps between consecutive exons.

Coordinates: the GTF is 1-based inclusive, BED (and the ENCODE peaks) is 0-based
half-open. Everything here is converted to and stored as 0-based half-open so it lines
up with the peak files without further thought.
"""

import gzip
from collections import defaultdict

import numpy as np

REGIONS = ("utr5", "utr3", "cds", "exon_nc", "intron")
MAIN_CHROMS = tuple(f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"])


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _attr(field, key):
    """Pull one value out of a GTF attribute column without a full parse."""
    i = field.find(key + ' "')
    if i == -1:
        return None
    j = field.index('"', i + len(key) + 2)
    return field[i + len(key) + 2:j]


def parse_transcripts(gtf, chroms=MAIN_CHROMS, want=("exon", "CDS", "UTR")):
    """transcript id -> dict with chrom, strand, and 0-based half-open feature lists."""
    keep = set(want)
    chrom_ok = set(chroms) if chroms else None
    tx = {}
    with _open(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if f[2] not in keep:
                continue
            if chrom_ok is not None and f[0] not in chrom_ok:
                continue
            tid = _attr(f[8], "transcript_id")
            if tid is None:
                continue
            rec = tx.get(tid)
            if rec is None:
                rec = tx[tid] = {"chrom": f[0], "strand": f[6],
                                 "exon": [], "CDS": [], "UTR": []}
            rec[f[2]].append((int(f[3]) - 1, int(f[4])))
    return tx


def split_utrs(utrs, cds, strand):
    """Assign each UTR to the 5' or 3' side using the CDS bounds and the strand."""
    if not cds:
        return [], []
    lo = min(s for s, _ in cds)
    hi = max(e for _, e in cds)
    five, three = [], []
    for s, e in utrs:
        if strand == "+":
            (five if e <= lo else three).append((s, e))
        else:
            (five if s >= hi else three).append((s, e))
    return five, three


def introns_of(exons):
    """Gaps between consecutive exons."""
    if len(exons) < 2:
        return []
    ex = sorted(exons)
    out = []
    prev = ex[0][1]
    for s, e in ex[1:]:
        if s > prev:
            out.append((prev, s))
        prev = max(prev, e)
    return out


def regions_of(rec):
    """region name -> list of 0-based half-open intervals for one transcript."""
    exons, cds, utrs = rec["exon"], rec["CDS"], rec["UTR"]
    utr5, utr3 = split_utrs(utrs, cds, rec["strand"])
    out = {
        "utr5": utr5,
        "utr3": utr3,
        "cds": sorted(cds),
        "exon_nc": sorted(exons) if not cds else [],
        "intron": introns_of(exons),
    }
    return {k: v for k, v in out.items() if v}


def merge_intervals(ivs):
    """Collapse a list of intervals into a sorted, disjoint list."""
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def build_index(gtf, chroms=MAIN_CHROMS):
    """region -> chrom -> (starts, ends): sorted, DISJOINT intervals as numpy arrays.

    Intervals are pooled across all transcripts and merged, which is what makes
    `classify` correct: with disjoint intervals, only the last one starting at or
    before a position can contain it, so one binary search suffices. Keeping them
    unmerged requires an unbounded backwards scan and is easy to get wrong.

    Strand is deliberately dropped. A window's strand comes from its peak, so the
    region's own strand is never needed for classification or negative matching.
    """
    tx = parse_transcripts(gtf, chroms=chroms)
    raw = {r: defaultdict(list) for r in REGIONS}
    for rec in tx.values():
        c = rec["chrom"]
        for region, ivs in regions_of(rec).items():
            for s, e in ivs:
                if e > s:
                    raw[region][c].append((s, e))
    del tx

    index = {}
    for region, per_chrom in raw.items():
        index[region] = {}
        for c, ivs in per_chrom.items():
            m = merge_intervals(ivs)
            index[region][c] = (np.array([s for s, _ in m], np.int64),
                                np.array([e for _, e in m], np.int64))
    return index


def classify(index, chrom, start, end, priority=REGIONS):
    """First region in `priority` order whose interval contains the window midpoint.

    Order matters: one position can sit in a UTR of one isoform and an intron of
    another, so the ranking (utr5, utr3, cds, exon_nc, intron) decides. Exonic classes
    win over intron because they are the more specific claim.
    """
    mid = (start + end) // 2
    for region in priority:
        per = index.get(region, {})
        if chrom not in per:
            continue
        s, e = per[chrom]
        i = int(np.searchsorted(s, mid, side="right")) - 1
        if i >= 0 and s[i] <= mid < e[i]:
            return region
    return None


def stats(index):
    """Interval count and total bases per region, for the EDA report."""
    rows = {}
    for region, per_chrom in index.items():
        n = sum(len(s) for s, _ in per_chrom.values())
        bases = sum(int((e - s).sum()) for s, e in per_chrom.values())
        rows[region] = {"intervals": n, "bases": bases, "chroms": len(per_chrom)}
    return rows


def build_gene_index(gtf, chroms=MAIN_CHROMS):
    """chrom -> (starts, ends, names) for gene features, sorted by start.

    Kept deliberately separate from `build_index`. The region index is already cached and
    every dataset on disk was built from it, so regenerating it would make those datasets
    stale for no reason. Gene identity is a new question, so it gets a new artefact.
    """
    by_chrom = defaultdict(list)
    with _open(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            if chroms and f[0] not in set(chroms):
                continue
            name = _attr(f[8], "gene_name") or _attr(f[8], "gene_id")
            by_chrom[f[0]].append((int(f[3]) - 1, int(f[4]), name))

    index = {}
    for c, rows in by_chrom.items():
        rows.sort()
        index[c] = (np.array([r[0] for r in rows], np.int64),
                    np.array([r[1] for r in rows], np.int64),
                    np.array([r[2] for r in rows], dtype=object))
    return index


def gene_at(gene_index, chrom, pos):
    """Name of the smallest gene containing `pos`, or None.

    Genes overlap -- nested genes, antisense pairs, readthrough transcripts -- so a
    position can fall in several. The smallest containing gene is the most specific
    claim, which is the convention used for annotation elsewhere in the project
    (`classify` prefers the more specific region for the same reason).
    """
    per = gene_index.get(chrom)
    if per is None:
        return None
    starts, ends, names = per
    # every gene starting at or before pos is a candidate; scan back while one can reach
    hi = int(np.searchsorted(starts, pos, side="right"))
    best, best_len = None, None
    for i in range(hi - 1, -1, -1):
        if ends[i] > pos:
            length = int(ends[i] - starts[i])
            if best_len is None or length < best_len:
                best, best_len = names[i], length
        # genes are sorted by start, but a long gene far back can still contain pos, so
        # bound the scan by the largest gene span rather than stopping at the first miss
        if starts[i] < pos - 2_500_000:
            break
    return best


def genes_of(gene_index, chrom, start, end):
    """Every gene overlapping [start, end)."""
    per = gene_index.get(chrom)
    if per is None:
        return []
    starts, ends, names = per
    hi = int(np.searchsorted(starts, end, side="left"))
    out = []
    for i in range(hi - 1, -1, -1):
        if ends[i] > start:
            out.append(names[i])
        if starts[i] < start - 2_500_000:
            break
    return out
