"""Build the positive windows: fixed-size sequences centred on eCLIP peaks.

A peak is a genomic interval, but the protein binds RNA. On the minus strand the RNA
is the reverse complement of the reference, so every window is strand-corrected before
it is handed to a model. Windows that are the wrong length or contain an unknown base
are dropped rather than padded, so no model ever sees a synthetic base.
"""

import gzip

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def to_rna(dna, strand):
    """Strand-correct and transcribe: what the protein actually sees."""
    return (revcomp(dna) if strand == "-" else dna).upper().replace("T", "U")


def gc_content(seq):
    s = seq.upper()
    n = len(s)
    return (s.count("G") + s.count("C")) / n if n else 0.0


def read_peaks(path):
    """Yield (chrom, start, end, strand) from a narrowPeak, skipping malformed lines."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            yield f[0], int(f[1]), int(f[2]), f[5].strip()


def window_bounds(start, end, size):
    """Half-open bounds of the size-nt window centred on a peak midpoint."""
    mid = (start + end) // 2
    half = size // 2
    return mid - half, mid + half + 1


def fetch(fasta, chrom, w0, w1):
    """Reference sequence for a window, or None if out of bounds."""
    if chrom not in fasta:
        return None
    if w0 < 0 or w1 > len(fasta[chrom]):
        return None
    return fasta[chrom][w0:w1].seq.upper()


def build_positives(peak_path, fasta, region_index, size, classify, drop_n=True,
                    chroms=None):
    """Positive windows for one protein.

    Returns (rows, dropped) where dropped counts the reason each peak was discarded,
    so the loss is reported rather than silently absorbed.
    """
    rows = []
    dropped = {"excluded_chrom": 0, "out_of_bounds": 0, "wrong_length": 0,
               "has_n": 0, "no_region": 0, "duplicate": 0}
    seen = set()
    allowed = set(chroms) if chroms else None

    for chrom, start, end, strand in read_peaks(peak_path):
        if allowed is not None and chrom not in allowed:
            dropped["excluded_chrom"] += 1
            continue
        w0, w1 = window_bounds(start, end, size)
        dna = fetch(fasta, chrom, w0, w1)
        if dna is None:
            dropped["out_of_bounds"] += 1
            continue
        if len(dna) != size:
            dropped["wrong_length"] += 1
            continue
        if drop_n and "N" in dna:
            dropped["has_n"] += 1
            continue
        region = classify(region_index, chrom, w0, w1)
        if region is None:
            dropped["no_region"] += 1
            continue
        key = (chrom, w0)
        if key in seen:
            dropped["duplicate"] += 1
            continue
        seen.add(key)
        rows.append({
            "chrom": chrom, "start": w0, "end": w1, "strand": strand,
            "region": region, "gc": round(gc_content(dna), 4),
            "seq_dna": dna, "seq_rna": to_rna(dna, strand),
        })
    return rows, dropped
