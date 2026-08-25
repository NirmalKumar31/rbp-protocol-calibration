"""Matched negatives: the design choice that makes the binding task meaningful.

Random genomic background would let a model score well by separating transcribed
sequence from empty regions, learning nothing about the protein. So each positive gets
one negative matched on:

  * the same region class (a 3'UTR positive gets a 3'UTR negative)
  * GC content within a tolerance
  * the same chromosome, so the pair cannot straddle a train/test split
  * at least `min_peak_distance` from ANY peak of that protein, so a "negative" is not
    an unlabelled binding site

With those held constant, a model that still separates the classes has learned
something specific to the protein rather than a shortcut.
"""

import numpy as np

from .windows import fetch, gc_content, to_rna


def merge(intervals):
    if not intervals:
        return []
    iv = sorted(intervals)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def exclusion_zones(peaks, margin):
    """chrom -> merged peak intervals padded by `margin` on both sides."""
    by_chrom = {}
    for chrom, start, end, _ in peaks:
        by_chrom.setdefault(chrom, []).append((max(0, start - margin), end + margin))
    return {c: merge(v) for c, v in by_chrom.items()}


def subtract(intervals, holes):
    """Interval difference: what is left of `intervals` after removing `holes`."""
    if not holes:
        return list(intervals)
    holes = merge(holes)
    out = []
    hs = [h[0] for h in holes]
    for s, e in intervals:
        cur = s
        i = int(np.searchsorted(hs, cur, side="right")) - 1
        i = max(i, 0)
        while i < len(holes) and holes[i][0] < e:
            h0, h1 = holes[i]
            if h1 <= cur:
                i += 1
                continue
            if h0 > cur:
                out.append((cur, min(h0, e)))
            cur = max(cur, h1)
            if cur >= e:
                break
            i += 1
        if cur < e:
            out.append((cur, e))
    return [(s, e) for s, e in out if e > s]


def available(region_index, region, exclusion, size):
    """chrom -> intervals of `region` that are clear of peaks and fit a whole window."""
    out = {}
    for chrom, (starts, ends) in region_index.get(region, {}).items():
        ivs = [(int(s), int(e)) for s, e in zip(starts, ends) if e - s >= size]
        if not ivs:
            continue
        free = subtract(ivs, exclusion.get(chrom, []))
        free = [(s, e) for s, e in free if e - s >= size]
        if free:
            out[chrom] = free
    return out


def _weighted_pick(rng, intervals, size):
    """Uniformly sample a window start across the usable span of `intervals`."""
    spans = np.array([e - s - size + 1 for s, e in intervals], dtype=np.int64)
    total = int(spans.sum())
    if total <= 0:
        return None
    k = int(rng.integers(0, total))
    idx = int(np.searchsorted(np.cumsum(spans), k, side="right"))
    s, _ = intervals[idx]
    offset = k - int(np.cumsum(spans)[idx - 1]) if idx else k
    return s + offset


def sample_negative(fasta, pool, chrom, target_gc, size, tolerance, rng,
                    tries=40, relax_after=25, drop_n=True):
    """One negative window matching `target_gc`, or None if the pool is too thin.

    GC is checked by actually reading the sequence, so a candidate is only accepted
    once it genuinely matches. After `relax_after` misses the tolerance widens, which
    keeps GC-extreme positives from being dropped outright.
    """
    ivs = pool.get(chrom)
    if not ivs:
        return None
    best = None
    for t in range(tries):
        w0 = _weighted_pick(rng, ivs, size)
        if w0 is None:
            return None
        dna = fetch(fasta, chrom, w0, w0 + size)
        if dna is None or len(dna) != size or (drop_n and "N" in dna):
            continue
        gc = gc_content(dna)
        gap = abs(gc - target_gc)
        if best is None or gap < best[0]:
            best = (gap, w0, dna, gc)
        tol = tolerance if t < relax_after else tolerance * 2
        if gap <= tol:
            return {"start": w0, "end": w0 + size, "gc": round(gc, 4), "seq_dna": dna}
    if best is not None and best[0] <= tolerance * 3:
        _, w0, dna, gc = best
        return {"start": w0, "end": w0 + size, "gc": round(gc, 4), "seq_dna": dna}
    return None


# ---------------------------------------------------------------------------------------
# Dinucleotide-matched negatives
# ---------------------------------------------------------------------------------------
#
# GC matching constrains G+C as a single total and leaves everything else free. Measured on
# our own data, that leaves a lot of room: TARDBP positives run GU at +1.94 log2 over their
# GC-matched negatives. A model can score 0.98 by detecting "GU-rich" without learning a
# motif, and the reported AUROC cannot tell the difference.
#
# The standard fix in this literature is to build negatives by dinucleotide-preserving
# SHUFFLE of the positives. We tested that and abandoned it for two independent reasons:
# TARDBP's motif survives the shuffle (UGUGU in 92.2% of positives, still 87.0% after)
# because the motif IS a dinucleotide repeat, and shuffled sequence is detectable per se at
# AUROC ~0.69 with no positives involved, so the arm has a floor that cannot be subtracted
# out.
#
# Matching instead of shuffling avoids both. We sample REAL genomic windows whose
# dinucleotide composition is close to the positive's. Both classes stay real sequence, so
# there is no shuffle artefact, and the motif is not reconstructed by the procedure because
# the procedure never touches the positive.
#
# The cost is that matching 16 frequencies is much harder than matching one number, so the
# achieved match quality has to be measured and reported rather than assumed.

ALPHABET = "ACGT"
DINUCS = [a + b for a in ALPHABET for b in ALPHABET]
_DI_INDEX = {d: i for i, d in enumerate(DINUCS)}


def dinuc_vector(seq, normalise=True):
    """16-vector of dinucleotide frequencies, or raw counts with normalise=False.

    THE COUNTS MODE EXISTS FOR REPRODUCIBILITY, NOT SPEED. A frequency is k/100 for a 101-nt
    window, and almost none of those are exactly representable in float64. The last-bit
    error that leaves behind is not the same on every CPU, and the nearest-neighbour search
    below compares distances that are frequently exactly tied -- so the tie is resolved by
    whichever way the rounding happened to fall. Measured on QKI: the same code, same numpy
    and same scipy chose a different negative for 8.3% of pairs on arm64 versus x86-64, of
    which 83% were exact ties and the rest were knock-on effects of the greedy assignment.
    Counts are integers, exactly representable, so their L1 distances are exact integers and
    identical on every machine.
    """
    v = np.zeros(16, dtype=np.float64)
    s = seq.upper()
    n = 0
    for i in range(len(s) - 1):
        j = _DI_INDEX.get(s[i:i + 2])
        if j is not None:
            v[j] += 1.0
            n += 1
    if not normalise:
        return v
    return v / n if n else v


def dinuc_matrix(seqs, normalise=True):
    return (np.vstack([dinuc_vector(s, normalise) for s in seqs]) if seqs
            else np.zeros((0, 16)))


def candidate_pool(fasta, pool, chrom, size, n_want, rng, drop_n=True, max_tries=None,
                   normalise=True):
    """Sample real windows from the allowed pool and return (starts, seqs, dinuc matrix).

    Built once per (region, chromosome) and then queried by nearest neighbour, rather than
    re-sampling per positive. Sampling 300 candidates for each of 24,401 positives would be
    seven million sequence reads for one protein; this is a few thousand.
    """
    ivs = pool.get(chrom)
    if not ivs:
        return np.zeros(0, dtype=np.int64), [], np.zeros((0, 16))
    max_tries = max_tries or n_want * 3
    starts, seqs, seen = [], [], set()
    for _ in range(max_tries):
        if len(starts) >= n_want:
            break
        w0 = _weighted_pick(rng, ivs, size)
        if w0 is None or w0 in seen:
            continue
        seen.add(w0)
        dna = fetch(fasta, chrom, w0, w0 + size)
        if dna is None or len(dna) != size or (drop_n and "N" in dna):
            continue
        starts.append(w0)
        seqs.append(dna)
    return np.array(starts, dtype=np.int64), seqs, dinuc_matrix(seqs, normalise)


def build_negatives_dinuc(positives, peaks, fasta, region_index, size,
                          min_peak_distance, seed=7, drop_n=True, pool_multiple=8,
                          pool_min=1500, max_l1=None):
    """One dinucleotide-matched negative per positive.

    Greedy nearest-neighbour assignment: for each positive, take the closest unused
    candidate by L1 distance over the 16 dinucleotide frequencies. Greedy rather than
    optimal because an exact assignment over tens of thousands of points buys very little
    here and costs a great deal; the achieved distances are reported so the quality of the
    approximation is visible.

    Returns (rows, dropped, distances) where `distances` is the achieved L1 per pair. L1 on
    frequency vectors runs 0 to 2, so 0.05 means the average dinucleotide frequency differs
    by about 0.3 percentage points.
    """
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(seed)
    excl = exclusion_zones(peaks, min_peak_distance)
    regions = {p["region"] for p in positives}
    pools = {r: available(region_index, r, excl, size) for r in regions}

    # group positives by (region, chromosome) so one candidate pool serves many of them
    buckets = {}
    for i, p in enumerate(positives):
        buckets.setdefault((p["region"], p["chrom"]), []).append(i)

    rows = [None] * len(positives)
    dists = np.full(len(positives), np.nan)
    dropped = {"no_pool": 0, "no_match": 0, "too_far": 0}
    # Scoped across ALL buckets, not per bucket. A genomic interval can be annotated as
    # both exon_nc and utr3 in different transcripts, so it appears in two region pools and
    # a per-bucket `used` set let the same window be chosen twice as a negative. Measured
    # cost before the fix: 20 duplicated rows in 1,797,048 (0.001%), so it changes no
    # result -- but a duplicated negative is double-counted and should not happen.
    used_global = set()

    for (region, chrom), idxs in buckets.items():
        pool = pools.get(region, {})
        if chrom not in pool:
            dropped["no_pool"] += len(idxs)
            continue
        n_want = max(pool_min, pool_multiple * len(idxs))
        # Counts, not frequencies -- see dinuc_vector. Everything from here to the reported
        # distance is integer arithmetic in float64, which is exact.
        starts, seqs, cand = candidate_pool(fasta, pool, chrom, size, n_want, rng,
                                            drop_n=drop_n, normalise=False)
        if len(starts) == 0:
            dropped["no_pool"] += len(idxs)
            continue

        target = dinuc_matrix([positives[i]["seq_dna"] for i in idxs], normalise=False)
        # A 17TH COLUMN THAT MAKES TIES IMPOSSIBLE.
        #
        # Integer counts removed the floating-point noise, but not the ties themselves --
        # hundreds of candidate windows can sit at exactly the same L1 distance. cKDTree
        # returns only the k nearest, and when more than k of them tie, WHICH k it returns
        # is an implementation detail that differs between the arm64 and x86-64 builds of
        # scipy. Re-sorting the returned window cannot fix that: the tied candidate it
        # should have picked was never in the window.
        #
        # So each candidate carries a tiny tiebreak, strictly increasing in genomic order.
        # Integer distances differ by at least 1 and the tiebreak stays below 0.5, so it can
        # never reorder candidates that genuinely differ; among equals it orders by
        # position. No ties means the k nearest are uniquely determined, so every platform
        # returns the same set. Targets carry 0 in that column.
        rank = np.argsort(np.argsort(starts)) / max(len(starts), 1) * 0.5
        cand = np.hstack([cand, rank.reshape(-1, 1)])
        target = np.hstack([target, np.zeros((len(idxs), 1))])
        tree = cKDTree(cand)
        # query more neighbours than we need so a used candidate can be skipped
        k = min(len(starts), 40)
        dist, nbr = tree.query(target, k=k, p=1)
        # reshape, NOT atleast_2d. With k == 1 -- a pool holding a single usable window --
        # query returns one value per positive as a flat (n,) array, and atleast_2d turns
        # that into (1, n) rather than (n, 1), so every row after the first reads another
        # positive's neighbour. Pre-existing; it needs a one-candidate pool to show up.
        dist = dist.reshape(len(idxs), -1)
        nbr = nbr.reshape(len(idxs), -1)
        # Exact distances still leave genuine ties, and which of them cKDTree returns first
        # is an implementation detail. Re-order each row by (distance, genomic start) so the
        # winner is decided by the data rather than by the tree's traversal.
        n_di = size - 1

        for row, i in enumerate(idxs):
            ok = nbr[row] < len(starts)
            ci = nbr[row][ok].astype(np.int64)
            order = np.lexsort((starts[ci], dist[row][ok]))
            chosen = None
            for j in ci[order]:
                j = int(j)
                if (chrom, int(starts[j])) in used_global:
                    continue
                chosen = j
                break
            if chosen is None:
                dropped["no_match"] += 1
                continue
            # Reported in frequency units, as before, so the thresholds and every published
            # dinuc_l1 figure keep their meaning. The tiebreak column is excluded: it is a
            # sorting device, not part of the composition distance.
            l1 = float(np.abs(cand[chosen, :16] - target[row, :16]).sum()) / n_di
            if max_l1 is not None and l1 > max_l1:
                dropped["too_far"] += 1
                continue
            used_global.add((chrom, int(starts[chosen])))
            p = positives[i]
            dna = seqs[chosen]
            dists[i] = l1
            rows[i] = {
                "chrom": chrom, "start": int(starts[chosen]),
                "end": int(starts[chosen]) + size, "strand": p["strand"],
                "region": region, "gc": round(gc_content(dna), 4),
                "seq_dna": dna, "seq_rna": to_rna(dna, p["strand"]),
            }
    return rows, dropped, dists


def build_negatives(positives, peaks, fasta, region_index, size, tolerance,
                    min_peak_distance, seed=7, drop_n=True):
    """One matched negative per positive. Returns (rows, dropped_counter)."""
    rng = np.random.default_rng(seed)
    excl = exclusion_zones(peaks, min_peak_distance)
    pools = {r: available(region_index, r, excl, size)
             for r in {p["region"] for p in positives}}

    rows, dropped = [], {"no_pool": 0, "no_match": 0}
    used = set()
    for p in positives:
        pool = pools.get(p["region"], {})
        if p["chrom"] not in pool:
            dropped["no_pool"] += 1
            rows.append(None)
            continue
        neg = None
        for _ in range(3):
            cand = sample_negative(fasta, pool, p["chrom"], p["gc"], size,
                                   tolerance, rng, drop_n=drop_n)
            if cand is None:
                break
            if (p["chrom"], cand["start"]) in used:
                continue
            neg = cand
            break
        if neg is None:
            dropped["no_match"] += 1
            rows.append(None)
            continue
        used.add((p["chrom"], neg["start"]))
        rows.append({
            "chrom": p["chrom"], "start": neg["start"], "end": neg["end"],
            "strand": p["strand"], "region": p["region"], "gc": neg["gc"],
            "seq_dna": neg["seq_dna"], "seq_rna": to_rna(neg["seq_dna"], p["strand"]),
        })
    return rows, dropped
