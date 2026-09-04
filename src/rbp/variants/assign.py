"""Attach ClinVar variants to a protein's binding sites and build their scoring windows.

Three details here are easy to get wrong and each one silently corrupts the result:

  * STRAND. The protein binds RNA, so a variant near a minus-strand peak has to be read
    as the reverse complement. Taking the reference strand throughout would score half
    the panel against a sequence the protein never sees. The strand comes from the
    nearest peak.
  * THE REFERENCE BASE. ClinVar's REF must match the genome at that position. When it
    does not -- wrong assembly, an indel misparsed as a substitution, a position off by
    one -- the "alternate" window differs from the reference in a way that has nothing to
    do with the variant. These are counted and dropped, never silently accepted.
  * THE FOLD. A variant's chromosome belongs to exactly one CV fold, and the variant must
    be scored by the model that did not train on that chromosome. Assigned here so the
    scoring code cannot get it wrong.
"""

import numpy as np

from ..data.windows import to_rna, window_bounds


def peak_index(peak_path, chroms=None):
    """chrom -> (starts, ends, strands) sorted by start, for nearest-peak lookup."""
    from ..data.windows import read_peaks
    per = {}
    keep = set(chroms) if chroms else None
    for chrom, start, end, strand in read_peaks(peak_path):
        if keep is not None and chrom not in keep:
            continue
        per.setdefault(chrom, []).append((start, end, strand))
    out = {}
    for chrom, rows in per.items():
        rows.sort()
        out[chrom] = (np.array([r[0] for r in rows]),
                      np.array([r[1] for r in rows]),
                      np.array([r[2] for r in rows], dtype=object))
    return out


def nearest_peak(index, chrom, pos):
    """(distance, strand) to the closest peak, or (inf, None) if the chromosome is absent.

    Distance is zero inside a peak. Only the two candidates bracketing the position are
    checked, which is enough because the intervals are sorted by start and we take the
    minimum of "gap to the one before" and "gap to the one after".
    """
    ent = index.get(chrom)
    if ent is None:
        return np.inf, None
    starts, ends, strands = ent
    i = int(np.searchsorted(starts, pos, side="right"))
    best, strand = np.inf, None
    for j in (i - 1, i):
        if not 0 <= j < len(starts):
            continue
        d = 0 if starts[j] <= pos < ends[j] else min(abs(pos - ends[j] + 1),
                                                     abs(starts[j] - pos))
        if d < best:
            best, strand = d, strands[j]
    return best, strand


def assign(variants, index, margin, fold_map):
    """Keep variants within `margin` of a peak, tagging strand, distance and fold.

    `variants` is an iterable of dicts from clinvar.load.
    """
    out = []
    for v in variants:
        d, strand = nearest_peak(index, v["chrom"], v["pos"])
        if d > margin:
            continue
        fold = fold_map.get(v["chrom"])
        if fold is None:                 # excluded chromosome, e.g. chrY
            continue
        out.append({**v, "peak_distance": int(d), "strand": strand, "fold": int(fold)})
    return out


def windows_for(fasta, chrom, pos, ref, alt, size, shifts, strand):
    """(ref_seqs, alt_seqs) as RNA, one pair per shift. Returns (None, None) on mismatch.

    The variant is placed `shift` bases from the window centre rather than always at it,
    because the discriminative signal sits 12-26 nt 5' of the peak midpoint we anchor on
    (an eCLIP artefact: reverse transcription stops at the crosslink, so called peaks
    extend 3' of the real site). One centred window would systematically miss it.
    """
    if chrom not in fasta:
        return None, None
    L = len(fasta[chrom])
    base = fasta[chrom][pos:pos + 1].seq.upper()
    if base != ref.upper():
        return None, None                       # assembly or parsing mismatch

    refs, alts = [], []
    for sh in shifts:
        w0, w1 = window_bounds(pos + sh, pos + sh, size)
        if w0 < 0 or w1 > L:
            continue
        dna = fasta[chrom][w0:w1].seq.upper()
        if len(dna) != size or "N" in dna:
            continue
        i = pos - w0
        if not 0 <= i < size or dna[i] != ref.upper():
            continue
        mut = dna[:i] + alt.upper() + dna[i + 1:]
        refs.append(to_rna(dna, strand))
        alts.append(to_rna(mut, strand))
    if not refs:
        return None, None
    return refs, alts


def build_scoring_table(assigned, fasta, size, shifts):
    """Flatten variants into one row per (variant, shift), ready for batch scoring.

    Returns (table, dropped) where `dropped` counts why rows were lost, so the loss is
    reported rather than absorbed. Flattening lets every window for the whole panel be
    scored in one vectorised pass instead of per variant.
    """
    rows = []
    dropped = {"ref_mismatch": 0, "no_usable_window": 0}
    for v in assigned:
        refs, alts = windows_for(fasta, v["chrom"], v["pos"], v["ref"], v["alt"],
                                 size, shifts, v["strand"])
        if refs is None:
            base = (fasta[v["chrom"]][v["pos"]:v["pos"] + 1].seq.upper()
                    if v["chrom"] in fasta else "")
            key = "ref_mismatch" if base and base != v["ref"].upper() \
                else "no_usable_window"
            dropped[key] += 1
            continue
        for j, (r, a) in enumerate(zip(refs, alts)):
            rows.append({"vid": v["vid"], "label": v["label"], "chrom": v["chrom"],
                         "pos": v["pos"], "pos_vcf": v["pos_vcf"], "fold": v["fold"],
                         "strand": v["strand"], "peak_distance": v["peak_distance"],
                         "shift_idx": j, "seq_ref": r, "seq_alt": a})
    return rows, dropped


def collapse_delta(vids, per_window_delta, how="max_abs"):
    """One delta per variant from its per-shift deltas.

    `max_abs` takes the largest disruption over the shifts, which is the pre-registered
    choice: the true binding site's offset from the peak midpoint is unknown per variant,
    so the maximum asks "is there ANY register in which this variant matters". Averaging
    would dilute a real effect with windows that do not contain the site.
    """
    vids = np.asarray(vids)
    d = np.asarray(per_window_delta, dtype=float)
    order = np.argsort(vids, kind="mergesort")
    vids, d = vids[order], d[order]
    edges = np.flatnonzero(np.r_[True, vids[1:] != vids[:-1], True])
    keys, vals = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        chunk = d[a:b]
        chunk = chunk[np.isfinite(chunk)]
        if len(chunk) == 0:
            keys.append(vids[a])
            vals.append(np.nan)
            continue
        if how == "max_abs":
            vals.append(float(chunk[np.argmax(np.abs(chunk))]))
        elif how == "mean_abs":
            vals.append(float(np.mean(np.abs(chunk))))
        else:
            raise ValueError(f"unknown collapse: {how}")
        keys.append(vids[a])
    return keys, vals
