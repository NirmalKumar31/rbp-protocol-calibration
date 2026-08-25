"""Positive control: prove the conservation-controlled test can detect a real effect.

A null result only means something if the test could have found a signal that was there.
Power analysis argues that theoretically; this argues it empirically, by manufacturing
variants whose effect we know and checking the test recovers them.

For each bound window containing a canonical literature motif we build two single-base
mutants:

    disruptive  the motif core base is changed          (label 1)
    neutral     the SAME base substitution applied      (label 0)
                >= min_distance nt from any motif hit

Same window, same substitution type, differing only in whether it lands on the motif.
Motifs come from the literature rather than from the model, so the control is not
circular -- picking "whichever position the model reacts to most" would be testing the
model against itself.

Conservation is the real phyloP at the mutated genomic position, so the full controlled
regression runs exactly as it does on ClinVar.

Expected outcome: a large delta coefficient that survives the conservation control. If it
does not, the method lacks power and no null result is interpretable.
"""

import numpy as np
import pandas as pd


def find_all(seq, motif):
    """Every start index of `motif` in `seq`, including overlapping occurrences."""
    hits, i = [], seq.find(motif)
    while i != -1:
        hits.append(i)
        i = seq.find(motif, i + 1)
    return hits


def genomic_pos(start, end, strand, i, one_based=True):
    """Genomic position of index `i` in a strand-corrected window.

    The window sequence is already reverse-complemented for minus-strand peaks, so index
    0 is the 5' end in RNA terms on both strands. Mapping back to the genome therefore
    has to invert that for the minus strand, or every minus-strand conservation lookup
    lands at the wrong end of the window.
    """
    pos0 = start + i if strand == "+" else end - 1 - i
    return pos0 + 1 if one_based else pos0


def build_pairs(df, motif, offset, replacement, min_distance=25, seed=0):
    """One disruptive and one matched neutral mutant per usable window.

    A window is usable when it contains the motif, the motif base to change is not
    already the replacement, and some position >= min_distance from every motif hit
    carries the same reference base -- so the neutral mutation is the same substitution,
    just somewhere uninformative.
    """
    rng = np.random.default_rng(seed)
    rows = []
    dropped = {"no_motif": 0, "already_replacement": 0, "no_neutral_site": 0}

    for r in df.itertuples():
        seq = r.seq_rna
        hits = find_all(seq, motif)
        if not hits:
            dropped["no_motif"] += 1
            continue

        mi = hits[0] + offset
        ref_base = seq[mi]
        if ref_base == replacement:
            dropped["already_replacement"] += 1
            continue

        forbidden = set()
        for h in hits:
            forbidden.update(range(h - min_distance, h + len(motif) + min_distance))
        cands = [j for j, b in enumerate(seq) if b == ref_base and j not in forbidden]
        if not cands:
            dropped["no_neutral_site"] += 1
            continue
        nj = int(rng.choice(cands))

        for kind, label, idx, alt_seq in (
            ("disruptive", 1, mi, seq[:mi] + replacement + seq[mi + 1:]),
            ("neutral", 0, nj, seq[:nj] + replacement + seq[nj + 1:]),
        ):
            rows.append({
                "id": r.id, "chrom": r.chrom, "start": r.start, "end": r.end,
                "strand": r.strand, "kind": kind, "label": label,
                "mut_index": idx, "ref_base": ref_base, "alt_base": replacement,
                "pos_vcf": genomic_pos(r.start, r.end, r.strand, idx),
                "ref_seq": seq, "alt_seq": alt_seq,
            })
    return pd.DataFrame(rows), dropped


def score_pairs(pairs, score_fn, batch_size=64):
    """Add p_ref, p_alt and delta. Reference windows are scored once each.

    Each window appears twice (once per mutant) but its reference score is identical, so
    deduplicating halves the forward passes.
    """
    uniq = pairs.drop_duplicates("ref_seq")["ref_seq"].tolist()
    p_ref = dict(zip(uniq, score_fn(uniq, batch_size=batch_size)))
    p_alt = score_fn(pairs["alt_seq"].tolist(), batch_size=batch_size)

    out = pairs.copy()
    out["p_ref"] = [p_ref[s] for s in out.ref_seq]
    out["p_alt"] = p_alt
    out["delta"] = (out.p_ref - out.p_alt).abs()
    out["delta_signed"] = out.p_alt - out.p_ref
    return out


def effect_size(scored):
    """Cohen's d between disruptive and neutral delta, plus the group means.

    Reported alongside the regression coefficient because d is scale-free and comparable
    across proteins, while the coefficient is not when a group approaches separation.
    """
    d = scored.loc[scored.kind == "disruptive", "delta"]
    n = scored.loc[scored.kind == "neutral", "delta"]
    pooled = np.sqrt((d.var(ddof=1) + n.var(ddof=1)) / 2) if len(d) > 1 and len(n) > 1 else np.nan
    return {
        "n_pairs": int(len(d)),
        "delta_disruptive_mean": round(float(d.mean()), 5),
        "delta_neutral_mean": round(float(n.mean()), 5),
        "delta_disruptive_median": round(float(d.median()), 5),
        "delta_neutral_median": round(float(n.median()), 5),
        "cohens_d": round(float((d.mean() - n.mean()) / pooled), 4) if pooled else np.nan,
    }


def motifs_from_config(cfg):
    """protein -> (motif, offset, replacement) from the `positive_control:` block."""
    raw = cfg["positive_control"]["motifs"]
    return {p: (v[0], int(v[1]), v[2]) for p, v in raw.items()}
