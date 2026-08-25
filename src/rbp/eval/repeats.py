"""Is a protein's enriched sequence a REPEAT? Because that decides whether our
composition correction is clean or over-corrects.

THE PROBLEM THIS EXISTS TO MEASURE. Rebuilding negatives with dinucleotide matching drops
the median model AUROC by 0.10, and we read that as removing a composition artefact. But
the positive control -- where we plant a known motif disruption ourselves -- gets WEAKER on
the matched arm: median Cohen's d 1.567 -> 0.807, with 2 of 9 proteins flipping from passing
to failing. So some real signal is being destroyed too, and the 0.10 is an upper bound on the
confound rather than a clean estimate.

The pattern in that failure is specific. Retention of the positive-control effect:

    RBFOX2 (GCAUG)  0.88     non-repeat motifs survive
    PUM2   (UGUA)   0.78
    QKI    (ACUAA)  0.74
    PTBP1  (UCUU)   0.54     CU-rich, partially repetitive
    TARDBP (UGUGU)  0.20     a dinucleotide repeat -- collapses

TARDBP binds (UG)n. If you match negatives on dinucleotide composition, negatives with
TARDBP's dinucleotide profile necessarily CONTAIN UG-rich sequence -- that is, they contain
the motif. So the correction removes the signal along with the confound. This is the same
flaw that killed the dinucleotide-SHUFFLE control (docs 21 Part 2): matching inherits the
shuffle's repeat-motif problem while fixing its detectable-artefact problem.

WHY A PROXY IS NEEDED. The positive control requires a literature-derived motif, which
exists for 9 of our 131 proteins. Inventing motifs for the rest would make the control
circular. So we need to measure repeat-ness from the DATA, with no motif and no model, and
then check the measure against the 9 proteins where ground truth exists.

WHAT THIS MEASURES. For each dataset, find the k-mers most enriched in bound windows versus
their matched negatives -- a model-free log-ratio, so no coefficients and no fitting -- and
ask how periodic those k-mers are. A k-mer is periodic if it repeats a shorter unit: UGUGU
is (UG) repeated, UUUUU is (U) repeated, GCAUG repeats nothing.

The output is one number per dataset, and it is a stratification variable, not a result: it
says whether to trust the composition correction for that protein.
"""

import numpy as np

ALPHABET = "ACGU"


def minimal_period(s):
    """Smallest p such that s[i] == s[i+p] for every valid i.

    p == len(s) means no repetition. Computed directly rather than with a string-matching
    algorithm because k is 4-6 here and clarity is worth more than the constant factor.

        UUUUU -> 1   (a homopolymer repeats a single base)
        UGUGU -> 2   (repeats UG)
        GCAUG -> 5   (repeats nothing)
    """
    n = len(s)
    for p in range(1, n):
        if all(s[i] == s[i + p] for i in range(n - p)):
            return p
    return n


def repeat_score(s):
    """0 for a k-mer that repeats nothing, 1 for a homopolymer.

    Scaled so the measure is comparable across k:  1 - (period - 1) / (k - 1).

        UUUUU -> 1.00
        UGUGU -> 0.75
        AUUUA -> 0.00   (period 5; near-homopolymer by eye but not periodic)
        GCAUG -> 0.00
    """
    n = len(s)
    if n < 2:
        return 0.0
    return 1.0 - (minimal_period(s) - 1) / (n - 1)


def low_complexity(s):
    """Fraction of positions inside a run of 3 or more identical bases.

    `repeat_score` only fires on exact global periodicity, so AUUUA scores 0 despite being
    obviously low-complexity. This catches that separately, and the two are combined below.
    """
    n = len(s)
    if n < 3:
        return 0.0
    covered = np.zeros(n, dtype=bool)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and s[j + 1] == s[i]:
            j += 1
        if j - i + 1 >= 3:
            covered[i:j + 1] = True
        i = j + 1
    return float(covered.mean())


def kmer_enrichment(pos_seqs, neg_seqs, k=5, pseudo=1.0):
    """log2 enrichment of every k-mer in positives versus negatives.

    Model-free on purpose. An earlier analysis in this project found that top k-mers ranked
    by L2 logistic coefficients include artefacts -- CG-rich k-mers surfacing for CU-binding
    proteins -- because L2 on sparse counts spreads weight oddly. A frequency ratio has no
    such failure mode.
    """
    from .nested import _counts
    kp = _counts(list(pos_seqs), k).sum(axis=0).astype(float)
    kn = _counts(list(neg_seqs), k).sum(axis=0).astype(float)
    fp = (kp + pseudo) / (kp.sum() + pseudo * len(kp))
    fn = (kn + pseudo) / (kn.sum() + pseudo * len(kn))
    return np.log2(fp / fn)


def kmers_of(k):
    from itertools import product
    return ["".join(t) for t in product(ALPHABET, repeat=k)]


def repetitiveness(pos_seqs, neg_seqs, k=5, top=20):
    """How repetitive is the sequence this protein prefers? One number per dataset.

    Enrichment-weighted mean of the per-k-mer repeat measure over the `top` most enriched
    k-mers. Weighting by enrichment rather than taking the single best k-mer, because a
    protein's preference is usually a family of related k-mers and the top one alone is
    noisy.
    """
    enr = kmer_enrichment(pos_seqs, neg_seqs, k=k)
    names = kmers_of(k)
    order = np.argsort(enr)[::-1][:top]
    w = np.clip(enr[order], 0.0, None)
    if w.sum() <= 0:
        return {"repetitiveness": np.nan, "top_kmers": [], "max_enrichment": float(enr.max())}
    sel = [names[i] for i in order]
    rep = np.array([repeat_score(s) for s in sel])
    lc = np.array([low_complexity(s) for s in sel])
    combined = np.maximum(rep, lc)
    return {"repetitiveness": float(np.average(combined, weights=w)),
            "repeat_only": float(np.average(rep, weights=w)),
            "low_complexity_only": float(np.average(lc, weights=w)),
            "top_kmers": sel[:5],
            "max_enrichment": float(enr[order[0]])}


def from_dataset(df, seq_col="seq_rna", label_col="label", **kw):
    """Repetitiveness from a prepared dataset.

    Uses `seq_rna`, the strand-corrected sequence, because that is what the protein actually
    sees and what a motif is defined in. `seq_dna` would report the complement's periodicity
    for minus-strand windows -- which has the same period, but the k-mer identities would be
    wrong and unreadable against the literature.
    """
    pos = df.loc[df[label_col] == 1, seq_col].tolist()
    neg = df.loc[df[label_col] == 0, seq_col].tolist()
    return repetitiveness(pos, neg, **kw)
