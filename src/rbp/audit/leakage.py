"""Sequence-similarity leakage audit, with no external tools.

Chromosome-level splitting guarantees no shared locus. It does NOT prevent homologous
sequence -- repeats, paralogues, gene families -- sitting on different chromosomes and
landing in different splits. That residual leakage inflates test scores and looks exactly
like success.

The usual tool here is MMseqs2 or CD-HIT clustering. We use exact long-k-mer sharing
instead, which for fixed-length 101-nt windows is both sufficient and easier to interpret:

  * a 101-nt window contains 70 distinct 32-mers
  * two unrelated sequences sharing even ONE 32-mer has probability ~4^-32, i.e. never
  * so a shared 32-mer is proof of common origin, not coincidence

That makes "fraction of shared 32-mers" a direct, threshold-free measure of how much of a
test window also appears in training. Clustering tools answer the same question with more
machinery and a tunable identity cut-off that then has to be justified.

The k-mer length is the one parameter, and it trades sensitivity for specificity: shorter
k catches more diverged homology but starts admitting coincidence.
"""

from collections import defaultdict

import numpy as np


def kmer_set(seq, k):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def build_reference(seqs, k):
    """k-mer -> how many reference sequences contain it."""
    idx = defaultdict(int)
    for s in seqs:
        for km in kmer_set(s, k):
            idx[km] += 1
    return idx


def overlap_profile(query_seqs, reference_index, k):
    """For each query, the fraction of its k-mers that also occur in the reference."""
    out = np.zeros(len(query_seqs))
    for i, s in enumerate(query_seqs):
        kms = kmer_set(s, k)
        if not kms:
            continue
        out[i] = sum(1 for km in kms if km in reference_index) / len(kms)
    return out


def audit_protein(df, k=32, thresholds=(0.0, 0.25, 0.5, 0.9)):
    """How much of this protein's test set is echoed in its training set.

    Reported at several thresholds because the distribution matters: a window sharing 5%
    of its k-mers with training is a brief repeat, one sharing 90% is effectively a
    duplicate, and those deserve different treatment.
    """
    tr = df[df.split == "train"].seq_rna.tolist()
    te = df[df.split == "test"].seq_rna.tolist()
    if not tr or not te:
        return None

    ref = build_reference(tr, k)
    frac = overlap_profile(te, ref, k)

    out = {"n_train": len(tr), "n_test": len(te), "k": k,
           "mean_shared": round(float(frac.mean()), 5),
           "median_shared": round(float(np.median(frac)), 5),
           "p95_shared": round(float(np.percentile(frac, 95)), 5)}
    for t in thresholds:
        label = "any" if t == 0.0 else f"gt{int(t * 100)}"
        out[f"frac_test_{label}"] = round(float((frac > t).mean()), 5)
    return out, frac


def audit_panel(datasets, k=32, **kw):
    """Run the audit for every protein. Returns (table, per-protein raw profiles)."""
    import pandas as pd
    rows, profiles = [], {}
    for p, df in datasets.items():
        res = audit_protein(df, k=k, **kw)
        if res is None:
            continue
        stats, frac = res
        rows.append({"protein": p, **stats})
        profiles[p] = frac
    return (pd.DataFrame(rows).sort_values("frac_test_any", ascending=False)
              .reset_index(drop=True), profiles)


def strict_filter_effect(datasets, k=32, cutoff=0.5):
    """How many test windows a strict-similarity filter would remove.

    The reported headline should be the leakage that REMAINS under our split, plus what a
    stricter filter would cost. Quoting only the filtered number hides the exposure.
    """
    import pandas as pd
    rows = []
    for p, df in datasets.items():
        res = audit_protein(df, k=k)
        if res is None:
            continue
        _, frac = res
        n = len(frac)
        rows.append({"protein": p, "n_test": n,
                     "would_remove": int((frac > cutoff).sum()),
                     "frac_removed": round(float((frac > cutoff).mean()), 5),
                     "n_test_after": int((frac <= cutoff).sum())})
    return pd.DataFrame(rows).sort_values("frac_removed", ascending=False).reset_index(drop=True)
