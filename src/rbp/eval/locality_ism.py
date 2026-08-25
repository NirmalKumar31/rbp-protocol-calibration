"""Is a model's score driven by LOCAL sequence features, or by global composition?

This replaces the probe in locality.py, which does not work. Reading why is the fastest way
to understand why this one is shaped the way it is.

WHAT THE OLD PROBE DID. Pick the k-mer most enriched in bound windows. Mutate its centre
(disruptive) and compare against the same substitution placed far from any occurrence
(neutral). Report Cohen's d between the two score changes. A local model should collapse on
the disruptive mutation and shrug at the neutral one.

WHY IT FAILED, and it is two separate faults:

  1. A BAG-OF-K-MERS MODEL IS LOCAL BY CONSTRUCTION. It represents composition AS weights on
     k-mers, so mutating a high-weight k-mer always moves the score more than mutating
     elsewhere -- whether the underlying signal is a motif or pure global composition.

  2. THE DISRUPTIVE SITE IS SELECTED AS THE MOST ENRICHED K-MER, which is exactly where a
     fitted model concentrates its weight. The comparison is biased before any data is seen.

Fault 1 is fixed by using a model that COULD be non-local -- a CNN or a transformer. Fault 2
is not: choosing the site by enrichment biases the contrast for any model. tests/unit/
test_locality.py pins the consequence -- a pure global-composition signal, with no local
feature anywhere, scores Cohen's d about 1.8 where a valid probe must score ~0.

THE MEASURE HERE. Do not choose a site at all. Mutate EVERY position to EVERY alternative
base, take the mean |change in score| per position, and ask how CONCENTRATED that profile is.

    global / compositional signal -> changing any base changes composition by about as much
                                     as changing any other, so the profile is FLAT
    local / motif signal          -> only the positions inside the motif matter, so the
                                     profile is SPIKY

Concentration is measured with a Gini coefficient: 0 is perfectly flat, 1 is all sensitivity
at a single position. No k-mer is chosen, no site is picked, and the statistic is a property
of the whole profile rather than of two hand-picked points. That removes fault 2 by
construction rather than by argument.

WHY THE MEAN OVER ALTERNATIVE BASES MATTERS. Substituting to a fixed base would confound
position with base identity: under a GC-driven model, a G position loses GC when mutated to A
while an A position does not, so the profile would look spiky for a purely global signal.
Averaging over all three alternatives makes every position's perturbation comparable.

WHAT IT DOES NOT MEASURE. Correctness. A model can depend sharply on a local feature that is
not the real motif; this says the dependence is local, not that it is right. That is a
limitation to state, not to hide.

`score_fn` is any callable taking a list of sequences and returning one score each, so the
same probe runs against a k-mer logistic model, a CNN, or SpliceBERT, and the three are
directly comparable.
"""

import numpy as np

ALPHABET = "ACGU"


def gini(x):
    """Gini coefficient of a non-negative vector. 0 = perfectly flat, ->1 = concentrated.

    Chosen over entropy or a top-k fraction because it is scale-free -- a model whose scores
    move by 0.01 and one whose scores move by 10 give the same answer if the SHAPE of the
    profile matches. Locality is a question about shape.
    """
    x = np.asarray(x, dtype=float)
    x = np.abs(x)
    if x.sum() <= 0:
        return 0.0
    x = np.sort(x)
    n = len(x)
    idx = np.arange(1, n + 1)
    return float((2.0 * (idx * x).sum()) / (n * x.sum()) - (n + 1.0) / n)


def ism_profile(score_fn, seq, batch=256, alphabet=ALPHABET):
    """Mean |delta score| per position, averaged over all alternative bases.

    Returns an array of length len(seq). One reference score plus 3*L mutants, scored in
    batches so a neural score_fn is not called once per mutant.
    """
    L = len(seq)
    variants, owner = [], []
    for i, ref_base in enumerate(seq):
        for b in alphabet:
            if b == ref_base:
                continue
            variants.append(seq[:i] + b + seq[i + 1:])
            owner.append(i)

    scores = []
    for s in range(0, len(variants), batch):
        scores.append(np.asarray(score_fn(variants[s:s + batch]), dtype=float))
    scores = np.concatenate(scores) if scores else np.zeros(0)
    ref = float(np.asarray(score_fn([seq]), dtype=float)[0])

    prof = np.zeros(L)
    counts = np.zeros(L)
    for pos, sc in zip(owner, scores, strict=True):
        prof[pos] += abs(sc - ref)
        counts[pos] += 1
    return prof / np.maximum(counts, 1)


def locality(score_fn, seqs, max_windows=40, seed=7, batch=256):
    """Concentration of positional sensitivity, averaged over windows.

    Per-window Gini is averaged rather than pooling all positions across windows, because a
    motif sits at a DIFFERENT offset in each window -- pooling would smear a sharp per-window
    profile into a flat aggregate one and understate locality for exactly the models this is
    meant to detect.
    """
    rng = np.random.default_rng(seed)
    seqs = list(seqs)
    if len(seqs) > max_windows:
        seqs = [seqs[i] for i in rng.choice(len(seqs), max_windows, replace=False)]
    if not seqs:
        return None

    ginis, profiles = [], []
    for s in seqs:
        p = ism_profile(score_fn, s, batch=batch)
        if p.sum() <= 0:                     # a model that ignores the sequence entirely
            continue
        ginis.append(gini(p))
        profiles.append(p / p.sum())

    if not ginis:
        return None
    P = np.vstack(profiles)
    # Share of sensitivity in each window's most sensitive 10% of positions. Reported
    # alongside Gini because it is easier to state in a paper: "x% of the model's
    # sensitivity sits in 10% of positions".
    k = max(1, int(round(0.10 * P.shape[1])))
    top10 = float(np.mean([np.sort(r)[::-1][:k].sum() for r in P]))
    return {
        "gini": float(np.mean(ginis)),
        "gini_sd": float(np.std(ginis)),
        "top10_frac": top10,
        "n_windows": len(ginis),
        "length": int(P.shape[1]),
    }


def kmer_score_fn(model, vec):
    """Adapter: a fitted sklearn model over a k-mer vectoriser."""
    def f(seqs):
        return model.decision_function(vec.transform(list(seqs)))
    return f


def torch_score_fn(handle, device, batch=128):
    """Adapter: a trained rbp.models handle. Eval mode, no grad."""
    import torch

    handle.model.eval()

    def f(seqs):
        out = []
        seqs = list(seqs)
        with torch.no_grad():
            for s in range(0, len(seqs), batch):
                chunk = seqs[s:s + batch]
                logits = handle.forward(handle.batch(chunk, device))
                out.append(logits.float().cpu().numpy().ravel())
        return np.concatenate(out)
    return f
