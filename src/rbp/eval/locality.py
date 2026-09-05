"""Is a model's score driven by LOCAL sequence, or by global composition?

WHY THIS EXISTS. Rebuilding negatives with dinucleotide matching drops the median model
AUROC by 0.10. We want to read that as removing a composition artefact, but the
literature-motif positive control gets weaker on the matched arm (median Cohen's d 1.567 ->
0.807), so some real signal is destroyed too. To report the correction honestly we need to
know, PER DATASET, whether it was clean -- and the literature control covers 9 of 131
proteins.

WHAT FAILED FIRST. I tried a proxy: measure how repetitive a protein's enriched k-mers are,
on the theory that dinucleotide matching necessarily puts a repeat motif into the negatives.
It did not validate -- correlation with positive-control retention came out +0.298, the wrong
sign. And with 8 ground-truth proteins the standard error of a correlation is about 0.45, so
that exercise could not have validated anything either way.

THE MEASURE HERE. Skip the motif. Ask the question the composition story actually turns on:

    does the score depend on WHERE you change a base, or only on the fact that you changed one?

A composition-driven model barely moves when one position in 101 is altered, wherever it is.
A model reading a local feature collapses when you hit that feature and shrugs when you do
not. So: take the k-mer most enriched in bound windows, mutate its centre, and compare
against the SAME substitution placed far from any occurrence. The difference is locality.

WHY THIS IS NOT CIRCULAR. The k-mer is chosen from the training folds only, and the windows
are scored by the fold model that never saw them. More importantly, the comparison is
disruptive-versus-neutral WITHIN one model, so "the model recognises its own training signal"
cancels: both mutants are the same substitution under the same model, and only the position
differs. What survives is positional dependence, which is the quantity in question.

It measures locality, not correctness. A model could depend sharply on a local feature that
is not the real motif. That is a limitation to state, not to hide -- but it is the right
measure for asking whether a composition correction removed composition or removed signal.
"""

import numpy as np

from .repeats import kmer_enrichment, kmers_of


def top_kmer(pos_seqs, neg_seqs, k=5):
    """The single most enriched k-mer in bound windows. Model-free."""
    enr = kmer_enrichment(pos_seqs, neg_seqs, k=k)
    names = kmers_of(k)
    i = int(np.argmax(enr))
    return names[i], float(enr[i])


def _occurrences(seq, kmer):
    out, i = [], seq.find(kmer)
    while i != -1:
        out.append(i)
        i = seq.find(kmer, i + 1)
    return out


def build_pairs(seqs, kmer, min_distance=25, replacement=None, rng=None):
    """One disruptive and one matched neutral mutant per usable window.

    Both mutants make the SAME base substitution; only the position differs. That is what
    isolates locality from "mutating anything changes the score".

    A window is usable when it contains the k-mer, and some position at least
    `min_distance` from every occurrence carries the same reference base as the one being
    changed -- so the neutral mutant is genuinely the same edit somewhere uninformative.
    """
    rng = rng or np.random.default_rng(0)
    centre = len(kmer) // 2
    rows = []
    for idx, s in enumerate(seqs):
        hits = _occurrences(s, kmer)
        if not hits:
            continue
        di = hits[0] + centre
        ref = s[di]
        alt = replacement or next(c for c in "ACGU" if c != ref)
        far = [j for j in range(len(s))
               if s[j] == ref and all(abs(j - (h + centre)) >= min_distance for h in hits)]
        if not far:
            continue
        ni = int(rng.choice(far))
        rows.append({"i": idx,
                     "disruptive": s[:di] + alt + s[di + 1:],
                     "neutral": s[:ni] + alt + s[ni + 1:],
                     "ref": s})
    return rows


def locality(df, k=5, kmer_k=5, min_distance=25, max_windows=400, seed=7,
             seq_col="seq_rna", label_col="label", fold_col="fold"):
    """Cohen's d between disruptive and neutral score changes, pooled over folds.

    The k-mer is selected and the model is fit on training folds only; the windows scored
    come from the held-out fold. Returns None when no fold yields enough usable pairs.
    """
    from sklearn.linear_model import LogisticRegression

    from ..data.splits import fold_roles
    from .baseline import kmer_matrix

    y = df[label_col].to_numpy()
    folds = df[fold_col].to_numpy()
    seqs = df[seq_col].tolist()
    n_folds = len(np.unique(folds))
    rng = np.random.default_rng(seed)

    X, vec = kmer_matrix(seqs, k)
    dis_all, neu_all, chosen = [], [], []

    for f in range(n_folds):
        te_f, _, tr_f = fold_roles(f, n_folds)
        te = folds == te_f
        tr = np.isin(folds, tr_f) | (folds == (te_f + 1) % n_folds)
        if len(np.unique(y[tr])) < 2 or te.sum() < 20:
            continue

        tr_pos = [s for s, m, lab in zip(seqs, tr, y) if m and lab == 1]
        tr_neg = [s for s, m, lab in zip(seqs, tr, y) if m and lab == 0]
        if len(tr_pos) < 20 or len(tr_neg) < 20:
            continue
        km, _ = top_kmer(tr_pos, tr_neg, k=kmer_k)
        chosen.append(km)

        te_pos = [s for s, m, lab in zip(seqs, te, y) if m and lab == 1]
        pairs = build_pairs(te_pos, km, min_distance=min_distance, rng=rng)
        if len(pairs) < 10:
            continue
        if len(pairs) > max_windows:
            pairs = [pairs[i] for i in rng.choice(len(pairs), max_windows, replace=False)]

        model = LogisticRegression(max_iter=3000, C=0.01).fit(X[tr], y[tr])
        refs = [p["ref"] for p in pairs]
        s_ref = model.decision_function(vec.transform(refs))
        s_dis = model.decision_function(vec.transform([p["disruptive"] for p in pairs]))
        s_neu = model.decision_function(vec.transform([p["neutral"] for p in pairs]))
        dis_all.append(np.abs(s_ref - s_dis))
        neu_all.append(np.abs(s_ref - s_neu))

    if not dis_all:
        return None
    d = np.concatenate(dis_all)
    n = np.concatenate(neu_all)
    pooled = np.sqrt((d.var(ddof=1) + n.var(ddof=1)) / 2) if len(d) > 1 and len(n) > 1 else np.nan
    return {"n_pairs": int(len(d)), "kmers": chosen,
            "delta_disruptive": float(d.mean()), "delta_neutral": float(n.mean()),
            "cohens_d": float((d.mean() - n.mean()) / pooled) if pooled else np.nan}
