"""Dataset summaries and figures.

Computation is kept separate from plotting so the numbers can be tested and reused in
the paper's tables without importing matplotlib.
"""

from pathlib import Path

import numpy as np
import pandas as pd

REGION_ORDER = ("utr5", "utr3", "cds", "exon_nc", "intron")
SPLIT_ORDER = ("train", "val", "test")


def load(processed_dir, proteins):
    """protein -> dataframe, from data/processed/<PROTEIN>/dataset.tsv."""
    out = {}
    for p in proteins:
        f = Path(processed_dir) / p / "dataset.tsv"
        if f.exists():
            out[p] = pd.read_csv(f, sep="\t")
    return out


def panel_summary(datasets):
    """One row per protein: sizes, split shares, GC match quality, region mix."""
    rows = []
    for p, df in datasets.items():
        pos = df[df.label == 1].reset_index(drop=True)
        negs = df[df.label == 0].reset_index(drop=True)
        gap = np.abs(pos.gc.values - negs.gc.values)
        r = {"protein": p, "pairs": len(pos), "rows": len(df)}
        for s in SPLIT_ORDER:
            r[f"{s}_pairs"] = int((pos.split == s).sum())
            r[f"{s}_frac"] = round(float((df.split == s).mean()), 3)
        r["gc_pos_mean"] = round(float(pos.gc.mean()), 4)
        r["gc_neg_mean"] = round(float(negs.gc.mean()), 4)
        r["gc_gap_mean"] = round(float(gap.mean()), 4)
        r["gc_gap_p95"] = round(float(np.percentile(gap, 95)), 4)
        r["gc_within_tol"] = round(float((gap <= 0.05).mean()), 3)
        vc = pos.region.value_counts(normalize=True)
        for reg in REGION_ORDER:
            r[f"frac_{reg}"] = round(float(vc.get(reg, 0.0)), 3)
        r["dominant_region"] = pos.region.value_counts().index[0]
        rows.append(r)
    return pd.DataFrame(rows).sort_values("protein").reset_index(drop=True)


def integrity(datasets):
    """Per-protein invariants, so the EDA re-states what the gate checked."""
    rows = []
    for p, df in datasets.items():
        pos = df[df.label == 1].reset_index(drop=True)
        negs = df[df.label == 0].reset_index(drop=True)
        rows.append({
            "protein": p,
            "balanced": len(pos) == len(negs),
            "region_match": round(float((pos.region.values == negs.region.values).mean()), 4),
            "chrom_match": round(float((pos.chrom.values == negs.chrom.values).mean()), 4),
            "split_match": round(float((pos.split.values == negs.split.values).mean()), 4),
            "all_101nt": bool((df.seq_rna.str.len() == 101).all()),
            "no_t_in_rna": bool(~df.seq_rna.str.contains("T").any()),
            "no_n": bool(~df.seq_dna.str.contains("N").any()),
            "dup_windows": int(df.duplicated(["chrom", "start"]).sum()),
        })
    return pd.DataFrame(rows).sort_values("protein").reset_index(drop=True)


def nucleotide_composition(datasets):
    """Base frequencies for positives and negatives, a shortcut-detection check."""
    rows = []
    for p, df in datasets.items():
        for label, name in ((1, "positive"), (0, "negative")):
            seqs = df[df.label == label].seq_rna
            joined = "".join(seqs.tolist())
            n = len(joined) or 1
            rows.append({"protein": p, "class": name,
                         **{b: round(joined.count(b) / n, 4) for b in "ACGU"}})
    return pd.DataFrame(rows)


def motif_enrichment(datasets, motifs):
    """Known-motif frequency in positives vs matched negatives.

    A motif that is enriched in positives but not in their GC- and region-matched
    negatives is evidence the task is about binding rather than composition.
    """
    rows = []
    for p, motif in motifs.items():
        if p not in datasets:
            continue
        df = datasets[p]
        pos = df[df.label == 1].seq_rna
        neg = df[df.label == 0].seq_rna
        fp = float(pos.str.contains(motif).mean())
        fn = float(neg.str.contains(motif).mean())
        rows.append({"protein": p, "motif": motif,
                     "positives": round(fp, 4), "negatives": round(fn, 4),
                     "enrichment": round(fp / fn, 2) if fn else np.inf})
    return pd.DataFrame(rows).sort_values("enrichment", ascending=False)


# --------------------------------------------------------------------------------------
# Analyses aimed at model design rather than description.
# --------------------------------------------------------------------------------------

def kmer_features(seqs, k):
    """Count matrix over all 4^k k-mers. Sparse, so k=6 stays cheap."""
    from sklearn.feature_extraction.text import CountVectorizer
    vec = CountVectorizer(analyzer="char", ngram_range=(k, k), lowercase=False)
    X = vec.fit_transform(seqs)
    return X, np.array(vec.get_feature_names_out())


def kmer_baseline(df, k=5, seed=7, top=12):
    """Logistic regression on k-mer counts, trained on train and scored on test.

    This is the honest floor for the whole project. If a bag of k-mers already reaches
    0.90 on a protein, a transformer beating it by 0.02 is a small gain and the paper
    should say so. It also reveals which k-mers carry the signal, which is a free
    interpretability check against known motifs.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    tr, te = df[df.split == "train"], df[df.split == "test"]
    if len(te) < 40 or te.label.nunique() < 2:
        return None
    X, names = kmer_features(pd.concat([tr, te]).seq_rna.tolist(), k)
    ntr = len(tr)
    model = LogisticRegression(max_iter=3000, C=1.0).fit(X[:ntr], tr.label.values)
    p = model.predict_proba(X[ntr:])[:, 1]
    auroc = float(roc_auc_score(te.label.values, p))
    coef = model.coef_[0]
    order = np.argsort(coef)
    return {
        "k": k, "auroc": round(auroc, 4),
        "n_train": ntr, "n_test": len(te),
        "top_positive": list(names[order[::-1][:top]]),
        "top_negative": list(names[order[:top]]),
    }


def kmer_sweep(datasets, ks=(3, 4, 5, 6)):
    """Baseline AUROC per protein for several k, to see what motif width matters."""
    rows = []
    for p, df in datasets.items():
        r = {"protein": p}
        for k in ks:
            out = kmer_baseline(df, k=k)
            r[f"k{k}"] = out["auroc"] if out else np.nan
        r["best_k"] = int(max(ks, key=lambda k: (r[f"k{k}"] if r[f"k{k}"] == r[f"k{k}"] else -1)))
        r["best_auroc"] = r[f"k{r['best_k']}"]
        rows.append(r)
    return pd.DataFrame(rows).sort_values("best_auroc", ascending=False).reset_index(drop=True)


def positional_profile(df, motif, bins=21):
    """Where in the window the motif sits, positives vs negatives.

    If the signal is sharply centred, a model can lean on position. If it is flat,
    the architecture needs position invariance, which is what global pooling provides.
    """
    size = len(df.seq_rna.iloc[0])
    edges = np.linspace(0, size, bins + 1)
    out = {}
    for label, name in ((1, "positive"), (0, "negative")):
        starts = []
        for s in df[df.label == label].seq_rna:
            i = s.find(motif)
            while i != -1:
                starts.append(i)
                i = s.find(motif, i + 1)
        hist, _ = np.histogram(starts, bins=edges)
        out[name] = hist / (hist.sum() or 1)
    out["centres"] = (edges[:-1] + edges[1:]) / 2
    return out


def peak_width_stats(peak_paths, window):
    """Peak widths against the fixed window size.

    Windows much narrower than a peak truncate it; much wider and the model is mostly
    reading flanking context. Either way it is a modelling assumption worth stating.
    """
    from ..data.windows import read_peaks
    rows = []
    for p, path in peak_paths.items():
        w = np.array([e - s for _, s, e, _ in read_peaks(path)])
        if not len(w):
            continue
        rows.append({"protein": p, "n": len(w), "median": int(np.median(w)),
                     "p10": int(np.percentile(w, 10)), "p90": int(np.percentile(w, 90)),
                     "max": int(w.max()),
                     "frac_wider_than_window": round(float((w > window).mean()), 3)})
    return pd.DataFrame(rows).sort_values("median").reset_index(drop=True)


def split_shift(datasets):
    """Is the test split distributionally like train? A shift changes what AUROC means."""
    rows = []
    for p, df in datasets.items():
        tr, te = df[df.split == "train"], df[df.split == "test"]
        r = {"protein": p,
             "gc_train": round(float(tr.gc.mean()), 4),
             "gc_test": round(float(te.gc.mean()), 4),
             "gc_shift": round(float(te.gc.mean() - tr.gc.mean()), 4)}
        a = tr[tr.label == 1].region.value_counts(normalize=True)
        b = te[te.label == 1].region.value_counts(normalize=True)
        regions = set(a.index) | set(b.index)
        r["region_l1_shift"] = round(
            float(sum(abs(a.get(x, 0) - b.get(x, 0)) for x in regions) / 2), 3)
        rows.append(r)
    return pd.DataFrame(rows).sort_values("region_l1_shift", ascending=False).reset_index(drop=True)


def redundancy(datasets):
    """Exact duplicate sequences within a protein, and how much they inflate the data.

    Duplicates make the effective sample size smaller than the row count, which matters
    when reading per-protein confidence intervals.
    """
    rows = []
    for p, df in datasets.items():
        n = len(df)
        uniq = df.seq_rna.nunique()
        pos = df[df.label == 1]
        rows.append({"protein": p, "rows": n, "unique_seqs": uniq,
                     "dup_frac": round(1 - uniq / n, 4),
                     "pos_dup_frac": round(1 - pos.seq_rna.nunique() / len(pos), 4)})
    return pd.DataFrame(rows).sort_values("dup_frac", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------------------
# General descriptive EDA.
# --------------------------------------------------------------------------------------

BASES = ("A", "C", "G", "U")


def _seq_matrix(seqs):
    """Sequences as a 2-D uint8 array of characters, for vectorised per-position work."""
    arr = np.frombuffer("".join(seqs).encode(), dtype=np.uint8)
    return arr.reshape(len(seqs), -1)


def positional_base_profile(df, label=1):
    """Base frequency at each position in the window. Shape (size, 4)."""
    seqs = df[df.label == label].seq_rna.tolist()
    if not seqs:
        return None
    m = _seq_matrix(seqs)
    return np.stack([(m == ord(b)).mean(axis=0) for b in BASES], axis=1)


def positional_signal(datasets):
    """Per-position |positive - negative| base-composition gap, summed over bases.

    A sharp central peak means the discriminative signal sits where the peak was
    centred, so position carries information. A flat profile means the model needs to
    be position invariant, which is what global pooling gives it.
    """
    rows = {}
    for p, df in datasets.items():
        a = positional_base_profile(df, 1)
        b = positional_base_profile(df, 0)
        if a is None or b is None:
            continue
        rows[p] = np.abs(a - b).sum(axis=1)
    return rows


def dinucleotide_enrichment(datasets):
    """log2(positive / negative) frequency for all 16 dinucleotides.

    GC matching constrains G+C but not how those bases are arranged, so this is where
    residual compositional shortcuts show up.
    """
    dins = [x + y for x in BASES for y in BASES]
    rows = []
    for p, df in datasets.items():
        out = {"protein": p}
        for label, name in ((1, "pos"), (0, "neg")):
            joined = "".join(df[df.label == label].seq_rna.tolist())
            tot = max(len(joined) - 1, 1)
            for d in dins:
                c = joined.count(d)
                out[f"{name}_{d}"] = c / tot
        for d in dins:
            a, b = out[f"pos_{d}"], out[f"neg_{d}"]
            out[f"lr_{d}"] = float(np.log2(a / b)) if a > 0 and b > 0 else np.nan
        rows.append(out)
    df = pd.DataFrame(rows)
    return df[["protein"] + [f"lr_{d}" for d in dins]]


def _max_run(s):
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best


def complexity(datasets, run_len=5):
    """Low-complexity content: homopolymer runs and per-window base entropy.

    Repeats are a classic genomics confound. If positives are systematically more
    repetitive than their matched negatives, a model can win on that alone.
    """
    rows = []
    for p, df in datasets.items():
        for label, name in ((1, "positive"), (0, "negative")):
            seqs = df[df.label == label].seq_rna.tolist()
            runs = np.array([_max_run(s) for s in seqs])
            m = _seq_matrix(seqs)
            freq = np.stack([(m == ord(b)).mean(axis=1) for b in BASES], axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                ent = -(freq * np.log2(np.where(freq > 0, freq, 1))).sum(axis=1)
            rows.append({"protein": p, "class": name,
                         "max_run_median": int(np.median(runs)),
                         "max_run_p95": int(np.percentile(runs, 95)),
                         "frac_with_run_ge": round(float((runs >= run_len).mean()), 3),
                         "entropy_mean": round(float(ent.mean()), 4)})
    return pd.DataFrame(rows)


def read_peaks_full(path):
    """narrowPeak with the quality columns we otherwise ignore."""
    import gzip
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            rows.append({"chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                         "strand": f[5], "signal": float(f[6]), "pval": float(f[7])})
    return pd.DataFrame(rows)


def peak_quality(peak_paths):
    """Fold-enrichment and significance per protein.

    ENCODE peaks are already IDR-filtered, but the spread tells us whether a
    confidence-stratified analysis is worth doing later.
    """
    rows = []
    for p, path in peak_paths.items():
        d = read_peaks_full(path)
        if d.empty:
            continue
        rows.append({"protein": p, "n": len(d),
                     "signal_median": round(float(d.signal.median()), 2),
                     "signal_p90": round(float(d.signal.quantile(0.9)), 2),
                     "pval_median": round(float(d.pval.median()), 2),
                     "pval_p90": round(float(d.pval.quantile(0.9)), 2),
                     "frac_plus": round(float((d.strand == "+").mean()), 3)})
    return pd.DataFrame(rows).sort_values("signal_median", ascending=False).reset_index(drop=True)


def chromosome_distribution(peak_paths, fai_path):
    """Peaks per chromosome, raw and per megabase, to spot genomic hot spots."""
    sizes = {}
    for line in open(fai_path):
        f = line.split("\t")
        sizes[f[0]] = int(f[1])
    rows = []
    for p, path in peak_paths.items():
        d = read_peaks_full(path)
        vc = d.chrom.value_counts()
        for c, n in vc.items():
            if c in sizes:
                rows.append({"protein": p, "chrom": c, "peaks": int(n),
                             "per_mb": round(n / (sizes[c] / 1e6), 3)})
    return pd.DataFrame(rows)


def gc_by_region(datasets):
    """GC content per region class, positives only.

    Regions differ a lot in GC, which is why negatives are matched within region
    rather than globally.
    """
    rows = []
    for p, df in datasets.items():
        pos = df[df.label == 1]
        for reg, g in pos.groupby("region"):
            rows.append({"protein": p, "region": reg, "n": len(g),
                         "gc_mean": round(float(g.gc.mean()), 4),
                         "gc_sd": round(float(g.gc.std()), 4)})
    return pd.DataFrame(rows)


def protein_clusters(cobinding_matrix, k=4):
    """Group proteins by their co-binding profile, so families show up explicitly."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    m = cobinding_matrix.to_numpy(dtype=float).copy()
    m = (m + m.T) / 2
    np.fill_diagonal(m, 1.0)
    dist = 1 - m
    np.fill_diagonal(dist, 0.0)
    dist = np.clip(dist, 0, None)
    Z = linkage(squareform(dist, checks=False), method="average")
    labels = fcluster(Z, k, criterion="maxclust")
    return pd.DataFrame({"protein": cobinding_matrix.index, "cluster": labels}
                        ).sort_values(["cluster", "protein"]).reset_index(drop=True), Z


def _gini(counts):
    """Gini coefficient of a count vector: 0 = perfectly even, 1 = all in one bin."""
    x = np.sort(np.asarray(counts, dtype=float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return np.nan
    return float((2.0 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def gene_concentration(datasets, gene_index, gene_at):
    """How concentrated each protein's binding is across genes.

    This bounds the EFFECTIVE sample size. Windows in the same transcript share
    sequence, expression and local context, so they are not independent observations.
    If 8,800 windows sit in 300 genes, a confidence interval computed as though n=8,800
    is optimistic. Reporting the concentration makes that visible rather than assumed.
    """
    rows = []
    for p, df in datasets.items():
        pos = df[df.label == 1]
        genes = [gene_at(gene_index, c, (s + e) // 2)
                 for c, s, e in zip(pos.chrom, pos.start, pos.end)]
        named = [g for g in genes if g]
        vc = pd.Series(named).value_counts()
        n = len(pos)
        rows.append({
            "protein": p,
            "n_windows": n,
            "n_genes": int(vc.size),
            "frac_unassigned": round(1 - len(named) / n, 4) if n else np.nan,
            "windows_per_gene_median": float(vc.median()) if vc.size else np.nan,
            "windows_per_gene_max": int(vc.iloc[0]) if vc.size else 0,
            "busiest_gene": vc.index[0] if vc.size else None,
            "top10_gene_share": round(float(vc.head(10).sum() / len(named)), 4) if named else np.nan,
            "gini": round(_gini(vc.values), 4) if vc.size else np.nan,
            # windows are between fully independent (n_windows) and fully dependent
            # (n_genes); both bounds are reported rather than picking one
            "eff_n_lower": int(vc.size),
            "eff_n_upper": n,
        })
    return pd.DataFrame(rows).sort_values("gini", ascending=False).reset_index(drop=True)
