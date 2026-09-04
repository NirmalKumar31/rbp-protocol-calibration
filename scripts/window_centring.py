"""B16: the window centre is a free parameter. What does moving it do to the measurement?

    python scripts/window_centring.py --n 10
    python scripts/window_centring.py --from-cache

THE PARAMETER. Every window in this study is 101 nt centred on the MIDPOINT of its narrowPeak
interval (`windows.window_bounds`). Peaks are not points, their widths vary, and the
crosslinking event that produced the read pileup is not at the centre of the called interval.
B14 measured that the discriminative signal sits a median 23 nt off the window centre, so the
centring is demonstrably not aligned with the signal.

WHAT WE CANNOT DO, AND IT IS WORTH STATING RATHER THAN OMITTING. The obvious comparison is
summit-centred against midpoint-centred windows. narrowPeak reserves column 10 for the
point-source summit, and in these ENCODE eCLIP files it is -1 in every row of every file: no
summit is provided. `windows.read_peaks` never reads that column, and it is right not to. So a
summit-centred arm cannot be built from this data at all, and the sensitivity that CAN be tested
is to the choice among the centres the interval does provide.

THREE CENTRES, all defensible from the same BED. The published midpoint; the peak's 5' boundary
in transcript orientation, which is the mechanistically motivated alternative because the eCLIP
crosslink sits at the 5' end of the read pileup; and the midpoint displaced 25 nt downstream,
which is a pure arbitrariness probe with no rationale at all and is therefore the useful bound.

Every centring rebuilds its own positives from the genome and its own matched negatives in both
composition-matched arms, because a negative is matched to a positive and cannot be carried
across. Folds come from the chromosome, so the fold design is identical throughout.

THE MIDPOINT ARM IS NOT A REPRODUCTION CONTROL, which is what it was written as. Both matchers
SAMPLE their candidate windows from an RNG, so rebuilding the published centring gives a fresh
draw of the same construction rather than the same negatives. That turns it into a more useful
measurement than the one intended: how much each arm's contribution moves between draws. The
answer is 5.29e-05 for the GC arm and 7.90e-03 for the dinucleotide arm, and the arm
constraining fifteen degrees of freedom is the more variable by two orders of magnitude,
because matching sixteen frequencies well depends on which candidates were sampled while a
five-point GC tolerance is met by a large fraction of any pool.

4-mer only: both neural models trained on the published windows.
"""

import argparse
import pickle
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.data import annotation as ann  # noqa: E402
from rbp.data import negatives as neg  # noqa: E402
from rbp.data import windows as win  # noqa: E402
from rbp.data.encode import peak_path  # noqa: E402
from rbp.eval.baseline import oof_scores as kmer_oof  # noqa: E402
from rbp.eval.delong import delong_test  # noqa: E402
from rbp.eval.nested import _oof_scores, composition_features  # noqa: E402
from rbp.stats import standardise  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

TABLES = ROOT / "results" / "tables"
DATA_ROOT = ROOT.parent / "rna-binding-proteins"
CENTRES = ("midpoint", "five_prime", "shift25")
ARMS = ("gc", "dn")


def log(m):
    print(m, flush=True)


def centre_of(start, end, strand, mode):
    """Where the window sits, for one peak, under one centring rule."""
    if mode == "midpoint":
        return (start + end) // 2
    if mode == "five_prime":
        # Transcript orientation: the 5' boundary is `start` on the plus strand and `end` on
        # the minus strand. Getting this backwards would silently make the minus-strand half
        # of the panel a 3'-centred arm, which is a different experiment.
        return start if strand != "-" else end - 1
    if mode == "shift25":
        mid = (start + end) // 2
        return mid + 25 if strand != "-" else mid - 25
    raise ValueError(mode)


def positives_at(peak_file, fasta, index, size, mode, keep_chroms):
    """Positive windows under one centring rule, mirroring windows.build_positives."""
    half = size // 2
    rows, seen = [], set()
    dropped = {"chrom": 0, "bounds": 0, "n": 0, "region": 0, "dup": 0}
    for chrom, start, end, strand in win.read_peaks(peak_file):
        if chrom not in keep_chroms:
            dropped["chrom"] += 1
            continue
        c = centre_of(start, end, strand, mode)
        w0, w1 = c - half, c + half + 1
        dna = win.fetch(fasta, chrom, w0, w1)
        if dna is None or len(dna) != size:
            dropped["bounds"] += 1
            continue
        if "N" in dna.upper():
            dropped["n"] += 1
            continue
        region = ann.classify(index, chrom, w0, w1)
        if region is None:
            dropped["region"] += 1
            continue
        key = (chrom, w0, strand)
        if key in seen:
            dropped["dup"] += 1
            continue
        seen.add(key)
        rows.append({"chrom": chrom, "start": w0, "end": w1, "strand": strand,
                     "region": region, "gc": round(win.gc_content(dna), 4),
                     "seq_dna": dna, "seq_rna": win.to_rna(dna, strand)})
    return rows, dropped


def score(pos, negs, folds):
    keep = [i for i, m in enumerate(negs) if m is not None]
    if len(keep) < 200:
        return None
    seqs = [pos[i]["seq_rna"] for i in keep] + [negs[i]["seq_rna"] for i in keep]
    y = np.array([1] * len(keep) + [0] * len(keep))
    fo = np.concatenate([folds[keep], folds[keep]])
    X, _ = composition_features(seqs)
    s_comp = _oof_scores(X, y, fo)
    sc, _, _ = kmer_oof(seqs, y, fo, k=4)
    s_full = _oof_scores(np.column_stack([X, standardise(sc)]), y, fo)
    good = np.isfinite(s_comp) & np.isfinite(s_full)
    r = delong_test(s_full[good], s_comp[good], y[good])
    return {"pairs": len(keep), "comp": float(r["auc_b"]), "gain": float(r["diff"])}


def build(store, index_path, fasta_path, want):
    from pyfaidx import Fasta

    cfg = cfgmod.load()
    size = cfg.windows["size"]
    mpd = cfg.negatives["min_peak_distance"]
    tol = cfg.negatives["gc_tolerance"]
    drop = set(cfg.encode.get("exclude_chroms", []))
    keep_chroms = {c for c in ann.MAIN_CHROMS if c not in drop}
    log("loading region index and genome ...")
    index = pickle.loads(Path(index_path).read_bytes())
    fasta = Fasta(str(fasta_path))
    folds_map = pd.read_csv(ROOT / "config" / "folds.tsv", sep="\t").set_index("chrom").fold

    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv").dropna(subset=["gain_dn"])
    pub = pub.sort_values("n_dn")
    pick = np.linspace(0, len(pub) - 1, want).round().astype(int)
    datasets = list(pub.iloc[pick].dataset)

    rows = []
    for i, ds in enumerate(datasets, 1):
        protein, cell = ds.split(":")
        try:
            pf = peak_path(DATA_ROOT, protein, cell)
        except (FileNotFoundError, RuntimeError) as e:
            log(f"  {ds}: {e}")
            continue
        peaks = [p for p in win.read_peaks(pf) if p[0] in keep_chroms]
        for mode in CENTRES:
            pos, dr = positives_at(pf, fasta, index, size, mode, keep_chroms)
            if len(pos) < 200:
                continue
            folds = np.array([folds_map.get(p["chrom"], -1) for p in pos])
            # The pools depend on the peaks and the annotation, not on the centring, but the
            # positives do, so they are recomputed per centring and shared across the two arms.
            excl = neg.exclusion_zones(peaks, mpd)
            pools = {r: neg.available(index, r, excl, size)
                     for r in {p["region"] for p in pos}}
            for arm in ARMS:
                t0 = time.time()
                if arm == "dn":
                    negs, _, dist = neg.build_negatives_dinuc(
                        pos, peaks, fasta, index, size, min_peak_distance=mpd,
                        seed=cfg.seed, drop_n=True, pools=pools)
                    l1 = float(np.nanmean(dist))
                else:
                    negs, _ = neg.build_negatives(
                        pos, peaks, fasta, index, size, tolerance=tol,
                        min_peak_distance=mpd, seed=cfg.seed, drop_n=True)
                    l1 = float("nan")
                s = score(pos, negs, folds)
                if s is None:
                    continue
                rows.append({"dataset": ds, "protein": protein, "cell": cell,
                             "centre": mode, "arm": arm, "n_positives": len(pos),
                             "dropped_bounds": dr["bounds"], "dropped_n": dr["n"],
                             "dropped_region": dr["region"], "mean_l1": l1,
                             "seconds": round(time.time() - t0, 1), **s})
                log(f"[{i:3d}/{len(datasets)}] {ds:18s} {mode:11s} {arm:3s} "
                    f"pairs {s['pairs']:6d}  comp {s['comp']:.4f}  gain {s['gain']:+.4f}  "
                    f"{round(time.time() - t0, 1)}s")
    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("nothing built; refusing to overwrite the committed table")
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--store", default=str(ROOT.parent / "rbp-store"))
    p.add_argument("--index", default=str(DATA_ROOT / "data/interim/regions.pkl"))
    p.add_argument("--fasta",
                   default=str(DATA_ROOT / "data/raw/GRCh38.primary_assembly.genome.fa"))
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()
    warnings.filterwarnings("ignore")

    per = TABLES / "window_centring_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
        if t.empty:
            sys.exit(f"{per} is empty; regenerate it without --from-cache")
    else:
        t = build(a.store, a.index, a.fasta, a.n)
        t.to_csv(per, index=False)

    out = []
    ds = sorted(t.dataset.unique())
    log(f"\n=== B16: window-centring sensitivity, {len(ds)} datasets, {len(t)} rebuilds ===\n")
    out.append({"check": "datasets rebuilt under every centring", "value": len(ds),
                "n": len(t)})
    out.append({"check": "window sets rebuilt", "value": len(t), "n": len(t)})

    # THE FACT THAT SHAPES THE SECTION: no summit is available, so the literal comparison a
    # reviewer asks for cannot be run. Committed as a number so the claim is checkable.
    n_files, n_summits = 0, 0
    for d in ds:
        protein, cell = d.split(":")
        try:
            pf = peak_path(DATA_ROOT, protein, cell)
        except (FileNotFoundError, RuntimeError):
            continue
        n_files += 1
        import gzip
        with gzip.open(pf, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 10 and f[9].strip() not in ("-1", ""):
                    n_summits += 1
    out.append({"check": "narrowPeak files checked for a point-source summit",
                "value": n_files, "n": len(t)})
    out.append({"check": "peaks carrying a summit in column 10", "value": n_summits,
                "n": len(t), "note": "-1 everywhere, so no summit-centred arm is possible"})
    log(f"  {n_summits} of the peaks in {n_files} narrowPeak files carry a column-10 summit, "
        f"so a summit-centred arm cannot be built from this data")

    # THE CONTROL. The published centring rebuilt through this path must return the published
    # per-dataset gain, or the three centrings are being compared with each other and not with
    # the paper.
    # REPORTED PER ARM, and the result is the reverse of what the code comments imply. BOTH
    # matchers sample their candidate windows from an RNG, so neither rebuild is bit-identical
    # to the published one; the dinucleotide matcher is deterministic only in the ASSIGNMENT
    # step, conditional on the pool it happened to draw. So this is not a reproduction check,
    # it is a measurement of each arm's run-to-run variability under a fresh draw of the same
    # construction.
    #
    # And the arm that constrains FIFTEEN degrees of freedom is the more variable of the two,
    # by two orders of magnitude. That is the right way round on reflection: matching sixteen
    # frequencies well depends heavily on which candidates were sampled, while a GC tolerance
    # of five points is satisfied by a large fraction of any pool, so the GC arm's achieved
    # composition barely moves between draws.
    pub = pd.read_csv(TABLES / "three_arm_per_dataset.csv").set_index("dataset")
    for arm, col in (("gc", "gain_gc"), ("dn", "gain_dn")):
        s = t[(t.centre == "midpoint") & (t.arm == arm)].set_index("dataset")
        if s.empty:
            continue
        w = float((s.gain - pub.loc[s.index, col]).abs().max())
        kind = "assignment deterministic, pool sampled" if arm == "dn" \
            else "candidate window sampled"
        out.append({"check": f"max |fresh-draw gain - published gain|, {arm} arm",
                    "value": w, "n": len(s), "note": kind})
        log(f"  a fresh draw of the midpoint construction lands within {w:.2e} of the "
            f"published per-dataset gain, {arm} arm ({kind})")

    log(f"\n  {'centre':12s} {'arm':4s} {'positives':>10s} {'pairs':>8s} "
        f"{'composition':>12s} {'contribution':>13s}")
    cells = {}
    for mode in CENTRES:
        for arm in ARMS:
            s = t[(t.centre == mode) & (t.arm == arm)]
            if s.empty:
                continue
            cells[(mode, arm)] = s
            for name, col in (("positives", "n_positives"), ("pairs matched", "pairs"),
                              ("composition AUROC", "comp"), ("4-mer contribution", "gain")):
                out.append({"check": f"{name}, {mode} centring, {arm} arm",
                            "value": float(s[col].mean()), "n": len(s)})
            log(f"  {mode:12s} {arm:4s} {s.n_positives.mean():10.0f} {s.pairs.mean():8.0f} "
                f"{s.comp.mean():12.4f} {s.gain.mean():+13.4f}")

    # THE CONTRAST UNDER EACH CENTRING, which is the paper's quantity. Datasets are matched
    # across arms within a centring, so the contrast is a paired difference as published.
    log("")
    contrasts = {}
    for mode in CENTRES:
        g, dn = cells.get((mode, "gc")), cells.get((mode, "dn"))
        if g is None or dn is None:
            continue
        m = g.set_index("dataset").join(dn.set_index("dataset"), rsuffix="_dn", how="inner")
        c = float((m.gain_dn - m.gain).mean())
        contrasts[mode] = c
        out.append({"check": f"two-arm contrast, {mode} centring", "value": c, "n": len(m)})
        log(f"  two-arm contrast, {mode:11s} {c:+.4f}   on {len(m)} shared datasets")
    if len(contrasts) > 1:
        rng_c = max(contrasts.values()) - min(contrasts.values())
        out.append({"check": "range of the two-arm contrast over window centrings",
                    "value": rng_c, "n": len(t)})
        out.append({"check": "smallest two-arm contrast over window centrings",
                    "value": float(min(contrasts.values())), "n": len(t)})
        log(f"\n  over the {len(contrasts)} centrings the contrast ranges {rng_c:.4f}, "
            f"from {min(contrasts.values()):+.4f} to {max(contrasts.values()):+.4f}")

    # HOW MANY POSITIVES EACH CENTRING KEEPS. Moving the centre moves the window off the end
    # of a contig or into an unannotated stretch for some peaks, so the arms are not on
    # identical row sets and the size of that difference has to be reported rather than
    # assumed negligible.
    base = float(cells[("midpoint", "dn")].n_positives.mean())
    for mode in CENTRES:
        if (mode, "dn") not in cells:
            continue
        v = float(cells[(mode, "dn")].n_positives.mean())
        out.append({"check": f"positives retained relative to the midpoint centring, {mode}",
                    "value": v / base, "n": len(t)})
    log("  the centrings are NOT on identical rows: moving the centre changes which peaks "
        "\n  survive the bounds, N and annotation filters, and the retention above says by "
        "how much.")

    pd.DataFrame(out).to_csv(TABLES / "window_centring.csv", index=False)
    log("\nwrote window_centring.csv and window_centring_per_dataset.csv")


if __name__ == "__main__":
    main()
