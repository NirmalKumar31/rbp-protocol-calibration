"""Stage 5a. Define THE study panel once, write it to GCS, and never decide again.

WHY THIS FILE EXISTS. In the original build the 95-dataset panel was an emergent property
of a command-line flag (`--every 2`) typed during one sweep. Nothing recorded it, so the
project ended up with four different counts in circulation -- 189, 187, 95, 94 -- each
correct for a different question and none written down. That is the single largest source of
confusion in the whole study.

Here the panel is an ARTEFACT, not a flag. It is computed once, uploaded, and every
downstream stage reads it. If a stage wants to know what the study runs on, it asks this
file rather than re-deriving it.

    python scripts/select_panel.py --every 2      # write manifest/study_panel.tsv
    python scripts/select_panel.py --show         # print what is already there

WHY SYSTEMATIC SAMPLING AND NOT A THRESHOLD. AUROC correlates with dataset size at r = +0.53
to +0.67 across every model class. So keeping "the biggest N" would confound the panel with
the very quantity being measured: the subset would look better than the population for a
reason that has nothing to do with the science. Sorting by pair count and keeping every Nth
is unbiased in size BY CONSTRUCTION -- the sample spans the full range, including the
smallest and largest datasets.

WHY THE DINUCLEOTIDE ARM DEFINES IT. The two arms do not contain the same datasets: matching
is a search, it succeeds to different degrees, and a dataset can clear the min_pairs floor
in one arm and miss it in the other. The dinucleotide arm is the primary one (it is the
harder control and the one the headline is reported against), so it defines membership, and
the GC arm is intersected against it later. Choosing the intersection up front would silently
shrink the panel for a reason unrelated to the study.
"""

import argparse
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402

from rbp.utils import cloud as cloudcfg  # noqa: E402

PANEL_KEY = "manifest/study_panel.tsv"
CELLS = ("K562", "HepG2")


def log(m):
    print(m, flush=True)


def read_arm_panels(bucket, arm):
    """Every dataset the panel stage kept for one arm, both cell lines."""
    rows = []
    for cell in CELLS:
        blob = bucket.blob(f"panel/{arm}/panel_final_{cell}_{arm}.tsv")
        if not blob.exists():
            raise SystemExit(
                f"missing gs://{bucket.name}/panel/{arm}/panel_final_{cell}_{arm}.tsv -- "
                f"stage 4 (build panel) has not run for arm={arm}")
        d = pd.read_csv(io.StringIO(blob.download_as_text()), sep="\t")
        rows.append(d)
    return pd.concat(rows, ignore_index=True)


def select(bucket, every, primary="dinuc"):
    full = read_arm_panels(bucket, primary)
    full["dataset"] = full.protein + ":" + full.cell_line
    # DETERMINISTIC SORT, AND THE TIEBREAKER IS NOT DECORATION.
    #
    # pandas sort_values defaults to quicksort, which is NOT stable, and three pairs of
    # datasets in this panel share a pair count exactly (539, 3640, 7988). Sorting on `pairs`
    # alone therefore leaves it arbitrary which member of a tied pair lands on an even index,
    # and `[::every]` then keeps a different one run to run. Measured: two of the three tie
    # pairs flipped between the original study and its reproduction, changing panel membership
    # by two datasets out of ninety-five.
    #
    # The science was unaffected -- tied datasets have identical size by definition, so the
    # panel's size distribution does not move -- but a panel that is not reproducible defeats
    # the entire purpose of writing it down as an artefact. Sorting on (pairs, dataset) with a
    # stable algorithm makes the selection a function of the data and nothing else.
    full = (full.drop_duplicates("dataset")
                .sort_values(["pairs", "dataset"], kind="mergesort")
                .reset_index(drop=True))

    picked = full.iloc[::every].reset_index(drop=True) if every > 1 else full
    picked = picked[["dataset", "protein", "cell_line", "pairs"]].copy()
    picked.insert(0, "idx", range(len(picked)))

    log(f"\n{primary} arm has {len(full)} datasets, {full.protein.nunique()} proteins")
    log(f"every={every} -> {len(picked)} datasets, {picked.protein.nunique()} proteins")
    log(f"pairs range kept {picked.pairs.min():,}-{picked.pairs.max():,} "
        f"(full panel {full.pairs.min():,}-{full.pairs.max():,})")

    # The check that the sample is not size-biased. If the kept set does not reach into both
    # tails of the full distribution, systematic sampling has been broken somewhere.
    lo_ok = picked.pairs.min() <= full.pairs.quantile(0.05)
    hi_ok = picked.pairs.max() >= full.pairs.quantile(0.95)
    if not (lo_ok and hi_ok):
        raise SystemExit(
            f"panel is size-biased: kept range does not span the full distribution "
            f"(reaches 5th pct: {lo_ok}, 95th pct: {hi_ok}). Refusing to write it.")
    log(f"size-unbiased: spans the 5th and 95th percentile of the full panel  OK")

    # How much of the panel also exists in the other arm, reported now so the R1 count is
    # never a surprise later.
    other = "gc" if primary == "dinuc" else "dinuc"
    try:
        oth = read_arm_panels(bucket, other)
        oth["dataset"] = oth.protein + ":" + oth.cell_line
        both = set(picked.dataset) & set(oth.dataset)
        log(f"of these, {len(both)} also clear the floor in the {other} arm "
            f"-> that is the n for the cost-of-matching result")
        picked["in_both_arms"] = picked.dataset.isin(both)
    except SystemExit:
        log(f"({other} arm not built yet; in_both_arms left blank)")
        picked["in_both_arms"] = pd.NA
    return picked


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--every", type=int, default=2,
                   help="keep every Nth dataset by pair rank. 1 = the whole panel")
    p.add_argument("--primary", default="dinuc", choices=["dinuc", "gc"])
    p.add_argument("--show", action="store_true", help="print the existing panel and exit")
    p.add_argument("--force", action="store_true", help="overwrite an existing panel")
    a = p.parse_args()

    bucket = cloudcfg.bucket()
    blob = bucket.blob(PANEL_KEY)

    if a.show:
        if not blob.exists():
            raise SystemExit(f"no panel at gs://{bucket.name}/{PANEL_KEY}")
        d = pd.read_csv(io.StringIO(blob.download_as_text()), sep="\t")
        log(f"{len(d)} datasets, {d.protein.nunique()} proteins, "
            f"pairs {d.pairs.min():,}-{d.pairs.max():,}")
        log(d.head(10).to_string(index=False))
        return

    # THE PANEL IS WRITTEN ONCE. Silently redefining it mid-study is how half the results
    # end up describing a different set of datasets than the other half.
    if blob.exists() and not a.force:
        d = pd.read_csv(io.StringIO(blob.download_as_text()), sep="\t")
        log(f"panel already exists: {len(d)} datasets. Refusing to redefine it.")
        log("Pass --force only if you intend every downstream result to be invalidated.")
        return

    picked = select(bucket, a.every, a.primary)
    blob.upload_from_string(picked.to_csv(sep="\t", index=False),
                            content_type="text/tab-separated-values")
    log(f"\nwrote gs://{bucket.name}/{PANEL_KEY}  ({len(picked)} datasets)")
    log("Every downstream stage reads this file. Nothing else decides the panel.")


if __name__ == "__main__":
    main()
