"""Stage 3e: the hard gate. Nothing trains until every dataset passes.

Exits non-zero on any violation. A warning that scrolls past is a warning that gets
ignored, and a silently corrupted dataset costs far more GPU time than a failed check.

    python scripts/validate.py                 # whole panel
    python scripts/validate.py --protein PUM2
"""

import argparse
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from rbp.data import annotation as ann  # noqa: E402
from rbp.data import splits  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Checks:
    def __init__(self):
        self.rows = []

    def add(self, name, ok, detail=""):
        self.rows.append((name, bool(ok), detail))
        return ok

    @property
    def failed(self):
        return [r for r in self.rows if not r[1]]


def digest(df):
    """Content hash, so a later silent change to preprocessing is detectable."""
    h = hashlib.sha256()
    for s in df.sort_values(["chrom", "start", "label"]).itertuples(index=False):
        h.update(f"{s.chrom}:{s.start}:{s.label}:{s.seq_rna}".encode())
    return h.hexdigest()[:16]


def validate_one(protein, cfg, path):
    c = Checks()
    df = pd.read_csv(path, sep="\t")
    size = cfg.windows["size"]
    pos = df[df.label == 1].reset_index(drop=True)
    negs = df[df.label == 0].reset_index(drop=True)

    c.add("non-empty", len(df) > 0, f"{len(df)} rows")
    c.add("balanced 1:1", len(pos) == len(negs), f"{len(pos)} pos / {len(negs)} neg")

    for sp in ("train", "val", "test"):
        s = df[df.split == sp]
        n1, n0 = int((s.label == 1).sum()), int((s.label == 0).sum())
        c.add(f"balanced in {sp}", n1 == n0, f"{n1}/{n0}")

    if len(pos) == len(negs):
        c.add("pair region match", (pos.region.values == negs.region.values).all())
        c.add("pair chromosome match", (pos.chrom.values == negs.chrom.values).all())
        c.add("pair split match", (pos.split.values == negs.split.values).all())
        gap = np.abs(pos.gc.values - negs.gc.values)
        tol = cfg.negatives["gc_tolerance"]
        c.add("GC gap median within tolerance", float(np.median(gap)) <= tol,
              f"median {np.median(gap):.4f} (tol {tol})")
        c.add("GC gap never extreme", float(gap.max()) <= tol * 4,
              f"max {gap.max():.4f}")

    c.add("no chromosome in two splits", splits.check_disjoint(
        df[["chrom", "split"]].to_dict("records")) == {})
    c.add("no duplicate windows", df.duplicated(["chrom", "start"]).sum() == 0,
          f"{int(df.duplicated(['chrom', 'start']).sum())} dupes")
    c.add("all sequences exact length", (df.seq_rna.str.len() == size).all())
    c.add("dna length matches", (df.seq_dna.str.len() == size).all())
    c.add("rna has no T", not df.seq_rna.str.contains("T").any())
    c.add("dna has no N", not df.seq_dna.str.contains("N").any())
    c.add("rna alphabet is ACGU", not df.seq_rna.str.contains("[^ACGU]").any())
    c.add("regions are known", set(df.region) <= set(ann.REGIONS),
          f"{sorted(set(df.region))}")

    excluded = set(cfg.encode.get("exclude_chroms", []))
    c.add("excluded chromosomes absent", not (set(df.chrom) & excluded))

    cfg_test, cfg_val = set(cfg.split["test"]), set(cfg.split["val"])
    c.add("split assignment matches config",
          all((r.chrom in cfg_test) == (r.split == "test")
              and (r.chrom in cfg_val) == (r.split == "val")
              for r in df.itertuples()))

    test_pairs = int(((df.split == "test") & (df.label == 1)).sum())
    c.add("enough test pairs", test_pairs >= cfg.panel["min_test_pairs"],
          f"{test_pairs} (need {cfg.panel['min_test_pairs']})")

    return c, digest(df), test_pairs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--protein", default=None)
    p.add_argument("--cell", default=None, help="default: encode.cell_line")
    p.add_argument("--arm", default=None, choices=sorted(panelmod.ARMS),
                   help="default: negatives.primary_arm")
    p.add_argument("--processed", default=None)
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--all-on-disk", action="store_true",
                   help="validate every dataset present, not just the panel")
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    cell = a.cell or cfg.encode["cell_line"]
    arm = panelmod.arm_of(cfg, a.arm)

    proc = Path(a.processed) if a.processed else panelmod.data_dir(cell, arm)
    panel = panelmod.read_panel(cell, arm)
    if a.protein:
        names = [a.protein]
    elif a.all_on_disk or not panel:
        names = sorted(d.name for d in proc.iterdir() if (d / "dataset.tsv").exists())
    else:
        names = [n for n, _ in panel]
        excl = panelmod.excluded_path(cell, arm)
        if excl.exists():
            n = len(excl.read_text().strip().splitlines()) - 1
            print(f"panel: {len(names)} proteins in {cell} {arm} ({n} excluded, "
                  f"see {excl.relative_to(ROOT)})")

    print(f"validating {len(names)} datasets against {cfg.windows['size']}-nt spec\n")
    print(f"{'protein':9} {'checks':>8} {'testpr':>7}  {'hash':16} status")
    failures = {}
    hashes = {}
    for name in names:
        c, h, tp = validate_one(name, cfg, proc / name / "dataset.tsv")
        hashes[name] = h
        bad = c.failed
        status = "PASS" if not bad else f"FAIL ({len(bad)})"
        print(f"{name:9} {len(c.rows):8d} {tp:7d}  {h:16} {status}")
        if bad:
            failures[name] = bad
        if a.verbose:
            for n, ok, d in c.rows:
                print(f"    {'ok ' if ok else 'FAIL'} {n}{'  ' + d if d else ''}")

    (ROOT / "results/tables").mkdir(parents=True, exist_ok=True)
    with open(ROOT / "results/tables/dataset_hashes.tsv", "w") as fh:
        fh.write("protein\tsha256_16\n")
        for k in sorted(hashes):
            fh.write(f"{k}\t{hashes[k]}\n")

    if failures:
        print("\nFAILURES")
        for prot, bad in failures.items():
            for n, _, d in bad:
                print(f"  {prot}: {n}  {d}")
        print(f"\nGATE CLOSED - {len(failures)} dataset(s) failed. Not safe to train.")
        sys.exit(1)

    print(f"\nGATE OPEN - all {len(names)} datasets pass. Hashes written to "
          f"results/tables/dataset_hashes.tsv")


if __name__ == "__main__":
    main()
