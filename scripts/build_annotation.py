"""Stage 3a: parse the GTF into region intervals, cache them, and report the mix.

    python scripts/build_annotation.py --config config/params.yaml
"""
import argparse, gzip, pickle, sys, time
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rbp.data import annotation as ann  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--gtf", default=str(ROOT / "data/raw/gencode.v45.annotation.gtf.gz"))
    p.add_argument("--out", default=str(ROOT / "data/interim/regions.pkl"))
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists() and not a.force:
        print(f"cached: {out.relative_to(ROOT)}")
        index = pickle.loads(out.read_bytes())
    else:
        t0 = time.time()
        print(f"parsing {Path(a.gtf).name} ...", flush=True)
        index = ann.build_index(a.gtf)
        out.write_bytes(pickle.dumps(index, protocol=5))
        print(f"  built in {time.time()-t0:.0f}s -> {out.relative_to(ROOT)} "
              f"({out.stat().st_size/1e6:.0f} MB)")

    print("\n=== region index ===")
    print(f"{'region':9} {'intervals':>10} {'bases':>14} {'chroms':>7}")
    for r, s in ann.stats(index).items():
        print(f"{r:9} {s['intervals']:10,} {s['bases']:14,} {s['chroms']:7d}")

    print("\n=== region of each protein's peak windows ===")
    size = cfg.windows["size"]
    half = size // 2
    rows = cfgmod.proteins()
    print(f"{'protein':9} {'utr3':>7} {'utr5':>6} {'cds':>7} {'exon_nc':>8} "
          f"{'intron':>7} {'none':>6}")
    totals = Counter()
    for r in rows:
        f = next((ROOT / "data/raw/peaks").glob(f"{r['protein']}.*.bed.gz"))
        c = Counter()
        with gzip.open(f, "rt") as fh:
            for line in fh:
                p_ = line.split("\t")
                mid = (int(p_[1]) + int(p_[2])) // 2
                c[ann.classify(index, p_[0], mid - half, mid + half + 1) or "none"] += 1
        n = sum(c.values()) or 1
        totals.update(c)
        print(f"{r['protein']:9} " + " ".join(
            f"{c[k]/n:6.1%}" for k in ("utr3", "utr5", "cds", "exon_nc", "intron", "none")))
    n = sum(totals.values())
    print(f"{'ALL':9} " + " ".join(
        f"{totals[k]/n:6.1%}" for k in ("utr3", "utr5", "cds", "exon_nc", "intron", "none")))


if __name__ == "__main__":
    main()
