"""Build the composition-control arm: negatives are shuffled copies of the positives.

    python scripts/prepare_shuffled.py --method dinuc     # -> data/processed_dinuc/
    python scripts/prepare_shuffled.py --method mono      # -> data/processed_mono/

Purpose. The primary arm matches negatives on region and GC, which the EDA showed is not
enough: dinucleotide arrangement still differs sharply (TARDBP GU +1.94 log2). Here each
negative is a shuffle of its own positive, so mononucleotide and (for `dinuc`)
dinucleotide frequencies are identical BY CONSTRUCTION. Whatever a model can still do is
motif recognition, not composition.

This is an ADDITIONAL arm. data/processed/ is untouched and remains the primary result.

Shuffled sequences have no genomic location, so `chrom`/`start` carry the source
positive's coordinates with a `shuf_` id prefix. They inherit the source's split, which
keeps the split assignment identical to the primary arm.
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd  # noqa: E402

from rbp.data import shuffle as sh  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
COLS = ["id", "label", "chrom", "start", "end", "strand", "region", "gc", "split",
        "seq_dna", "seq_rna"]


def gc(s):
    u = s.upper()
    return round((u.count("G") + u.count("C")) / len(u), 4) if u else 0.0


def build(protein, method, seed, src, out):
    df = pd.read_csv(src / protein / "dataset.tsv", sep="\t")
    pos = df[df.label == 1].reset_index(drop=True)
    seqs, dropped = sh.shuffled_negatives(pos.seq_rna.tolist(), seed=seed, method=method)

    rows, n = [], 0
    for i, r in pos.iterrows():
        s = seqs[i]
        if s is None:
            continue
        rows.append({**{c: r[c] for c in COLS if c != "id"},
                     "id": f"{protein}_pos_{n}", "label": 1})
        rows.append({"id": f"{protein}_shuf_{n}", "label": 0,
                     "chrom": r.chrom, "start": r.start, "end": r.end,
                     "strand": r.strand, "region": r.region, "gc": gc(s),
                     "split": r.split, "seq_dna": s.replace("U", "T"), "seq_rna": s})
        n += 1

    d = out / protein
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "dataset.tsv", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for x in rows:
            fh.write("\t".join(str(x[c]) for c in COLS) + "\n")
    return {"protein": protein, "pairs": n, "rows": len(rows),
            "test_pairs": sum(1 for x in rows if x["split"] == "test" and x["label"] == 1),
            **dropped}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--method", default="dinuc", choices=["dinuc", "mono"])
    p.add_argument("--protein", default=None)
    p.add_argument("--cell", default=None, help="default: encode.cell_line")
    p.add_argument("--arm", default="gc", choices=sorted(panelmod.ARMS),
                   help="which arm to shuffle FROM. gc, because the shuffle arm exists "
                        "to be compared against the original matching")
    p.add_argument("--src", default=None)
    a = p.parse_args()
    cfg = cfgmod.load()
    cell = a.cell or cfg.encode["cell_line"]
    src = Path(a.src) if a.src else panelmod.data_dir(cell, a.arm)
    # Keyed by cell line, or the second cell line silently overwrites the first.
    out = ROOT / f"data/processed_{a.method}" / cell
    names = ([a.protein] if a.protein else
             [n for n, pairs in panelmod.read_panel(cell, a.arm)
              if pairs >= cfg.cv["min_pairs"]])
    if not names:
        sys.exit(f"no panel for {cell} {a.arm}")

    print(f"{a.method} shuffle of {cell} {a.arm} -> {out.relative_to(ROOT)}\n")
    print(f"{'protein':9} {'pairs':>7} {'testpr':>7} {'failed':>7} {'similar':>8} {'sec':>6}")
    reps = []
    for nm in names:
        t0 = time.time()
        r = build(nm, a.method, cfg.seed, src, out)
        reps.append(r)
        print(f"{nm:9} {r['pairs']:7d} {r['test_pairs']:7d} {r['failed']:7d} "
              f"{r['too_similar']:8d} {time.time()-t0:6.1f}", flush=True)
    tot = sum(r["rows"] for r in reps)
    print(f"\n{len(reps)} proteins, {tot:,} rows, "
          f"{sum(r['failed'] for r in reps)} failed, "
          f"{sum(r['too_similar'] for r in reps)} rejected as too similar")


if __name__ == "__main__":
    main()
