"""Stages 3b-3d: positives, matched negatives, and chromosome splits for one protein
or all of them.

    python scripts/prepare.py                      # every protein in the panel
    python scripts/prepare.py --protein RBFOX2     # just one (array-friendly)

Writes data/processed/<PROTEIN>/dataset.tsv and a per-protein prep report. Positives
that cannot be matched to a negative are dropped so the classes stay balanced 1:1,
and every drop is counted in the report rather than passing unnoticed.
"""

import argparse
import os
import pickle
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
from pyfaidx import Fasta  # noqa: E402

from rbp.data import annotation as ann  # noqa: E402
from rbp.data import encode  # noqa: E402
from rbp.data import negatives as neg  # noqa: E402
from rbp.data import splits, windows as win  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
COLS = ["id", "label", "chrom", "start", "end", "strand", "region", "gc", "split",
        "fold", "seq_dna", "seq_rna"]


def load_fold_map(cfg):
    """chrom -> fold, from the frozen assignment solved by optimize_folds.py."""
    f = ROOT / cfg.cv["folds"]
    if not f.exists():
        raise SystemExit(f"missing {f} - run scripts/optimize_folds.py first")
    lines = f.read_text().strip().splitlines()[1:]
    return {c: int(k) for c, k in (ln.split("\t") for ln in lines)}


def existing_report(outdir, protein, k, match="gc"):
    """Parse a previous run's report if the dataset is already complete and current.

    Preparing the full panel is ~2 hours of CPU, and without this a restart -- a laptop
    sleeping, a cluster preemption, a crash on one protein -- redoes all of it. The check
    is not just "the file exists": a report written under a different fold count belongs
    to a different protocol and must be redone.
    """
    d = outdir / protein
    rep, ds = d / "prep_report.txt", d / "dataset.tsv"
    if not (rep.exists() and ds.exists()):
        return None
    out = {}
    for ln in rep.read_text().splitlines():
        if ":" not in ln:
            continue
        key, val = ln.split(":", 1)
        out[key.strip()] = val.strip()
    if "fold_proportions" not in out:
        return None                       # written before folds existed
    # A dataset built with a different negative design is a different dataset. Reusing a
    # GC-matched one for a dinucleotide-matched run would silently mix the two arms.
    if out.get("negative_match", "gc") != match:
        return None
    try:
        folds = eval(out["fold_proportions"], {"__builtins__": {}})
        if len(folds) != k:
            return None
        return {"protein": protein, "peaks": int(out["peaks"]),
                "pairs": int(out["pairs"]), "rows": int(out["rows"]),
                "seconds": 0.0, "fold_proportions": folds,
                "fold_problems": eval(out.get("fold_problems", "[]"),
                                      {"__builtins__": {}}),
                "negative_match": out.get("negative_match", "gc"),
                "cached": True}
    except (KeyError, ValueError, SyntaxError):
        return None


def prepare_one(protein, cfg, fasta, index, outdir, fold_map, match="gc"):
    t0 = time.time()
    size = cfg.windows["size"]
    pth = encode.peak_path(ROOT, protein, cfg.encode["cell_line"])

    drop_chroms = set(cfg.encode.get("exclude_chroms", []))
    keep = [c for c in ann.MAIN_CHROMS if c not in drop_chroms]
    pos, pos_dropped = win.build_positives(
        pth, fasta, index, size, ann.classify, drop_n=cfg.windows["drop_n"],
        chroms=keep)
    peaks = [p_ for p_ in win.read_peaks(pth) if p_[0] in set(keep)]

    if match == "dinuc":
        negs, neg_dropped, dists = neg.build_negatives_dinuc(
            pos, peaks, fasta, index, size,
            min_peak_distance=cfg.negatives["min_peak_distance"],
            seed=cfg.seed, drop_n=cfg.windows["drop_n"])
    else:
        negs, neg_dropped = neg.build_negatives(
            pos, peaks, fasta, index, size,
            tolerance=cfg.negatives["gc_tolerance"],
            min_peak_distance=cfg.negatives["min_peak_distance"],
            seed=cfg.seed, drop_n=cfg.windows["drop_n"])
        dists = None

    rows = []
    n = 0
    for p, m in zip(pos, negs):
        if m is None:
            continue
        for label, rec in ((1, p), (0, m)):
            rows.append({"id": f"{protein}_{'pos' if label else 'neg'}_{n}",
                         "label": label, **{k: rec[k] for k in
                                            ("chrom", "start", "end", "strand",
                                             "region", "gc", "seq_dna", "seq_rna")}})
        n += 1

    splits.assign(rows, cfg.split)
    leaked = splits.check_disjoint(rows)
    if leaked:
        raise SystemExit(f"FATAL {protein}: chromosomes in multiple splits: {leaked}")

    # The CV fold is stored once per row and the train/val/test roles are derived per
    # iteration, so one dataset file serves all k runs. Writing k copies would multiply
    # storage and let them drift apart.
    splits.assign_folds(rows, fold_map)
    # A tiny protein can genuinely have no peaks on some fold's chromosomes, which makes
    # the k-fold protocol unrunnable for it -- but every such dataset measured so far is
    # far below the inclusion threshold and would be dropped anyway. So this is recorded
    # as an exclusion reason and only escalated to fatal for datasets that DO clear the
    # threshold, where an empty fold would mean something changed and we need to know.
    # Aborting the whole panel run on AARS (32 windows) was the wrong behaviour.
    fold_problems = splits.check_folds(rows, cfg.cv["k"])

    d = outdir / protein
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "dataset.tsv", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")

    prop = splits.proportions(rows)
    report = {
        "protein": protein, "peaks": len(peaks), "positives": len(pos), "pairs": n,
        "test_pairs": sum(1 for x in rows if x["split"] == "test" and x["label"] == 1),
        "rows": len(rows), "seconds": round(time.time() - t0, 1),
        "pos_dropped": pos_dropped, "neg_dropped": neg_dropped,
        "split_proportions": {k: round(v, 3) for k, v in prop.items()},
        "fold_proportions": {f: round(v, 3) for f, v in
                             splits.fold_proportions(rows, cfg.cv["k"]).items()},
        "fold_problems": fold_problems,
        "negative_match": match,
        "regions": {r: sum(1 for x in rows if x["region"] == r and x["label"] == 1)
                    for r in ann.REGIONS},
    }
    if dists is not None:
        ok = dists[~np.isnan(dists)]
        # Match quality has to be reported, not assumed: matching 16 frequencies can fail
        # quietly and leave the "matched" negatives no better than GC-matched ones.
        report["dinuc_l1_median"] = round(float(np.median(ok)), 5) if len(ok) else None
        report["dinuc_l1_p90"] = round(float(np.percentile(ok, 90)), 5) if len(ok) else None

    (d / "prep_report.txt").write_text(
        "\n".join(f"{k}: {v}" for k, v in report.items()) + "\n")
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--protein", default=None)
    p.add_argument("--index", default=str(ROOT / "data/interim/regions.pkl"))
    p.add_argument("--outdir", default=None, help="default data/processed/<cell_line>")
    p.add_argument("--panel", default=None, help="panel TSV, default config/proteins.tsv")
    p.add_argument("--cell-line", default=None, help="override encode.cell_line")
    p.add_argument("--match", default="gc", choices=["gc", "dinuc"],
                   help="negative matching: gc (primary) or dinuc (stronger control)")
    p.add_argument("--force", action="store_true",
                   help="re-prepare datasets that are already complete")
    a = p.parse_args()
    cfg = cfgmod.load(a.config)
    if a.cell_line:
        cfg["encode"]["cell_line"] = a.cell_line
    cell = cfg.encode["cell_line"]

    names = ([a.protein] if a.protein
             else [r["protein"] for r in cfgmod.proteins(panel=a.panel)])
    if not a.protein and "SLURM_ARRAY_TASK_ID" in os.environ:
        names = [names[int(os.environ["SLURM_ARRAY_TASK_ID"])]]

    print("loading region index, fold map and genome ...", flush=True)
    index = pickle.loads(Path(a.index).read_bytes())
    fold_map = load_fold_map(cfg)
    fasta = Fasta(str(ROOT / "data/raw/GRCh38.primary_assembly.genome.fa"))
    # Default output follows the ARM, not just the cell line. It used to be
    # data/processed/<cell> unconditionally, so `--match dinuc` without an explicit
    # --outdir wrote dinucleotide-matched datasets on top of the GC-matched ones.
    outdir = Path(a.outdir) if a.outdir else panelmod.data_dir(cell, a.match)

    k = cfg.cv["k"]
    min_pairs = cfg.cv["min_pairs"]
    print(f"{cell}: {len(names)} proteins, k={k}, threshold {min_pairs} pairs\n")
    print(f"{'protein':10} {'peaks':>7} {'pairs':>7} {'fold spread':>13} {'sec':>6}")
    reports = []
    n_cached = 0
    for name in names:
        r = None if a.force else existing_report(outdir, name, k, a.match)
        if r is None:
            r = prepare_one(name, cfg, fasta, index, outdir, fold_map, a.match)
        else:
            n_cached += 1
        reports.append(r)
        fp = r["fold_proportions"]
        flag = "" if r["pairs"] >= min_pairs else f"  BELOW min_pairs={min_pairs}"
        if r.get("cached"):
            flag += "  cached"
        print(f"{name:10} {r['peaks']:7d} {r['pairs']:7d} "
              f"{min(fp.values()):.2f}-{max(fp.values()):.2f}"
              f"{'':>5} {r['seconds']:6.1f}{flag}", flush=True)
    if n_cached:
        print(f"\n  reused {n_cached} already-prepared datasets; "
              f"--force to redo them")

    if len(reports) > 1:
        def reason(r):
            if r["pairs"] < min_pairs:
                return f"pairs<{min_pairs}"
            if r["fold_problems"]:
                return "; ".join(r["fold_problems"])
            return None

        keep = [r for r in reports if reason(r) is None]
        drop = [(r, reason(r)) for r in reports if reason(r) is not None]
        # An empty fold in a dataset that clears the threshold means the fold map and the
        # data have gone out of step. Silently excluding it would shrink the panel without
        # anyone noticing why.
        broken = [r for r, why in drop if r["pairs"] >= min_pairs]
        print(f"\n{cell}: {len(keep)} of {len(reports)} clear min_pairs={min_pairs}")
        if drop:
            print("  dropped: " + ", ".join(f"{r['protein']}({r['pairs']})"
                                            for r, _ in drop))
        print(f"  total rows across kept datasets: {sum(r['rows'] for r in keep):,}")
        if broken:
            raise SystemExit(
                "FATAL: these datasets clear the threshold but have a broken fold "
                "assignment, so the fold map no longer matches the data: "
                + ", ".join(f"{r['protein']} {r['fold_problems']}" for r in broken))

        # The final panel is decided by the data, so record it explicitly rather than
        # leaving it implicit in which datasets happen to exist on disk. Keyed by cell
        # line AND by negative arm: a single panel_final.tsv would have the second cell
        # line overwrite the first and silently halve the study, and for a while the
        # filename omitted the arm, so the dinuc run overwrote the GC run's pair counts
        # and nothing said so. One helper writes both files so the two axes cannot drift.
        pf, xf = panelmod.write_panel(
            cell, a.match,
            [(r["protein"], r["pairs"]) for r in keep],
            [(r["protein"], r["pairs"], why) for r, why in drop])
        print(f"  wrote {pf.relative_to(ROOT)} ({len(keep)}) and "
              f"{xf.relative_to(ROOT)} ({len(drop)})")


if __name__ == "__main__":
    main()
