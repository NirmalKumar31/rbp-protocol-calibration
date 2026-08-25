"""Stage 1: discover the protein panel from the live ENCODE API.

Writes config/proteins.tsv (the single source of truth for every later stage) and
config/panel_report.txt recording what was skipped and why.

    python scripts/build_panel.py --config config/params.yaml
    python scripts/build_panel.py --limit 16          # take the first N by name
    python scripts/build_panel.py --only TARDBP,FUS   # restrict to named proteins
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rbp.data import encode  # noqa: E402
from rbp.utils import config as cfgmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", default=None)
    p.add_argument("--limit", type=int, default=None, help="keep only the first N proteins")
    p.add_argument("--only", default=None, help="comma-separated protein names to keep")
    p.add_argument("--all", action="store_true", help="ignore panel.candidates, keep every protein")
    p.add_argument("--cell-line", default=None, help="override encode.cell_line")
    p.add_argument("--out", default=str(ROOT / "config" / "proteins.tsv"))
    a = p.parse_args()

    cfg = cfgmod.load(a.config)
    if a.cell_line:
        cfg["encode"]["cell_line"] = a.cell_line
    print(f"querying ENCODE: {cfg.encode['assay']} / {cfg.encode['cell_line']} / "
          f"{cfg.encode['assembly']}", flush=True)

    rows, dupes, skipped = encode.build_panel(cfg)
    print(f"found {len(rows)} proteins with '{cfg.encode['output_type']}'")
    print(f"  {len(dupes)} additional experiments for proteins already covered")
    print(f"  {len(skipped)} experiments skipped")

    all_found = len(rows)
    want = None
    if a.only:
        want = {s.strip().upper() for s in a.only.split(",")}
    elif cfg.get("panel", {}).get("mode") == "subset" and not a.all:
        want = {p.upper() for p in cfg["panel"]["candidates"]}
    if want:
        rows = [r for r in rows if r["protein"].upper() in want]
        missing = want - {r["protein"].upper() for r in rows}
        print(f"  panel subset: {len(rows)}/{len(want)} candidates available in ENCODE")
        if missing:
            print(f"  WARNING requested but unavailable: {sorted(missing)}")
    if a.limit:
        rows = rows[:a.limit]
    print(f"  ({all_found} proteins exist in total; the cloud phase can use all of them)")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["protein", "accession", "cell_line", "experiment", "n_replicates"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in cols) + "\n")
    print(f"\nwrote {out} with {len(rows)} proteins")
    print("  " + ", ".join(r["protein"] for r in rows))

    rep = out.parent / "panel_report.txt"
    with open(rep, "w") as fh:
        fh.write(f"assay={cfg.encode['assay']} cell_line={cfg.encode['cell_line']} "
                 f"assembly={cfg.encode['assembly']}\n")
        fh.write(f"output_type={cfg.encode['output_type']}\n\n")
        fh.write(f"SELECTED ({len(rows)})\n")
        for r in rows:
            fh.write(f"  {r['protein']}\t{r['accession']}\t{r['experiment']}\n")
        fh.write(f"\nDUPLICATE EXPERIMENTS ({len(dupes)})\n")
        for r in dupes:
            fh.write(f"  {r['protein']}\t{r['accession']}\t{r['experiment']}\n")
        fh.write(f"\nSKIPPED ({len(skipped)})\n")
        for acc, why in skipped:
            fh.write(f"  {acc}\t{why}\n")
    print(f"wrote {rep}")


if __name__ == "__main__":
    main()
