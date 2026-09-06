"""Were the positives filtered? ENCODE already did it, and the Limitations used to deny it.

    python scripts/peak_thresholds.py               # downloads the 95 narrowPeak files
    python scripts/peak_thresholds.py --from-cache

An earlier draft conceded that "the positive set includes weak peaks that a study applying
ENCODE's standard thresholds would exclude". That is false, and it was the kind of concession
that invites a rejection for a flaw the study does not have. ENCODE's released eCLIP `peaks`
files are already thresholded, so there were no weak peaks available to exclude. This reads
every file in the panel and records the minimum signal and significance actually present.

The real remaining limitation is peak WIDTH, which the pipeline discards: every positive is a
fixed 101 nt window on the midpoint, so a 20 nt peak and a 2 kb peak are treated identically.
"""

import argparse
import gzip
import io
import sys
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
URL = "https://www.encodeproject.org/files/{a}/@@download/{a}.bed.gz"
# Van Nostrand et al. 2020: log2 fold-enrichment >= 3 and p <= 0.001.
FC_MIN = 3.0
P_MIN = 3.0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--from-cache", action="store_true")
    a = p.parse_args()

    per = TABLES / "peak_thresholds_per_dataset.csv"
    if a.from_cache:
        t = pd.read_csv(per)
    else:
        panel = pd.read_csv(TABLES / "supplementary_table_s1.csv")
        rows = []
        for i, r in enumerate(panel.itertuples(), 1):
            try:
                raw = urllib.request.urlopen(URL.format(a=r.accession), timeout=120).read()
            except Exception as e:
                sys.exit(f"could not fetch {r.accession}: {e}")
            b = pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep="\t", header=None,
                            usecols=[3, 6, 7], names=["name", "log2fc", "neglog10p"])
            rows.append({"dataset": r.dataset, "accession": r.accession, "peaks": len(b),
                         "min_log2fc": float(b.log2fc.min()),
                         "min_neglog10p": float(b.neglog10p.min()),
                         "idr_named": int(b.name.str.contains("IDR").sum())})
            if i % 20 == 0:
                print(f"  [{i}/{len(panel)}]", flush=True)
        t = pd.DataFrame(rows)
        t.to_csv(per, index=False)

    n_peaks = int(t.peaks.sum())
    idr = int(t.idr_named.sum())
    out = [{"check": "peak files audited", "value": len(t), "n": len(t), "note": ""},
           {"check": "peaks in the panel", "value": n_peaks, "n": len(t), "note": ""},
           {"check": "minimum log2 fold-enrichment over the panel",
            "value": float(t.min_log2fc.min()), "n": len(t),
            "note": f"threshold {FC_MIN}"},
           {"check": "minimum -log10 p over the panel",
            "value": float(t.min_neglog10p.min()), "n": len(t), "note": f"threshold {P_MIN}"},
           {"check": "fraction of peaks carrying ENCODE's IDR label",
            "value": idr / n_peaks, "n": len(t), "note": f"{idr} peaks"},
           {"check": "datasets with no IDR-named peaks",
            "value": int((t.idr_named == 0).sum()), "n": len(t), "note": ""},
           {"check": "files below either threshold",
            "value": int(((t.min_log2fc < FC_MIN) | (t.min_neglog10p < P_MIN)).sum()),
            "n": len(t), "note": "zero means ENCODE filtered upstream"}]
    pd.DataFrame(out).to_csv(TABLES / "peak_thresholds.csv", index=False)
    print(f"\n=== {len(t)} files, {n_peaks} peaks ===")
    print(f"  minimum log2 fold-enrichment: {t.min_log2fc.min():.4f}  (threshold {FC_MIN})")
    print(f"  minimum -log10 p:             {t.min_neglog10p.min():.4f}  (threshold {P_MIN})")
    print(f"  IDR-named:                    {idr}/{n_peaks} ({100 * idr / n_peaks:.1f}%), "
          f"{int((t.idr_named == 0).sum())} datasets with none")
    print(f"  files below either threshold: "
          f"{int(((t.min_log2fc < FC_MIN) | (t.min_neglog10p < P_MIN)).sum())}")
    print("\nwrote peak_thresholds.csv and peak_thresholds_per_dataset.csv")


if __name__ == "__main__":
    main()
