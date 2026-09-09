"""The record of the external-benchmark search. A curated finding, emitted as a table.

    python scripts/external_search.py

Why a script for something nobody computed. docs/EXTERNAL_BENCHMARK_PROTOCOL.md promises that
every candidate considered is recorded with the criterion it failed on, and a promise like that
is worth what the artefact behind it is worth. These findings came from reading papers and data
releases, not from a computation, so this file IS the record: versioned, diffable, and something
provenance.py can name as a producing script instead of classifying the table as an orphan.

Nothing here is inferred. Each row is a claim about a public dataset that a reader can check
against the cited source, and the outcome column applies the eligibility criteria that were
committed at e76a80c before any of this was looked at.

THE HEADLINE OF THE SEARCH is that the protocol's own premise was wrong. It asserted that the
Horlacher et al. 2023 benchmark releases one negative-set construction and therefore could not
test this paper's central claim. It releases two, over the same positives, for all 223 ENCODE
datasets, in a deposit this repository already had unpacked on disk. See the protocol's addendum.
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rbp.utils.log import log  # noqa: E402

OUT = ROOT / "results" / "tables" / "external_search.csv"
FIELDS = ["candidate", "constructions", "datasets", "outcome", "failed_criterion", "note"]

CANDIDATES = [
    dict(
        candidate="Horlacher et al. 2023, Benchmark-RBP "
                  "(Brief Bioinform 24(5) bbad307; Zenodo 10.5281/zenodo.10600977)",
        constructions="2 (negative-1 uniform transcript positions; "
                      "negative-2 other RBPs' crosslink sites)",
        datasets="223 ENCODE eCLIP, of which 135 are outside our study panel",
        outcome="QUALIFIES on the disjoint subset",
        failed_criterion="",
        note="The protocol asserted this benchmark releases ONE construction and therefore "
             "could not test Claim A. That was wrong: it releases both, per fold, over the same "
             "positives. Analysed in scripts/external_replication.py on the 135 datasets our "
             "panel does not contain, which is what criterion 1 requires",
    ),
    dict(
        candidate="Horlacher et al. 2023, the 45-dataset intersection already published as our "
                  "external arm",
        constructions="2",
        datasets="45",
        outcome="EXCLUDED",
        failed_criterion="1, independence",
        note="All 45 are also in our 94-dataset panel by design, because horlacher_arm.py "
             "deliberately restricts to the overlap so the measurement is on the same proteins. "
             "That makes it an independent CONSTRUCTION and not an independent sample. Retained "
             "and published as what it is",
    ),
    dict(
        candidate="iONMF (Strazar et al. 2016), as redistributed in the Horlacher benchmark",
        constructions="2, in the Horlacher reprocessing",
        datasets="31 mixed-protocol CLIP",
        outcome="EXCLUDED",
        failed_criterion="1, independence of construction lineage",
        note="Its two negative sets are Horlacher's rather than iONMF's, so adding it beside "
             "the ENCODE subset would count the same construction twice instead of adding a "
             "second one",
    ),
    dict(
        candidate="Mukherjee et al. PAR-CLIP, as redistributed in the Horlacher benchmark",
        constructions="2, in the Horlacher reprocessing",
        datasets="59 PAR-CLIP",
        outcome="CANDIDATE, not analysed",
        failed_criterion="",
        note="A genuinely different assay and the strongest remaining candidate for widening "
             "scope beyond eCLIP. Not analysed: a different crosslinking chemistry has "
             "different positional bias, so the 101-nt window convention and the composition "
             "baseline would both need re-justification rather than transporting unchanged, and "
             "the protocol fixes the estimand as transported UNCHANGED. Recorded as the next "
             "experiment, not as a result",
    ),
    dict(
        candidate="RNAcompete and RNA Bind-n-Seq in vitro affinity data",
        constructions="n/a",
        datasets="many",
        outcome="EXCLUDED",
        failed_criterion="exclusion list, in vitro",
        note="The negative is a synthetic library rather than a genomic window, so a "
             "composition baseline over it measures something different",
    ),
    dict(
        candidate="GraphProt, DeepBind, iDeepS, DeepCLIP, PrismNet and Pysster training sets",
        constructions="1 each",
        datasets="varies",
        outcome="EXCLUDED",
        failed_criterion="2, needs two constructions over the same positives",
        note="Each ships the single negative construction its authors chose. That is the state "
             "of practice this paper calibrates, and it is why the search was expected to come "
             "back empty",
    ),
    dict(
        candidate="Benchmarks whose only negatives are dinucleotide shuffles of the positives",
        constructions="1",
        datasets="varies",
        outcome="EXCLUDED",
        failed_criterion="exclusion list, degenerate construction",
        note="Our own shuffled arm pins the composition baseline at exactly 0.5000 on all 94 "
             "datasets, so this is a degenerate case rather than a second protocol",
    ),
]


def main():
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(CANDIDATES)
    n_q = sum(1 for r in CANDIDATES if r["outcome"].startswith("QUALIFIES"))
    log(f"  {len(CANDIDATES)} candidates considered, {n_q} qualifying "
        f"-> {OUT.relative_to(ROOT)}")
    for r in CANDIDATES:
        log(f"    {r['outcome']:26s} {r['candidate'][:60]}")


if __name__ == "__main__":
    main()
