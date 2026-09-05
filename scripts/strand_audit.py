"""Are the negatives real RNA? And does the answer break the model-dependence claim?

THE BUG THIS AUDITS, which is in this pipeline and in the published benchmarks it copies.
`negatives.py:328` gives each sampled negative the POSITIVE's strand, and `to_rna` reverse
complements on that basis. A negative's genomic location is chosen to match composition, not
strand, so whether it lands on a gene transcribed in the assigned direction is a coin flip.
Measured here: only ~55% of negatives carry the strand their own gene is on. The other ~45%
are antisense sequence that no transcript in the cell ever produces, while 100% of positives
are true sense RNA.

WHY IT THREATENS THE HEADLINE. Sense/antisense is a NON-COMPOSITIONAL cue. Dinucleotide
matching is done on forward DNA and revcomp is applied to both members of a pair, so the
composition match survives into RNA space -- but the directional cue does not. A model
pretrained on pre-mRNA can read directionality (splice motifs, polypyrimidine tracts, polyA
signals all have a direction); a mononucleotide+dinucleotide model largely cannot. So "the
composition share is lower for SpliceBERT" has an alternative reading: SpliceBERT is not
finding more motif, it is noticing which negatives are backwards.

WHY THAT READING LOSES, and this is the whole point of the script. The claim is a CONTRAST
between two shares, not a level. The artifact has to act DIFFERENTIALLY on the two models to
create a contrast. It does not:

  - frac_sense vs k-mer share          rho +0.41, p=0.009
  - frac_sense vs SpliceBERT share     rho +0.52, p=0.001
  - frac_sense vs THE CONTRAST         rho +0.23, p=0.16   <- not significant

and the model-free version, splitting the panel at the median sense fraction:

  antisense-rich half   k-mer 0.578  SpliceBERT 0.329  contrast +0.249
  antisense-poor half   k-mer 0.785  SpliceBERT 0.520  contrast +0.264

The artifact moves both shares hard and in the same direction. The contrast barely moves, and
what movement there is goes the WRONG way for the objection. So the bug is real, it inflates
every absolute number in R1, and it does not manufacture the model-dependence.

WHAT THIS SCRIPT DOES NOT SETTLE. The sense fraction has a narrow spread (sd 0.038, range
0.43-0.62), so a correlation across datasets has limited power to detect a differential
effect. The stratified contrast is the more trustworthy half of the evidence because it is a
direct comparison rather than a regression on a barely-varying regressor. Fixing the sampler
and re-running is the real answer; this bounds the damage in the meantime.
"""

import argparse
import glob
import gzip
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from rbp.utils.log import log

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "tables"



def gene_index(gtf):
    """chrom -> (starts, ends, strands), sorted by start, for interval lookup."""
    genes = defaultdict(list)
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as f:
        for ln in f:
            if ln[0] == "#":
                continue
            p = ln.split("\t", 8)
            if p[2] != "gene":
                continue
            genes[p[0]].append((int(p[3]) - 1, int(p[4]), p[6]))
    idx = {}
    for c, v in genes.items():
        v.sort()
        idx[c] = (np.array([g[0] for g in v]), np.array([g[1] for g in v]),
                  np.array([g[2] for g in v]))
    log(f"gene index: {sum(len(v) for v in genes.values()):,} genes, {len(idx)} contigs")
    return idx


def own_strands(idx, chrom, start, end, back=400):
    """Strands of every gene overlapping this window.

    `back` bounds the linear scan: genes are sorted by start, so an overlapping gene must
    start before `end`, and only the last few hundred can also end after `start` unless a
    single gene is megabases long. 400 covers GENCODE.
    """
    if chrom not in idx:
        return set()
    s, e, st = idx[chrom]
    i = np.searchsorted(s, end, side="right")
    lo = max(0, i - back)
    hit = (s[lo:i] < end) & (e[lo:i] > start)
    return set(st[lo:i][hit])


def audit(idx, files):
    rows = []
    for f in sorted(files):
        d = pd.read_csv(f, sep="\t", usecols=["label", "chrom", "start", "end", "strand"])
        neg = d[d.label == 0]
        ok = amb = none = 0
        for c, s, e, a in zip(neg.chrom, neg.start, neg.end, neg.strand):
            ss = own_strands(idx, c, s, e)
            if not ss:
                none += 1
            elif len(ss) > 1:
                amb += 1                 # gene on each strand; direction is genuinely unclear
            elif a in ss:
                ok += 1
        n = len(neg)
        det = n - none - amb
        cell, prot = os.path.basename(f)[:-4].split("_", 1)
        rows.append({"dataset": f"{prot}:{cell}", "n_neg": n, "unambiguous": det,
                     "frac_sense": ok / det if det else np.nan,
                     "frac_no_gene": none / n, "frac_ambiguous": amb / n})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--datasets", required=True,
                   help="glob of processed dataset.tsv files, named {cell}_{protein}.tsv")
    a = p.parse_args()

    r = audit(gene_index(a.gtf), glob.glob(a.datasets))
    if r.empty:
        raise SystemExit(f"no files matched {a.datasets}")

    log(f"\n{len(r)} datasets")
    log(f"negatives whose ASSIGNED strand matches their own gene: mean {r.frac_sense.mean():.1%}"
        f" (range {r.frac_sense.min():.1%}-{r.frac_sense.max():.1%}, sd {r.frac_sense.std():.3f})")
    log(f"negatives in no annotated gene:        {r.frac_no_gene.mean():.1%}")
    log(f"negatives in genes on BOTH strands:    {r.frac_ambiguous.mean():.1%}")

    fm = pd.read_csv(TABLES / "matched_four_models.csv")
    m = r.merge(fm[["dataset", "composition_auroc", "kmer_auroc", "splicebert"]], on="dataset")
    m["share_kmer"] = (m.composition_auroc - 0.5) / (m.kmer_auroc - 0.5)
    m["share_splicebert"] = (m.composition_auroc - 0.5) / (m.splicebert - 0.5)
    m["share_contrast"] = m.share_kmer - m.share_splicebert
    m = m.replace([np.inf, -np.inf], np.nan).dropna(subset=["share_contrast"])

    out = []
    for c, lab in (("share_contrast", "THE CLAIM: k-mer share minus SpliceBERT share"),
                   ("share_kmer", "k-mer share"),
                   ("share_splicebert", "SpliceBERT share"),
                   ("composition_auroc", "composition AUROC")):
        rho, pv = spearmanr(m.frac_sense, m[c])
        out.append({"check": f"frac_sense vs {lab}", "value": rho, "p": pv, "n": len(m)})
        log(f"  frac_sense vs {lab:46} rho {rho:+.4f} p={pv:.3f}")

    med = m.frac_sense.median()

    def sh(x, col):
        return (x.composition_auroc.mean() - 0.5) / (x[col].mean() - 0.5)

    log("")
    for lab, sub in (("antisense-rich half", m[m.frac_sense <= med]),
                     ("antisense-poor half", m[m.frac_sense > med])):
        k, s = sh(sub, "kmer_auroc"), sh(sub, "splicebert")
        log(f"  {lab:22} n={len(sub):2d}  k-mer {k:.4f}  SpliceBERT {s:.4f}  contrast {k - s:+.4f}")
        out.append({"check": f"share contrast, {lab}", "value": k - s, "p": np.nan,
                    "n": len(sub)})

    r.to_csv(TABLES / "strand_audit.csv", index=False)
    pd.DataFrame(out).to_csv(TABLES / "strand_audit_summary.csv", index=False)
    log("\nwrote strand_audit.csv, strand_audit_summary.csv")


if __name__ == "__main__":
    main()
