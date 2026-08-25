"""Load ClinVar variants, reconciled to the coordinate conventions used everywhere else.

ClinVar's VCF names chromosomes `1`, `2`, `MT`; the genome, GENCODE annotation and
ENCODE peaks all use `chr1`, `chr2`, `chrM`. Joining them without normalising matches
nothing and raises no error, so normalisation happens here, once, at the point of
loading, and is covered by tests.

The VCF is also 1-based while every interval in this project is 0-based half-open, so
positions are converted on the way in.
"""

import gzip
import re

# ClinVar reports a molecular consequence as MC=SO:0001627|intron_variant
_MC = re.compile(r"MC=([^;]+)")
_INFO = {
    "CLNSIG": re.compile(r"(?:^|;)CLNSIG=([^;]+)"),
    "CLNVC": re.compile(r"(?:^|;)CLNVC=([^;]+)"),
    "GENEINFO": re.compile(r"(?:^|;)GENEINFO=([^;]+)"),
}


def normalise_chrom(chrom):
    """ClinVar-style chromosome name -> the `chr`-prefixed form used project-wide."""
    c = str(chrom).strip()
    if c.startswith("chr"):
        base = c[3:]
    else:
        base = c
    if base in ("MT", "M"):
        return "chrM"
    return f"chr{base}"


def consequences(info):
    """Every molecular consequence term on the record, lowercased."""
    m = _MC.search(info)
    if not m:
        return []
    out = []
    for part in m.group(1).split(","):
        term = part.split("|")[-1].strip().lower()
        if term:
            out.append(term)
    return out


def _field(info, key):
    m = _INFO[key].search(info)
    return m.group(1) if m else None


def is_snv(ref, alt, clnvc=None):
    """Single-nucleotide substitution only: one base to one different base."""
    if clnvc and clnvc.lower() != "single_nucleotide_variant":
        return False
    return (len(ref) == 1 and len(alt) == 1
            and ref in "ACGT" and alt in "ACGT" and ref != alt)


def load(vcf, pathogenic, benign, noncoding=None, chroms=None, snv_only=True):
    """Yield dicts for variants whose significance is an exact match to the label sets.

    Positions come out 0-based. `noncoding`, when given, keeps only records carrying at
    least one of those consequence terms. Ambiguous significance strings such as
    "Pathogenic/Likely_pathogenic" or "Conflicting_classifications_of_pathogenicity"
    are excluded by construction, because only exact members of the label sets pass.
    """
    want_sig = {s: 1 for s in pathogenic} | {s: 0 for s in benign}
    nc = {t.lower() for t in noncoding} if noncoding else None
    keep_chrom = set(chroms) if chroms else None

    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 8:
                continue
            info = f[7]
            sig = _field(info, "CLNSIG")
            if sig is None or sig not in want_sig:
                continue
            ref, alt = f[3], f[4]
            if snv_only and not is_snv(ref, alt, _field(info, "CLNVC")):
                continue
            terms = consequences(info)
            if nc is not None and not (set(terms) & nc):
                continue
            chrom = normalise_chrom(f[0])
            if keep_chrom is not None and chrom not in keep_chrom:
                continue
            pos1 = int(f[1])
            yield {
                "vid": f"{chrom}:{pos1}:{ref}:{alt}",
                "chrom": chrom,
                "pos": pos1 - 1,          # 0-based, matching every interval here
                "pos_vcf": pos1,
                "ref": ref,
                "alt": alt,
                "clnsig": sig,
                "label": want_sig[sig],
                "consequence": terms[0] if terms else None,
                "gene": (_field(info, "GENEINFO") or "").split(":")[0] or None,
            }
