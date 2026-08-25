"""ClinVar loading: the chromosome-naming reconciliation and label strictness."""

import gzip

import pytest

from rbp.variants import clinvar as cv


class TestNormaliseChrom:
    """The bug guard: ClinVar says "1", the genome/GTF/peaks all say "chr1"."""

    @pytest.mark.parametrize("raw,expected", [
        ("1", "chr1"), ("22", "chr22"), ("X", "chrX"), ("Y", "chrY"),
        ("MT", "chrM"), ("M", "chrM"),
        ("chr1", "chr1"), ("chrX", "chrX"), (" 7 ", "chr7"), (7, "chr7"),
    ])
    def test_values(self, raw, expected):
        assert cv.normalise_chrom(raw) == expected

    def test_idempotent(self):
        once = cv.normalise_chrom("1")
        assert cv.normalise_chrom(once) == once


class TestIsSnv:
    def test_accepts_substitution(self):
        assert cv.is_snv("G", "A")

    @pytest.mark.parametrize("ref,alt", [("AG", "A"), ("A", "AG"), ("A", "A"),
                                         ("N", "A"), ("A", "N"), ("", "A")])
    def test_rejects_non_snv(self, ref, alt):
        assert not cv.is_snv(ref, alt)

    def test_respects_clnvc(self):
        assert not cv.is_snv("G", "A", clnvc="Deletion")
        assert cv.is_snv("G", "A", clnvc="single_nucleotide_variant")


class TestConsequences:
    def test_single(self):
        assert cv.consequences("MC=SO:0001627|intron_variant") == ["intron_variant"]

    def test_multiple(self):
        got = cv.consequences("AF=0.1;MC=SO:1|intron_variant,SO:2|3_prime_UTR_variant;X=1")
        assert got == ["intron_variant", "3_prime_utr_variant"]

    def test_absent(self):
        assert cv.consequences("AF=0.1") == []


HEADER = "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _vcf(tmp_path, rows):
    p = tmp_path / "clinvar.vcf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(HEADER)
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return p


def test_load_normalises_and_filters(tmp_path):
    rows = [
        # kept: pathogenic SNV, intronic, chromosome gets the chr prefix
        ("1", "100", "1", "G", "A",  ".", ".",
         "CLNSIG=Pathogenic;CLNVC=single_nucleotide_variant;MC=SO:1|intron_variant;GENEINFO=FOO:1"),
        # kept: benign SNV in a 3'UTR
        ("X", "200", "2", "C", "T", ".", ".",
         "CLNSIG=Benign;CLNVC=single_nucleotide_variant;MC=SO:2|3_prime_UTR_variant"),
        # dropped: ambiguous significance must not slip through
        ("2", "300", "3", "G", "A", ".", ".",
         "CLNSIG=Pathogenic/Likely_pathogenic;CLNVC=single_nucleotide_variant;MC=SO:1|intron_variant"),
        ("2", "301", "4", "G", "A", ".", ".",
         "CLNSIG=Conflicting_classifications_of_pathogenicity;CLNVC=single_nucleotide_variant;MC=SO:1|intron_variant"),
        ("2", "302", "5", "G", "A", ".", ".",
         "CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;MC=SO:1|intron_variant"),
        # dropped: an indel
        ("3", "400", "6", "AG", "A", ".", ".",
         "CLNSIG=Pathogenic;CLNVC=Deletion;MC=SO:1|intron_variant"),
        # dropped: coding consequence when a noncoding filter is applied
        ("4", "500", "7", "G", "A", ".", ".",
         "CLNSIG=Pathogenic;CLNVC=single_nucleotide_variant;MC=SO:9|missense_variant"),
    ]
    vcf = _vcf(tmp_path, rows)
    out = list(cv.load(vcf, ["Pathogenic"], ["Benign"],
                       noncoding=["intron_variant", "3_prime_UTR_variant"]))

    assert [v["chrom"] for v in out] == ["chr1", "chrX"]     # normalised
    assert [v["label"] for v in out] == [1, 0]
    assert [v["pos"] for v in out] == [99, 199]              # 0-based
    assert [v["pos_vcf"] for v in out] == [100, 200]
    assert out[0]["vid"] == "chr1:100:G:A"
    assert out[0]["gene"] == "FOO"


def test_load_chromosome_filter_uses_normalised_names(tmp_path):
    rows = [("1", "100", "1", "G", "A", ".", ".",
             "CLNSIG=Pathogenic;CLNVC=single_nucleotide_variant;MC=SO:1|intron_variant")]
    vcf = _vcf(tmp_path, rows)
    assert len(list(cv.load(vcf, ["Pathogenic"], ["Benign"], chroms=["chr1"]))) == 1
    assert len(list(cv.load(vcf, ["Pathogenic"], ["Benign"], chroms=["1"]))) == 0
