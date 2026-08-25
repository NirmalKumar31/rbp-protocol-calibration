"""Unit tests for the tricky parts of GTF parsing: UTR sides, introns, coordinates."""

import numpy as np
import pytest

from rbp.data import annotation as ann


def test_attr_extraction():
    field = 'gene_id "ENSG1"; transcript_id "ENST9"; gene_name "FOO";'
    assert ann._attr(field, "transcript_id") == "ENST9"
    assert ann._attr(field, "gene_id") == "ENSG1"
    assert ann._attr(field, "missing") is None


class TestSplitUtrs:
    """On the plus strand the 5' UTR sits before the CDS; on the minus strand, after."""

    cds = [(1000, 2000)]

    def test_plus_strand(self):
        five, three = ann.split_utrs([(500, 1000), (2000, 2500)], self.cds, "+")
        assert five == [(500, 1000)]
        assert three == [(2000, 2500)]

    def test_minus_strand_is_reversed(self):
        five, three = ann.split_utrs([(500, 1000), (2000, 2500)], self.cds, "-")
        assert five == [(2000, 2500)]
        assert three == [(500, 1000)]

    def test_noncoding_transcript_has_neither(self):
        assert ann.split_utrs([(500, 1000)], [], "+") == ([], [])

    def test_multi_exon_cds_uses_outer_bounds(self):
        cds = [(1000, 1200), (1800, 2000)]
        five, three = ann.split_utrs([(500, 1000), (2000, 2500)], cds, "+")
        assert five == [(500, 1000)] and three == [(2000, 2500)]


class TestIntrons:
    def test_gaps_between_exons(self):
        assert ann.introns_of([(0, 100), (200, 300), (500, 600)]) == [(100, 200), (300, 500)]

    def test_single_exon_has_none(self):
        assert ann.introns_of([(0, 100)]) == []

    def test_unsorted_input(self):
        assert ann.introns_of([(500, 600), (0, 100)]) == [(100, 500)]

    def test_adjacent_exons_make_no_intron(self):
        assert ann.introns_of([(0, 100), (100, 200)]) == []

    def test_overlapping_exons_do_not_make_negative_introns(self):
        # alternative isoform exons can overlap; must not emit a reversed interval
        out = ann.introns_of([(0, 200), (100, 300), (500, 600)])
        assert out == [(300, 500)]
        assert all(e > s for s, e in out)


class TestRegionsOf:
    def test_coding_transcript(self):
        rec = {"chrom": "chr1", "strand": "+",
               "exon": [(500, 1000), (1000, 2000), (2000, 2500)],
               "CDS": [(1000, 2000)],
               "UTR": [(500, 1000), (2000, 2500)]}
        r = ann.regions_of(rec)
        assert r["utr5"] == [(500, 1000)]
        assert r["utr3"] == [(2000, 2500)]
        assert r["cds"] == [(1000, 2000)]
        assert "exon_nc" not in r          # coding, so exons are not labelled noncoding

    def test_noncoding_transcript(self):
        rec = {"chrom": "chr1", "strand": "+",
               "exon": [(0, 100), (200, 300)], "CDS": [], "UTR": []}
        r = ann.regions_of(rec)
        assert r["exon_nc"] == [(0, 100), (200, 300)]
        assert r["intron"] == [(100, 200)]
        assert "utr5" not in r and "cds" not in r


def _index(intervals):
    """Build a minimal index for classify() tests, merged as build_index would."""
    out = {}
    for region, ivs in intervals.items():
        m = ann.merge_intervals(ivs)
        s = np.array([a for a, _ in m], np.int64)
        e = np.array([b for _, b in m], np.int64)
        out[region] = {"chr1": (s, e)}
    return out


class TestClassify:
    def test_midpoint_inside_interval(self):
        idx = _index({"intron": [(0, 1000)]})
        assert ann.classify(idx, "chr1", 400, 501) == "intron"

    def test_outside_returns_none(self):
        idx = _index({"intron": [(0, 100)]})
        assert ann.classify(idx, "chr1", 500, 601) is None

    def test_unknown_chromosome_returns_none(self):
        idx = _index({"intron": [(0, 100)]})
        assert ann.classify(idx, "chrZ", 0, 101) is None

    def test_priority_prefers_exonic_over_intron(self):
        # same midpoint claimed by a UTR of one isoform and an intron of another
        idx = _index({"utr3": [(0, 1000)], "intron": [(0, 1000)]})
        assert ann.classify(idx, "chr1", 400, 501) == "utr3"

    def test_finds_interval_when_not_the_last_one_before_midpoint(self):
        # a long interval followed by short ones that start before the midpoint
        idx = _index({"intron": [(0, 10_000), (100, 200), (300, 400)]})
        assert ann.classify(idx, "chr1", 4950, 5051) == "intron"


@pytest.mark.parametrize("gtf_start,gtf_end,bed_start,bed_end", [(1, 100, 0, 100),
                                                                 (500, 500, 499, 500)])
def test_gtf_to_bed_coordinate_convention(gtf_start, gtf_end, bed_start, bed_end):
    """GTF is 1-based inclusive; we store 0-based half-open."""
    assert (gtf_start - 1, gtf_end) == (bed_start, bed_end)


class TestMergeIntervals:
    def test_overlapping_collapse(self):
        assert ann.merge_intervals([(0, 100), (50, 200)]) == [(0, 200)]

    def test_touching_collapse(self):
        assert ann.merge_intervals([(0, 100), (100, 200)]) == [(0, 200)]

    def test_disjoint_preserved_and_sorted(self):
        assert ann.merge_intervals([(300, 400), (0, 100)]) == [(0, 100), (300, 400)]

    def test_nested_absorbed(self):
        assert ann.merge_intervals([(0, 1000), (100, 200)]) == [(0, 1000)]

    def test_empty(self):
        assert ann.merge_intervals([]) == []
