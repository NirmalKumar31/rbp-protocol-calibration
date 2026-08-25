"""Tests for attaching variants to binding sites and building their scoring windows.

Strand, the reference base and the fold are the three things that corrupt the variant
analysis without raising anything. Each gets its own constructed case here.
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp.data.windows import revcomp  # noqa: E402
from rbp.variants import assign  # noqa: E402

SIZE = 101


class FakeSeq:
    def __init__(self, s):
        self.seq = s


class FakeChrom:
    """Minimal stand-in for a pyfaidx chromosome: slicing and len()."""

    def __init__(self, s):
        self.s = s

    def __getitem__(self, sl):
        return FakeSeq(self.s[sl])

    def __len__(self):
        return len(self.s)


class FakeFasta(dict):
    def __init__(self, mapping):
        super().__init__({k: FakeChrom(v) for k, v in mapping.items()})


@pytest.fixture
def genome():
    rng = np.random.default_rng(5)
    return FakeFasta({"chr1": "".join(rng.choice(list("ACGT"), size=2000))})


def write_peaks(tmp_path, rows):
    p = tmp_path / "peaks.bed.gz"
    with gzip.open(p, "wt") as fh:
        for chrom, s, e, strand in rows:
            fh.write(f"{chrom}\t{s}\t{e}\tpk\t0\t{strand}\n")
    return p


# =======================================================================================
# Peak proximity
# =======================================================================================

class TestNearestPeak:
    @pytest.fixture
    def index(self, tmp_path):
        return assign.peak_index(write_peaks(tmp_path, [
            ("chr1", 1000, 1100, "+"),
            ("chr1", 2000, 2050, "-"),
            ("chr2", 500, 600, "+"),
        ]))

    def test_inside_a_peak_is_distance_zero(self, index):
        assert assign.nearest_peak(index, "chr1", 1050)[0] == 0

    def test_strand_comes_from_the_nearest_peak(self, index):
        assert assign.nearest_peak(index, "chr1", 1050)[1] == "+"
        assert assign.nearest_peak(index, "chr1", 2020)[1] == "-"

    def test_distance_to_the_left_of_a_peak(self, index):
        assert assign.nearest_peak(index, "chr1", 990)[0] == 10

    def test_distance_to_the_right_of_a_peak(self, index):
        assert assign.nearest_peak(index, "chr1", 1105)[0] == 6

    def test_picks_the_closer_of_two_peaks(self, index):
        """Between the two chr1 peaks, 1200 is nearer the first one."""
        assert assign.nearest_peak(index, "chr1", 1200)[1] == "+"
        assert assign.nearest_peak(index, "chr1", 1900)[1] == "-"

    def test_absent_chromosome_is_infinitely_far(self, index):
        d, strand = assign.nearest_peak(index, "chrZ", 100)
        assert d == np.inf and strand is None

    def test_before_the_first_peak(self, index):
        assert np.isfinite(assign.nearest_peak(index, "chr1", 0)[0])

    def test_after_the_last_peak(self, index):
        assert np.isfinite(assign.nearest_peak(index, "chr1", 999999)[0])


class TestAssign:
    @pytest.fixture
    def index(self, tmp_path):
        return assign.peak_index(write_peaks(tmp_path, [("chr1", 1000, 1100, "-")]))

    def variant(self, pos, chrom="chr1"):
        return {"vid": f"{chrom}:{pos}", "chrom": chrom, "pos": pos, "pos_vcf": pos + 1,
                "ref": "A", "alt": "G", "label": 1}

    def test_keeps_a_variant_inside_the_margin(self, index):
        out = assign.assign([self.variant(1110)], index, 25, {"chr1": 2})
        assert len(out) == 1 and out[0]["peak_distance"] == 11

    def test_drops_a_variant_beyond_the_margin(self, index):
        assert assign.assign([self.variant(1200)], index, 25, {"chr1": 2}) == []

    def test_tags_the_peak_strand(self, index):
        assert assign.assign([self.variant(1050)], index, 25, {"chr1": 2})[0]["strand"] == "-"

    def test_tags_the_fold_from_the_chromosome(self, index):
        assert assign.assign([self.variant(1050)], index, 25, {"chr1": 3})[0]["fold"] == 3

    def test_drops_variants_on_excluded_chromosomes(self, index):
        """chrY is excluded from the study, so its fold map entry is absent."""
        assert assign.assign([self.variant(1050)], index, 25, {"chr2": 0}) == []


# =======================================================================================
# Window construction
# =======================================================================================

class TestWindowsFor:
    def test_reference_window_matches_the_genome(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        refs, _ = assign.windows_for(genome, "chr1", pos, ref, "G" if ref != "G" else "C",
                                     SIZE, [0], "+")
        assert refs[0] == genome["chr1"][pos - 50:pos + 51].seq.replace("T", "U")

    def test_alternate_differs_from_reference_at_exactly_one_position(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        alt = "G" if ref != "G" else "C"
        refs, alts = assign.windows_for(genome, "chr1", pos, ref, alt, SIZE, [0], "+")
        assert sum(a != b for a, b in zip(refs[0], alts[0])) == 1

    def test_reference_base_mismatch_is_refused(self, genome):
        """A wrong REF means the 'alternate' window is wrong for reasons unrelated to
        the variant, so it must not be scored."""
        pos = 500
        wrong = "A" if genome["chr1"][pos:pos + 1].seq != "A" else "C"
        assert assign.windows_for(genome, "chr1", pos, wrong, "G", SIZE, [0], "+") \
            == (None, None)

    def test_minus_strand_is_reverse_complemented(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        alt = "G" if ref != "G" else "C"
        plus, _ = assign.windows_for(genome, "chr1", pos, ref, alt, SIZE, [0], "+")
        minus, _ = assign.windows_for(genome, "chr1", pos, ref, alt, SIZE, [0], "-")
        assert minus[0] == revcomp(plus[0].replace("U", "T")).replace("T", "U")

    def test_minus_strand_alt_still_differs_at_one_position(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        alt = "G" if ref != "G" else "C"
        refs, alts = assign.windows_for(genome, "chr1", pos, ref, alt, SIZE, [0], "-")
        assert sum(a != b for a, b in zip(refs[0], alts[0])) == 1

    def test_one_window_per_shift(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        refs, alts = assign.windows_for(genome, "chr1", pos, ref, "G" if ref != "G" else "C",
                                        SIZE, [-40, -20, 0, 20, 40], "+")
        assert len(refs) == len(alts) == 5

    def test_every_window_is_the_configured_length(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        refs, _ = assign.windows_for(genome, "chr1", pos, ref, "G" if ref != "G" else "C",
                                     SIZE, [-40, 0, 40], "+")
        assert all(len(s) == SIZE for s in refs)

    def test_shifted_windows_still_contain_the_variant(self, genome):
        """A shift that moved the variant outside its own window would score nothing."""
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        alt = "G" if ref != "G" else "C"
        refs, alts = assign.windows_for(genome, "chr1", pos, ref, alt, SIZE,
                                        [-40, -20, 0, 20, 40], "+")
        for r, a in zip(refs, alts):
            assert sum(x != y for x, y in zip(r, a)) == 1

    def test_positions_near_the_chromosome_edge_are_dropped(self, genome):
        ref = genome["chr1"][2:3].seq
        assert assign.windows_for(genome, "chr1", 2, ref, "G", SIZE, [0], "+") \
            == (None, None)

    def test_absent_chromosome_returns_nothing(self, genome):
        assert assign.windows_for(genome, "chrZ", 500, "A", "G", SIZE, [0], "+") \
            == (None, None)

    def test_no_thymine_survives(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        refs, alts = assign.windows_for(genome, "chr1", pos, ref,
                                        "G" if ref != "G" else "C", SIZE, [0], "+")
        assert "T" not in refs[0] and "T" not in alts[0]


class TestBuildScoringTable:
    def test_flattens_to_one_row_per_shift(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        v = [{"vid": "v1", "chrom": "chr1", "pos": pos, "pos_vcf": pos + 1, "ref": ref,
              "alt": "G" if ref != "G" else "C", "label": 1, "fold": 0, "strand": "+",
              "peak_distance": 3}]
        rows, dropped = assign.build_scoring_table(v, genome, SIZE, [-20, 0, 20])
        assert len(rows) == 3 and sum(dropped.values()) == 0

    def test_counts_reference_mismatches(self, genome):
        pos = 500
        wrong = "A" if genome["chr1"][pos:pos + 1].seq != "A" else "C"
        v = [{"vid": "v1", "chrom": "chr1", "pos": pos, "pos_vcf": pos + 1, "ref": wrong,
              "alt": "G", "label": 1, "fold": 0, "strand": "+", "peak_distance": 3}]
        rows, dropped = assign.build_scoring_table(v, genome, SIZE, [0])
        assert rows == [] and dropped["ref_mismatch"] == 1

    def test_carries_the_fold_through(self, genome):
        pos = 500
        ref = genome["chr1"][pos:pos + 1].seq
        v = [{"vid": "v1", "chrom": "chr1", "pos": pos, "pos_vcf": pos + 1, "ref": ref,
              "alt": "G" if ref != "G" else "C", "label": 0, "fold": 4, "strand": "+",
              "peak_distance": 0}]
        rows, _ = assign.build_scoring_table(v, genome, SIZE, [0])
        assert rows[0]["fold"] == 4


class TestCollapseDelta:
    def test_max_abs_takes_the_largest_disruption(self):
        keys, vals = assign.collapse_delta(["a"] * 3, [0.1, -0.9, 0.4])
        assert keys == ["a"] and vals == [-0.9]

    def test_keeps_the_sign_of_the_chosen_shift(self):
        _, vals = assign.collapse_delta(["a", "a"], [0.2, -0.5])
        assert vals[0] == -0.5

    def test_separates_variants(self):
        keys, vals = assign.collapse_delta(["a", "b", "a"], [0.1, 0.7, 0.9])
        assert dict(zip(keys, vals)) == {"a": 0.9, "b": 0.7}

    def test_ignores_nan_shifts(self):
        _, vals = assign.collapse_delta(["a"] * 3, [np.nan, 0.3, np.nan])
        assert vals == [0.3]

    def test_all_nan_gives_nan(self):
        _, vals = assign.collapse_delta(["a", "a"], [np.nan, np.nan])
        assert np.isnan(vals[0])

    def test_mean_abs_is_available(self):
        _, vals = assign.collapse_delta(["a"] * 2, [0.2, -0.4], how="mean_abs")
        assert vals[0] == pytest.approx(0.3)

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError):
            assign.collapse_delta(["a"], [0.1], how="median")

    def test_output_is_one_row_per_variant(self):
        rng = np.random.default_rng(3)
        vids = rng.choice(["a", "b", "c", "d"], size=40).tolist()
        keys, vals = assign.collapse_delta(vids, rng.normal(size=40))
        assert len(keys) == len(set(vids)) == len(vals)
