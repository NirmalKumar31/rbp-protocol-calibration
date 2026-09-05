"""Panel membership: which datasets are in the study, for a given negative arm.

These exist because the arm was omitted from the filename for a while, so the dinucleotide
run silently overwrote the GC run's pair counts and the analysis scripts then used those
counts to filter GC data.
"""

import pytest

from rbp.utils import panel


@pytest.fixture
def repo(tmp_path):
    (tmp_path / "config").mkdir()
    for _arm, sub in panel.ARMS.items():
        for cell in ("K562", "HepG2"):
            (tmp_path / sub / cell).mkdir(parents=True)
    return tmp_path


def _dataset(repo, cell, arm, protein):
    d = repo / panel.ARMS[arm] / cell / protein
    d.mkdir(parents=True, exist_ok=True)
    (d / "dataset.tsv").write_text("id\tlabel\n")
    return d


class TestArms:
    def test_the_two_arms_map_to_different_directories(self):
        assert panel.ARMS["gc"] != panel.ARMS["dinuc"]

    def test_unknown_arm_is_rejected_loudly(self):
        with pytest.raises(ValueError, match="unknown negative arm"):
            panel.check_arm("gcc")

    @pytest.mark.parametrize("fn", [panel.panel_path, panel.excluded_path, panel.data_dir])
    def test_every_path_helper_validates_the_arm(self, fn):
        with pytest.raises(ValueError):
            fn("K562", "not-an-arm")


class TestPaths:
    def test_panel_filename_carries_both_cell_and_arm(self, repo):
        p = panel.panel_path("K562", "dinuc", root=repo)
        assert p.name == "panel_final_K562_dinuc.tsv"

    def test_the_two_arms_never_share_a_panel_filename(self, repo):
        """The actual bug: one filename for two arms means the second run wins."""
        a = panel.panel_path("K562", "gc", root=repo)
        b = panel.panel_path("K562", "dinuc", root=repo)
        assert a != b

    def test_the_two_cells_never_share_a_panel_filename(self, repo):
        a = panel.panel_path("K562", "gc", root=repo)
        b = panel.panel_path("HepG2", "gc", root=repo)
        assert a != b

    def test_data_dir_follows_the_arm(self, repo):
        assert panel.data_dir("K562", "gc", root=repo) != \
               panel.data_dir("K562", "dinuc", root=repo)


class TestReadWrite:
    def test_round_trip(self, repo):
        panel.write_panel("K562", "gc", [("QKI", 2189), ("FUS", 900)],
                          [("AARS", 32, "pairs<400")], root=repo)
        assert panel.read_panel("K562", "gc", root=repo) == [("FUS", 900), ("QKI", 2189)]

    def test_missing_panel_reads_as_empty_not_an_error(self, repo):
        assert panel.read_panel("K562", "dinuc", root=repo) == []

    def test_exclusions_record_a_reason(self, repo):
        _, xf = panel.write_panel("K562", "gc", [], [("AARS", 32, "pairs<400")], root=repo)
        assert "pairs<400" in xf.read_text()

    def test_writing_one_arm_does_not_touch_the_other(self, repo):
        """This is the regression. Before the fix, the second write clobbered the first."""
        panel.write_panel("K562", "gc", [("AATF", 1043)], [], root=repo)
        panel.write_panel("K562", "dinuc", [("AATF", 1081)], [], root=repo)
        assert panel.read_panel("K562", "gc", root=repo) == [("AATF", 1043)]
        assert panel.read_panel("K562", "dinuc", root=repo) == [("AATF", 1081)]


class TestDatasets:
    def test_resolves_to_the_arms_own_directory(self, repo):
        _dataset(repo, "K562", "gc", "QKI")
        _dataset(repo, "K562", "dinuc", "QKI")
        panel.write_panel("K562", "gc", [("QKI", 2189)], [], root=repo)
        panel.write_panel("K562", "dinuc", [("QKI", 2189)], [], root=repo)
        gc, _ = panel.datasets(["K562"], "gc", root=repo)
        di, _ = panel.datasets(["K562"], "dinuc", root=repo)
        assert gc["QKI:K562"] != di["QKI:K562"]
        assert panel.ARMS["gc"] in str(gc["QKI:K562"])
        assert panel.ARMS["dinuc"] in str(di["QKI:K562"])

    def test_min_pairs_filters(self, repo):
        for p in ("BIG", "SMALL"):
            _dataset(repo, "K562", "gc", p)
        panel.write_panel("K562", "gc", [("BIG", 900), ("SMALL", 100)], [], root=repo)
        out, _ = panel.datasets(["K562"], "gc", min_pairs=400, root=repo)
        assert set(out) == {"BIG:K562"}

    def test_a_panel_entry_with_no_file_is_reported_not_silently_dropped(self, repo):
        panel.write_panel("K562", "gc", [("GHOST", 900)], [], root=repo)
        out, missing = panel.datasets(["K562"], "gc", root=repo)
        assert out == {}
        assert missing == ["GHOST:K562"]

    def test_keys_are_protein_colon_cell(self, repo):
        _dataset(repo, "HepG2", "gc", "QKI")
        panel.write_panel("HepG2", "gc", [("QKI", 900)], [], root=repo)
        out, _ = panel.datasets(["HepG2"], "gc", root=repo)
        assert list(out) == ["QKI:HepG2"]

    def test_both_cells_combine_without_collision(self, repo):
        for cell in ("K562", "HepG2"):
            _dataset(repo, cell, "gc", "QKI")
            panel.write_panel(cell, "gc", [("QKI", 900)], [], root=repo)
        out, _ = panel.datasets(["K562", "HepG2"], "gc", root=repo)
        assert set(out) == {"QKI:K562", "QKI:HepG2"}
