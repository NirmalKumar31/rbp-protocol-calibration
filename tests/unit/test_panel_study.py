"""The study definition: which datasets, which arm, which cell lines.

These pin the fix for a bug that was live for four days. Four scripts each answered "what
is the panel?" independently, all four read `config/panel_final.tsv` -- a 17-protein file
keyed by neither cell line nor negative arm, left over from the development panel -- and
all four then loaded from the GC directory. A sweep driven by that would have trained 85
runs on the wrong data in one cell line and reported them as the study.

The permanent fix is that the question has exactly one answer, `panel.study(cfg)`, and the
file that caused it no longer exists. Both are asserted here.
"""

import ast
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rbp.utils import config as cfgmod  # noqa: E402
from rbp.utils import panel as panelmod  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIRS = ("src", "scripts")


@pytest.fixture(scope="module")
def cfg():
    return cfgmod.load()


class TestTheStaleFileIsGone:
    """A deleted file comes back the moment one script still writes it."""

    @pytest.mark.parametrize("name", ["panel_final.tsv", "panel_excluded.tsv"])
    def test_file_does_not_exist(self, name):
        assert not (ROOT / "config" / name).exists(), (
            f"config/{name} is the unkeyed panel. Anything reading it gets one cell line "
            "and one arm by accident. Use rbp.utils.panel instead.")

    @pytest.mark.parametrize("name", ["panel_final.tsv", "panel_excluded.tsv"])
    def test_no_source_file_reads_it(self, name):
        """Prose may mention it as history. Code may not name it as a path.

        Docstrings are excluded deliberately: the modules that used to read this file
        explain why they no longer do, and that explanation is worth keeping. What is
        forbidden is a string literal that some line of code could open.
        """
        offenders = []
        for d in SOURCE_DIRS:
            for p in (ROOT / d).rglob("*.py"):
                tree = ast.parse(p.read_text(), filename=str(p))
                docs = {id(ast.get_docstring(n, clean=False))
                        for n in ast.walk(tree)
                        if isinstance(n, (ast.Module, ast.ClassDef, ast.FunctionDef,
                                          ast.AsyncFunctionDef))}
                for node in ast.walk(tree):
                    if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                            and name in node.value and id(node.value) not in docs):
                        offenders.append(f"{p.relative_to(ROOT)}:{node.lineno}")
        assert not offenders, f"still referencing config/{name}: {offenders}"


class TestArmAndCells:
    def test_primary_arm_comes_from_config(self, cfg):
        assert panelmod.arm_of(cfg) == cfg["negatives"]["primary_arm"]

    def test_explicit_arm_wins(self, cfg):
        assert panelmod.arm_of(cfg, "gc") == "gc"

    def test_unknown_arm_is_rejected(self, cfg):
        with pytest.raises(ValueError):
            panelmod.arm_of(cfg, "shuffled")

    def test_both_cell_lines_are_in_the_study(self, cfg):
        assert set(panelmod.cells_of(cfg)) == {"K562", "HepG2"}

    def test_singular_cell_line_is_a_subset_of_the_plural(self, cfg):
        """encode.cell_line still exists for download and census, which run one at a
        time. It must not drift out of the study's set."""
        assert cfg.encode["cell_line"] in panelmod.cells_of(cfg)

    def test_arm_binds_both_the_panel_file_and_the_data_directory(self):
        """The original bug in one assertion: these two were chosen independently."""
        for arm in panelmod.ARMS:
            assert arm in str(panelmod.panel_path("K562", arm))
            assert panelmod.ARMS[arm] in str(panelmod.data_dir("K562", arm))


class TestStudy:
    def test_keys_are_protein_and_cell(self, cfg):
        paths, _ = panelmod.study(cfg, require_files=False)
        assert paths, "no datasets found; the panel files are missing"
        for key in paths:
            protein, cell = key.split(":")
            assert protein and cell in panelmod.cells_of(cfg)

    def test_every_path_is_in_the_primary_arm_directory(self, cfg):
        paths, _ = panelmod.study(cfg, require_files=False)
        arm_dir = panelmod.ARMS[panelmod.arm_of(cfg)]
        assert all(arm_dir in str(p) for p in paths.values())

    def test_the_two_arms_select_different_directories(self, cfg):
        gc, _ = panelmod.study(cfg, "gc", require_files=False)
        dinuc, _ = panelmod.study(cfg, "dinuc", require_files=False)
        common = set(gc) & set(dinuc)
        assert common, "the arms share no datasets, which cannot be right"
        assert all(gc[k] != dinuc[k] for k in common)

    def test_min_pairs_is_enforced(self, cfg):
        paths, _ = panelmod.study(cfg, require_files=False)
        floor = cfg.cv["min_pairs"]
        arm = panelmod.arm_of(cfg)
        for cell in panelmod.cells_of(cfg):
            for protein, pairs in panelmod.read_panel(cell, arm):
                if pairs < floor:
                    assert f"{protein}:{cell}" not in paths

    def test_the_panel_is_the_size_the_paper_claims(self, cfg):
        """189 datasets across two cell lines. Drafts said 187 for a while; this is the
        assertion that stops the number drifting again."""
        paths, _ = panelmod.study(cfg, require_files=False)
        assert len(paths) == 189

    def test_the_sweep_is_the_size_the_plan_claims(self, cfg):
        """189 datasets x 5 models x 5 folds. If any of the three moves, the cost
        estimate and the overnight schedule move with it."""
        paths, _ = panelmod.study(cfg, require_files=False)
        assert len(paths) * len(cfg["models"]) * cfg.cv["k"] == 4725

    def test_require_files_is_what_separates_intent_from_disk(self, cfg):
        """The container has config/ but not data/. Both answers are legitimate; what
        would not be is the same call meaning different things in the two places."""
        intended, _ = panelmod.study(cfg, require_files=False)
        present, missing = panelmod.study(cfg)
        assert set(present) <= set(intended)
        assert set(intended) - set(present) <= set(missing)
