"""The split-word gate must catch the two defects it was written for, and not fire on prose.

A checker that only ever returns "clean" is indistinguishable from one that looks at nothing,
and both of this project's real instances of this defect were found by hand after the gate that
would have caught them did not exist. So the detector is exercised against the actual broken
text from history, not only against the fixed tree.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import tex_line_breaks as t  # noqa: E402

MAN = ROOT / "manuscript"

# The container image carries src, scripts, config and tests only, so there is no manuscript/
# to check there. The detector tests below build their own files and run everywhere; only the
# whole-tree test needs the real sources, and it says so rather than asserting an empty list.
needs_manuscript = pytest.mark.skipif(
    not MAN.exists(),
    reason="no manuscript/ here: this is the container file set")


def _write(tmp_path, body):
    p = tmp_path / "s.tex"
    p.write_text(body)
    return p


def test_it_catches_the_introduction_defect(tmp_path):
    p = _write(tmp_path, "The magnitude of this contribution varied several-\nfold across it.\n")
    assert [o[0] for o in t.offenders(p)] == [1]


def test_it_catches_the_results_defect(tmp_path):
    p = _write(tmp_path, "the first is a protocol-\nspecific offset after a curve.\n")
    assert [o[0] for o in t.offenders(p)] == [1]


def test_the_joined_form_is_clean(tmp_path):
    p = _write(tmp_path, "The magnitude of this contribution varied\nseveral-fold across it.\n")
    assert t.offenders(p) == []


def test_a_minus_sign_ending_a_display_equation_is_not_a_split_word(tmp_path):
    """Methods ends a line with `) -` inside an equation; a space precedes the hyphen."""
    p = _write(tmp_path, "\\mathrm{AUROC}(\\text{composition} + \\text{score}) -\n"
                         "\\mathrm{AUROC}(\\text{composition alone}).\n")
    assert t.offenders(p) == []


def test_a_comment_rule_is_not_a_split_word(tmp_path):
    p = _write(tmp_path, "% ------------------------------------------\nsome prose follows.\n")
    assert t.offenders(p) == []


def test_a_dash_before_a_capital_or_a_macro_is_left_alone(tmp_path):
    p = _write(tmp_path, "the well-known result-\n\\citet{horlacher2023} showed it.\n")
    assert t.offenders(p) == []


@needs_manuscript
def test_it_catches_an_unlabelled_subsection(tmp_path):
    p = _write(tmp_path, "\\subsection{Standard reparameterisations}\n\nSome prose.\n")
    assert [o[0] for o in t.unlabelled(p)] == [1]


def test_a_labelled_subsection_is_clean(tmp_path):
    p = _write(tmp_path, "\\subsection{Standard reparameterisations}\n\\label{sec:r}\n\nProse.\n")
    assert t.unlabelled(p) == []


def test_a_title_spanning_two_lines_still_finds_its_label(tmp_path):
    """sec:families' title spans two source lines; a naive next-line check misreports it."""
    p = _write(tmp_path, "\\subsection{Composition explains more variation than protocol\n"
                         "labels for a $k$-mer model}\n\\label{sec:families}\n\nProse.\n")
    assert t.unlabelled(p) == []


@needs_manuscript
def test_every_shipped_subsection_is_labelled():
    srcs = sorted(MAN.glob("*.tex")) + sorted((MAN / "sections").glob("*.tex"))
    assert srcs, "no manuscript sources found; this test checked nothing"
    assert {p.name: t.unlabelled(p) for p in srcs if t.unlabelled(p)} == {}


@needs_manuscript
def test_the_shipped_manuscript_is_clean():
    srcs = sorted(MAN.glob("*.tex")) + sorted((MAN / "sections").glob("*.tex"))
    assert srcs, "no manuscript sources found; this test checked nothing"
    assert {p.name: t.offenders(p) for p in srcs if t.offenders(p)} == {}
